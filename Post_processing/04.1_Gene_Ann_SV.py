#!/usr/bin/env python3
import argparse
import pysam
import pandas as pd
import pyranges as pr
import re

def parse_bnd_alt(alt_string):
    """Extracts chrom and pos from a BND string (e.g., N]9:20441690])"""
    match = re.search(r'([^\[\]:]+):(\d+)', alt_string)
    if match:
        return match.group(1), int(match.group(2))
    return None, None

def get_annotated_df(vcf_rows, gtf_pr):
    """Cross-references VCF coordinates with the GTF"""
    if not vcf_rows:
        return pd.DataFrame()

    df_vcf = pd.DataFrame(vcf_rows)
    # PyRanges requires Start and End. For an exact point: [pos-1, pos]
    pr_vcf = pr.PyRanges(df_vcf)

    # Intersection with the GTF
    annotated = pr_vcf.join(gtf_pr, how="left")
    return annotated.df

def main(vcf_file, gtf_file, min_support, output):

    # -------------------------------
    # 1. Load GTF
    # -------------------------------
    print(f"Loading GTF...")
    gtf = pr.read_gtf(gtf_file)
    gtf_genes = gtf[gtf.Feature == 'gene'][['Chromosome', 'Start', 'End', 'gene_name']]
    # Normalize GTF chromosomes to pure strings (1, 2, X...)
    gtf_genes.Chromosome = gtf_genes.Chromosome.astype(str).str.replace("chr", "")

    # -------------------------------
    # 2. Process VCF
    # -------------------------------
    vcf_in = pysam.VariantFile(vcf_file)
    bnd_points = []
    other_data = []

    print("Reading VCF...")
    for rec in vcf_in:
        # FILTER 1: Only the first version of the variant (avoid Lumpy's _2)
        if rec.id.endswith("_2"):
            continue

        # FILTER 2: Minimum support
        pe = 0
        sr = 0
        ## Check Manta Format (PR/SR in sample fields)
        sample_data = rec.samples[0] # Pysam allows indexing by position too
        if "PR" in sample_data:
            # Manta PR is (ref_support, alt_support)
            vals = sample_data["PR"]
            pe = vals[1] if isinstance(vals, (tuple, list)) else 0
        if "SR" in sample_data:
            vals = sample_data["SR"]
            sr = vals[1] if isinstance(vals, (tuple, list)) else vals

        ## Check Lumpy Format (PE/SR in info fields)
        if "PE" in rec.info:
            val = rec.info["PE"]
            pe = sum(val) if isinstance(val, (tuple, list)) else val
        if "SR" in rec.info:
            val = rec.info["SR"]
            sr = sum(val) if isinstance(val, (tuple, list)) else val

        total_support = pe + sr
        if total_support < min_support:
            continue

        chrom_a = str(rec.chrom).replace("chr", "")
        pos_a = rec.pos  # Original position (1-based)
        svtype = rec.info.get("SVTYPE", "NA")

        if svtype == "BND":
            alt_str = str(rec.alts[0])
            chrom_b, pos_b = parse_bnd_alt(alt_str)
            if chrom_b:
                chrom_b = chrom_b.replace("chr", "")
                coord_str = f"{chrom_a}:{pos_a}-{chrom_b}:{pos_b}"

                # Register both points so PyRanges searches for the gene at each location
                bnd_points.append({"Chromosome": chrom_a, "Start": pos_a-1, "End": pos_a, "VAR_ID": rec.id, "Side": "A", "Coord": coord_str, "Support": total_support})
                bnd_points.append({"Chromosome": chrom_b, "Start": pos_b-1, "End": pos_b, "VAR_ID": rec.id, "Side": "B", "Coord": coord_str, "Support": total_support})
        else:
            coord_str = f"{chrom_a}:{pos_a}"
            other_data.append({
                "Chromosome": chrom_a, "Start": pos_a-1, "End": pos_a,
                "VAR_ID": rec.id, "Support": total_support, "SVType": svtype, "Coord": coord_str
            })

    # -------------------------------
    # 3. Annotate
    # -------------------------------
    print("Annotating genes...")
    res_bnd = get_annotated_df(bnd_points, gtf_genes)
    res_other = get_annotated_df(other_data, gtf_genes)

    final_rows = []

    # Reassemble BND fusions
    if not res_bnd.empty:
        for var_id, group in res_bnd.groupby("VAR_ID"):
            # Obtain gene names (or Intergenic if no hit)
            gene_a = group[group.Side == "A"]["gene_name"].fillna("Intergenic").iloc[0]
            gene_b = group[group.Side == "B"]["gene_name"].fillna("Intergenic").iloc[0]

            # Save only if at least one side overlaps a gene
            if gene_a != "Intergenic" or gene_b != "Intergenic":
                final_rows.append({
                    "Gene": f"{gene_a}::{gene_b}",
                    "SVType": "BND",
                    "Support": group["Support"].iloc[0],
                    "Coordinates": group["Coord"].iloc[0]
                })

    # Reassemble other SVs (DEL, DUP, INV)
    if not res_other.empty:
        for _, row in res_other.iterrows():
            raw_gene = str(row.get("gene_name", "Intergenic"))

            # Clean up null values or -1 results from the join
            if raw_gene in ["nan", "-1", "None", ""]:
                gene = "Intergenic"
            else:
                gene = raw_gene
            if gene != "Intergenic":
                final_rows.append({
                    "Gene": gene,
                    "SVType": row["SVType"],
                    "Support": row["Support"],
                    "Coordinates": row["Coord"]
                })

    # -------------------------------
    # 4. Output
    # -------------------------------
    df = pd.DataFrame(final_rows)
    if not df.empty:
        # Sort by support from highest to lowest
        df = df.sort_values(by="Support", ascending=False)
        df.to_csv(output, sep="\t", index=False)
        print(f"Table saved to {output}. Found {len(df)} variants.")
    else:
        print("No variants found meeting the criteria.")

if __name__ == "__main__":
    # Parameters from Snakemake
    if "snakemake" in locals() or "snakemake" in globals():
        main(
        vcf_file = snakemake.input.vcf,
        gtf_file = snakemake.input.gtf,
        min_support = snakemake.params.get("min_support", 10),
        output = snakemake.output.tsv
        )
    # Manual parameter parsing
    else:
        parser = argparse.ArgumentParser(description="Annotate LUMPY VCF")
        parser.add_argument("--vcf_file", required=True)
        parser.add_argument("--gtf_file", required=True)
        parser.add_argument("--min_support", type=int, default=10)
        parser.add_argument("--output", required=True)
        args = parser.parse_args()

        main(args.vcf_file, args.gtf_file, args.min_support, args.output)