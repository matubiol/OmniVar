#!/usr/bin/env python3
import pandas as pd
import pyranges as pr
import argparse
import sys
import numpy as np

def main(call, bed, cyto, output, min_probes=3, format_type='gatk', sample_name=None):

    # 1. Load and Standardize Input
    if format_type == 'gatk':
        # GATK Logic: Original approach that worked
        df = pd.read_csv(call, sep=None, engine='python')
        df['Chromosome'] = df['Chromosome'].astype(str).str.replace('chr', '', case=False)
        df = df[df['Num_Probes'] >= min_probes].copy()

        # --- GENE ANNOTATION (GATK Specific) ---
        if bed:
            bed_df = pd.read_csv(bed, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Gene'])
            bed_df['Chromosome'] = bed_df['Chromosome'].astype(str).str.replace('chr', '', case=False)

            pr_df = pr.PyRanges(df)
            pr_bed = pr.PyRanges(bed_df)
            genes_joined = pr_df.join(pr_bed).df

            genes_grouped = (
                genes_joined.groupby(['Sample', 'Chromosome', 'Start', 'End'], observed=True)['Gene']
                .apply(lambda x: sorted(list(x.unique())) if x.notnull().any() else [])
                .reset_index()
            )
            df = pd.merge(df, genes_grouped, on=['Sample', 'Chromosome', 'Start', 'End'], how='left')
        else:
            df['Gene'] = [[] for _ in range(len(df))]

    elif format_type == 'cnvkit':
        # CNVkit Logic
        df = pd.read_csv(call, sep='\t')
        df = df.rename(columns={'chromosome': 'Chromosome', 'start': 'Start', 'end': 'End', 'probes': 'Num_Probes', 'gene': 'Gene'})
        df['Chromosome'] = df['Chromosome'].astype(str).str.replace('chr', '', case=False)
        df = df[df['Num_Probes'] >= min_probes].copy()

        # Inject Sample name if provided
        if sample_name:
            df['Sample'] = sample_name
        elif 'sample' not in df.columns:
            df['Sample'] = 'Unknown'

        # Derive Call (Gain/Loss/Neutral)
        df['Call'] = np.where(df['cn'] > 2, 'Gain', np.where(df['cn'] < 2, 'Loss', 'Neutral'))

        # Convert comma-separated genes string to list
        df['Gene'] = df['Gene'].apply(lambda x: sorted(str(x).split(',')) if pd.notnull(x) and x != '' else [])

    # Ensure Gene is a list for the final merge (safeguard)
    df['Gene'] = df['Gene'].apply(lambda x: x if isinstance(x, list) else [])

    # 2. Cytoband Annotation (Common to both)
    cyto_df = pd.read_csv(cyto, sep='\t', header=None,
                          names=['Chromosome', 'Start', 'End', 'Band', 'Stain'])
    cyto_df['Chromosome'] = cyto_df['Chromosome'].astype(str).str.replace('chr', '', case=False)

    pr_df = pr.PyRanges(df)
    pr_cyto = pr.PyRanges(cyto_df)
    cyto_joined = pr_df.join(pr_cyto).df

    def get_arm(bands):
        arms = sorted(list(set([str(b)[0] for b in bands if pd.notnull(b) and str(b) != ''])))
        return "".join(arms)

    cyto_grouped = (
        cyto_joined.groupby(['Sample', 'Chromosome', 'Start', 'End'], observed=True)
        .agg(
            Cytoband=('Band', lambda x: f"{x.iloc[0]}-{x.iloc[-1]}" if len(x) > 1 else str(x.iloc[0])),
            Arm=('Band', get_arm)
        )
        .reset_index()
    )

    # 3. Final Merge
    final_df = pd.merge(df, cyto_grouped, on=['Sample', 'Chromosome', 'Start', 'End'], how='left')
    final_df['Arm'] = final_df['Arm'].fillna('')
    final_df['Cytoband'] = final_df['Cytoband'].fillna('')

    # Custom Sorting (1-22, X, Y)
    def chrom_key(chrom):
        val = str(chrom).upper().replace('X', '23').replace('Y', '24').replace('MT', '25').replace('M', '25')
        try:
            return int(val)
        except ValueError:
            return 99

    final_df['sort_key'] = final_df['Chromosome'].apply(chrom_key)
    final_df = final_df.sort_values(by=['Sample', 'sort_key', 'Start']).drop(columns=['sort_key'])

    # 4. Save Output
    final_df.to_csv(output, sep='\t', index=False)
    print(f"Done! Output saved to: {output}")

if __name__ == "__main__":
    # Parameters from Snakemake
    if "snakemake" in locals() or "snakemake" in globals():
        main(
            call = snakemake.input.call,
            bed = snakemake.input.get("bed"),
            cyto = snakemake.input.cyto,
            output = snakemake.output.tsv,
            min_probes = snakemake.params.get("min_probes", 3),
            format_type = snakemake.params.format,
            sample_name = snakemake.params.sample
        )
    # Manual parameter parsing
    else:
        parser = argparse.ArgumentParser(description="Annotate CNV-call with BED (optional) + cytobands. Supports GATK and CNVkit.")
        parser.add_argument("--call", required=True, help="Input CNV file")
        parser.add_argument("--format", choices=['gatk', 'cnvkit'], required=True, help="Input format")
        parser.add_argument("--bed", help="BED file for gene annotation (required for GATK)")
        parser.add_argument("--cyto", required=True, help="Cytoband file")
        parser.add_argument("--min_probes", type=int, default=3)
        parser.add_argument("--sample", help="Sample name (required for CNVkit)")
        parser.add_argument("--output", required=True)

        args = parser.parse_args()

        main(args.call, args.bed, args.cyto, args.output, args.min_probes, args.format, args.sample)
