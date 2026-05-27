import pandas as pd
import sys
import argparse

def classify_acmg(row):
    """
    Classification logic updated to avoid filtering critical FLT3/splicing variants.
    """
    # 1. ClinVar Priority
    clinvar = str(row['CLNSIG']).lower()
    if clinvar and clinvar != 'nan' and clinvar != '.':
        if 'pathogenic' in clinvar: return "Patogénica"
        if 'likely_pathogenic' in clinvar: return "Probablemente patogénica"
        if 'uncertain_significance' in clinvar or 'conflicting' in clinvar: return "Significado incierto"
        if 'benign' in clinvar: return "Benigna"

    # 2. SnpEff Impact & Effect Logic
    effect = str(row['ANN[0].EFFECT']).lower()
    impact = str(row['ANN[0].IMPACT']).lower()
    gene = str(row['ANN[0].GENE']).upper()

    # High Impact: Frameshift, stop, splice_site (canonical)
    if any(x in effect for x in ['frameshift', 'stop_gained', 'splice_site_variant']):
        return "Patogénica"

    # RESCUE: Splice regions, duplications or insertions (Critical for FLT3/ITD)
    # If it's a splice region or a complex variant in a key gene, we mark as VUS to keep it
    if 'splice_region' in effect or 'insertion' in effect or 'duplication' in effect:
        return "Significado incierto"

    # Missense
    if 'missense_variant' in effect or impact == 'moderate':
        return "Significado incierto"

    return "Benigna"

def main(input, min_depth, output):

    # Read the annotated TSV file (already parsed from VCF)
    df = pd.read_csv(input, sep='\t')

    # --- FILTERS ---
    # 1. Technical Filter: Alt Depth (AD[1]) >= min_depth
    df['GEN[0].AD[1]'] = pd.to_numeric(df['GEN[0].AD[1]'], errors='coerce')
    df = df[df['GEN[0].AD[1]'] >= min_depth]

    # 2. Protein Effect Filter:
    # Remove synonymous, but keep 'splice' or HIGH/MODERATE impact
    is_synonymous = df['ANN[0].EFFECT'].str.contains('synonymous', case=False, na=False)
    is_splice = df['ANN[0].EFFECT'].str.contains('splice', case=False, na=False)

    # Filter: remove if synonymous AND NOT splice-related
    df = df[~(is_synonymous & ~is_splice)]

    # Remove MODIFIER impact unless it's a splice region
    df = df[(df['ANN[0].IMPACT'].str.upper() != 'MODIFIER') | (is_splice)]

    # --- CLASSIFICATION ---
    df['Clasificación'] = df.apply(classify_acmg, axis=1)

    # Remove if it's strictly classified as "Benigna"
    df = df[df['Clasificación'] != "Benigna"]

    # --- FORMATTING ---
    if 'ANN[0].RANK' in df.columns:
        df['ANN[0].RANK'] = df['ANN[0].RANK'].astype(str).str.split('/').str[0]

    rename_map = {
        'ANN[0].GENE': 'Gen',
        'ANN[0].FEATUREID': 'Transcrito Ref.',
        'ANN[0].RANK': 'Exón',
        'ANN[0].HGVS_C': 'c.DNA',
        'ANN[0].HGVS_P': 'Proteína',
        'GEN[0].AD[1]': 'Alt_depth',
        'GEN[0].AF': 'VAF',
        'VARTYPE': 'Tipo'
    }
    df = df.rename(columns=rename_map)

    expected_columns = ['Gen', 'Transcrito Ref.', 'Exón', 'c.DNA', 'Proteína', 'Alt_depth', 'VAF', 'Tipo', 'Clasificación']
    final_cols = [c for c in expected_columns if c in df.columns]

    df_final = df[final_cols]

    df_final.to_csv(output, sep='\t', index=False)
    print(f"Report generated: {output} ({len(df_final)} variants)", file=sys.stderr)

if __name__ == "__main__":
    # Parameters from Snakemake
    if "snakemake" in locals() or "snakemake" in globals():
        main(
        input = snakemake.input.vars,
        min_depth = snakemake.params.get("min_depth", 8),
        output = snakemake.output.report
        )
    # Manual parameter parsing
    else:
        parser = argparse.ArgumentParser(description="Format SnpEff VCF")
        parser.add_argument("--vars", required=True)
        parser.add_argument("--min_depth", type=int, default=8)
        parser.add_argument("--report", required=True)
        args = parser.parse_args()

        main(args.vars, args.min_depth, args.report)