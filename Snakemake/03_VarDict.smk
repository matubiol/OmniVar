"""
PIPELINE: VarDict Small Variation Analysis
DESCRIPTION:
This Snakemake workflow performs somatic variant calling using VarDict, followed by annotation with SnpEff and ClinVar.
It includes configurable parameters for VarDict to allow users to adjust sensitivity and specificity based on their needs.
"""

# =============================================================================
# PIPELINE RULES
# =============================================================================
# Variant calling
rule vardict:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_recal.bam",
        ref = REFS['genome'],
        gnomad = REFS['gnomad'],
        pon = REFS['pon'],
        intervals = config['proj']['intervals']
    output:
        vcf = f"{VAR_DIR}/01_Variants/{{sample}}.vcf"
    log:
        f"{VAR_DIR}/01_Variants/Logs/{{sample}}.log",
    benchmark:
        f"{VAR_DIR}/01_Variants/Logs/{{sample}}.benchmark"
    conda:
        ENVS['vardict']
    resources:
        mem_mb = 48000
    threads: 8
    shell:
        """
        vardict-java \
            -G {input.ref} \
            -f 0.01 \
            -N {wildcards.sample} \
            -b {input.bam} \
            --verbose -th {threads} \
            -c 1 -S 2 -E 3 -g 4 {input.intervals} 2> {log} | \
            teststrandbias.R 2>> {log} | \
            var2vcf_valid.pl \
            -A -N {wildcards.sample} -E -f 0.01 > {output.vcf} 2>> {log}
        """

# Annotate variants with SnpEff
rule annotate_vardict:
    input:
        vcf = f"{VAR_DIR}/01_Variants/{{sample}}.vcf",
        clinvar = REFS['clinvar'],
        intervals = config['proj']['intervals']
    params:
        db = config["snpeff_db"],
        java_opts = "-Xmx16G"
    output:
        vcf_eff = temp(f"{VAR_DIR}/02_Annotated/{{sample}}_ann.vcf"),
        vcf_clinvar = temp(f"{VAR_DIR}/02_Annotated/{{sample}}_ann_clinvar.vcf"),
        tsv = temp(f"{VAR_DIR}/02_Annotated/{{sample}}_ann_clinvar.tsv")
    envmodules:
        "apps/java"
    resources:
        mem_mb = 18000
    shell:
        """
        set -o pipefail

        # Annotate variants with SnpEff
        snpeff \
            {params.java_opts} \
            -verbose \
            -noStats \
            -filterInterval {input.intervals} \
            {params.db} \
            {input.vcf} | \
        snpsift varType - | \
        tee {output.vcf_eff} | \
        snpsift annotate \
            {input.clinvar} \
            -info CLNSIG \
            - | \
        tee {output.vcf_clinvar} | \
        snpsift extractFields \
            - \
            "CHROM" "POS" "ANN[0].GENE" "ANN[0].FEATUREID" "ANN[0].RANK" "ANN[0].HGVS_C" \
            "ANN[0].HGVS_P" "GEN[0].AD[1]" "GEN[0].AF" "ANN[0].EFFECT" "ANN[0].IMPACT" CLNSIG "VARTYPE" \
            > {output.tsv}
        """
    
# Format SnpEff vcf to create a report
rule report_vardict:
    input:
        vars = f"{VAR_DIR}/02_Annotated/{{sample}}_ann_clinvar.tsv"
    params:
        min_depth = config.get("vc_min_depth", 8)
    output:
        report = f"{VAR_DIR}/02_Annotated/{{sample}}_report.tsv"
    conda:
        ENVS['python']
    script:
        SRC['snpeff_report']