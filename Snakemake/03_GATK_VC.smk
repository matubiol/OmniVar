"""
PIPELINE: GATK Small Variation Analysis
DESCRIPTION:
This Snakemake workflow performs somatic variant calling using GATK Mutect2, followed by annotation with SnpEff and ClinVar.
It includes configurable parameters for Mutect2 to allow users to adjust sensitivity and specificity based on their needs.
"""

# =============================================================================
# PIPELINE RULES
# =============================================================================
# Dynamically resolve Mutect2 execution parameters based on configuration modes.
def get_mutect_params(wildcards):    
    mode = config.get("mutect_mode", "moderate")
    common = "--pair-hmm-implementation FASTEST_AVAILABLE --genotype-pon-sites true"
    
    # Define extra flags as a clean continuous string (no commas at line ends)
    extra = (
        "--use-pdhmm true --use-pdhmm-overlap-optimization true "
        "--disable-read-filter SoftClippedReadFilter "
        "--disable-adaptive-pruning true --max-mnp-distance 10"
    )

    if mode == "default":
        return common

    elif mode == "moderate":
        moderated_flags = [
            common, 
            extra,
            "--initial-tumor-lod 0.5",
            "--tumor-lod-to-emit 0.5",
            "--kmer-size 10",
            "--kmer-size 25"
        ]
        # Join all elements into a single clean space-separated string
        return " ".join(moderated_flags)

    elif mode == "relaxed":
        relaxed_flags = [
            common, 
            extra,
            "--initial-tumor-lod 0.05",
            "--tumor-lod-to-emit 0.05",
            "--kmer-size 10",
            "--kmer-size 15",
            "--kmer-size 25",
            "--active-probability-threshold 0.002"
        ]
        return " ".join(relaxed_flags)
        
    return common

# Run Mutect2 for somatic variant calling
rule mutect2:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_recal.bam",
        ref = REFS['genome'],
        gnomad = REFS['gnomad'],
        pon = REFS['pon'],
        intervals = config['proj']['intervals']
    params:
        mode_args = get_mutect_params,
        java_opts = "-Xmx45G"
    output:
        vcf_temp = temp(f"{GATK_DIR}/01_Variants/{{sample}}_temp.vcf"),
        f1r2 = f"{GATK_DIR}/01_Variants/{{sample}}_f1r2.tar.gz",
        vcf_norm = protected(f"{GATK_DIR}/01_Variants/{{sample}}.vcf")
    log:
        f"{GATK_DIR}/01_Variants/Logs/{{sample}}.log",
    benchmark:
        f"{GATK_DIR}/01_Variants/Logs/{{sample}}.benchmark"
    container:
        APPS['gatk']
    resources:
        mem_mb = 48000
    threads: 8
    shell:
        """
        # Run Mutect2 with specified parameters
        gatk Mutect2 \
            --java-options "{params.java_opts}" \
            --native-pair-hmm-threads {threads} \
            -R {input.ref} \
            -I {input.bam} \
            -tumor {wildcards.sample} \
            -L {input.intervals} \
            --germline-resource {input.gnomad} \
            --panel-of-normals {input.pon} \
            -O {output.vcf_temp} \
            --f1r2-tar-gz {output.f1r2} \
            {params.mode_args} \
            --verbosity INFO &> {log}

        # Normalize VCF to left-align and split multi-allelics
        gatk LeftAlignAndTrimVariants \
            -R {input.ref} \
            -V {output.vcf_temp} \
            -O {output.vcf_norm} \
            --split-multi-allelics true
        """

# Annotate variants with SnpEff
rule annotate_gatk:
    input:
        vcf = f"{GATK_DIR}/01_Variants/{{sample}}.vcf",
        clinvar = REFS['clinvar'],
        intervals = config['proj']['intervals']
    params:
        db = config['snpeff_db'],
        java_opts = "-Xmx16G"
    output:
        vcf_eff = temp(f"{GATK_DIR}/02_Annotated/{{sample}}_ann.vcf"),
        vcf_clinvar = temp(f"{GATK_DIR}/02_Annotated/{{sample}}_ann_clinvar.vcf"),
        tsv = temp(f"{GATK_DIR}/02_Annotated/{{sample}}_ann_clinvar.tsv")
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
rule report_gatk:
    input:
        vars = f"{GATK_DIR}/02_Annotated/{{sample}}_ann_clinvar.tsv"
    params:
        min_depth = config.get("vc_min_depth", 8)
    output:
        report = f"{GATK_DIR}/02_Annotated/{{sample}}_report.tsv"
    conda:
        ENVS['python']
    script:
        SRC['snpeff_report']