"""
PIPELINE: Raw Reads Quality Control
DESCRIPTION:
This Snakemake workflow performs quality control on raw DNA-seq reads, including FastQC for quality assessment, Fastp for trimming and filtering, and MultiQC for aggregating QC reports.
"""

# =============================================================================
# GLOBAL CONFIGURATION
# =============================================================================
# Dynamically find FASTQ files using sample ID and read type ('R1' or 'R2')
def get_raw(wildcards, read):
    files = list(Path(RAW_DIR).glob(f"{wildcards.sample}_*_{read}*.fastq.gz"))
    if not files:
        files = list(Path(RAW_DIR).glob(f"{wildcards.sample}_{read}.fastq.gz"))
    return str(files[0]) if files else []

# =============================================================================
# PIPELINE RULES
# =============================================================================
# FastQC for raw reads
rule fastqc:
    input:
        lambda wildcards: get_raw(wildcards, wildcards.read)
    output:        
        html = f"{RES_DIR}/01_RawReads_QC/01_FastQC/{{sample}}_{{read}}_fastqc.html",
        zip  = f"{RES_DIR}/01_RawReads_QC/01_FastQC/{{sample}}_{{read}}_fastqc.zip"
    threads: 4
    resources:
        mem_mb = 4096
    wrapper:
        "v5.7.0/bio/fastqc"

# Fastp for trimming and quality filtering
rule fastp:
    input:
        r1 = lambda wildcards: get_raw(wildcards, "R1"),
        r2 = lambda wildcards: get_raw(wildcards, "R2")
    params:
        len_min = 35
    output:
        r1 = f"{QC_DIR}/{{sample}}_trimmed_R1.fastq.gz",
        r2 = f"{QC_DIR}/{{sample}}_trimmed_R2.fastq.gz",
        html = f"{QC_DIR}/{{sample}}_fastp_report.html",
        json = f"{QC_DIR}/{{sample}}_fastp_report.json"
    threads: 8
    envmodules:
        "apps/fastp"
    shell:
        """
        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --length_required {params.len_min} \
            --cut_right \
            --thread {threads} \
            --html {output.html} \
            --json {output.json}
        """

# MultiQC to aggregate FastQC and Fastp reports
rule raw_reads_multiqc:
    input:
        fastqc_files = expand(f"{RES_DIR}/01_RawReads_QC/01_FastQC/{{sample}}_{{read}}_fastqc.zip", sample=SAMPLES, read=["R1", "R2"]),
        fastp_html   = expand(f"{RES_DIR}/01_RawReads_QC/02_Fastp/{{sample}}_fastp_report.html", sample=SAMPLES),        
        fastp_json   = expand(f"{RES_DIR}/01_RawReads_QC/02_Fastp/{{sample}}_fastp_report.json", sample=SAMPLES)
    output:
        html = f"{RES_DIR}/01_RawReads_QC/03_MultiQC/MultiQC_report.html"
    threads: 8
    wrapper:
        "v5.7.0/bio/multiqc"