"""
PIPELINE: Alignment and QC of DNA-seq data
DESCRIPTION:
This Snakemake workflow performs read alignment, duplicate marking, base quality score recalibration, and collects alignment QC metrics for DNA-seq data.
It uses BWA-MEM2 for alignment and GATK for post-alignment processing.
"""

# =============================================================================
# PIPELINE RULES
# =============================================================================
# Reads alignment
rule bwa_align:
    input:
        ref = REFS['genome'],
        r1 = f"{QC_DIR}/{{sample}}_trimmed_R1.fastq.gz",
        r2 = f"{QC_DIR}/{{sample}}_trimmed_R2.fastq.gz"
    output:
        temp(f"{RES_DIR}/02_Alignment/Alignments/01_RawBAM/{{sample}}.bam")
    log:
        f"{RES_DIR}/02_Alignment/Logs/{{sample}}.log",
    benchmark:
        f"{RES_DIR}/02_Alignment/Logs/{{sample}}.benchmark"
    threads: 24
    envmodules:
        "apps/samtools"
    resources:
        mem_mb = 64000
    shell:
        """
        set -o pipefail
        bwa-mem2 mem -t {threads} -M -Y \
            -R "@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA" \
            {input.ref} {input.r1} {input.r2} 2> {log} \
        | samtools view -@ {threads} -b -o {output} 2>> {log}
        """

# Mark duplicates
rule mark_duplicates:
    input:
        f"{RES_DIR}/02_Alignment/Alignments/01_RawBAM/{{sample}}.bam"
    output:
        dedup = temp(f"{RES_DIR}/02_Alignment/Alignments/02_DedupBAM/{{sample}}_dedup.bam"),
        index = temp(f"{RES_DIR}/02_Alignment/Alignments/02_DedupBAM/{{sample}}_dedup.bam.bai"),
        stats = f"{RES_DIR}/02_Alignment/QC/{{sample}}_dedup_stats.txt"
    container:
        APPS['gatk']
    threads: 24
    resources:
        mem_mb = 36000
    shell:
        """
        gatk MarkDuplicatesSpark \
            -I {input} \
            -O {output.dedup} \
            -M {output.stats} \
            --conf 'spark.executor.cores={threads}' \
            --conf "spark.executor.memory=32g" \
            --conf "spark.driver.maxResultSize=4g" \
            --create-output-bam-index true
        """

# Generates a recalibration table
rule base_recalibrator:
    input:
        dedup = f"{RES_DIR}/02_Alignment/Alignments/02_DedupBAM/{{sample}}_dedup.bam",
        index = f"{RES_DIR}/02_Alignment/Alignments/02_DedupBAM/{{sample}}_dedup.bam.bai",
        bed = config['proj']['intervals'],
        ref = REFS['genome'],     
        dbsnp = REFS['dbsnp'],
        mills = REFS['mills'],
        high_conf = REFS['high_conf']
    output:
        table = f"{RES_DIR}/02_Alignment/Alignments/{{sample}}_recal_data.tsv"
    container:
        APPS['gatk']
    threads: 24
    resources:
        mem_mb = 36000
    shell:
        """
        gatk BaseRecalibratorSpark \
            -R {input.ref} \
            -I {input.dedup} \
            -L {input.bed} \
            --known-sites {input.dbsnp} \
            --known-sites {input.mills} \
            --known-sites {input.high_conf} \
            -O {output.table} \
            -- \
            --conf 'spark.executor.cores={threads}' \
            --conf "spark.executor.memory=32g" \
            --conf "spark.driver.maxResultSize=4g"
        """

# Apply recalibration
rule apply_bqsr:
    input:
        dedup = f"{RES_DIR}/02_Alignment/Alignments/02_DedupBAM/{{sample}}_dedup.bam",
        ref = REFS['genome'],
        table = f"{RES_DIR}/02_Alignment/Alignments/{{sample}}_recal_data.tsv"
    output:
        recal = protected(f"{RES_DIR}/02_Alignment/Alignments/{{sample}}_recal.bam")
    container:
        APPS['gatk']
    threads: 4
    resources:
        mem_mb = 8000
    shell:
        """
        gatk ApplyBQSR \
            -R {input.ref} \
            -I {input.dedup} \
            --bqsr-recal-file {input.table} \
            -O {output.recal}
        """

# Collect alignment QC stats
rule collect_metrics:
    input:
        recal = f"{RES_DIR}/02_Alignment/Alignments/{{sample}}_recal.bam",
        ref = REFS['genome']
    output:
        stats = f"{RES_DIR}/02_Alignment/QC/{{sample}}_recal_stats.txt"
    container:
        APPS['gatk']
    threads: 2
    resources:
        mem_mb = 4000
    shell:
        """
        gatk CollectAlignmentSummaryMetrics \
            -I {input.recal} \
            -R {input.ref} \
            -O {output.stats}
        """

# Aggregate QC results with MultiQC
rule multiqc_align:
    input:
        expand(f"{RES_DIR}/02_Alignment/QC/{{sample}}_dedup_stats.txt", sample=SAMPLES),
        expand(f"{RES_DIR}/02_Alignment/QC/{{sample}}_recal_stats.txt", sample=SAMPLES)
    params:
        work_dir = f"{RES_DIR}/02_Alignment/QC"
    output:
        report = f"{RES_DIR}/02_Alignment/QC/multiqc_report.html"
    envmodules:
        "apps/multiqc"
    shell:
        """
        multiqc {params.work_dir} \
            -o {params.work_dir} \
            -n multiqc_report.html \
            -f
        """
