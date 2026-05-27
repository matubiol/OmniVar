"""
PIPELINE: Lumpy SV Detection & Clinical Annotation
DESCRIPTION:
This Snakemake workflow performs structural variant detection using Lumpy, followed by gene annotation and clinical annotation using the CIViC database.
It includes configurable parameters for Lumpy to allow users to adjust sensitivity and specificity based on their needs.
"""

# =============================================================================
# PIPELINE RULES
# =============================================================================
# Extract split reads and discordant pairs.
rule extract_reads:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_recal.bam"
    output:
        split = temp(f"{LUMPY_DIR}/01_SV_Calling/{{sample}}_splitters.bam"),
        disc  = temp(f"{LUMPY_DIR}/01_SV_Calling/{{sample}}_discordants.bam")
    conda:
        ENVS['lumpy']
    threads: 8
    resources:
        mem_mb=20000
    shell:
        """
        set -o pipefail
        # Extract Split Reads
        samtools view -h -F 1294 {input.bam} \
            | extractSplitReads_BwaMem -i stdin \
            | samtools view -Sb - \
            | samtools sort -@ {threads} -o {output.split}

        # Extract Discordant Pairs
        samtools view -b -F 1294 {input.bam} \
            | samtools sort -@ {threads} -o {output.disc}
        """

# Run LUMPY structural variant caller.
rule run_lumpy:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_recal.bam",
        split = f"{LUMPY_DIR}/01_SV_Calling/{{sample}}_splitters.bam",
        disc = f"{LUMPY_DIR}/01_SV_Calling/{{sample}}_discordants.bam",
        ref = REFS['genome']
    output:
        vcf = f"{LUMPY_DIR}/01_SV_Calling/{{sample}}.vcf"
    log:
        f"{LUMPY_DIR}/01_SV_Calling/Logs/{{sample}}.log",
    benchmark:
        f"{LUMPY_DIR}/01_SV_Calling/Logs/{{sample}}.benchmark"
    conda:
        ENVS['lumpy']
    resources:
        mem_mb=20000
    shell:
        """
        lumpyexpress \
            -B {input.bam} \
            -R {input.ref} \
            -S {input.split} \
            -D {input.disc} \
            -o {output.vcf} 2> {log}
        """

# Primary gene annotation.
rule ann_genes_lumpy:
    input:
        vcf = f"{LUMPY_DIR}/01_SV_Calling/{{sample}}.vcf",
        gtf = REFS['gtf']
    params:
        min_support = config.get("sv_min_support", 3)
    output:
        tsv = temp(f"{LUMPY_DIR}/02_Annotated/{{sample}}_ann_gene.tsv")
    conda:
        ENVS['python']
    script:
        SRC['SV_gene_ann']

# Clinical annotation using CIViC database.
rule ann_civic_lumpy:
    input:
        tsv = f"{LUMPY_DIR}/02_Annotated/{{sample}}_ann_gene.tsv",
        civic = REFS['civic']
    output:
        xlsx = f"{LUMPY_DIR}/02_Annotated/{{sample}}_ann_civic.xlsx"
    envmodules:
        "apps/R"
    script:
        SRC['SV_CIVIC_ann']