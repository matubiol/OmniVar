"""
PIPELINE: Whole Genome Bisulfite Sequencing Analysis using BSBOLT
DESCRIPTION:
This Snakemake workflow performs structural variant detection using Manta, followed by gene annotation and clinical annotation using the CIViC database.
It includes configurable parameters for Manta to allow users to adjust sensitivity and specificity based on their needs.
"""

# =============================================================================
# PIPELINE RULES
# =============================================================================
# configure and run Manta
rule run_manta:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_recal.bam",
        db = REFS['genome']
    params:
        dir_name = f"{MANTA_DIR}/01_SV_Calling/{{sample}}",
        script_dir = f"{MANTA_DIR}/01_SV_Calling/{{sample}}/runWorkflow.py"
    output:
        vcf = f"{MANTA_DIR}/01_SV_Calling/{{sample}}/results/variants/tumorSV.vcf.gz"
    log:
        f"{MANTA_DIR}/01_SV_Calling/{{sample}}/Slurm_logs/{{sample}}.log",
    benchmark:
        f"{MANTA_DIR}/01_SV_Calling/{{sample}}/Slurm_logs/{{sample}}.benchmark"
    container:
        APPS['manta']
    threads: 16
    resources:
        mem_mb=32000
    shell:
        """
        # Configure Manta for somatic SV calling
        configManta.py \
            --tumorBam {input.bam} \
            --referenceFasta {input.db} \
            --exome \
            --runDir {params.dir_name}
        # Run Manta workflow with specified parameters for improved sensitivity
        {params.script_dir} \
            -m local \
            -j {threads} 2> {log}
        """

# Primary gene annotation.
rule ann_genes_manta:
    input:
        vcf = f"{MANTA_DIR}/01_SV_Calling/{{sample}}/results/variants/tumorSV.vcf.gz",
        gtf = REFS['gtf']
    params:
        min_support = config.get("sv_min_support", 3)
    output:
        tsv = temp(f"{MANTA_DIR}/02_Annotated/{{sample}}_ann_gene.tsv")
    conda:
        ENVS['python']
    script:
        SRC['SV_gene_ann']

# Clinical annotation using CIViC database.
rule ann_civic_manta:
    input:
        tsv = f"{MANTA_DIR}/02_Annotated/{{sample}}_ann_gene.tsv",
        civic = REFS['civic']
    output:
        xlsx = f"{MANTA_DIR}/02_Annotated/{{sample}}_ann_civic.xlsx"
    envmodules:
        "apps/R"
    script:
        SRC['SV_CIVIC_ann']
