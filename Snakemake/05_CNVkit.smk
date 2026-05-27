"""
PIPELINE: CNVkit Copy Number Variation Analysis
DESCRIPTION:
This Snakemake workflow performs copy number variation (CNV) analysis using CNVkit, including panel of normals creation, CNV calling, and annotation with gene information.
It uses VCF information from small variant callers to improve CNV calling accuracy and includes configurable parameters for CNVkit to allow users to adjust sensitivity and specificity based on their needs.
"""

# =============================================================================
# GLOBAL CONFIGURATION
# =============================================================================
# Extract sample names from BAM file names in the input directory
NORMALS = config["pon_samples"]

# Determine which variant caller's VCFs to use for CNVkit calling
caller = config["callers"]["small"]
if caller in ["gatk", "both"]:
    VCF_DIR = f"{GATK_DIR}/01_Variants"
else:
    VCF_DIR = f"{VAR_DIR}/01_Variants"

# =============================================================================
# PIPELINE RULES
# =============================================================================
# Create panel of normals from normal samples
rule create_pon_cnvkit:
    input:
        bams = [f"{ALIGN_DIR}/{normal}_recal.bam" for normal in NORMALS],
        targets = config['proj']['intervals'],
        fasta = REFS['genome'],
        gtf = REFS['gtf']
    params:
        out_dir = f"{CNVkit_DIR}/00_PanelOfNormal"
    output:
        f"{CNVkit_DIR}/00_PanelOfNormal/reference.cnn"
    container:
        APPS['cnvkit']
    threads: 16
    resources:
        mem_mb=32000
    shell:
        """
        cnvkit.py batch {input.bams} \
            --normal \
            --targets {input.targets} \
            --fasta {input.fasta} \
            --annotate {input.gtf} \
            --output-dir {params.out_dir} \
            --processes {threads}
        """

# Run CNVkit batch analysis on all samples
rule cnvkit_batch:
    input:
        bams = [f"{ALIGN_DIR}/{sample}_recal.bam" for sample in SAMPLES],
        targets = config['proj']['intervals'],
        fasta = REFS['genome'],
        gtf = REFS['gtf'],
        ref = f"{CNVkit_DIR}/00_PanelOfNormal/reference.cnn"
    params:
        out_dir = f"{CNVkit_DIR}/01_Batch"
    output:
        cns=expand(f"{CNVkit_DIR}/01_Batch/{{sample}}.cns", sample=SAMPLES),
        cnr=expand(f"{CNVkit_DIR}/01_Batch/{{sample}}.cnr", sample=SAMPLES)
    container:
        APPS['cnvkit']
    threads: 16
    resources:
        mem_mb=32000
    shell:
        """
        cnvkit.py batch {input.bams} \
            --reference {input.ref} \
            --output-dir {params.out_dir} \
            --diagram \
            --scatter \
            --processes {threads}
        """

# Call CNVs for each sample using VCF information
rule cnvkit_call:
    input:
        cns=f"{CNVkit_DIR}/01_Batch/{{sample}}.cns",
        vcf=f"{VCF_DIR}/{{sample}}.vcf"
    output:
        f"{CNVkit_DIR}/02_CNV_Calling/{{sample}}.call.cns"
    container:
        APPS['cnvkit']
    shell:
        """
        cnvkit.py call {input.cns} \
            --vcf {input.vcf} \
            --purity 0.25 \
            --ploidy 2 \
            --filter ampdel \
            -o {output}
        """

# Generate scatter plot visualization for each sample
rule cnvkit_scatter:
    input:
        cnr=f"{CNVkit_DIR}/01_Batch/{{sample}}.cnr",
        cns=f"{CNVkit_DIR}/02_CNV_Calling/{{sample}}.call.cns"
    output:
        f"{CNVkit_DIR}/02_CNV_Calling/{{sample}}_scatter.pdf"
    container:
        APPS['cnvkit']
    shell:
        """
        cnvkit.py scatter {input.cnr} -s {input.cns} -o {output}
        """

# Generate diagram visualization for each sample
rule cnvkit_diagram:
    input:
        cnr=f"{CNVkit_DIR}/01_Batch/{{sample}}.cnr",
        cns=f"{CNVkit_DIR}/02_CNV_Calling/{{sample}}.call.cns"
    output:
        f"{CNVkit_DIR}/02_CNV_Calling/{{sample}}_diagram.pdf"
    container:
        APPS['cnvkit']
    shell:
        """
        cnvkit.py diagram {input.cnr} -s {input.cns} -o {output}
        """

# Annotate coordinates with chromosome arms
rule annotate_coord_cnvkit:
    input:
        call = f"{CNVkit_DIR}/02_CNV_Calling/{{sample}}.call.cns",
        cyto = REFS['cyto']
    params:
        format = "cnvkit",
        min_probes = 0,
        sample = "{sample}"
    output:
        tsv = f"{CNVkit_DIR}/03_Annotated/{{sample}}.ann.tsv"
    conda:
        ENVS['python']
    script:
        SRC['CNV_gene_ann']