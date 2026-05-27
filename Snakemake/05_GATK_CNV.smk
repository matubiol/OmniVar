"""
PIPELINE: GATK Copy Number Variation Analysis
USAGE:
    snakemake \
        -s /shared/users/72181407R/proj/2026_Test_Hemato_Pipelines/src/06_GATK_CNV.smk \
        --directory /shared/users/72181407R/proj/2026_Test_Hemato_Pipelines/results/04_VariantCalling/03_CNVs/GATK \
        --profile slurm
"""

# =============================================================================
# GLOBAL CONFIGURATION
# =============================================================================
# Extract sample names from BAM file names in the input directory
NORMALS = config["pon_samples"]

# =============================================================================
# PIPELINE RULES
# =============================================================================
# Process intervals for CNV analysis
rule process_int:
    input:
        bed = config['proj']['intervals'],
        ref = REFS['genome']
    output:
        f"{PROJ_DIR}/data/Ref/CNV/targets.interval_list"
    container:
        APPS['gatk']
    shell:
        """
        gatk PreprocessIntervals \
            -L {input.bed} \
            -R {input.ref} \
            --bin-length 0 \
            -imr OVERLAPPING_ONLY \
            -O {output}
        """

# Collect read counts for normal samples
rule read_count_normals:
    input:
        bam = f"{ALIGN_DIR}/{{normal}}_recal.bam",
        ref = REFS['genome'],
        intervals = config['proj']['intervals']
    output:
        f"{GATK_CNV_DIR}/00_PanelOfNormal/{{normal}}.hdf5"
    container:
        APPS['gatk']
    shell:
        """
        gatk CollectReadCounts \
            -I {input.bam} \
            -R {input.ref} \
            -L {input.intervals} \
            -imr OVERLAPPING_ONLY \
            --format HDF5 \
            -O {output}
        """

# Create panel of normals from normal samples
rule create_pon_gatk:
    input:
        counts = expand(f"{GATK_CNV_DIR}/00_PanelOfNormal/{{normal}}.hdf5", normal=NORMALS)
    output:
        f"{GATK_CNV_DIR}/00_PanelOfNormal/panel_of_normals.hdf5"
    container:
        APPS['gatk']
    shell:
        """
        gatk CreateReadCountPanelOfNormals \
            -I {input.counts[0]} \
            -I {input.counts[1]} \
            -O {output}
        """

# Collect read counts for all samples
rule read_count_all:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_recal.bam",
        ref = REFS['genome'],
        intervals = config['proj']['intervals']
    output:
        temp(f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.counts.hdf5")
    container:
        APPS['gatk']
    shell:
        """
        gatk CollectReadCounts \
            -I {input.bam} \
            -R {input.ref} \
            -L {input.intervals} \
            -imr OVERLAPPING_ONLY \
            --format HDF5 \
            -O {output}
        """

rule allelic_counts:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_recal.bam",
        ref = REFS['genome'],
        gnomad = REFS['gnomad'],
    output:
        f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.allelicCounts.tsv"
    container:
        APPS['gatk']
    shell:
        """
        gatk CollectAllelicCounts \
            -I {input.bam} \
            -R {input.ref} \
            -L {input.gnomad} \
            -O {output}
        """

# Denoise read counts for all samples using panel of normals
rule denoise_counts:
    input:
        counts = f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.counts.hdf5",
        pon = f"{GATK_CNV_DIR}/00_PanelOfNormal/panel_of_normals.hdf5"
    output:
        scr = f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.standardizedCR.tsv",
        dcr = f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.denoisedCR.tsv"
    container:
        APPS['gatk']
    shell:
        """
        gatk DenoiseReadCounts \
            -I {input.counts} \
            --count-panel-of-normals {input.pon} \
            --standardized-copy-ratios {output.scr} \
            --denoised-copy-ratios {output.dcr}
        """

# Plot CNV profiles for all samples
rule plot_denoised:
    input:
        scr = f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.standardizedCR.tsv",
        dcr = f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.denoisedCR.tsv",
        dic = REFS['genome'].with_suffix(".dict")
    output:
        f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.denoised.png",
        temp([
            f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.{suffix}"
            for suffix in [
                "deltaMAD.txt", "denoisedMAD.txt",
                "scaledDeltaMAD.txt", "standardizedMAD.txt"
            ]
        ])
    container:
        APPS['gatk']
    shell:
        """
        gatk PlotDenoisedCopyRatios \
            --standardized-copy-ratios {input.scr} \
            --denoised-copy-ratios {input.dcr} \
            --sequence-dictionary {input.dic} \
            --minimum-contig-length 46709983 \
            --output-prefix {wildcards.sample} \
            --output {GATK_CNV_DIR}/01_Denoising
        """

# Segment copy ratio profiles for all samples
rule model_segments:
    input:
        dcr = f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.denoisedCR.tsv",
        ac = f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.allelicCounts.tsv"
    output:
        [f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.{suffix}"
        for suffix in [
            "cr.seg", "cr.igv.seg",
            "af.igv.seg", "hets.tsv", "modelFinal.seg"
            ]
        ],
        temp([f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.{suffix}"
            for suffix in [
                "modelBegin.seg",
                "modelBegin.cr.param", "modelFinal.cr.param",
                "modelBegin.af.param", "modelFinal.af.param"
            ]])
    container:
        APPS['gatk']
    shell:
        """
        gatk ModelSegments \
            --denoised-copy-ratios {input.dcr} \
            --allelic-counts {input.ac} \
            --output-prefix {wildcards.sample} \
            -O {GATK_CNV_DIR}/02_Segmentation
        """

# Plot segmented CNV profiles for all samples
rule plot_segmented:
    input:
        dcr = f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.denoisedCR.tsv",
        ac = f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.hets.tsv",
        seg = f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.modelFinal.seg",
        dic = REFS['genome'].with_suffix(".dict")
    output:
        f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.modeled.png"
    container:
        APPS['gatk']
    shell:
        """
        gatk PlotModeledSegments \
            --denoised-copy-ratios {input.dcr} \
            --allelic-counts {input.ac} \
            --segments {input.seg} \
            --sequence-dictionary {input.dic} \
            --minimum-contig-length 46709983 \
            --output-prefix {wildcards.sample} \
            --output {GATK_CNV_DIR}/02_Segmentation
        """

# Call CNVs for all samples based on segmented profiles
rule call_cnvs:
    input:
        f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.cr.seg"
    output:
        call = f"{GATK_CNV_DIR}/03_CNV_Calling/{{sample}}.called.seg",
        igv = f"{GATK_CNV_DIR}/03_CNV_Calling/{{sample}}.called.igv.seg"
    container:
        APPS['gatk']
    shell:
        """
        gatk CallCopyRatioSegments \
            -I {input} \
            --neutral-segment-copy-ratio-lower-bound 0.95 \
            --neutral-segment-copy-ratio-upper-bound 1.05 \
            -O {output.call}
        """

# Annotate coordinates with genes and chromosome arms
rule annotate_coord_gatk:
    input:
        call = f"{GATK_CNV_DIR}/03_CNV_Calling/{{sample}}.called.seg",
        cyto = REFS['cyto'],
        bed = config['proj']['intervals']
    params:
        format = "gatk",
        min_probes = 0,
        sample = "{sample}"
    output:
        tsv = f"{GATK_CNV_DIR}/04_Annotated/{{sample}}.ann.tsv"
    conda:
        ENVS['python']
    script:
        SRC['CNV_gene_ann']
