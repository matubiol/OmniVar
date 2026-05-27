"""
PIPELINE: Hemato Somatic Variant & Structural/CNV Calling Pipeline (2026)
DESCRIPTION:
    A comprehensive high-throughput pipeline designed for Hemato-Oncology
    NGS panels. It integrates raw data QC, alignment,
    small variant calling (GATK Mutect2, VarDict), structural variant
    detection (Manta, Lumpy), and copy number variation profiling
    (CNVkit, GATK CNV).
USAGE:
    snakemake -s src/Snakefile.smk --configfile src/config.yml --profile slurm
AUTHOR:
    Matías Fernández Huarte (matubiol@gmail.com)
"""

# =============================================================================
# GLOBAL CONFIGURATION
# =============================================================================
from pathlib import Path

# Load global configuration file
configfile: "config.yaml"

# Set up paths based on config
PROJ_DIR = config['proj']['proj_dir']
RAW_DIR = config['proj']["raw_dir"]
RES_DIR = f"{PROJ_DIR}/results"
REF_DIR = Path(config['proj']['ref_dir'])
OPT_DIR = Path(config['proj']['opt_dir'])
SRC_DIR = Path(config['proj']['src_dir'])

# Resolve absolute paths for references, containers, environments, and post-processing scripts
REFS = {key: REF_DIR / path for key, path in config["refs"].items()}
APPS = {key: str(OPT_DIR / path) for key, path in config["containers"].items()}
ENVS = {key: OPT_DIR / path for key, path in config["envs"].items()}
SMK = {key: SRC_DIR / path for key, path in config["smk"].items()}
SRC = {key: SRC_DIR / path for key, path in config["post"].items()}

# Extract sample names from file names in the input directory
raw_files = [p.name for p in Path(RAW_DIR).glob("*_R1*.fastq.gz")]
SAMPLES = sorted(list(set([re.split(r'_', f)[0] for f in raw_files])))

print(f"Discovered samples: {SAMPLES}")

# Set paths
QC_DIR = f"{RES_DIR}/01_RawReads_QC/02_Fastp"
ALIGN_DIR = f"{RES_DIR}/02_Alignment/Alignments"
GATK_DIR = f"{RES_DIR}/03_VariantCalling/01_SmallVariants/GATK"
VAR_DIR = f"{RES_DIR}/03_VariantCalling/01_SmallVariants/VarDict"
LUMPY_DIR = f"{RES_DIR}/03_VariantCalling/02_StructuralVariants/Lumpy"
MANTA_DIR = f"{RES_DIR}/03_VariantCalling/02_StructuralVariants/Manta"
GATK_CNV_DIR = f"{RES_DIR}/03_VariantCalling/03_CNVs/GATK"
CNVkit_DIR = f"{RES_DIR}/03_VariantCalling/03_CNVs/CNVkit"

# =============================================================================
# PIPELINE CONFIGURATION & VALIDATION
# =============================================================================
# Set type of variant calls to run based on config
VALID_TYPES = {"small", "sv", "cnv"}

# Allow user to specify a single type as a string, or multiple types as a list
calls_input = config["calls"]
if isinstance(calls_input, str):
    calls_input = [calls_input]

# Validate that the input is either 'all' or a list of valid types
if "all" in calls_input:
    CALL_TYPES = ["small", "sv", "cnv"]
else:
    # Validate that all specified call types are valid
    invalid_calls = set(calls_input) - VALID_TYPES
    if invalid_calls:
        raise ValueError(
            f"Invalid call type(s) specified: {invalid_calls}. "
            f"Choose from 'small', 'sv', 'cnv', or 'all'."
        )    
    # Assign the validated list of call types
    CALL_TYPES = calls_input

# Set small variant caller(s) to run based on config
if config["callers"]["small"] == "gatk":
    SNP_CALLERS = ["gatk"]
elif config["callers"]["small"] == "vardict":
    SNP_CALLERS = ["vardict"]
elif config["callers"]["small"] == "both":
    SNP_CALLERS = ["gatk", "vardict"]
else:
    raise ValueError("Invalid variant caller specified in config. Choose 'gatk', 'vardict', or 'both'.")

# Set structural variant caller(s) to run based on config
if config["callers"]["sv"] == "manta":
    SV_CALLERS = ["manta"]
elif config["callers"]["sv"] == "lumpy":
    SV_CALLERS = ["lumpy"]
elif config["callers"]["sv"] == "both":
    SV_CALLERS = ["manta", "lumpy"]
else:
    raise ValueError("Invalid SV caller specified in config. Choose 'manta', 'lumpy', or 'both'.")

# Set CNV caller(s) to run based on config
if config["callers"]["cnv"] == "cnvkit":
    CNV_CALLERS = ["cnvkit"]
elif config["callers"]["cnv"] == "gatk_cnv":
    CNV_CALLERS = ["gatk_cnv"]
elif config["callers"]["cnv"] == "both":
    CNV_CALLERS = ["cnvkit", "gatk_cnv"]
else:
    raise ValueError("Invalid CNV caller specified in config. Choose 'cnvkit', 'gatk_cnv', or 'both'.")

# =============================================================================
# DYNAMIC INPUT BUILDING FOR TARGET RULE
# =============================================================================
# Always include core QC and alignment reports
final_outputs = [
    f"{RES_DIR}/01_RawReads_QC/03_MultiQC/MultiQC_report.html",
    f"{RES_DIR}/02_Alignment/QC/multiqc_report.html"
]

# Append small variant workflow outputs if requested
if "small" in CALL_TYPES or "all" in CALL_TYPES:
    if "gatk" in SNP_CALLERS:
        final_outputs.extend(expand(f"{GATK_DIR}/02_Annotated/{{sample}}_report.tsv", sample=SAMPLES))
    if "vardict" in SNP_CALLERS:
        final_outputs.extend(expand(f"{VAR_DIR}/02_Annotated/{{sample}}_report.tsv", sample=SAMPLES))

# Append structural variant workflow outputs if requested
if "sv" in CALL_TYPES or "all" in CALL_TYPES:
    if "manta" in SV_CALLERS:
        final_outputs.extend(expand(f"{MANTA_DIR}/02_Annotated/{{sample}}_ann_civic.xlsx", sample=SAMPLES))
    if "lumpy" in SV_CALLERS:
        final_outputs.extend(expand(f"{LUMPY_DIR}/02_Annotated/{{sample}}_ann_civic.xlsx", sample=SAMPLES))

# Append copy number variant workflow outputs if requested
if "cnv" in CALL_TYPES or "all" in CALL_TYPES:
    if "cnvkit" in CNV_CALLERS:
        final_outputs.extend(
            expand([            
                f"{CNVkit_DIR}/02_CNV_Calling/{{sample}}_scatter.pdf",
                f"{CNVkit_DIR}/02_CNV_Calling/{{sample}}_diagram.pdf",
                f"{CNVkit_DIR}/03_Annotated/{{sample}}.ann.tsv"
            ], sample=SAMPLES)
        )
    if "gatk_cnv" in CNV_CALLERS:  
        final_outputs.extend(
            expand([
                f"{GATK_CNV_DIR}/01_Denoising/{{sample}}.denoised.png",
                f"{GATK_CNV_DIR}/02_Segmentation/{{sample}}.modeled.png",
                f"{GATK_CNV_DIR}/03_CNV_Calling/{{sample}}.called.seg",
                f"{GATK_CNV_DIR}/04_Annotated/{{sample}}.ann.tsv"
            ], sample=SAMPLES)
        )

# =============================================================================
# TARGET RULE - Defines final expected outputs for the entire pipeline
# =============================================================================
rule all:
    input:
        final_outputs

# =============================================================================
# MODULE INCLUSIONS
# =============================================================================
include: SMK['qc']
include: SMK['align']
include: SMK['vc_gatk']
include: SMK['vc_vardict']
include: SMK['sv_manta']
include: SMK['sv_lumpy']
include: SMK['cnv_cnvkit']
include: SMK['cnv_gatk']