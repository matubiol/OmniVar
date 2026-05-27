<div align="center">
  <h1> Hemato-Oncology pipeline </h1>
  <p><strong>A Comprehensive High-Throughput Snakemake Pipeline for Hemato-Oncology NGS Panels</strong></p>

  <p>
    <img src="https://img.shields.io/badge/Workflow-Snakemake-emerald?style=flat-square&logo=snakemake" alt="Snakemake">
    <img src="https://img.shields.io/badge/Language-Python%20%7C%20R%20%7C%20Bash-blue?style=flat-square&logo=python" alt="Languages">
    <img src="https://img.shields.io/badge/Genome-GRCh37%20%2F%20hg19-orange?style=flat-square" alt="Genome Assembly">
    <img src="https://img.shields.io/badge/License-MIT-yellow?style=flat-square" alt="License">
    <img src="https://img.shields.io/badge/Platform-HPC%20%7C%20Slurm-lightgrey?style=flat-square" alt="Platform">
  </p>
</div>

<hr />

<h2> Overview</h2>
<p>
  <code>Hemato_Somatic_Variants</code> is a production-ready, high-throughput bioinformatic workflow designed specifically for the analysis of Next-Generation Sequencing (NGS) data from hematological malignancies. Somatic variant calling in liquid tumors introduces unique challenges, such as low tumor purity, subclonal mutations, and highly recurrent structural fusions. 
</p>
<p>
  This pipeline automates the entire process from reference infrastructure setup to fully annotated, clinical-grade variant reports by integrating state-of-the-art tools for small variants (SNVs/Indels), structural variations (SVs), and copy number alterations (CNVs).
</p>

<h2> Key Features & Modular Architecture</h2>

<table>
  <thead>
    <tr>
      <th>Module</th>
      <th>Tools Integrated</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><strong>0. Reference Build</strong></td>
      <td><code>BWA-MEM2</code>, <code>GATK</code>, <code>BCFtools</code>, <code>SnpEff</code></td>
      <td>Automated script to build indexes, download dbsnp/gnomAD, filter targeted panel beds, and pre-process GATK interval lists.</td>
    </tr>
    <tr>
      <td><strong>1. Quality Control & Alignment</strong></td>
      <td><code>FastQC</code>, <code>Fastp</code>, <code>BWA-MEM</code>, <code>Samtools</code></td>
      <td>Raw reads QC, adapter trimming, high-precision GRCh37 genome alignment, and post-alignment processing (BQSR, deduplication, sorting).</td>
    </tr>
    <tr>
      <td><strong>2. Small Somatic Variants</strong></td>
      <td><code>GATK Mutect2</code>, <code>VarDict</code></td>
      <td>Dual-caller strategy optimizing sensitivity and specificity for low-VAF (Variant Allele Frequency) somatic mutations typical in hemato-oncology samples.</td>
    </tr>
    <tr>
      <td><strong>3. Structural Variants & Fusions</strong></td>
      <td><code>Manta</code>, <code>Lumpy</code></td>
      <td>Comprehensive detection of structural anomalies, chromosomal rearrangements, and key oncogenic gene fusions (e.g., <i>BCR::ABL1</i>, <i>KMT2A</i> rearrangements).</td>
    </tr>
    <tr>
      <td><strong>4. Copy Number Profiling (CNV)</strong></td>
      <td><code>CNVkit</code>, <code>GATK CNV</code></td>
      <td>Accurate profiling of copy number alterations (amplifications and deletions) tailored for targeted gene panels using preprocessed target intervals.</td>
    </tr>
    <tr>
      <td><strong>5. Automated Clinical Reporting</strong></td>
      <td><code>SnpEff</code>, <code>ClinVar</code>, <code>CIViC Database</code></td>
      <td>Custom downstream annotation filtering out synonymous variants, applying ACMG-like prioritizing tiers, and mapping fusions to actionable clinical evidence.</td>
    </tr>
  </tbody>
</table>

<hr />

<h2> Prerequisites & Environment Setup</h2>
<p>The pipeline relies on <strong>Conda / Mamba</strong> for automated package management and environment isolation, alongside support for HPC cluster environments via <strong>Environment Modules (Lmod/Tcl)</strong>.</p>

<details>
  <summary><b>Click to expand installation steps</b></summary>
  
  <h3>1. Clone the Repository</h3>
  <pre><code>git clone https://github.com/your-username/Hemato_Somatic_Variants.git
cd Hemato_Somatic_Variants</code></pre>

  <h3>2. Install Conda/Mamba & Snakemake</h3>
  <p>Ensure you have Mamba installed. Then create the base execution environment:</p>
  <pre><code>mamba create -c conda-forge -c bioconda -n snakemake snakemake pandas openpyxl
conda activate snakemake</code></pre>
</details>

<hr />

<h2> Reference & Data Preparation</h2>
<p>
  Before running the Snakemake workflow, you must download and index the GRCh37 human reference genome, resource bundles (dbSNP, Mills, gold standard indels), variant calling files (gnomAD panels), and clinical databases (CIViC, ClinVar).
</p>
<p>
  A comprehensive reference builder script is included in the repository. To set up the infrastructure:
</p>

<ol>
  <li>Open <code>scripts/00_BuildRef_GRCh37.sh</code> and configure your desired output directory (<code>OUTDIR</code>) along with your targeted panel regions (<code>intervals.bed</code>).</li>
  <li>Run the script using an environment with access to cluster modules:</li>
</ol>

<pre><code>bash scripts/00_BuildRef_GRCh37.sh</code></pre>

<p>
  This script will automatically index the reference genome, build the sequence dictionaries, intersect and format the GTF annotations with gene names instead of Ensembl IDs, filter gnomAD to your targeted regions via <code>bcftools</code>, pre-process GATK CNV target intervals, and cache the <strong>CIViC Database</strong>.
</p>

<hr />

<h2> Configuration</h2>
<p>Edit the <code>config/config.yaml</code> file to point to the resource directory generated by your build script:</p>

<pre><code># Reference and Databases (GRCh37)
reference_genome: "/path/to/resources/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
gtf_annotation: "/path/to/resources/Homo_sapiens.GRCh37.87.chr.mod.gtf"
civic_database: "/path/to/resources/Annotation/01-Apr-2026-AcceptedAndSubmittedClinicalEvidenceSummaries.tsv"
gnomad_vcf: "/path/to/resources/VariantCalling/gnomad.oncokit.vcf.gz"

# Variant Filtering Thresholds
min_depth: 8
sv_min_support: 3</code></pre>

<p>Configure your sample list in <code>config/samples.tsv</code> using a standard tab-separated format containing your sample identifiers and raw FASTQ paths.</p>

<hr />

<h2> Running the Pipeline</h2>

<h3>Local Execution</h3>
<p>To run the pipeline locally using 8 cores:</p>
<pre><code>snakemake --cores 8 --use-conda</code></pre>

<h3>HPC Cluster Execution (Slurm)</h3>
<p>The workflow is natively optimized for cluster scheduling. To submit jobs to a Slurm partition using Environment Modules:</p>
<pre><code>snakemake --cluster "sbatch --partition={resources.slurm_partition} --account={resources.slurm_account} --mem={resources.mem_mb} --cpus-per-task={threads}" --jobs 10 --use-conda --use-envmodules</code></pre>

<hr />

<h2> Outputs & Automated Reports</h2>
<p>Results are systematically structured inside the <code>results/</code> directory. The core actionable outputs of the pipeline include:</p>

<ul>
  <li>
    <strong>Small Variants Tier Report (<code>results/.../01_SmallVariants/GATK/02_Annotated/{sample}_report.tsv</code>):</strong> 
    Filtered somatic variants classified into pathogenic tiers based on ClinVar status, variant impacts, and molecular effects. It includes an automated rescue protocol for critical hematological insertion/duplication drivers (e.g., <i>FLT3-ITD</i>) and splice regions.
  </li>
  <li>
    <strong>Clinical SV/Fusion Annotation (<code>results/.../02_StructuralVariants/Manta/02_Annotated/{sample}_ann_civic.xlsx</code>):</strong> 
    Automated Excel summary mapping structural variants and fusions against the <strong>CIViC database</strong>, pre-filtered for a curated panel of 29 critical hemato-oncology genes (including <code>ABL1</code>, <code>BCR</code>, <code>JAK2</code>, <code>KMT2A</code>, <code>RUNX1</code>, <code>NPM1</code>, and <code>TP53</code>).
  </li>
</ul>

<hr />

<h2> License</h2>
<p>This project is licensed under the <strong>MIT License</strong> - see the <a href="LICENSE">LICENSE</a> file for details. It grants free permission to use, modify, distribute, and commercially exploit the software, provided the original copyright notice is preserved.</p>

<hr />

<div align="center">
  <p>Developed for Hemato-Oncology Research Groups. For issues, bug reports, or feature requests, please open a GitHub Issue.</p>
</div>
