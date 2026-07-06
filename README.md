<div align="center">
  <h1> OmniVar </h1>
  <p><strong>A Comprehensive High-Throughput Snakemake Pipeline for Somatic NGS Panels</strong></p>

  <p>
    <img src="https://img.shields.io/badge/Workflow-Snakemake-emerald?style=flat-square&logo=snakemake" alt="Snakemake">
    <img src="https://img.shields.io/badge/Language-Python%20%7C%20R%20%7C%20Bash-blue?style=flat-square&logo=python" alt="Languages">
    <img src="https://img.shields.io/badge/Genome-GRCh37%20%2F%20hg19-orange?style=flat-square" alt="Genome Assembly">
    <img src="https://img.shields.io/badge/License-MIT-yellow?style=flat-square" alt="License">
    <img src="https://img.shields.io/badge/Platform-HPC%20%7C%20Slurm-lightgrey?style=flat-square" alt="Platform">
  </p>
</div>

<hr />

<hr />

<h2> Overview</h2>
<p>
  <code>OmniVar</code> is a high-throughput bioinformatic workflow designed for the analysis of Next-Generation Sequencing (NGS) data in clinical oncology. Whether applied to liquid biopsies or solid tumor tissue, somatic variant calling introduces unique challenges such as low tumor purity, subclonal mutations, and complex structural rearrangements.
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
      <td><strong>0. Reference Preparation</strong></td>
      <td></td>
      <td>Harmonizes the chromosome naming convention of the intervals bed file with that of the reference fasta (e.g. "chr1" vs "1").</td>
    </tr>
    <tr>
      <td><strong>1. Quality Control </strong></td>
      <td><code>FastQC</code>, <code>Fastp</code>, <code>MultiQC</code></td>
      <td>Quality control on raw DNA-seq reads, including FastQC for quality assessment, Fastp for trimming and filtering, and MultiQC for aggregating QC reports.</td>
    </tr>
    <tr>
      <td><strong>2. Alignment </strong></td>
      <td><code>BWA-MEM2</code>, <code>GATK</code>, <code>MultiQC</code></td>
      <td>Aligns the processed reads to the reference genome using BWA-MEM2 and performs subsequent processing with GATK.</td>
    </tr>
    <tr>
      <td><strong>3. Small Somatic Variants</strong></td>
      <td><code>GATK Mutect2</code>, <code>VarDict</code>, <code>SnpEff</code></td>
      <td>Somatic variant calling using GATK Mutect2 and VarDict, followed by annotation with SnpEff and ClinVar. It includes configurable parameters for Mutect2 to allow users to adjust sensitivity and specificity based on their needs.</td>
    </tr>
    <tr>
      <td><strong>4. Structural Variants & Fusions</strong></td>
      <td><code>Manta</code>, <code>Lumpy</code></td>
      <td>Detection of genomic rearrangements, translocations, and pathologically relevant gene fusions, followed by gene annotation and clinical annotation using the CIViC database.</td>
    </tr>
    <tr>
      <td><strong>5. Copy Number Profiling</strong></td>
      <td><code>CNVkit</code>, <code>GATK</code></td>
      <td>Copy number profiling using CNVkit and GATK CNV, including panel of normals creation, CNV calling, and annotation with gene information. .</td>
    </tr>
  </tbody>
</table>

<hr />

<h2> Prerequisites & Environment Setup</h2>
<p>The pipeline relies on <strong>Conda / Mamba</strong> for automated package management and environment isolation, alongside support for HPC cluster environments via <strong>Environment Modules</strong>.</p>

<details>
  <summary><b>Click to expand installation steps</b></summary>
  
  <h3>1. Clone the Repository</h3>
  <pre><code>git clone https://git.iconcologia.net/72181407R/OmniVar.git
cd OmniVar</code></pre>

  <h3>2. Install Conda/Mamba & Snakemake</h3>
  <p>Ensure you have Mamba installed. Then create the base execution environment:</p>
  <pre><code>mamba create -c conda-forge -c bioconda -n snakemake snakemake pandas openpyxl
mamba activate snakemake</code></pre>
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
  <li>Open <code>src/BuildRef_GRCh37.sh</code> and configure your desired output directory (<code>OUTDIR</code>).</li>
  <li>Run the script using an environment with access to cluster modules:</li>
</ol>

<pre><code>bash src/BuildRef_GRCh37.sh</code></pre>

<p>
  This script will automatically index the reference genome, build the sequence dictionaries, intersect and format the GTF annotations with gene names instead of Ensembl IDs, filter gnomAD to your targeted regions via <code>bcftools</code>, pre-process GATK CNV target intervals, and cache the <strong>CIViC Database</strong>.
</p>

<hr />

<h2> Configuration</h2>
<p>Edit the <code>config.yml</code> file to point to the resource directory generated by your build script:</p>
<p>Configure your sample list in <code>config/samples.tsv</code> using a standard tab-separated format containing your sample identifiers and raw FASTQ paths.</p>

<hr />

<h2> Running the Pipeline</h2>

<h3>Execution</h3>
<p>To run the pipeline locally using 8 cores:</p>
<pre><code>snakemake -s Snakefile.smk --configfile config.yml --profile slurm</code></pre>
<p>Slurm profiles can be configured in <code>config/slurm/</code> to specify partition, account, memory, and other cluster-specific parameters. Other options can be specified as needed. Find more information in the <a href="https://snakemake.readthedocs.io/en/stable/executing/cli.html#defining-global-profiles">Snakemake documentation</a>.</p>

<hr />

<h2> License</h2>
<p>This project is licensed under the <strong>MIT License</strong> - see the <a href="LICENSE">LICENSE</a> file for details. It grants free permission to use, modify, distribute, and commercially exploit the software, provided the original copyright notice is preserved.</p>

<hr />

<div align="center">
  <p>Developed for Hemato-Oncology Research Groups. For issues, bug reports, or feature requests, please open a GitHub Issue.</p>
</div>
