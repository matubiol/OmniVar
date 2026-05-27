#!/bin/bash -l

# =========================
# Load modules
# =========================
ml purge
ml load apps/samtools
ml load apps/bcftools
ml load apps/java/jdk-21-corretto

# =========================
# Set working paths
# =========================
OUTDIR="/path/to/resources/"
mkdir -p ${OUTDIR}

# =========================
# GRCh37 reference FASTA and index for BWA and GATK
# =========================
# Get reference genome
REF_URL="https://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"
wget -c ${REF_URL} -P ${OUTDIR}

# Index reference genome
bwa-mem2 index ${OUTDIR}/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
samtools faidx ${OUTDIR}/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
gatk CreateSequenceDictionary \
   -R ${OUTDIR}/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
   -O ${OUTDIR}/Homo_sapiens.GRCh37.75.dna.primary_assembly.dict

# Get reference genome annotations
GFF_URL="https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz"
wget -c ${GFF_URL} -P ${OUTDIR}

# Filter annotation for probes regions only
bedtools intersect \
   -a ${OUTDIR}/Homo_sapiens.GRCh37.87.chr.gtf.gz \
   -b ${OUTDIR}/intervals.bed \
   -wa > ${OUTDIR}/Homo_sapiens.GRCh37.87.intervals.gtf.gz

# Modify the GTF file to include gene names instead of gene ensemble IDs
sed -r 's/gene_id "[^"]+";(.*)gene_name "([^"]+)";/gene_id "\2";\1gene_name "\2";/' \
   ${OUTDIR}/Homo_sapiens.GRCh37.87.intervals.gtf.gz > ${OUTDIR}/Homo_sapiens.GRCh37.87.chr.mod.gtf

# =========================
# BaseRecalibration resources
# =========================
URLS_BR=(
   "https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz"
   "https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz.tbi"
   "https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
   "https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Mills_and_1000G_gold_standard.indels.b37.vcf.gz.tbi"
   "https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/1000G_phase1.snps.high_confidence.b37.vcf.gz"
   "https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/1000G_phase1.snps.high_confidence.b37.vcf.gz.tbi"
)

for URL in "${URLS_BR[@]}"; do
   wget -c ${URL} -P ${OUTDIR}/BaseRecalibration
done

# =========================
# Variant calling resources
# =========================
URLS_VC=(
   "https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf"
   "https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf"
   "https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx"
)

for URL in "${URLS_VC[@]}"; do
   wget -c ${URL} -P ${OUTDIR}/VariantCalling
done

# Filter the gnomAD VCF to include only variants in the target regions
bgzip ${OUTDIR}/VariantCalling/af-only-gnomad.raw.sites.vcf
tabix -p vcf ${OUTDIR}/VariantCalling/af-only-gnomad.raw.sites.vcf.gz

bcftools view \
    -R ${OUTDIR}/Haematology_Oncokkit_HS2_GRCh37_GenesOnly.bed.gz \
    ${OUTDIR}/VariantCalling/af-only-gnomad.raw.sites.vcf.gz \
    -Oz \
    -o ${OUTDIR}/VariantCalling/gnomad.oncokit.vcf.gz
tabix -p vcf ${OUTDIR}/VariantCalling/gnomad.oncokit.vcf.gz

# =========================
# SnpEff (annotation) resources
# =========================
snpeff download -v GRCh37.p13
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi

# =========================
# Copy number variants calling, interval file preprocessing
# =========================
gatk BedToIntervalList \
    -I ${OUTDIR}/intervals.bed \
    -O ${OUTDIR}/CopyNumberVariant/targets.interval_list \
    -SD ${OUTDIR}/Homo_sapiens.GRCh37.75.dna.primary_assembly.dict

gatk PreprocessIntervals \
    -L ${OUTDIR}/CopyNumberVariant/targets.interval_list \
    -R ${OUTDIR}/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    --bin-length 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O ${OUTDIR}/CopyNumberVariant/targets.preprocessed.interval_list

# =========================
# CIViC (annotation) cache
# =========================
URL_CIViC="https://civicdb.org/downloads/01-Apr-2026/01-Apr-2026-AcceptedAndSubmittedClinicalEvidenceSummaries.tsv"
wget ${URL_CIViC} -P ${OUTDIR}/Annotation

echo "=== Setup complete ==="
echo "Files in ${OUTDIR}:"
ls -l ${OUTDIR}