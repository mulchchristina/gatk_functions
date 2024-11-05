#!/bin/bash
#
#SBATCH --job-name=AnnotateSVdataSubsets
#SBATCH --ntasks=1
#SBATCH --get-user-env
#SBATCH --output=AnnotateSVdataSubsets.log
#SBATCH --error=AnnotateSVdataSubsets.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=64000
#SBATCH --partition=exacloud
#SBATCH --time=0-8

# Pipeline Overview:
# This script subsets, annotates, and extracts structural variants (SVs) from a VCF file using GATK, based on multiple filters:
# - Uses a BED file for initial filtering to include only specific regions
# - Filters based on impact and excess heterozygosity
# - Extracts summary metrics and produces detailed log outputs for each subset created

set -e  # Exit if a command fails
set -x  # Print commands and their arguments as they execute

# Paths
export PATH=/home/exacloud/gscratch/prime-seq/bin:$PATH
SV_VCF="/path/to/your/input.vcf.gz"
BED_FILE="/path/to/your/regions.bed"
GATK="/path/to/GenomeAnalysisTK4.jar"
JAVA="/path/to/java"
REFERENCE_GENOME="/path/to/128_Mmul_10.fasta"

# Intermediate Files
BASENAME="AnnotatedSVSubsets"
LOG_FILE="${BASENAME}_run.log"

# Function for processing and annotating SV data with GATK and generating summary outputs
runVariantsToTable() {
    INPUT=$1
    TBL=$(echo "$INPUT" | sed 's/vcf\.gz/table\.txt/g')
    
    # Convert VCF to Table
    $JAVA -jar $GATK VariantsToTable \
        -R $REFERENCE_GENOME \
        -V $INPUT \
        -O $TBL \
        -F CHROM -F POS -F END -F ID -F REF -F ALT \
        -F SVTYPE -F AF -F NCALLED -F HET -F HOM-REF -F HOM-VAR \
        -F SVLEN -F IMPACT -F OG -F ExcHet -F FILTER -F CLNSIG
    
    # Clean up long sequence names
    awk -v OFS='\t' '{ if (length($5) > 50) $5="<LONG>"; if (length($6) > 50) $6="<LONG>"; print $0 }' $TBL > tmp.txt
    mv tmp.txt $TBL
    gzip $TBL

    # Log summary statistics
    echo "Summary for $TBL:" >> "$LOG_FILE"
    zcat $TBL | awk 'BEGIN { OFS="\t" } NR>1 { count[$7]++ } END { for (type in count) print type, count[type] }' >> "$LOG_FILE"
}

# Initial BED file annotation
AnnotatedVCF="${BASENAME}.bedannot.vcf.gz"
if [ ! -e ${AnnotatedVCF}.done ]; then
    $JAVA -jar $GATK SelectVariants \
        -R $REFERENCE_GENOME \
        -V $SV_VCF \
        -L $BED_FILE \
        -O $AnnotatedVCF

    echo "Initial BED annotation complete." >> "$LOG_FILE"
    runVariantsToTable "$AnnotatedVCF"
    touch ${AnnotatedVCF}.done
fi

# Subsetting by impact, excess heterozygosity, and SV length
HighImpactVCF="${BASENAME}.highImpact.vcf.gz"
if [ ! -e ${HighImpactVCF}.done ]; then
    $JAVA -jar $GATK SelectVariants \
        -R $REFERENCE_GENOME \
        -V $AnnotatedVCF \
        -O $HighImpactVCF \
        -select "IMPACT == 'HIGH'"

    echo "High-impact subset complete." >> "$LOG_FILE"
    runVariantsToTable "$HighImpactVCF"
    touch ${HighImpactVCF}.done
fi

LowExcessHetVCF="${BASENAME}.lowExHet.vcf.gz"
if [ ! -e ${LowExcessHetVCF}.done ]; then
    $JAVA -jar $GATK SelectVariants \
        -R $REFERENCE_GENOME \
        -V $AnnotatedVCF \
        -O $LowExcessHetVCF \
        -select "ExcHet < 0.05"

    echo "Low Excess Het subset complete." >> "$LOG_FILE"
    runVariantsToTable "$LowExcessHetVCF"
    touch ${LowExcessHetVCF}.done
fi

LongSVVCF="${BASENAME}.longSV.vcf.gz"
if [ ! -e ${LongSVVCF}.done ]; then
    $JAVA -jar $GATK SelectVariants \
        -R $REFERENCE_GENOME \
        -V $AnnotatedVCF \
        -O $LongSVVCF \
        -select "SVLEN > 125000"

    echo "Long SV length subset complete." >> "$LOG_FILE"
    runVariantsToTable "$LongSVVCF"
    touch ${LongSVVCF}.done
fi

# Log completion and summary of subsets
echo "Pipeline completed successfully. Subset logs and summaries:" >> "$LOG_FILE"
cat "$LOG_FILE"
