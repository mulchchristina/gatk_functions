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

# Pipeline script for looking at specific aspects of an structural variant data set, with the input being a bed file to get the subsets of interest and a vcf to add annotation data to make the variants in each subset relevant to downstream plotting 
# This script subsets, annotates, and extracts structural variants (SVs) from a VCF file using GATK:
# - Inputs are a vcf file with or without prior annotations and a BED file for filtering to include only specific regions
# - Filters based on impact and excess heterozygosity
# - for each subset( ie something like variants with excess heterozygosity) creates a log and a summary of the subset

#the importance of the subsets has to do with functional aspects we are interested in adding downstream. We want to remove technical bias at allele frequency .5 by removing those with exceptional excess heterozygosity from the  whole( downstream filtering includes taking the ids with very low excess het an filtering them)
#impact is of interest because we want to contrast high-impact variants to those with medium, low, or without an impact score( when the impact is added by snpEff upstream only those in coding regions have scores)
#SV length the longest structural variants in the data are of interest to us bc they may affect multiple genes or coding regions


set -e  # forces abort if a step fails to produce a done file
set -x  # prints the steps as they run

# Paths to files of interest and GATK toolkit, current java is also needed for GATK 
export PATH=/home/exacloud/gscratch/prime-seq/bin:$PATH
SV_VCF="/path/to/input.vcf.gz"
BED_FILE="/path/to/regions.bed"
GATK="/path/to/GenomeAnalysisTK4.jar"
JAVA="/path/to/java"
REFERENCE_GENOME="/path/to/128_Mmul_10.fasta"

# each step has output files with the basename for file tracking for this project
BASENAME="AnnotatedSVSubsets"
LOG_FILE="${BASENAME}_run.log"

# this step uses variants to table to make an easily referenced version of variants SV data with GATK and generating summary outputs
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
    
    # tidy long sequence names
    awk -v OFS='\t' '{ if (length($5) > 50) $5="<LONG>"; if (length($6) > 50) $6="<LONG>"; print $0 }' $TBL > tmp.txt
    mv tmp.txt $TBL
    gzip $TBL

    # add summary statistics to log files
    echo "Summary for $TBL:" >> "$LOG_FILE"
    zcat $TBL | awk 'BEGIN { OFS="\t" } NR>1 { count[$7]++ } END { for (type in count) print type, count[type] }' >> "$LOG_FILE"
}

# Annotate the bed file with GATK
AnnotatedVCF="${BASENAME}.bedannot.vcf.gz"
if [ ! -e ${AnnotatedVCF}.done ]; then
    $JAVA -jar $GATK SelectVariants \
        -R $REFERENCE_GENOME \
        -V $SV_VCF \
        -L $BED_FILE \
        -O $AnnotatedVCF

    echo " input bed files annotated with GATK" >> "$LOG_FILE"
    runVariantsToTable "$AnnotatedVCF"
    touch ${AnnotatedVCF}.done
fi

# Subsetting by impact, excess heterozygosity, and SV length

# high impact subset
HighImpactVCF="${BASENAME}.highImpact.vcf.gz"
if [ ! -e ${HighImpactVCF}.done ]; then
    $JAVA -jar $GATK SelectVariants \
        -R $REFERENCE_GENOME \
        -V $AnnotatedVCF \
        -O $HighImpactVCF \
        -select "IMPACT == 'HIGH'"

    echo "high impact variants subset" >> "$LOG_FILE"
    runVariantsToTable "$HighImpactVCF"
    touch ${HighImpactVCF}.done
fi

# very low excess het pulled out for filtering

LowExcessHetVCF="${BASENAME}.lowExHet.vcf.gz"
if [ ! -e ${LowExcessHetVCF}.done ]; then
    $JAVA -jar $GATK SelectVariants \
        -R $REFERENCE_GENOME \
        -V $AnnotatedVCF \
        -O $LowExcessHetVCF \
        -select "ExcHet < 0.05"

    echo "excess het variants subset" >> "$LOG_FILE"
    runVariantsToTable "$LowExcessHetVCF"
    touch ${LowExcessHetVCF}.done
fi

#very long SVs over 125k bp pulled out, the range for svs in this data is 50 to ~150k
LongSVVCF="${BASENAME}.longSV.vcf.gz"
if [ ! -e ${LongSVVCF}.done ]; then
    $JAVA -jar $GATK SelectVariants \
        -R $REFERENCE_GENOME \
        -V $AnnotatedVCF \
        -O $LongSVVCF \
        -select "SVLEN > 125000"

    echo "Long SVs subset " >> "$LOG_FILE"
    runVariantsToTable "$LongSVVCF"
    touch ${LongSVVCF}.done
fi

#subsets of SVs can now be added to a table in R and plotted to explore these subsets of interest

# Log completion and summary of subsets
echo "Pipeline completed successfully. Subset logs and summaries:" >> "$LOG_FILE"
cat "$LOG_FILE"
