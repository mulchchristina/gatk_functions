#!/bin/bash
#
#SBATCH --job-name=DoSvSubsetting
#SBATCH --ntasks=1
#SBATCH --get-user-env
#SBATCH --output=DoSvSubsetting.txt
#SBATCH --error=DoSvSubsetting.txt
#SBATCH --cpus-per-task=2
#SBATCH --mem=64000
#SBATCH --partition=exacloud
#SBATCH --time=0-8

# This should run VTT

set -e
set -x

export PATH=/home/exacloud/gscratch/prime-seq/bin:$PATH

# Output ID: 577018
SV_VCF=/home/groups/OnprcColonyData/ColonyData/322/@files/sequenceOutputs/fixploidy.vcf.gz

GATK=/home/exacloud/gscratch/prime-seq/bin/GenomeAnalysisTK4.jar
JAVA=/home/groups/prime-seq/exacloud/java/current/bin/java

MMUL10=/home/groups/prime-seq/production/Shared/@files/.referenceLibraries/128/128_Mmul_10.fasta
GTF=/home/groups/prime-seq/production/Shared/@files/.referenceLibraries/128/tracks/NCBI.Mmul_10.103.translated.gtf

BASENAME=PacBioPara

ARG_FILE=VariantsToTable.args.list
echo "-F CHROM" > $ARG_FILE
echo "-F POS" >> $ARG_FILE
echo "-F END" >> $ARG_FILE
echo "-F ID" >> $ARG_FILE
echo "-F REF" >> $ARG_FILE
echo "-F ALT" >> $ARG_FILE
echo "-F SVTYPE" >> $ARG_FILE
echo "-F AF" >> $ARG_FILE
echo "-F NCALLED" >> $ARG_FILE
echo "-F HET" >> $ARG_FILE
echo "-F HOM-REF" >> $ARG_FILE
echo "-F HOM-VAR" >> $ARG_FILE
echo "-F ME" >> $ARG_FILE
echo "-F SVLEN" >> $ARG_FILE
echo "-F IMPACT" >> $ARG_FILE
echo "-F OG" >> $ARG_FILE
echo "-F HIG" >> $ARG_FILE
echo "-F VE" >> $ARG_FILE
echo "-F ExcHet" >> $ARG_FILE
echo "-F ExcessHet" >> $ARG_FILE
echo "-F FILTER" >> $ARG_FILE
echo "-F CLNSIG" >> $ARG_FILE
echo "-F FILTER" >> $ARG_FILE
echo "-F FILTER" >> $ARG_FILE

runVariantsToTable() {
	INPUT=$1
	TBL=$(echo "$INPUT" | sed 's/vcf\.gz/table\.txt/g')
	
	$JAVA -jar $GATK VariantsToTable \
		-R $MMUL10 \
		-V $INPUT \
		-O $TBL \
		--arguments_file $ARG_FILE
		
	cat $TBL | grep -v 'BND'  | awk -v OFS='\t' ' { if( length($5) > 50) $5="<LONG>"; if ( length($6) > 50) $6="<LONG>"; print $0 } ' > tmp.txt
	rm $TBL
	mv tmp.txt $TBL
	gzip $TBL
}
