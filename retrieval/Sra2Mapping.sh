#!/bin/bash

HS_GENOME_INDEX="/scratch/dragon/amd/thimmamp/References/hisat2_humangenome/grch38/genome"
MS_GENOME_INDEX="/scratch/dragon/amd/thimmamp/References/hisat2_mousegenome/grcm38/genome"

export PRJ=$PWD
mkdir -p $PRJ/fastq;
mkdir -p $PRJ/mappedFiles;
mkdir -p $PRJ/logs ;
mylogs=$PRJ/logs;
#fastqlocation=$PRJ/fastq;
#mappedlocation=$PRJ/mappedFiles

#for sample in `cat sample.txt`;
for sample in `cat $1`;
do

### To get the prefix of the file 

prefix=`basename -s .sra $sample`

### To get the directory 
location=${sample%/*};
#location=$PRJ/fastq
echo "sample=$sample prefix=$x location=$location"

## Submit SRA to FASTQ
module load sratoolkit/2.9.0
#SRR_CMD="fastq-dump.2 --split-files -O $fastqlocation $sample";
SRR_CMD="fastq-dump.2 --split-files -O $PRJ/fastq $sample";
SRR_Job="sbatch --partition=batch --job-name=SRA_to_Fastq.$prefix --time=2:00:00 --output=$mylogs/$prefix.SRA_2_Fastq-%J.out --error=$mylogs/$prefix.SRA_2_Fastq-%J.err --nodes=1 --cpus-per-task=8";
SRR_ID=$(${SRR_Job} --parsable --wrap="${SRR_CMD}");
echo "SRR Job submitted (\" ${SRR_CMD} is executing \") and this job id is " ${SRR_ID};


## Identify the reference 

#location=$PRJ/mappedFiles
ref=`cat /scratch/dragon/amd/thimmamp/testFiles/lookup.txt | grep $prefix | awk -F '\t' '{print $2}'`
lib=`cat /scratch/dragon/amd/thimmamp/testFiles/lookup.txt | grep $prefix | awk -F '\t' '{print $3}'`
#
### Execute the Mapping job only when First job is successful
module load hisat2/2.1.0
if [ "$ref" == "Mus musculus" ]; then
  myref=$MS_GENOME_INDEX;
else
  myref=$HS_GENOME_INDEX;
fi

if [ -f "$PRJ/fastq/${prefix}_2.fastq" ]; then
   HISAT_CMD="hisat2 -p 16 -x $myref -1 $PRJ/fastq/${prefix}_1.fastq -2 $PRJ/fastq/${prefix}_2.fastq -S $PRJ/mappedFiles/$prefix.sam";
else   
   HISAT_CMD="hisat2 -p 16 -x $myref -U $PRJ/fastq/${prefix}_1.fastq -S $PRJ/mappedFiles/$prefix.sam";
fi

#HISAT_CMD="hisat2 -p 16 -x $myref -U $PRJ/fastq/${prefix}_1.fastq -S $PRJ/mappedFiles/$prefix.sam";
#
#if [ "$lib" == "SINGLE" ]; then
#HISAT_CMD="hisat2 -p 16 -x $myref -U $fastqlocation/${prefix}_1.fastq -S $mappedlocation/$prefix.sam";
#else
#HISAT_CMD="hisat2 -p 16 -x $myref -1 $fastqlocation/${prefix}_1.fastq -2 $fastqlocation/${prefix}_2.fastq -S $mappedlocation/$prefix.sam";
#fi
echo "${HISAT_CMD}"
HISAT_Job="sbatch --partition=batch --job-name=HISAT.$prefix --time=8:00:00 --output=$mylogs/$prefix.HISAT-%J.out --error=$mylogs/$prefix.HISAT-%J.err --nodes=1 --cpus-per-task=16 --mem=115gb";
HISAT_ID=$(${HISAT_Job} --parsable --dependency=afterok:${SRR_ID} --wrap="${HISAT_CMD}");
#HISAT_ID=$(${HISAT_Job} --parsable --wrap="${HISAT_CMD}");
echo " HISAT Job (\" ${HISAT_CMD} \")  was submitted (Job_ID=${HISAT_ID}) and it will execute when the SRR Job_ID=${SRR_ID} is successful"
done 
