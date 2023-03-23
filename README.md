# Sequence capture of internal probes (SCIP) workflow tutorial

_n.b._ End-to-end script to analyze paired-end DNA sequencing data from SCIP experiments

Dependecies:

[bwa](https://sourceforge.net/projects/bio-bwa/files/)
[picard](http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics)  
[samtools](http://www.htslib.org/)
[delly](https://github.com/dellytools/delly)
[bcftools](http://www.htslib.org/)
[bedtools](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)

## Step 1. Create folder structure to store the output
```
mkdir $(pwd)/bam_non_dup_picard
mkdir $(pwd)/bam_non_dup_picard/delly
mkdir $(pwd)/bam_non_dup_picard/discordants
```

## Step 1. Map and sort reads to the merged genomes using [bwa](https://sourceforge.net/projects/bio-bwa/files/)

General option
```
DIR_FASTQS="/PATH/TO/FASTQS"
FILES=(SPACE SEPARATED LIST OF GZIPPED FASTQ BASENAMES I.E. WITHOUT _1.fastq.gz or _2.fastq.gz)
DIR_GENOMES="/PATH/TO/GENOMES"
```
RF_PICR=/scratch/gnnxm/Repositories/Pharma/Genomes/CriGri-PICR/CriGri-PICR_Original/ncbi-genomes-2022-11-16/GCF_003668045.1_CriGri-PICR_genomic.fna
RF_PICR_VEC39079=/scratch/gnnxm/Repositories/Pharma/Genomes/CriGri-PICR/Merged_CriGri-PICR_VEC_39079/Merged_VEC_39079_GCF_003668045.1_CriGri-PICR_genomic.fna
RF_PICRH=/scratch/gnnxm/Repositories/Pharma/Genomes/CriGri-PICRH/CriGri-PICRH-1.0_Original/ncbi-genomes-2022-04-08/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna
RF_PICRH_VEC33310=/scratch/gnnxm/Repositories/Pharma/Genomes/CriGri-PICRH/Merged_CriGri-PICRH-1.0_VEC_33310/Merged_CriGri-PICRH-1.0_VEC_33310/Merged_VEC_33310_GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna
RF_PICRH_VEC39079=/scratch_raid/gnnxm/LSC_Location/Repositories/Pharma/Genomes/CriGri-PICRH/Merged_CriGri-PICRH-1.0_VEC_39079/Merged_CriGri-PICRH-1.0_VEC_39079/Merged_VEC_39079_GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna
RF_CHO_Bayer_hifiasm=/scratch_raid/gnnxm/LSC_Location/Repositories/Pharma/Genomes/CHO_Bayer/assem/asm.fasta
RF_CHO_Bayer_ragoo=/scratch_raid/gnnxm/LSC_Location/Repositories/Pharma/Genomes/CHO_Bayer/scaffolding/asm_ragoo.fasta
```


for i in "${FILES[@]}"
  do
    java -jar /usr/local/bin/Trimmomatic-0.39/trimmomatic-0.39.jar \
    PE \
    -threads 12 \
    -phred33 \
    ${DIR}/${i}_1.fastq.gz \
    ${DIR}/${i}_2.fastq.gz \
    ${DIR}/${i}_1_P.fastq.gz \
    ${DIR}/${i}_1_U.fastq.gz \
    ${DIR}/${i}_2_P.fastq.gz \
    ${DIR}/${i}_2_U.fastq.gz \
    ILLUMINACLIP:${DIR}/adapters.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:30
  done
```



The above will perform the following:

Remove adapters (`ILLUMINACLIP:adapters.fa:2:30:10`)
Note the additional `:2` in front of `keepBothReads` this is the minimum adapter length in palindrome mode
Remove leading low quality or N bases (below quality 3) (`LEADING:3`)
Remove trailing low quality or N bases (below quality 3) (`TRAILING:3`)
Drop reads below the 30 bases long (`MINLEN:30`)

## Step 2. Build index and map reads using [bowtie](http://bowtie-bio.sourceforge.net/manual.shtml)

Index genome assembly fasta:

```
INDEX="SHORT IDENTIFIER FOR REFERENCE GENOME ASSEMBLY"
REF="REFERENCE GENOME ASSEMBLY FASTA FILE NAME"
DIR="/PATH/TO/REFERENCE/FASTA"
cd ${DIR}
bowtie-build ${REF} ${INDEX}
```

Map ATAC-seq reads to indexed genome assembly:

```
INDEX="/PATH/TO/INDEX"
DIR="/PATH/TO/TRIMMED/FASTQS"
FILES=(SPACE SEPARATED LIST OF TRIMMED FASTQ BASENAMES I.E. WITHOUT _1_P.fastq or _2_P.fastq)

for i in "${FILES[@]}"
  do
    gunzip ${DIR}/${i}_1_P.fastq.gz
    gunzip ${DIR}/${i}_2_P.fastq.gz
    bowtie ${INDEX} -1 ${DIR}/${i}_1_P.fastq -2 ${DIR}/${i}_2_P.fastq -S ${DIR}/${i}.sam -t -p 24 -v 2 -X 1000 --best --strata -m 1
    gzip ${DIR}/${i}_1_P.fastq
    gzip ${DIR}/${i}_2_P.fastq
    samtools sort -O 'bam' -o ${DIR}/${i}_sorted.bam -T tmp -@ 24 ${DIR}/${i}.sam
    rm ${DIR}/${i}.sam
  done
```

In addtion to mapping the ATAC-seq reads, the above will perform the following:

`gunzip` fastqs prior to running `bowtie`
`gzip` fastqs after mapping
Convert output sam file to bam, sort bam, and delete sam file (`samtools` and `rm`)

## Step 3. Remove PCR and optical duplicates from bam using [picard](http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics)

```
DIR="/PATH/TO/SORTED/BAMS"
FILES=(SPACE SEPARATED LIST OF BAM BASENAMES I.E. WITHOUT _sorted.bam)

for i in "${FILES[@]}"
  do
    java -Xmx48G -jar picard.jar MarkDuplicates INPUT=${DIR}/${i}_sorted.bam OUTPUT=${DIR}/${i}_rmdup.bam METRICS_FILE=${DIR}/${i}_METRICS_FILE.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
  done
```

## Step 4. Merge bam technical replicates

```
DIR="/PATH/TO/REMOVED/DUPLICATE/BAMS"
LINE="LINE NAME E.G. B73"
TIS="NAME OF TISSUE E.G. leaf"

samtools merge ${DIR}/${LINE}_${TIS}.bam -@ 12 ${DIR}/FILE1_rmdup.bam ${DIR}/FILE2_rmdup.bam ${DIR}/FILE3_rmdup.bam
```

The above requires the user to specify the removed duplicate bam files. Repeat code for additional tissues.

## Step 5. Call peaks using [macs2](https://github.com/macs3-project/MACS)

```
DIR="/PATH/TO/MERGED/BAM"
FILES=(SPACE SEPARATED LIST OF MERGED REMOVED BAM BASENAMES I.E. WITHOUT .bam)

for i in "${FILES[@]}"
  do
    macs2 callpeak \
    --verbose 3 \
    --treatment ${DIR}/${i}.bam \
    --format BAMPE \
    --name ${i} \
    --outdir ${DIR} \
    --bdg \
    --qvalue 0.01 \
    --gsize 2.4e9 \
    --keep-dup all \
    --tempdir ${DIR}/temp
  done
```

The above is uses a false-discovery rate (`--qvalue`) of 0.01 estimated genome size of 2.4 Gbp for _Zea mays_ (`--gsize`).

## Assess Fraction of Reads in Peaks (FRiP) score

"_The fraction of reads in called peak regions (FRiP score) should be >0.3, though values greater than 0.2 are acceptable_". [ENCODE](https://www.encodeproject.org/atac-seq/)

```
DIR="/PATH/TO/MERGED/BAM/AND/NARROWPEAKS"
BAM="BASENAME OF BAM I.E. WITHOUT .bam"
PEAK="BASENAME OF NARROWPEAK I.E. WITHOUT .narrowPeak"

total_reads=$(samtools view -@ 24 -c ${BAM}.bam)
reads_in_peaks=$(bedtools sort -i ${PEAK}.narrowPeak \
| bedtools merge -i stdin \
| bedtools intersect -u -nonamecheck -a ${BAM}.bam -b stdin -ubam \
| samtools view -@ 24 -c)
FRiP=$(awk "BEGIN {print "${reads_in_peaks}"/"${total_reads}"}")
echo "tissue=${i}, total_reads=${total_reads}, reads_in_peaks=${reads_in_peaks}, FRiP=${FRiP}" >> ${BAM}_frip.txt
```
