# Sequence capture of internal probes (SCIP) workflow tutorial

_n.b._ End-to-end script to analyze paired-end DNA sequencing data from SCIP experiments

Dependecies:

[bwa](https://sourceforge.net/projects/bio-bwa/files/)
[picard](http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics/)  
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

## Step 2. Map and sort reads to the merged genomes using [bwa](https://sourceforge.net/projects/bio-bwa/files/)

General option
```
DIR_FASTQS="/PATH/TO/FASTQS"
FILES=(SPACE SEPARATED LIST OF GZIPPED FASTQ BASENAMES I.E. WITHOUT _1.fastq.gz or _2.fastq.gz)
DIR_GENOMES="/PATH/TO/GENOMES"
```
Specific option with precise file paths 
```
RF_PICR=/scratch/gnnxm/Repositories/Pharma/Genomes/CriGri-PICR/CriGri-PICR_Original/ncbi-genomes-2022-11-16/GCF_003668045.1_CriGri-PICR_genomic.fna
RF_PICR_VEC39079=/scratch/gnnxm/Repositories/Pharma/Genomes/CriGri-PICR/Merged_CriGri-PICR_VEC_39079/Merged_VEC_39079_GCF_003668045.1_CriGri-PICR_genomic.fna
RF_PICRH=/scratch/gnnxm/Repositories/Pharma/Genomes/CriGri-PICRH/CriGri-PICRH-1.0_Original/ncbi-genomes-2022-04-08/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna
RF_PICRH_VEC33310=/scratch/gnnxm/Repositories/Pharma/Genomes/CriGri-PICRH/Merged_CriGri-PICRH-1.0_VEC_33310/Merged_CriGri-PICRH-1.0_VEC_33310/Merged_VEC_33310_GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna
RF_PICRH_VEC39079=/scratch_raid/gnnxm/LSC_Location/Repositories/Pharma/Genomes/CriGri-PICRH/Merged_CriGri-PICRH-1.0_VEC_39079/Merged_CriGri-PICRH-1.0_VEC_39079/Merged_VEC_39079_GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna
RF_CHO_Bayer_hifiasm=/scratch_raid/gnnxm/LSC_Location/Repositories/Pharma/Genomes/CHO_Bayer/assem/asm.fasta
RF_CHO_Bayer_ragoo=/scratch_raid/gnnxm/LSC_Location/Repositories/Pharma/Genomes/CHO_Bayer/scaffolding/asm_ragoo.fasta
```

```
for i in "${FILES[@]}"
  do
    bwa mem -t 20 $RF_PICRH_VEC39079 $MATE1 $MATE2 | samtools sort -@10 -o  $(pwd)/$BASENAME.sorted.bam
  done
```

The above will perform the following:

Aligned the paired-end reads to the merged reference genome from the options above using bwa-mem
Sort the resulting bam file based on position

## Step 3. Mark duplicates on the sorted bam files using [picard](http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics) 

Detects and tags PCR and optical duplicates in the sorted bam file
```
for i in *.bam
do
        BASENAME=${i%%.bam}                           
        echo $BASENAME
        picard MarkDuplicates VALIDATION_STRINGENCY=LENIENT I=${BASENAME}.bam O=$(pwd)/bam_non_dup_picard/${BASENAME}_nodups.bam M=$(pwd)/bam_non_dup_picard/${BASENAME}_marked_dup_metrics.txt TMP_DIR= /scratch_raid/gnnxm/Tempo/

done
```

## Step 4. Create a stats summary of the bam file and estimes the average insert size for each library

```
rm *.bam
cd $(pwd)/bam_non_dup_picard/

for i in *.bam
do
        echo $i
        samtools stats $i | grep "insert size average:" | cut -f3 >> bam_average_insert_size.tsv
done
```

## Step 5. Calling interchromosomal translocations is performed using [delly](https://github.com/dellytools/delly)

```
for i in *.bam
do
        BASENAME=${i%%.bam}
        echo $BASENAME
        T=${BASENAME}.bam
        echo -e  $T "\t"
        samtools index $i
        delly call --svtype BND \ 
                    -o $(pwd)/delly/$BASENAME.bcf \
                    -g $RF_PICRH_VEC39079 $T

done
```

The above calls all the interchromosomal translocations in the sample, however, we just care about the translocations where the transfected vector is involved, (e.g. VEC39079) and these will be filtered in step 8
 
## Step 6. Fetch all discordants reads from the bam file using [samtools](http://www.htslib.org/)

```
for i in *.bam
do
        BASENAME=${i%%.bam}
        echo $BASENAME
        echo $i
        samtools view -b -F 1294 $i > $(pwd)/discordants/$BASENAME.discordants.bam
        samtools index $(pwd)/discordants/$BASENAME.discordants.bam
done
```

The above selects all the discordant reads from the sample and indexes the file to be uploaded as tracks on any genome browser (e.g. IGV)



