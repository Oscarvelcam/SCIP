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

## Step 7. Create a variant call file (VCF) from the binary output of delly using [bcftools](http://www.htslib.org/)  

```
cd  $(pwd)/delly
for i in *.bcf
do
        BASENAME=${i%%.bcf}
        echo $BASENAME
        bcftools convert $BASENAME.bcf -o $BASENAME.vcf
done
```
The above selects all the discordant reads from the sample and indexes the file to be uploaded as tracks on any genome browser (e.g. IGV)

## Step 8. Filter the list of putative integration sites based on quality and type of discordant reads supporting it 

```
for i in *.vcf
do
        BASENAME=${i%%.vcf}
        echo $BASENAME
        less $BASENAME.vcf| grep "VEC39079" | grep -v "IMPRECISE" | grep -v "LowQual" | grep -v "##" | cut -f8 | cut -d ";" -f5,6 > $BASENAME.bed
done
```

```
Rscript SCIP_2_Ideogram.R path2input/file.bed output_file.txt
```

```
─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.2.2 (2022-10-31)
 os       macOS Ventura 13.2
 system   aarch64, darwin22.1.0
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Europe/Berlin
 date     2023-03-23
 rstudio  2022.07.1+554 Spotted Wakerobin (desktop)
 pandoc   NA

─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 package     * version   date (UTC) lib source
 assertthat    0.2.1     2019-03-21 [1] CRAN (R 4.2.1)
 base64enc     0.1-3     2015-07-28 [1] CRAN (R 4.2.1)
 cli           3.4.1     2022-09-23 [1] CRAN (R 4.2.1)
 colorspace    2.0-3     2022-02-21 [1] CRAN (R 4.2.1)
 DBI           1.1.3     2022-06-18 [1] CRAN (R 4.2.1)
 dplyr         1.0.10    2022-09-01 [1] CRAN (R 4.2.1)
 fansi         1.0.3     2022-03-24 [1] CRAN (R 4.2.1)
 generics      0.1.3     2022-07-05 [1] CRAN (R 4.2.1)
 ggplot2       3.4.0     2022-11-04 [1] CRAN (R 4.2.1)
 glue          1.6.2     2022-02-24 [1] CRAN (R 4.2.1)
 grImport2     0.2-0     2019-11-25 [1] CRAN (R 4.2.1)
 gtable        0.3.1     2022-09-01 [1] CRAN (R 4.2.1)
 jpeg          0.1-9     2021-07-24 [1] CRAN (R 4.2.1)
 lifecycle     1.0.3     2022-10-07 [1] CRAN (R 4.2.1)
 magrittr      2.0.3     2022-03-30 [1] CRAN (R 4.2.1)
 munsell       0.5.0     2018-06-12 [1] CRAN (R 4.2.1)
 pillar        1.8.1     2022-08-19 [1] CRAN (R 4.2.1)
 pkgconfig     2.0.3     2019-09-22 [1] CRAN (R 4.2.1)
 png           0.1-7     2013-12-03 [1] CRAN (R 4.2.1)
 purrr         0.3.5     2022-10-06 [1] CRAN (R 4.2.1)
 R6            2.5.1     2021-08-19 [1] CRAN (R 4.2.1)
 RIdeogram   * 0.2.2     2020-01-20 [1] CRAN (R 4.2.2)
 rlang         1.0.6     2022-09-24 [1] CRAN (R 4.2.1)
 rstudioapi    0.14      2022-08-22 [1] CRAN (R 4.2.1)
 rsvg          2.3.2     2022-10-18 [1] CRAN (R 4.2.1)
 scales        1.2.1     2022-08-20 [1] CRAN (R 4.2.1)
 sessioninfo   1.2.2     2021-12-06 [1] CRAN (R 4.2.1)
 tibble        3.1.8     2022-07-22 [1] CRAN (R 4.2.1)
 tidyr         1.2.1     2022-09-08 [1] CRAN (R 4.2.1)
 tidyselect    1.2.0     2022-10-10 [1] CRAN (R 4.2.1)
 utf8          1.2.2     2021-07-24 [1] CRAN (R 4.2.1)
 vctrs         0.5.0     2022-10-22 [1] CRAN (R 4.2.1)
 XML           3.99-0.12 2022-10-28 [1] CRAN (R 4.2.1)

 [1] /opt/homebrew/lib/R/4.2/site-library
 [2] /opt/homebrew/Cellar/r/4.2.2_1/lib/R/library

───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

```
