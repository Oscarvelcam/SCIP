install.packages("RIdeogram")
require(RIdeogram)
#> Loading required package: RIdeogram

# Load the files necessary for the display of the ideogram data 
data(gene_density, package="RIdeogram")
data(Random_RNAs_500, package="RIdeogram")


# This function is quite important since it helps us create an karyotype file for our own species
setwd("/Users/oscar_vc/Documents/bin/CriGri-PICRH-1.0/genome_assemblies_genome_gff/ncbi-genomes-2022-04-08/")

#Create a karyotype file for our species with 3 columns only since we don't know the centromere position
karyotype = "/Users/oscar_vc/Documents/bin/CriGri-PICRH-1.0/genome_assemblies_genome_gff/ncbi-genomes-2022-04-08/CriGri_PICRH_1_0_karyotype.tsv"

gff <- read.table("/Users/oscar_vc/Documents/bin/CriGri-PICRH-1.0/genome_assemblies_genome_gff/ncbi-genomes-2022-04-08/GCF_003668045.3_CriGri-PICRH-1.0_genomic_NOCOMMENTS.gff",
                  stringsAsFactors = F,
                  header = F,
                  comment.char = "#",
                  sep = '\t',
                  quote = ""
)

karyotype <- read.table("/Users/oscar_vc/Documents/bin/CriGri-PICRH-1.0/genome_assemblies_genome_gff/ncbi-genomes-2022-04-08/CriGri_PICRH_1_0_karyotype.tsv",
                        sep = "\t",
                        header = T,
                        stringsAsFactors = F
)

# The original file for the 39 insertions is located in the following direction 
# /data/gnnxm/AM_scMultiome_4992/FASTQ_scMultiome/2_4992Bayer_V261-44ofTOP48_Multiome_CriGri_2/outs
# This list was obtained by running the delly_run.sh on top of the aggregated file "atac_possorted_bam.bam"
# wc /data/gnnxm/AM_scMultiome_4992/FASTQ_scMultiome/2_4992Bayer_V261-44ofTOP48_Multiome_CriGri_2/outs/atac_possorted_bam.bed

#T_insertions <- read.table("/Users/oscar_vc/Documents/bin/CriGri-PICRH-1.0/genome_assemblies_genome_gff/ncbi-genomes-2022-04-08/71_BND_merged.txt",
#                        sep = "\t",
#                        header = T,
#                        stringsAsFactors = F
#)


T_insertions <- read.table("/Users/oscar_vc/Documents/LLL/scMultiomics/ActivMotif/ARC_Cellranger_Sep09/1_4992Bayer_V217-Pool1_Multiome_CriGri_2/Integrations_1.tsv",
                           sep = "\t",
                           header = T,
                           stringsAsFactors = F
)


feature = "gene"
window = 1000000
gff <- subset(gff, gff$V1 %in% karyotype$Chr & gff$V3 == feature)

list_chr <- vector("list", length(names(table(gff$V1))))
names(list_chr) <- names(table(gff$V1))

for (i in 1:(length(list_chr))){
  list_chr[[i]] <- as.data.frame(table(cut(subset(gff, gff$V1 == names(list_chr[i]))$V4,
                                           breaks = c(seq(0, subset(karyotype, karyotype$Chr == names(list_chr[i]))[1,3], window),
                                                      subset(karyotype, karyotype$Chr == names(list_chr[i]))[1,3]))))
  list_chr[[i]] <- tidyr::separate(list_chr[[i]], 1, into = c("Start","End"), sep = ",")
  list_chr[[i]]$Start <- gsub('\\(', '', list_chr[[i]]$Start)
  list_chr[[i]]$End <- gsub('\\]', '', list_chr[[i]]$End)
  list_chr[[i]]$Start <- as.numeric(list_chr[[i]]$Start)
  list_chr[[i]]$End <- as.numeric(list_chr[[i]]$End)
  list_chr[[i]]$Start <- list_chr[[i]]$Start + 1
  list_chr[[i]][nrow(list_chr[[i]]),2] <- subset(karyotype, karyotype$Chr == names(list_chr[i]))[1,3]
  list_chr[[i]]$Chr <- names(list_chr[i])
  list_chr[[i]] <- cbind(list_chr[[i]][,4], list_chr[[i]][,1:3])
  colnames(list_chr[[i]]) <- c("Chr", "Start", "End", "Value")
  list_chr[[i]]$Chr <- as.character(list_chr[[i]]$Chr)
}

l <- data.frame()
for (i in 1:(length(list_chr))){
  df.now <- list_chr[[i]]
  l <- rbind(l, df.now)
}

gene_density <- data.frame(l)
karyotype <- karyotype[1:12,]
#gene_density <- gene_density[1:10,]
ideogram(karyotype = karyotype, overlaid = gene_density, label = T_insertions, label_type = "marker")

convertSVG("chromosome.svg", device = "png")
