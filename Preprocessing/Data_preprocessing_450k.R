############### prepare phenotype data ###############
library("GEOquery")
gse=getGEO(filename="/ix/ksoyeon/YQ/data/schizophrenia/GSE152027_series_matrix.txt")
pheno <- as.data.frame(cbind(id = as.character(map(strsplit(gse$title, split = " "), 1)), 
                             pheno = gse$`status:ch1`, 
                             age = gse$`ageatbloodcollection:ch1`, 
                             gender = gse$`gender:ch1`))
write.table(pheno, "/ix/ksoyeon/YQ/data/schizophrenia/phenotype.txt", quote = F, row.names = F)

############### prepare methylation data (quality control) ###############
library(data.table)
raw <- fread("/ix/ksoyeon/YQ/data/schizophrenia/GSE152027_IOP_raw_Signal.csv")
rownames(raw) <- raw$V1
cpg.names <- raw$V1
raw <- raw[,-1]
rownames(raw) <- cpg.names
ind <- seq(1, 2400, by = 3)
data <- as.data.frame(raw)

intensity <- data[,-ind]
rownames(intensity) <- cpg.names

##### low quality probes
detp <- data[,ind]
failed <- detp > 0.01
failed.probe <- rowMeans(failed) > 0.1
sum(failed.probe) # 625
failed.sample <- apply(failed, 2, mean) > 0.1
sum(failed.sample) > 0 # FALSE

##### snps info from GEO
info <- read.csv("/ix/ksoyeon/YQ/data/PBMC_inner_city/GPL13534_HumanMethylation450_15017482_v.1.1.copy.csv", header = T)
snp.probe <- info$Name[info$Probe_SNPs != "" | info$Probe_SNPs_10 != ""] # 89,678

##### BOWTIE2 mapping of 450k probes
multi.map <- read.table("/ix/ksoyeon/YQ/data/HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt", header = F)
sum(rownames(raw) %in% multi.map$V1 == T) # 33,457
multi.map.probe <- rownames(raw)[rownames(raw) %in% multi.map$V1 == T]

##### Chen et al., identified non-specific probes across the 450k
cross.react <- read.csv("/ix/ksoyeon/YQ/data/48639-non-specific-probes-Illumina450k.csv")
sum(rownames(raw) %in% cross.react$TargetID == T) # 29,233
cross.react.probe <-  rownames(raw)[rownames(raw) %in% cross.react$TargetID == T]

##### probes on sex chr
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
sum(rownames(raw) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")]) # 11,648
sex.chr.probe <- rownames(raw)[rownames(raw) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")]]

##### filter by all steps
bad.probe <- unique(c(failed.probe, snp.probe, multi.map.probe, cross.react.probe, sex.chr.probe))

filter.bad <- rownames(intensity) %in% bad.probe
intensity.filtered <- intensity[!filter.bad,]
dim(intensity.filtered) # 357,598 * 1600
head(intensity.filtered)

############### prepare methylation data (gene level) ###############
##### get beta value from intensity
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotation <- ann450k[match(rownames(intensity.filtered),ann450k$Name),
                      c(1:4,24:ncol(ann450k))]
probe.without.anno <- which(is.na(rownames(annotation) == T)) # no such probes
intensity.filtered.with.anno <- intensity.filtered
library(dplyr)
intensity.merge <- cbind(intensity.filtered.with.anno,annotation)
library(stringr)
library(purrr)
intensity.merge$gene <- as.character(map(strsplit(intensity.merge$UCSC_RefGene_Name, split = ";"), 1))
intensity.merge <- intensity.merge[intensity.merge$gene != "NULL",] # 275,754 * 1600
gene.cnt <- as.data.frame(table(intensity.merge$gene))
summary(gene.cnt$Freq)
library(ggplot2)
ggplot(gene.cnt, aes(x = Freq)) +
  geom_bar() + 
  xlab("Number of probes per gene") +
  theme_bw()

gene.list <- unique(intensity.merge$gene) # 18,936
gene.average.beta <- as.data.frame(matrix(data = NA, nrow = length(gene.list), ncol = 800))
rownames(gene.average.beta) <- gene.list
colnames(gene.average.beta) <- as.character(paste0(map(strsplit(colnames(intensity.filtered[,seq(2,1600,by = 2)]), split = "_"),1),
                                                   "_",
                                                   map(strsplit(colnames(intensity.filtered[,seq(2,1600,by = 2)]), split = "_"),2)))
i <- 1
for(i in 1:length(gene.list)){
  print(i)
  gene.average.beta[i,1:800] <- colMeans(intensity.merge[which(intensity.merge$gene == gene.list[i]),seq(1,1599,by = 2)])/(colMeans(intensity.merge[which(intensity.merge$gene == gene.list[i]),seq(1,1599,by = 2)])+colMeans(intensity.merge[which(intensity.merge$gene == gene.list[i]),seq(2,1600,by = 2)])+1)
}
write.table(gene.average.beta, "/ix/ksoyeon/YQ/data/schizophrenia/gene.average.beta.by.intensity.txt", quote = F, row.names = T)

############### prepare methylation data (CpG site level) ###############
beta <- intensity.merge[,seq(1,1599,by = 2)]/(intensity.merge[,seq(1,1599,by = 2)]+intensity.merge[,seq(2,1600,by = 2)]+1)
colnames(beta) <- colnames(gene.average.beta)
write.table(beta, "/ix/ksoyeon/YQ/data/schizophrenia/beta.by.intensity.txt", quote = F, row.names = T)
beta <- read.table("/ix/ksoyeon/YQ/data/schizophrenia/beta.by.intensity.txt")
