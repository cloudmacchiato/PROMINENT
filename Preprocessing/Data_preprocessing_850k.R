library('GEOquery')
library('limma')
library('minfi')
library('IlluminaHumanMethylationEPICanno.ilm10b4.hg19')
library('sva')
library('IlluminaHumanMethylationEPICmanifest')
library('ENmix')
library('lumi')
library('knitr')
library('sesame')
library('dplyr')
library('BiocParallel')
library('stringr')

################### save raw data ###################
sesameDataCache()
sdfs <- searchIDATprefixes("/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis")
mft = sesameDataGet("MM285.address")$ordering
betas = openSesame("/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis", manifest = mft)
idat_dir <- "/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis"
sdfs <- openSesame(idat_dir, prep="QCDPB", func = NULL, BPPARAM = BiocParallel::MulticoreParam(20))
betas <- openSesame(idat_dir, prep="QCDPB", func = getBetas, BPPARAM = BiocParallel::MulticoreParam(20)) 
head(sdfs)
head(betas)
dim(betas)
save(sdfs, file = "/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/sdfs.RData")
save(betas, file = "/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/betas.RData")

sdfs = lapply(searchIDATprefixes(system.file(
  "extdata/", package = "sesameData")), readIDATpair)

################### prepare phenotype data ###################
setwd('/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/')
head(list.files("/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis", pattern = "idat"))
idatFiles <- list.files("/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)
#series_matrix <- getGEO(filename="/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/GSE175459-GPL21145_series_matrix.txt")
pheno <- as.data.frame(
  cbind(Sample_Name = series_matrix$geo_accession,
        Sample_Well = series_matrix$`methylation_well:ch1`,
        Sample_Plate = '',
        Sample_Group	= series_matrix$`case_control_status:ch1`,
        Pool_ID = '',
        Sentrix_ID	= str_split_i(series_matrix$`methylationid:ch1`,"_",1),
        Sentrix_Position	= str_split_i(series_matrix$`methylationid:ch1`,"_",-1),
        age = series_matrix$`age:ch1`,
        sex = series_matrix$`Sex:ch1`,
        race = series_matrix$`race:ch1`,
        ethnicity = series_matrix$`ethnicity:ch1`,
        smoking = series_matrix$`smoking_ever_never:ch1`))

write.csv(pheno, "/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/sample_sheet.csv",quote = F,row.names = F)

################### prepare methylation data (quality control) ###################
dataDirectory <- "/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis"
# list the files
list.files(dataDirectory, recursive = TRUE)
setwd('/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/')

# read in the sample sheet for the experiment
targets <- read.metharray.sheet(dataDirectory, pattern="sample_sheet.csv")
targets$Basename <- paste0(dataDirectory,"/",targets$Sample_Name,"_",targets$Slide,"_",targets$Array,"_noid")

cat('load idat files.\n')
# read in the raw data from the IDAT files; warnings can be ignored.
rgSet <- read.metharray.exp(targets=targets,force=TRUE)
save(rgSet, file="/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/rgSet.RData")
load("/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/rgSet.RData")

cat('detection p-values.\n')
detP <- detectionP(rgSet)
save(detP, file="/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/detP.RData")
load("/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/detP.RData")

MSet <- preprocessRaw(rgSet)

Meth <- getMeth(MSet)
Unmeth <- getUnmeth(MSet)

mSetSq <- preprocessQuantile(rgSet)
save(mSetSq, file="/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/mSetSq.RData")

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detectionP(rgSet)
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]

# remove any probes that have failed in one or more samples; this next line
# checks for each row of detP whether the number of values < 0.01 is equal
# to the number of samples (TRUE) or not (FALSE)
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
table(keep)
# keep
#FALSE   TRUE 
#40535 825324
# Subset the GenomicRatioSet
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt
# 799006 left

# remove probes that have been demonstrated to map to multiple places in the genome
# probes from Pidsley 2016 (EPIC)
epic.cross1 <- read.csv("/zfs1/hpark/YQ/sepsis_methylation/Problematic_Probes_Pidsley2016/13059_2016_1066_MOESM1_ESM.csv", head = T)
epic.cross2 <- read.csv("/zfs1/hpark/YQ/sepsis_methylation/Problematic_Probes_Pidsley2016/13059_2016_1066_MOESM2_ESM.csv", head = T)
epic.cross3 <- read.csv("/zfs1/hpark/YQ/sepsis_methylation/Problematic_Probes_Pidsley2016/13059_2016_1066_MOESM3_ESM.csv", head = T)
epic.cross4 <- read.csv("/zfs1/hpark/YQ/sepsis_methylation/Problematic_Probes_Pidsley2016/13059_2016_1066_MOESM4_ESM.csv", head = T)
epic.cross5 <- read.csv("/zfs1/hpark/YQ/sepsis_methylation/Problematic_Probes_Pidsley2016/13059_2016_1066_MOESM5_ESM.csv", head = T)
epic.cross6 <- read.csv("/zfs1/hpark/YQ/sepsis_methylation/Problematic_Probes_Pidsley2016/13059_2016_1066_MOESM6_ESM.csv", head = T)
epic.cross.probes <- unique(c(epic.cross1$...1,epic.cross2$PROBE,epic.cross3$PROBE,epic.cross4$PROBE,epic.cross5$PROBE,epic.cross6$PROBE))
keep <- !(featureNames(mSetSqFlt) %in% epic.cross.probes)
table(keep)
# keep
# FALSE   TRUE 
# 121842 677164 
mSetSqFlt <- mSetSqFlt[keep,]

# remove probes on Sex chromosome
library(missMethyl)
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
keep <- !(featureNames(mSetSqFlt) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])
table(keep)
#  FALSE   TRUE 
# 15618 661546 
mSetSqFlt <- mSetSqFlt[keep,]
save(mSetSqFlt, file="/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/mSetSqFlt.RData")

MethFlt <- Meth[rownames(mSetSqFlt),]
UnmethFlt <- Unmeth[rownames(mSetSqFlt),]
save(MethFlt, file="/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/MethFlt.RData")
save(UnmethFlt, file="/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/UnmethFlt.RData")
load("/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/MethFlt.RData")
load("/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/UnmethFlt.RData")

################### prepare methylation data (gene level) ###################
annotation <- ann850k[match(rownames(MethFlt),ann850k$Name),
                      c(1:4,22:ncol(ann850k))]
write.table(annotation, "/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/gene.annotation.txt", quote = F, row.names = T)

probe.without.anno <- which(is.na(rownames(annotation) == T)) # 0 
MethFlt.with.anno <- MethFlt
library(dplyr)
intensity.merge <- cbind(MethFlt.with.anno,annotation)
library(stringr)
library(purrr)
gene <- as.character(map(strsplit(annotation$UCSC_RefGene_Name, split = ";"), 1))
MethFlt <- MethFlt[gene != "NULL",] # 481378 * 547
UnmethFlt <- UnmethFlt[gene != "NULL",] # 481378 * 547
gene <- gene[gene != "NULL"] # 481378
gene.cnt <- as.data.frame(table(gene))
summary(gene.cnt$Freq)
library(ggplot2)
ggplot(gene.cnt, aes(x = Freq)) +
  geom_bar() + 
  xlab("Number of probes per gene") +
  theme_bw()

gene.list <- unique(gene)
gene.average.beta <- as.data.frame(matrix(data = NA, nrow = length(gene.list), ncol = 547))
rownames(gene.average.beta) <- gene.list
colnames(gene.average.beta) <- colnames(MethFlt)
i <- 1
for(i in 1:length(gene.list)){
  print(c(i,gene.list[i]))
  if(length(which(gene == gene.list[i])) > 1){
    gene.average.beta[i,1:547] <- colMeans(MethFlt[which(gene == gene.list[i]),])/(colMeans(MethFlt[which(gene == gene.list[i]),])+colMeans(UnmethFlt[which(gene == gene.list[i]),])+1)
  }
  else{
    gene.average.beta[i,1:547] <- mean(MethFlt[which(gene == gene.list[i]),])/(mean(MethFlt[which(gene == gene.list[i]),])+mean(UnmethFlt[which(gene == gene.list[i]),])+1)
  }
}
write.table(gene.average.beta, "/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/gene.average.beta.by.intensity.txt", quote = F, row.names = T) #24728
gene.average.beta <- read.table("/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/gene.average.beta.by.intensity.txt")

################### prepare methylation data (CpG site level) ###################
beta <- MethFlt/(MethFlt+UnmethFlt+1)
write.table(beta, "/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/beta.by.intensity.all.regions.txt", quote = F, row.names = T) # 688,195 * 100
write.table(beta, "/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/beta.by.intensity.txt", quote = F, row.names = T) # 481378
beta <- read.table("/ix/ksoyeon/YQ/data/Idiopathic_Pulmonary_Fibrosis/beta.by.intensity.txt")

