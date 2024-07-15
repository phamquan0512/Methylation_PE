# tao dia chi work direction (vd: wd1,wd2,...)
wd1 ="D:/QUAN/CAO_HOC/LUAN_VAN/METHYLATION/Script/DATA/GSE98224D"
wd2 ="D:/QUAN/CAO_HOC/LUAN_VAN/METHYLATION/Script/DATA/GSE44667D"
wd3 ="D:/QUAN/CAO_HOC/LUAN_VAN/METHYLATION/Script/DATA/GSE100197_D" #OK
setwd(wd1)
library(datasets)


# tai 1 so thu vien dac biet
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("missMethyl")

# khai bao thu vien
library(limma)
library(minfi)
library(RColorBrewer)
library(missMethyl)
library(matrixStats)
library(umap)
library(stringr)
library(config)
library(GEOquery)
library(rtracklayer)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(missMethyl)

args = commandArgs(trailingOnly = TRUE)
# Input args
data_dir = "D:/QUAN/CAO_HOC/LUAN_VAN/METHYLATION/Script/DATA"
grCol="classify:ch1"
configFile = "D:/QUAN/CAO_HOC/LUAN_VAN/METHYLATION/Script/01.Downloads_and_process_data/GSE162984.config"
gse="GSE98224"
type="raw"

### Load config
config <- config::get(file=configFile)
# Get variables from config file
Pvalue <- config$Pvalue
xReactiveProbes <- config$xReactProbes
C <- config$C
lambda <- config$lambda
fdr <- config$fdr
genome1 <- config$genome1
genome2 <- config$genome2
array <- config$array
anno_dir <- config$anno_dir

Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10000000000000000000000000000000000000000000)
options(timeout=1000000000000000000000000000000000000000000000000000000000000000000000)
set.seed(123)

raw <- paste0(data_dir, "/GSE98224D")
list.files(raw, recursive=TRUE)
targets <- read.table(paste0(raw, "/SampleSheet.tsv"), header=T, sep="\t", check.names=F, stringsAsFactors=F)
targets

# Make sure no space in group's names
targets$base <- str_replace_all(targets[,grCol], " ", ".")
targets$Basename <- paste0(raw, "/", targets$Slide, "_", targets$Array)
QC <- paste0(data_dir, "/QC")
if (!dir.exists(QC)) {dir.create(QC)}

if (type=="raw") {
### Loading data
# read in the raw data from the IDAT files
rgSet <- read.metharray.exp(targets=targets, force=T)
setwd(raw)
system("gunzip -f *.idat.gz")
rgSet <- read.metharray(basenames=targets$Basename, force=T)
rgSet

sampleNames(rgSet) <- targets$sample_acc
rgSet
allProbes <- nrow(rgSet)
save(rgSet, file=paste0(raw, "/rgSet.RData"))
### Quality control
print("QC step")
# Calculate p-values
print("detection P")
detP <- detectionP(rgSet)
head(detP)

# Normalization
print("Normalization using Quantile")
mSetSq <- preprocessQuantile(rgSet)

# Filter based on p-values
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
keep <- rowSums(detP < Pvalue) == ncol(mSetSq)
table(keep)
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt
gtPval <- data.frame(table(keep))[1,2]
ltPval <- data.frame(table(keep))[2,2]

# P_value
pdf(file = paste0(QC, "/Mean_detP.pdf"))
pal = brewer.pal(8,"Dark2")
barplot(colMeans(detP),col=pal[factor(targets$base)],las=2,cex.names=0.8,
        main="Mean detection p-values", xaxt='n')
abline(h=Pvalue, col="red")
dev.off()

sampleP=data.frame(sample_acc=rownames(data.frame(colMeans(detP))), detP=colMeans(detP))
write.table(sampleP, paste0(QC, "/sampleP.tsv"), sep="\t", quote=F, row.names=F)

}else {
load(paste0(raw, "/mSetSqFlt_processed.RData"))
sampleNames(mSetSqFlt) <- targets$sample_acc
gtPval <- 0
allProbes <- nrow(mSetSqFlt)
}

# Probes with snps at CpGs
snpInfo <- getSnpInfo(mSetSqFlt)
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt
snpProbes <- nrow(data.frame(snpInfo)[!is.na(data.frame(snpInfo)$Probe_rs),])

# Exclude cross reactive probes
xReactive <- readRDS(xReactiveProbes)
d <- c()
for(i in which((featureNames(mSetSqFlt) %in% xReactive))) {d <- c(d, i)}
CrosReactiveProbes <- featureNames(mSetSqFlt)[d]
write.table(CrosReactiveProbes, paste0(QC, "/CrosReactiveProbes.tsv"), sep="\t", quote=F, col.names=NA) 
# Remove Probes
keep <- !(featureNames(mSetSqFlt) %in% xReactive)
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]
mSetSqFlt
ReactiveProbes <- data.frame(table(keep))[1,2]

# After filter Probes
FurtherUse <- allProbes-(ReactiveProbes + snpProbes + gtPval)

# Filer pie chart
df <- t(data.frame(Pvalue.Filter=gtPval, SNP.Filter=snpProbes, CrosReactive.Filter=ReactiveProbes, FurtherUse=FurtherUse))
color <- brewer.pal(nrow(df), "Set2")
pdf(file = paste0(QC, "/Probe_filtering.pdf"))
pie_labels <- paste0(rownames(df), " = ", round(100 * df[,1]/sum(df[,1]), 2), "%")
par(mar=c(5.1, 4.1, 4.1, 6.5))
pie(df, main="Probe filtering", labels = pie_labels, col = color)
dev.off()

# Density
pdf(file = paste0(QC, "/Dens.Plot.pdf"))
pal = brewer.pal(8,"Dark2")
densityPlot(getBeta(mSetSqFlt), sampGroups = targets$base, main="densityPlot_normalized")
dev.off()

# Write processed data
methProbes <- paste0(data_dir, "/methProbes")
if (!dir.exists(methProbes)) {dir.create(methProbes)}
save(mSetSqFlt, file=paste0(methProbes, "/mSetSqFlt.RData"))

### Get Beta_value
bVals <- as.data.frame(getBeta(mSetSqFlt))
bigTable <- as.data.frame(apply(bVals, 2, function(x) round(x, 2)))
probeID <- row.names(bigTable)
bigTable$probeID <- probeID

# Coordination
# Genome1
probeID <- read.table(paste0(anno_dir, "/", genome1, "/", array, "_CpG.tsv"), header=T)
dt <- base::merge(bigTable, probeID, by="probeID")
bT.g1 <- dt[,c("CpG_chrm", "CpG_beg", "CpG_end", "probe_strand", "probeID", targets$sample_acc)]
bT.g1 <- bT.g1[order(bT.g1$CpG_chrm, bT.g1$CpG_beg),]
names(bT.g1)[1:5] <- c("chr", "start", "end", "strand", "probeID")
write.table(bT.g1, paste0(methProbes, "/bigTable.", genome1,  ".tsv"), sep="\t", col.names=T, quote=F, row.names=F)
# Genome2
probeID <- read.table(paste0(anno_dir, "/", genome2, "/", array, "_CpG.tsv"), header=T)
dt <- base::merge(bigTable, probeID, by="probeID")
bT.g2 <- dt[,c("CpG_chrm", "CpG_beg", "CpG_end", "probe_strand", "probeID", targets$sample_acc)]
bT.g2 <- bT.g2[order(bT.g2$CpG_chrm, bT.g2$CpG_beg),]
names(bT.g2)[1:5] <- c("chr", "start", "end", "strand", "probeID")
write.table(bT.g2, paste0(methProbes, "/bigTable.", genome2,  ".tsv"), sep="\t", col.names=T, quote=F, row.names=F)

savehistory("D:/QUAN/CAO_HOC/LUAN_VAN/METHYLATION/Script/DATA/IDAT2Bigtable.Rhistory")
