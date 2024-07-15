setwd("D:/QUAN/CAO_HOC/LUAN_VAN/METHYLATION/Script/DATA")
library(GEOquery)
library(curl)
library(minfi)
require(doParallel)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#tai package Illumina
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
BiocManager::install("GEOquery")



#args=commandArgs(trailingOnly = TRUE)
# Input args
data_dir="D:/QUAN/CAO_HOC/LUAN_VAN/METHYLATION/Script/DATA"
gse="GSE181034"
type="raw"
array="450K"

Sys.setenv("VROOM_CONNECTION_SIZE"=131072*100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000)
options(timeout=100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000)

### Download data
if (!dir.exists(data_dir)) {dir.create(data_dir)}
raw=paste0(data_dir, "/GSE181034D")
dir.create(raw)
gseInfo=getGEO(gse,GSEMatrix=TRUE, destdir=raw)
print("Be careful, columns name might different with different GSE. Check column detailedly before setting up 'type'")

#pData và phenoData cho phép truy cap vào data.
#pData se show nguyen ban du lieu ra
# phenoData show thong tin du lieu nhu so luong bien, so luong ID

#liet ke ten cac cot trong bang du lieu
colnames(pData(phenoData(gseInfo[[1]])))

#tmp là lenh de ghi nho tam thoi/ lenh gan
tmp=pData(phenoData(gseInfo[[1]]))[,c(grep(":ch1", colnames(pData(phenoData(gseInfo[[1]])))), 1, 8, grep("supplementary_file", colnames(pData(phenoData(gseInfo[[1]])))))]

# Remove white space 
tmp=as.data.frame(apply(tmp,2,function(x)gsub('\\s+', '_',x)))

# Array colums
for (i in 1:length(tmp$supplementary_file)){
	tmp$Array[i]=strsplit(strsplit(as.character(tmp$supplementary_file[i]), split="/")[[1]][9], split="_")[[1]][3]
}
# Slide columns
for (i in 1:length(tmp$supplementary_file)){
	tmp$Slide[i]=gsub("_R.*", "", strsplit(as.character(tmp$supplementary_file[i]), split="/")[[1]][9])
}
tmp$sample_acc=rownames(tmp)
head(tmp)

# Write out - lenh này de tao ra bang du lieu moi
write.table(tmp, paste0(raw, "/SampleSheet.tsv"), quote=F, sep="\t", dec=",", row.names=F)
write.table(pData(phenoData(gseInfo[[1]])), paste0(raw, "/SampleSheet_full.tsv"), quote=F, sep="\t", dec=",", row.names=F)

# Downloads
if(type=="raw") {
url=tmp[,grepl("supplementary_file",names(tmp))]
registerDoParallel(8, .libPaths("D:/QUAN/R-4.3.1/library/"))
foreach(i=1:length(url[,1])) %dopar% {
	down1=as.character(url[i,1])
	name1=strsplit(as.character(url[i,1]), "/")[[1]][9]
	name2=strsplit(as.character(url[i,2]), "/")[[1]][9]
	down2=as.character(url[i,2])
	download.file(down1, paste0(raw, "/", name1),quiet=F)
	download.file(down2, paste0(raw, "/", name2),quiet=F)
                                        }
} else {
		print("You are using processed dataset - without idat files")
		if (array == "450k"){

		mSetSqFlt=getGenomicRatioSetFromGEO(gse, path=raw)
		save(mSetSqFlt, file=paste0(raw, "/mSetSqFlt_processed.RData"))
	} else {

		EPIC_getGenomicRatioSetFromGEO = function (GSE = NULL, path = NULL, what = "Beta", mergeManifest = FALSE) {
		gset <- getGEO(gse)
		gset <- gset[[1]]

		ann <- paste0("IlluminaHumanMethylation", array, "anno.ilm10b2.hg19")
		if (!require(ann, character.only = TRUE)) {
		    stop(sprintf("cannot load annotation package %s", ann))
		}

		object <- get(ann)
		gr <- getLocations(object = object, mergeManifest = F, orderByLocation = TRUE)
		locusNames <- names(gr)
		sampleNames(gset) <- gset$title
		common <- intersect(locusNames, fData(gset)$ ID)
		ind1 <- match(common, fData(gset)$ ID)
		ind2 <- match(common, locusNames)
		preprocessing <- c(rg.norm = paste0("See GEO ", gse, " for details"))
		out <- GenomicRatioSet(gr = gr[ind2, ], Beta = exprs(gset)[ind1, , drop = FALSE], colData = as(pData(gset), "DataFrame"), annotation = c(array = paste0("IlluminaHumanMethylation", array), annotation="ilm10b2.hg19"), preprocessMethod = preprocessing)
		
		out
		}

		mSetSqFlt=EPIC_getGenomicRatioSetFromGEO(gse, path=raw)
		save(mSetSqFlt, file=paste0(raw, "/mSetSqFlt_processed.RData"))
	}

}
  
### Final check
listF=list.files(raw, pattern=".idat.gz", full.names=T)
for (idat in listF[1:3]) {
  print(idat)
	check=system(paste0('gzip -t ', idat))
	if(check != 0 ) {
		while(TRUE){
		print(paste0(idat, " was corrupted. Redownload idat..."))
		name=strsplit(idat, "/")[[1]][9]
		unlist(apply(url, 2, function(x) grep(name, x, fixed = TRUE, value=T)))
		down=unlist(apply(url, 2, function(x) grep(name, x, fixed = TRUE, value=T)))
		download.file(down, paste0(raw, "/", name))
		check=system(paste0('gzip -t ', paste0(raw, "/", name)))

		if(check==0) break
		}
	}	

}
