
options(java.parameters = "-Xmx12000m")

library(recount3)
library(Biobase)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(openxlsx)

dir.create("~/ADC/data", recursive=T)

ap<-available_projects(organism = "human")
gtex<-subset(ap, file_source == "gtex" & project_type == "data_sources")
gtex<-gtex[gtex$project!="STUDY_NA",]
tcga<-subset(ap, file_source == "tcga" & project_type == "data_sources")

#************************************
# download GTEx transcriptomic data
#************************************

gtexList.tpm<-vector("list", length=length(gtex$project))
names(gtexList.tpm)<-gtex$project

gtex2<-gtex[, c("project", "organism", "project_home")]

for(i in 1:nrow(gtex2)){
  print(paste(i, "/", nrow(gtex2)))
  x <- create_rse(project_info = gtex2[i,],
                  type = "gene",
                  annotation = "gencode_v29")
  
  tpm<-getTPM(x, length_var = "bp_length")
  tpm.coll<-aggregate(tpm,by=list(Symbol=rowData(x)$gene_name),sum)
  rownames(tpm.coll)<-tpm.coll$Symbol
  tpm.coll<-as.matrix(tpm.coll[,-1])
  
  pdata<-as.data.frame(colData(x))
  fdata<-data.frame(Symbol=rownames(raw.coll), stringsAsFactors = F)
  rownames(fdata)<-fdata$Symbol
  
  gtexList.tpm[[i]]<-ExpressionSet(assayData=as.matrix(tpm.coll),
                                   phenoData=new("AnnotatedDataFrame",pdata),
                                   featureData=new("AnnotatedDataFrame",fdata))
}

# generate TPM ExpressionSet object

pdatacol<-vector("list", length=length(gtexList.tpm))
names(pdatacol)<-names(gtexList.tpm)

pdatacol<-lapply(gtexList.tpm, function(x)colnames(pData(x)))
pdatacol<-Reduce(intersect, pdatacol)

expMat<-exprs(gtexList.tpm[[1]])
pdata<-pData(gtexList.tpm[[1]])
fdata<-data.frame(Symbol=rownames(expMat), stringsAsFactors = F)
entrez<-unlist(mapIds(org.Hs.eg.db, unique(fdata$Symbol),"ENTREZID", "SYMBOL"))
entrez<-data.frame(Symbol=names(entrez), Entrez=entrez, stringsAsFactors = F)
fdata<-merge(fdata, entrez, by="Symbol", all.x=T)
rownames(fdata)<-fdata$Symbol
fdata<-fdata[rownames(expMat),]

for(i in 2:length(gtexList.tpm)){
  print(i)
  expMat<-cbind(expMat, exprs(gtexList.tpm[[i]]))
  pdataTemp<-pData(gtexList.tpm[[i]])
  pdata<-rbind(pdata, pdataTemp)
}

gtex_tpm<-ExpressionSet(assayData=as.matrix(expMat),
                          phenoData=new("AnnotatedDataFrame",pdata),
                          featureData=new("AnnotatedDataFrame",fdata))
save(gtex_tpm, file="~/ADC/data/eSet_GTEx_recount3_tpm.RData")

#************************************
# download GTEx proteomic data
#************************************

download.file(url = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867420310783-mmc3.xlsx",
              destfile="~/ADC/data/TableS2_Jiang_et_al.xlsx")

df<-read.xlsx(xlsxFile = "~/ADC/data/TableS2_Jiang_et_al.xlsx", sheet = "D protein relative abundance", startRow = 4)
write.table(df, file="~/ADC/data/TableS2D_protein_relative_abundance_Jiang_et_al.txt",sep="\t", row.names=F, quote=F)

df<-read.xlsx(xlsxFile = "~/ADC/data/TableS2_Jiang_et_al.xlsx", sheet = "G protein TS score", startRow = 3)
write.table(df[,c("ensembl_id","hgnc_symbol")], file="~/ADC/data/geneID_conversion.txt",sep="\t", row.names=F, quote=F)


#************************************
# download TCGA transcriptomic data
#************************************

tcga.tpm<-vector("list", length=length(tcga$project))
names(tcga.tpm)<-tcga$project

tcga2<-tcga[, c("project", "organism", "project_home")]

for(i in 1:nrow(tcga2)){
  print(paste(i, "/", nrow(tcga2)))
  x <- create_rse(project_info = tcga2[i,],
                  type = "gene",
                  annotation = "gencode_v29")
  
  tpm<-getTPM(x, length_var = "bp_length")
  tpm.coll<-aggregate(tpm,by=list(Symbol=rowData(x)$gene_name),sum)
  rownames(tpm.coll)<-tpm.coll$Symbol
  tpm.coll<-as.matrix(tpm.coll[,-1])
  
  pdata<-as.data.frame(colData(x))
  fdata<-data.frame(Symbol=rownames(raw.coll), stringsAsFactors = F)
  rownames(fdata)<-fdata$Symbol
  
  tcga.tpm[[i]]<-ExpressionSet(assayData=as.matrix(tpm.coll),
                                   phenoData=new("AnnotatedDataFrame",pdata),
                                   featureData=new("AnnotatedDataFrame",fdata))
}

# generate TPM ExpressionSet object

pdatacol<-vector("list", length=length(tcga.tpm))
names(pdatacol)<-names(tcga.tpm)

pdatacol<-lapply(tcga.tpm, function(x)colnames(pData(x)))
pdatacol<-Reduce(intersect, pdatacol)

expMat<-exprs(tcga.tpm[[1]])
pdata<-pData(tcga.tpm[[1]])
fdata<-data.frame(Symbol=rownames(expMat), stringsAsFactors = F)
entrez<-unlist(mapIds(org.Hs.eg.db, unique(fdata$Symbol),"ENTREZID", "SYMBOL"))
entrez<-data.frame(Symbol=names(entrez), Entrez=entrez, stringsAsFactors = F)
fdata<-merge(fdata, entrez, by="Symbol", all.x=T)
rownames(fdata)<-fdata$Symbol
fdata<-fdata[rownames(expMat),]

for(i in 2:length(tcga.tpm)){
  print(i)
  expMat<-cbind(expMat, exprs(tcga.tpm[[i]]))
  pdataTemp<-pData(tcga.tpm[[i]])
  pdata<-rbind(pdata, pdataTemp)
}

pancan_tpm<-ExpressionSet(assayData=as.matrix(expMat),
                           phenoData=new("AnnotatedDataFrame",pdata),
                           featureData=new("AnnotatedDataFrame",fdata))

pancan_tpm$tcga.gdc_cases.samples.sample_type[pancan_tpm$tcga.gdc_cases.samples.sample_type=="Primary solid Tumor"] <- "Primary Tumor"
pancan_tpm$tcga.gdc_cases.project.name[pancan_tpm$tcga.gdc_cases.project.name=="ACUTE MYELOID LEUKEMIA"]<-"Acute Myeloid Leukemia"
pancan_tpm$tcga.gdc_cases.project.name[pancan_tpm$tcga.gdc_cases.project.name=="BREAST INVASIVE CARCINOMA"]<-"Breast Invasive Carcinoma"
pancan_tpm$tcga.gdc_cases.project.name[pancan_tpm$tcga.gdc_cases.project.name=="KIDNEY RENAL CLEAR CELL CARCINOMA"]<-"Kidney Renal Clear Cell Carcinoma"

# remove duplicate samples

sum(table(pancan_tpm$tcga.tcga_barcode)>1)
dup<-names(which(table(pancan_tpm$tcga.tcga_barcode)>1))

pancan_dup<-pancan_tpm[,pancan_tpm$tcga.tcga_barcode%in%dup]
toex<-data.frame(pData(pancan_dup)[,c("rail_id","tcga.tcga_barcode")], Mean=0)
toex$Mean<-colMeans(exprs(pancan_dup))
toKeep<-NULL
for(i in unique(toex$tcga.tcga_barcode)){
  temp<-toex[toex$tcga.tcga_barcode==i,]
  highD<-temp$rail_id[temp$Mean==max(temp$Mean)]
  if(length(highD)>1){
    highD<-max(highD)
    toKeep<-c(toKeep,highD)
  } else {
    toKeep<-c(toKeep,highD)
  }
}
toExclude<-pancan_dup$rail_id[!pancan_dup$rail_id%in%toKeep]

pancan_tpm<-pancan_tpm[,!pancan_tpm$rail_id%in%toExclude]
save(pancan_tpm, file="~/ADC/data/rnaseq/eSet_TCGA_recount3_tpm.RData")

#************************************
# download TCGA proteomic data
#************************************

download.file(url = "http://api.gdc.cancer.gov/data/fcbb373e-28d4-4818-92f3-601ede3da5e1",
              destfile="~/ADC/data/TCGA-RPPA-pancan-clean.txt")

download.file(url = "https://api.gdc.cancer.gov/v0/data/62647302-b4d3-4a81-a7c0-d141f5dbd300",
              destfile="~/ADC/data/TCGA_antibodies_descriptions.gencode.v36.tsv")

#***************************************************
# download Supplementary Tables from Bosi et al
#***************************************************

download.file(url = "https://ars.els-cdn.com/content/image/1-s2.0-S0959804923006810-mmc2.xlsx",
              destfile="~/ADC/data/1-s2.0-S0959804923006810-mmc2.xlsx")

