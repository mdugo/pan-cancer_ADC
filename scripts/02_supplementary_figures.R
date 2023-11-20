

library(Biobase)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(ggnewscale)
library(ggridges)
library(ggcorrplot)
library(rcartocolor)
library(RColorBrewer)
library(cowplot)
library(openxlsx)
library(viridis)

#***********************************************************************
# Supplementary Figure 1
#***********************************************************************

adc<-read.xlsx("~/ADC/data/1-s2.0-S0959804923006810-mmc2.xlsx", sheet = 5, startRow = 2)
adc<-adc[!is.na(adc$Gene_Symbol),]
adc<-adc[adc$Drug!="Multiple_ADCs",]
tmp1<-adc[grep(";",adc$Gene_Symbol),]
tmp1$Drug<-gsub("\\+.*","",tmp1$Drug)
tmp1$Gene_Symbol<-gsub(";.*","",tmp1$Gene_Symbol)
tmp1$class_payload<-gsub(";.*","",tmp1$class_payload)
tmp2<-adc[grep(";",adc$Gene_Symbol),]
tmp2$Drug<-gsub(".*\\+","",tmp2$Drug)
tmp2$Gene_Symbol<-gsub(".*;","",tmp2$Gene_Symbol)
tmp2$class_payload<-gsub(".*;","",tmp2$class_payload)
adc<-rbind(adc, tmp1, tmp2)
adc<-adc[-grep(";", adc$Gene_Symbol),]

df<-unique(adc[,c("Drug","Gene_Symbol","class_payload")])
df2<-as.data.frame(table(df$Gene_Symbol, df$class_payload))
df2<-df2[df2$Freq>0,]
ntot<-sort(tapply(df2$Freq, df2$Var1, sum), decreasing = T)
df2$Var1<-factor(df2$Var1, levels=names(ntot))
df2$Var2<-gsub("_", " ", df2$Var2)
g1<-ggplot(data=df2, aes(x=Var1, y=Freq, fill=Var2)) +
  geom_bar(stat="identity") +
  theme_pubr(legend="top") +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=10)) +
  xlab("") +
  ylab("Number of ADC drugs") +
  scale_y_continuous(breaks = seq(0,20,5), limits = c(0,20)) +
  scale_fill_manual(name="", values=c(brewer.pal(9,"Set1"),brewer.pal(8,"Accent")), na.translate = F) +
  guides(fill=guide_legend(ncol=8))


df<-as.data.frame(table(adc2$Gene_Symbol, adc2$Tumor_Site))
df<-df[df$Freq>0,]
colnames(df)<-c("Gene_Symbol","Tumor_Site","N_trials")
totTrials<-aggregate(df$N_trials, by=list(Gene_Symbol=df$Gene_Symbol),sum)
colnames(totTrials)[2]<-"Tot_trials"
df<-merge(df,totTrials,by="Gene_Symbol")
df<-df[order(df$Tot_trials, decreasing=T),]
df$Gene_Symbol<-factor(df$Gene_Symbol, levels=unique(df$Gene_Symbol))

g2<-ggplot(data=df, aes(x=Gene_Symbol, y=N_trials, fill=Tumor_Site)) +
  geom_bar(stat="identity") +
  theme_pubr(legend="top") +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=10)) +
  xlab("") +
  ylab("Number of clinical trials") +
  scale_fill_manual(name="", values=c(carto_pal(12,"Safe")[1:12],carto_pal(12,"Prism")[1:11]), na.translate = F) +
  guides(fill=guide_legend(ncol=9))

pdf("~/ADC/results/supp_figure_1.pdf", width=13, height=11)
plot_grid(g1, g2, ncol=1, labels = c("A","B"), label_size = 18)
dev.off()

#**********************************************
# Supplementary Figure 2
#**********************************************

pairing<-read.xlsx("~/ADC/data/1-s2.0-S0959804923006810-mmc2.xlsx", sheet = 1, startRow = 2)
colnames(pairing)<-c("Cancer","TCGA_study","Primary_Tumor","Solid_Tissue_Normal","Reference_Tissue")
load("~/ADC/data/eSet_TCGA_recount3_tpm.RData")
load("~/ADC/data/eSet_GTEx_recount3_tpm.RData")
nSamplesGTEx<-data.frame(table(gtex_tpm$study))
colnames(nSamplesGTEx)<-c("Reference_Tissue","nSamplesGTEx")
pairing<-merge(pairing,nSamplesGTEx,by="Reference_Tissue",all.x=T)
pairing<-pairing[pairing$Solid_Tissue_Normal>=30 & pairing$nSamplesGTEx>=30,]
subsetting<-unique(c(pairing$Cancer,pairing$Reference_Tissue))

pancan_tpm<-pancan_tpm[,pancan_tpm$tcga.gdc_cases.project.name %in% subsetting]
gtex_tpm<-gtex_tpm[,gtex_tpm$study %in% subsetting]
pairing<-pairing[pairing$TCGA_study%in%pancan_tpm$study & pairing$Reference_Tissue%in%gtex_tpm$study,]

# combine phenodata

pdata_tcga<-pData(pancan_tpm)[,c(1,3,4,6,11)]
colnames(pdata_tcga)[c(3:5)]<-c("sample_id", "project", "sample_type")
pdata_tcga$dataset<-"TCGA"
pdata_gtex<-pData(gtex_tpm)[,c(1,3,9)]
pdata_gtex$project<-pdata_gtex$study
pdata_gtex$sample_type<-"Normal"
pdata_gtex$dataset<-"GTEx"
colnames(pdata_gtex)[3]<-"sample_id"
pdata_gtex<-pdata_gtex[,colnames(pdata_tcga)]
pdata<-rbind(pdata_tcga, pdata_gtex)

# combine expression data

expMat<-cbind(exprs(pancan_tpm), exprs(gtex_tpm))
rm(list=c("pancan_tpm","gtex_tpm", "pdata_gtex", "pdata_tcga"))
gc()
expMat<-log2(expMat+1)

myplots<-vector("list", length=nrow(pairing))
for(i in 1:nrow(pairing)){
  print(i)
  expMat.temp<-expMat[,pdata$study%in%c(pairing$TCGA_study[i],pairing$Reference_Tissue[i])]
  pdata.temp<-pdata[pdata$study%in%c(pairing$TCGA_study[i],pairing$Reference_Tissue[i]),]
  stdv<-sort(rowSds(expMat.temp), decreasing=T)
  expMat.temp<-expMat.temp[names(stdv)[1:2000],]
  pca<-prcomp(x=t(expMat.temp))
  percentVar <- round(100 * base::summary(pca)$importance[2,1:2],1)
  df<-data.frame(PC1=pca$x[,1], PC2=pca$x[,2],sample_type=pdata.temp$sample_type,dataset=pdata.temp$dataset)
  df$sample_type[!df$sample_type%in%c("Normal","Solid Tissue Normal")]<-"Tumor"
  df$sample_type[df$sample_type=="Normal"]<-"Healthy tissue"
  df$sample_type[df$sample_type=="Solid Tissue Normal"]<-"Tumor-paired normal"
  
  gx<-ggplot(data = df, aes(x = PC1, y = PC2, color = sample_type, shape=dataset)) + 
    geom_point(size=2) + 
    theme_pubr(legend="bottom") +
    scale_color_manual(name="",values=carto_pal(12,"Vivid")[1:3]) +
    scale_shape_manual(name="", values=c(4,16)) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    labs(title = paste("Cancer:",pairing$TCGA_study[i], ",\nHealthy tissue:", gsub("_"," ",pairing$Reference_Tissue[i])))
  myplots[[i]]<-gx
}

pdf("~/ADC/results/supp_figure_2.pdf", width=12, height=16)
ggpubr::ggarrange(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]],
                  myplots[[5]], myplots[[6]], myplots[[7]], myplots[[8]],
                  myplots[[9]], myplots[[10]], myplots[[11]], myplots[[12]],
                  common.legend = T,
                  legend = "bottom",
                  align = "hv",
                  ncol = 3,
                  nrow=4)
dev.off()

#****************************************************************************
# Supplementary figure 3
#****************************************************************************

#### TCGA

adc<-read.xlsx("~/ADC/data/1-s2.0-S0959804923006810-mmc2.xlsx", sheet = 5, startRow = 2)
adc<-adc[!is.na(adc$Gene_Symbol),]
adc<-adc[-grep(";", adc$Gene_Symbol),]

load("~/ADC/data/eSet_TCGA_recount3_tpm_GDCmatched_noDup.RData")
tumors <- pancan_tpm[,pancan_tpm$tcga.gdc_cases.samples.sample_type %in%c("Primary Tumor")]
tumors<-tumors[,!tumors$tcga.cgc_sample_sample_type%in%""]
tumors<-tumors[,!is.na(tumors$tcga.cgc_sample_sample_type)]
tumors$sample_barcode<-substr(tumors$tcga.tcga_barcode,1,16)

rppa<-read.table("~/ADC/data/TCGA-RPPA-pancan-clean.txt", sep="\t", header=T, as.is=T, row.names=1)
rppa<-rppa[,-1]
rppa<-t(rppa)
rownames(rppa)<-gsub("^X","",rownames(rppa))
rppaAnnot<-read.table("~/ADC/data/TCGA_antibodies_descriptions.gencode.v36.tsv", sep="\t", header=T, as.is=T, row.names=2)
cg<-intersect(rownames(rppa),rownames(rppaAnnot))
rppaAnnot<-rppaAnnot[cg,]
rppa<-rppa[cg,]
identical(rownames(rppa), rownames(rppaAnnot))
rppaPheno<-data.frame(tcga_barcode=colnames(rppa),sample_barcode=substr(colnames(rppa),1,16), stringsAsFactors = F)
rownames(rppaPheno)<-rppaPheno$tcga_barcode
rppa<-ExpressionSet(assayData = as.matrix(rppa),
                    phenoData = new("AnnotatedDataFrame",rppaPheno),
                    featureData = new("AnnotatedDataFrame",rppaAnnot))
cs<-intersect(rppa$sample_barcode, tumors$sample_barcode)
tumors<-tumors[, match(cs, tumors$sample_barcode)]
rppa<-rppa[, match(cs, rppa$sample_barcode)]
identical(rppa$sample_barcode, tumors$sample_barcode)
rppa<-rppa[fData(rppa)$gene_name%in%adc$Gene_Symbol,]
rppa<-rppa[-grep("_p", rownames(rppa)),]
rownames(rppa)<-fData(rppa)$gene_name
tumors<-tumors[rownames(rppa),]
identical(rownames(rppa), rownames(tumors))
exprs(tumors)<-log2(exprs(tumors)+1)

correl<-data.frame(Symbol=rownames(tumors),spearman=0)
for(i in 1:nrow(tumors)){
  correl$spearman[i]<-cor(exprs(tumors)[i,], exprs(rppa)[i,], use="pairwise.complete.obs")
}
correl<-correl[order(correl$spearman, decreasing = T),]
correl$Symbol<-factor(correl$Symbol, levels=correl$Symbol)
correl$Dataset<-"TCGA"

#### GTEx

load("~/ADC/data/eSet_GTEx_recount3_tpm.RData")
prot<-read.table("~/ADC/data/TableS2D_protein_relative_abundance_Jiang_et_al.txt", sep="\t", header=T, as.is=T)
colnames(prot)[2]<-"ensembl_id"
ids<-read.table("~/ADC/data/geneID_conversion.txt", sep="\t", header=T, as.is=T)
prot<-merge(prot, ids, by="ensembl_id")
prot<-prot[prot$hgnc_symbol%in%adc$Gene_Symbol,]
rownames(prot)<-prot$hgnc_symbol
prot<-prot[,!colnames(prot)%in%c("ensembl_id","gene.id.full","hgnc_symbol")]
phenoProt<-data.frame(assayID=colnames(prot), stringsAsFactors = F)
phenoProt$subjectID<-sapply(strsplit(colnames(prot),"\\."), function(x)paste(x[1],x[2],sep="."))
phenoProt$sampleID<-sapply(strsplit(colnames(prot),"\\."), function(x)paste(x[1],x[2],x[3],sep="."))
rownames(phenoProt)<-phenoProt$assayID
prot<-ExpressionSet(assayData = as.matrix(prot),
                    phenoData = new("AnnotatedDataFrame",phenoProt))

normals <- gtex_tpm[rownames(gtex_tpm) %in% adc$Gene_Symbol,]
normals$sampleID<-sapply(strsplit(normals$gtex.sampid,"-"), function(x)paste(x[1],x[2],x[3],sep="."))
cg<-intersect(rownames(prot), rownames(normals))
cs<-intersect(prot$sampleID,normals$sampleID)
normals<-normals[cg,match(cs,normals$sampleID)]
prot<-prot[cg,match(cs,prot$sampleID)]
identical(rownames(prot), rownames(normals))
identical(prot$sampleID, normals$sampleID)
exprs(normals) <- log2(exprs(normals) + 1)

correl2<-data.frame(Symbol=rownames(normals),spearman=0)
for(i in 1:nrow(normals)){
  correl2$spearman[i]<-cor(exprs(normals)[i,], exprs(prot)[i,], use="pairwise.complete.obs", method="spearman")
}
correl2$Dataset<-"GTEx"
correl2<-correl2[order(correl2$spearman, decreasing = T),]
correl3<-rbind(correl2, correl)
correl3$ID<-correl3$Symbol
correl3$ID[correl3$Dataset=="TCGA"]<-paste0(" ",correl3$ID[correl3$Dataset=="TCGA"])
correl3$ID<-factor(correl3$ID, levels=correl3$ID)

pdf("~/ADC/results/supp_figure_3.pdf", width=13, height=4.5)
ggplot(data=correl3, aes(x=ID, y=spearman, fill=Dataset)) +
  geom_bar(stat="identity") +
  theme_pubr(border=T) + 
  xlab("") +
  ylab("Spearman's correlation") +
  theme(axis.text.x = element_text(angle=45, size=10, hjust=1),
        strip.text.x = element_text(size=14, face="bold")) +
  scale_y_continuous(breaks=seq(-0.4,0.9,0.1)) +
  #geom_hline(yintercept = median(correl2$spearman), linetype="dashed") +
  facet_grid(. ~ Dataset, scales = "free", space='free') +
  #scale_x_discrete(labels=gsub(" .*","",levels(correl3$ID))) +
  scale_fill_manual(name="", values=carto_pal(12,"Vivid")[1:2], guide="none")
dev.off()


#*********************************************************************************
# Supplementary Figure 4
#*********************************************************************************

load("~/ADC/data/eSet_TCGA_recount3_tpm_GDCmatched_noDup.RData")
load("~/ADC/data/eSet_GTEx_recount3_tpm.RData")

adc<-read.xlsx("~/ADC/data/ADC_target_list.xlsx", sheet = 1)
adc<-adc[adc$Tested.in.hematological.malignancies.only==0,]
adc$phase2_count<-sapply(adc$condition.phase2, function(x)length(unlist(strsplit(x, ","))))
adc$phase2_count[is.na(adc$condition.phase2)]<-0
adc$phase3_count<-sapply(adc$condition.phase3, function(x)length(unlist(strsplit(x, ","))))
adc$phase3_count[is.na(adc$condition.phase3)]<-0
adc$approved_count<-sapply(adc$approved_indications, function(x)length(unlist(strsplit(x, ","))))
adc$approved_count[is.na(adc$approved_indications)]<-0
tumors <- pancan_tpm[rownames(pancan_tpm) %in% adc$Gene_Symbol , pancan_tpm$tcga.gdc_cases.samples.sample_type %in%c("Primary Tumor")]
tumors<-tumors[,!tumors$tcga.cgc_sample_sample_type%in%""]
tumors<-tumors[,!is.na(tumors$tcga.cgc_sample_sample_type)]
tumors<-tumors[order(rownames(tumors)),]
dim(tumors)
exprs(tumors) <- log2(exprs(tumors) + 1)
nTumorSamples<-as.data.frame(table(tumors$study))
colnames(nTumorSamples)<-c("Cancer","N")

normals <- gtex_tpm[rownames(gtex_tpm) %in% adc$Gene_Symbol ,]
normals<-normals[order(rownames(normals)),]
dim(normals)
exprs(normals) <- log2(exprs(normals) + 1)
identical(rownames(tumors), rownames(normals))
adc<-adc[match(rownames(tumors), adc$Gene_Symbol),]
identical(rownames(tumors), adc$Gene_Symbol)

th<-apply(exprs(tumors), 1, quantile, probs = 0.80)
detection<-apply(exprs(tumors), 2, function(x)x>th)
detection<-apply(detection, 1, function(x)tapply(x, tumors$study, sum)/table(tumors$study)*100)
df.det<-data.frame(Cancer = rownames(detection), detection)
df.det = reshape::melt(df.det, id = c("Cancer"))
colnames(df.det)[3]<-"Perc_samples"
df.det$ID<-paste(df.det$Cancer, df.det$variable, sep="_")

phaseStudy<-detection
phaseStudy[phaseStudy>=0]<-NA
for(i in 1:ncol(phaseStudy)){
  phaseStudy[rownames(phaseStudy)%in%unlist(strsplit(adc$condition.phase2[i],",")),i]<-"Phase 2"
  phaseStudy[rownames(phaseStudy)%in%unlist(strsplit(adc$condition.phase3[i],",")),i]<-"Phase 3"
  phaseStudy[rownames(phaseStudy)%in%unlist(strsplit(adc$approved_indications[i],",")),i]<-"Approved"
}
df.phase<-data.frame(Cancer = rownames(phaseStudy), phaseStudy)
df.phase = reshape::melt(df.phase, id = c("Cancer"))
colnames(df.phase)[3]<-"Phase"
df.phase$ID<-paste(df.phase$Cancer, df.phase$variable, sep="_")

fc<-data.frame(Cancer=rep(unique(tumors$study), each=nrow(tumors)),
               Symbol=rep(rownames(tumors), length(unique(tumors$study))), logFC=0, p.value=0, Bonferroni=0)
for(i in 1:nrow(tumors)){
  print(i)
  for(j in unique(tumors$study)){
    statTest <- t.test(exprs(tumors)[i,tumors$study==j], exprs(tumors)[i,tumors$study!=j], alternative = "greater")
    fc$logFC[fc$Symbol==rownames(tumors)[i] & fc$Cancer==j] <- -diff(statTest$estimate)
    fc$p.value[fc$Symbol==rownames(tumors)[i] & fc$Cancer==j] <- statTest$p.value
  }
}
fc$Bonferroni<-p.adjust(fc$p.value, method="bonferroni")
fc$logBonferroni <- -log10(fc$Bonferroni)
fc$ID<-paste(fc$Cancer, fc$Symbol, sep="_")

fcMat<-reshape(data = fc[,1:3], idvar="Cancer", timevar="Symbol", direction="wide")
rownames(fcMat)<-fcMat$Cancer
colnames(fcMat)<-gsub("logFC\\.", "", colnames(fcMat))
fcMat<-as.matrix(fcMat[,-1])
pMat<-reshape(data = fc[,c(1,2,5)], idvar="Cancer", timevar="Symbol", direction="wide")
rownames(pMat)<-pMat$Cancer
colnames(pMat)<-gsub("Bonferroni\\.", "", colnames(pMat))
pMat<-as.matrix(pMat[,-1])
fcMat[fcMat<log2(2) | pMat >= 0.05]<-0

hcrow<-hclust(dist(fcMat, method="euclidean"), method="ward.D2")
hccol<-hclust(dist(t(fcMat), method="euclidean"), method="ward.D2")
fcMat<-fcMat[hcrow$order, hccol$order]
df<-data.frame(Cancer = rownames(fcMat), fcMat)
df = reshape::melt(df, id = c("Cancer"))
colnames(df)[2]<-"Symbol"
colnames(df)[3]<-"logFC"
df$ID<-paste(df$Cancer, df$Symbol, sep="_")
df$Index<-1:nrow(df)
df<-merge(df, fc[,-c(1:3)],by="ID",all=T)
df<-merge(df, df.det[,-c(1:2)],by="ID",all=T)
df<-merge(df, df.phase[,-c(1:2)],by="ID",all=T)
df<-merge(df, nTumorSamples, by="Cancer",all=T)
df<-df[order(df$Index),]
df$CancerN<-paste0(df$Cancer, " (N = ", df$N, ")")
df$CancerN<-factor(df$CancerN, levels=unique(df$CancerN))
df$Symbol<-factor(df$Symbol, levels=unique(df$Symbol))
df$Bonferroni_cat<-1
df$Bonferroni_cat[df$Bonferroni<0.05]<-2
df$Bonferroni_cat[df$logFC<0]<-1
df$Bonferroni_cat[df$Bonferroni<0.01]<-3
df$Bonferroni_cat[df$Bonferroni<0.001]<-4
df$Phase<-factor(df$Phase, levels=c("Phase 2","Phase 3","Approved"))
df$All_higher<-factor(ifelse(df$Perc_samples==100,"100%","< 100%"))
df$All_higher_num<-ifelse(df$Perc_samples==100,1,0.25)

pl1<-ggplot(data = df, aes(x = Symbol, y = CancerN)) +
  geom_tile(aes(fill = Phase),color="black") +
  scale_fill_manual(name="Clinical status", breaks=c("Phase 2","Phase 3","Approved"), values=c(carto_pal(7,"Emrld")[c(2,4,6)]),na.value="white") +
  theme_pubr() +
  theme(legend.key=element_blank(),
        axis.text.x = element_text(colour = "black", size = 9, face = "plain", angle = 45, vjust = 1, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "plain", size = 11), 
        legend.text = element_text(size = 10, face ="plain", colour ="black"), 
        legend.title = element_text(size = 12, face = "plain"), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right",
        plot.margin=unit(c(6, 6, 6, 6), "pt")) +
  xlab("") +
  ylab("") +
  new_scale_fill() +
  geom_point(aes(fill = logFC,  size = Perc_samples, color=All_higher, stroke=All_higher_num), shape = 21) +
  scale_size_continuous(name="% samples", range = c(-1,4), breaks=seq(20,100,20), limits = c(0,100)) +
  scale_color_manual(name="", values=c("black","red"), guide="none") +
  binned_scale(aesthetics = "fill",
               scale_name = "stepsn",
               name=expression(log['2']*FC),
               #palette = function(x)inferno(n=5),
               palette = function(x)c("grey90",carto_pal(7,"ag_Sunset")[5:1]),
               breaks = c(0,seq(1, ceiling(max(df$logFC)), 2)),
               labels=c("< 1", seq(1, ceiling(max(df$logFC)), 2)),
               limits = c(0, ceiling(max(df$logFC))),
               guide = guide_colorsteps(frame.colour = "black")
  ) +
  guides(size = guide_legend(override.aes = list(stroke = c(0.25,0.25,0.25,0.25,1), color=c("black","black","black","black","red"))))

pdf("~/ADC/results/supp_figure_4.pdf", width=13, height=6.5)
pl1
dev.off()

#*************************************************************
# Supplementary figure 5
#*************************************************************

comparisonBoxPlot<-function(tumors, normals, gene, show_th=TRUE, tumor_quantile_th=0.8, normal_quantile_th=0.8, pdf.path, plot.width = 12, plot.height = 5, save.pdf=FALSE){
  df.t<-data.frame(Exp=exprs(tumors)[gene,],
                   Dataset="Cancer",
                   Study=tumors$study)
  mediane<-sort(tapply(df.t$Exp, df.t$Study, median), decreasing=T)
  df.t$Study<-factor(df.t$Study, levels=names(mediane))
  df.n<-data.frame(Exp=exprs(normals)[gene,],
                   Dataset="Normals",
                   Study=normals$study)
  mediane<-sort(tapply(df.n$Exp, df.n$Study, median), decreasing=T)
  df.n$Study<-factor(df.n$Study, levels=names(mediane))
  df<-rbind(df.t,df.n)
  df$Fill<-c(as.numeric(df$Study[df$Dataset=="Cancer"]),as.numeric(df$Study[df$Dataset=="Normals"])-31)
  g1<-ggplot(data=df, aes(x=Study,y=Exp, fill=Fill)) +
    geom_boxplot(coef=NULL) +
    scale_fill_gradient(high="#266b6e", low="#c4e6c3") +
    facet_wrap(~Dataset, scales = "free") +
    theme_pubr(border=T, legend="none") +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=8),
          plot.margin = margin(t=8,l=2,r=5,b=2)) +
    xlab("") +
    ylab(paste(gene,"log2(TPM+1)")) +
    ylim(0,round(max(df$Exp))+1)
  if(show_th==TRUE){
    g1 <- g1 +
      geom_hline(yintercept = quantile(df$Exp[df$Dataset=="Cancer"], probs=tumor_quantile_th), linetype="dashed",color="red2") +
      geom_hline(yintercept = quantile(df$Exp[df$Dataset=="Normals"], probs=normal_quantile_th), linetype="dashed",color="royalblue") 
  }
  g2<-ggplot(data=df[df$Dataset=="Healthy tissue",], aes(x=Exp)) +
    geom_density(fill="grey") + 
    geom_hline(yintercept = quantile(df$Exp[df$Dataset=="Healthy tissue"], probs=normal_quantile_th), linetype="dashed",color="royalblue") +
    theme_pubr() +
    theme(axis.title.y=element_blank(),
          #axis.line.y = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.margin = margin(t=27,r=10,l=10,b=46)) +
    xlim(0,round(max(df$Exp))+1) +
    coord_flip()
  return(g1)
}

adc<-read.xlsx("~/ADC/data/1-s2.0-S0959804923006810-mmc2.xlsx", sheet = 5, startRow = 2)
adc<-adc[!is.na(adc$Gene_Symbol),]
adc<-adc[-grep(";", adc$Gene_Symbol),]

load("~/ADC/data/eSet_TCGA_recount3_tpm_GDCmatched_noDup.RData")
load("~/ADC/data/eSet_GTEx_recount3_tpm.RData")
tumors <- pancan_tpm[rownames(pancan_tpm) %in% adc$Gene_Symbol, pancan_tpm$tcga.gdc_cases.samples.sample_type %in%c("Primary Tumor")]
tumors<-tumors[,!tumors$tcga.cgc_sample_sample_type%in%""]
tumors<-tumors[,!is.na(tumors$tcga.cgc_sample_sample_type)]
exprs(tumors) <- log2(exprs(tumors) + 1)
normals <- gtex_tpm[rownames(gtex_tpm) %in% adc$Gene_Symbol ,]
exprs(normals) <- log2(exprs(normals) + 1)

g1<-comparisonBoxPlot(tumors=tumors, normals=normals, gene = "ERBB2", show_th = T) +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=12), axis.title = element_text(size=12), strip.text = element_text(size=15), plot.title = element_text(size=18, face="bold")) +
  ylab(expression(log['2']*"(TPM+1)")) +
  labs(title = "ERBB2")
g2<-comparisonBoxPlot(tumors=tumors, normals=normals, gene = "F3") +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=12), axis.title = element_text(size=12), strip.text = element_text(size=15), plot.title = element_text(size=18, face="bold")) +
  ylab(expression(log['2']*"(TPM+1)")) +
  labs(title = "F3")
g3<-comparisonBoxPlot(tumors=tumors, normals=normals, gene = "FOLR1") +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=12), axis.title = element_text(size=12), strip.text = element_text(size=15), plot.title = element_text(size=18, face="bold")) +
  ylab(expression(log['2']*"(TPM+1)")) +
  labs(title = "FOLR1")
g4<-comparisonBoxPlot(tumors=tumors, normals=normals, gene = "NECTIN4") +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=12), axis.title = element_text(size=12), strip.text = element_text(size=15), plot.title = element_text(size=18, face="bold")) +
  ylab(expression(log['2']*"(TPM+1)")) +
  labs(title = "NECTIN4")
g5<-comparisonBoxPlot(tumors=tumors, normals=normals, gene = "TACSTD2") +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=12), axis.title = element_text(size=12), strip.text = element_text(size=15), plot.title = element_text(size=18, face="bold")) +
  ylab(expression(log['2']*"(TPM+1)")) +
  labs(title = "TACSTD2")

pdf("~/ADC/results/supp_figure_5.pdf", width=13, height=25)
plot_grid(g1, g2, g3, g4, g5, ncol=1)
dev.off()


#*************************************************************
# Supplementary figure 6
#*************************************************************

adc<-read.xlsx("~/ADC/data/1-s2.0-S0959804923006810-mmc2.xlsx", sheet = 5, startRow = 2)
adc<-adc[!is.na(adc$Gene_Symbol),]
adc<-adc[-grep(";", adc$Gene_Symbol),]

load("~/ADC/data/eSet_TCGA_recount3_tpm_GDCmatched_noDup.RData")
load("~/ADC/data/eSet_GTEx_recount3_tpm.RData")
tumors <- pancan_tpm[rownames(pancan_tpm) %in% adc$Gene_Symbol, pancan_tpm$tcga.gdc_cases.samples.sample_type %in%c("Primary Tumor")]
tumors<-tumors[,!tumors$tcga.cgc_sample_sample_type%in%""]
tumors<-tumors[,!is.na(tumors$tcga.cgc_sample_sample_type)]
exprs(tumors) <- log2(exprs(tumors) + 1)
correl<-cor(t(exprs(tumors)), method="spearman")

col<- colorRampPalette(c("slateblue4","white","firebrick3"))(10)

pdf("~/ADC/results/supp_figure_6.pdf", width=14, height=14)
corrplot::corrplot(correl, outline = TRUE,
         order = 'hclust',
         hclust.method = 'ward.D2',
         diag=F,
         method="circle", 
         type="upper", 
         col=col,
         cl.ratio=0.1,
         cl.cex = 1.1,
         pch.cex=6,
         tl.cex=1,
         tl.col="black", 
         tl.srt=45)
dev.off()


#*************************************************************
# Supplementary Figure 7
#*************************************************************

load("~/ADC/data/eSet_GTEx_recount3_tpm.RData")

adc<-read.xlsx("~/ADC/data/1-s2.0-S0959804923006810-mmc2.xlsx", sheet = 5, startRow = 2)
adc<-adc[!is.na(adc$Gene_Symbol),]
adc<-adc[-grep(";", adc$Gene_Symbol),]

normals <- gtex_tpm[rownames(gtex_tpm) %in% adc$Gene_Symbol ,]
normals<-normals[order(rownames(normals)),]
dim(normals)
exprs(normals) <- log2(exprs(normals) + 1)
adc<-adc[match(rownames(normals), adc$Gene_Symbol),]
identical(rownames(normals), adc$Gene_Symbol)

th<-apply(exprs(normals), 1, quantile, probs = 0.80)
detection<-apply(exprs(normals), 2, function(x)x>th)
detection<-apply(detection, 1, function(x)tapply(x, normals$study, sum)/table(normals$study)*100)
df.det<-data.frame(Tissue = rownames(detection), detection)
df.det = reshape::melt(df.det, id = c("Tissue"))
colnames(df.det)[3]<-"Perc_samples"
df.det$ID<-paste(df.det$Tissue, df.det$variable, sep="_")

fc<-data.frame(Tissue=rep(unique(normals$study), each=nrow(normals)),
               Symbol=rep(rownames(normals), length(unique(normals$study))), logFC=0, p.value=0, Bonferroni=0)
for(i in 1:nrow(normals)){
  print(i)
  for(j in unique(normals$study)){
    statTest <- t.test(exprs(normals)[i,normals$study==j], exprs(normals)[i,normals$study!=j], alternative = "greater")
    fc$logFC[fc$Symbol==rownames(normals)[i] & fc$Tissue==j] <- -diff(statTest$estimate)
    fc$p.value[fc$Symbol==rownames(normals)[i] & fc$Tissue==j] <- statTest$p.value
  }
}
fc$Bonferroni<-p.adjust(fc$p.value, method="bonferroni")
fc$logBonferroni <- -log10(fc$Bonferroni)
fc$ID<-paste(fc$Tissue, fc$Symbol, sep="_")

fcMat<-reshape(data = fc[,1:3], idvar="Tissue", timevar="Symbol", direction="wide")
rownames(fcMat)<-fcMat$Tissue
colnames(fcMat)<-gsub("logFC\\.", "", colnames(fcMat))
fcMat<-as.matrix(fcMat[,-1])
pMat<-reshape(data = fc[,c(1,2,5)], idvar="Tissue", timevar="Symbol", direction="wide")
rownames(pMat)<-pMat$Tissue
colnames(pMat)<-gsub("Bonferroni\\.", "", colnames(pMat))
pMat<-as.matrix(pMat[,-1])
fcMat[fcMat<log2(2) | pMat >= 0.05]<-0

hcrow<-hclust(dist(fcMat, method="euclidean"), method="ward.D2")
hccol<-hclust(dist(t(fcMat), method="euclidean"), method="ward.D2")
fcMat<-fcMat[hcrow$order, hccol$order]
df<-data.frame(Tissue = rownames(fcMat), fcMat)
df = reshape::melt(df, id = c("Tissue"))
colnames(df)[2]<-"Symbol"
colnames(df)[3]<-"logFC"
df$ID<-paste(df$Tissue, df$Symbol, sep="_")
df$Index<-1:nrow(df)
df<-merge(df, fc[,-c(1:3)],by="ID",all=T)
df<-merge(df, df.det[,-c(1:2)],by="ID",all=T)
nTissueSamples<-as.data.frame(table(normals$study))
colnames(nTissueSamples)<-c("Tissue","N")
df<-merge(df, nTissueSamples, by="Tissue",all=T)
df<-df[order(df$Index),]
df$TissueN<-paste0(df$Tissue, " (N = ", df$N, ")")
df$TissueN<-factor(df$TissueN, levels=unique(df$TissueN))
df$Symbol<-factor(df$Symbol, levels=unique(df$Symbol))
df$Bonferroni_cat<-1
df$Bonferroni_cat[df$Bonferroni<0.05]<-2
df$Bonferroni_cat[df$logFC<0]<-1
df$Bonferroni_cat[df$Bonferroni<0.01]<-3
df$Bonferroni_cat[df$Bonferroni<0.001]<-4
df$All_higher<-factor(ifelse(df$Perc_samples==100,"100%","< 100%"))
df$All_higher_num<-ifelse(df$Perc_samples==100,1,0.25)

toxic<-read.xlsx("~/ADC/data/ADC_toxicities.xlsx", sheet = 1)
toxic<-toxic[,c(1,7:ncol(toxic))]
toxic<-aggregate(toxic[,-1], by=list(Symbol=toxic$Symbol), max)
toxic<-reshape2::melt(toxic, id="Symbol")
toxic$ID<-paste(toxic$variable, toxic$Symbol, sep="_")
colnames(toxic)[3]<-"Toxicity"
toxic<-toxic[toxic$Toxicity==1,]

df<-merge(df, toxic[3:4], by="ID", all.x=T)
df$Toxicity[!is.na(df$Toxicity)]<-"Frequency of adverse\nreaction > 9%"

pl1<-ggplot(data = df, aes(x = Symbol, y = TissueN)) +
  geom_tile(color="black", aes(fill=Toxicity)) +
  scale_fill_manual(name="", values = c("olivedrab2"), na.translate = F) +
  #scale_fill_manual(name="Clinical status", breaks=c("Phase 2","Phase 3","Approved"), values=c(carto_pal(7,"Emrld")[c(2,4,6)]),na.value="white") +
  theme_pubr() +
  theme(legend.key=element_blank(),
        axis.text.x = element_text(colour = "black", size = 9, face = "plain", angle = 45, vjust = 1, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "plain", size = 11), 
        legend.text = element_text(size = 11, face ="plain", colour ="black"), 
        legend.title = element_text(size = 12, face = "plain"), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right",
        plot.margin=unit(c(6, 6, 6, 6), "pt")) +
  xlab("") +
  ylab("") +
  new_scale_fill() +
  geom_point(aes(fill = logFC,  size = Perc_samples, color=All_higher, stroke=All_higher_num), shape = 21) +
  scale_size_continuous(name="% samples", range = c(-1,4), breaks=seq(20,100,20), limits = c(0,100)) +
  scale_color_manual(name="", values=c("black","red"), guide="none") +
  binned_scale(aesthetics = "fill",
               scale_name = "stepsn",
               name=expression(log['2']*FC),
               #palette = function(x)inferno(n=5),
               palette = function(x)c("grey90",carto_pal(7,"ag_Sunset")[5:1]),
               breaks = c(0,seq(1, ceiling(max(df$logFC)), 2)),
               limits = c(0, ceiling(max(df$logFC))),
               labels=c("< 1", seq(1, ceiling(max(df$logFC)), 2)),
               guide = guide_colorsteps(frame.colour = "black")
  ) +
  guides(size = guide_legend(override.aes = list(stroke = c(0.25,0.25,0.25,0.25,1), color=c("black","black","black","black","red"))))

pdf("~/ADC/results/supp_figure_7.pdf", width=13, height=5.9)
pl1
dev.off()


#*************************************************************
# Supplementary Figure 8
#*************************************************************

load("~/ADC/data/eSet_GTEx_recount3_tpm.RData")

resMark<-read.xlsx("~/ADC/data/1-s2.0-S0959804923006810-mmc2.xlsx", sheet = 8, startRow = 2)
resMark<-resMark[!is.na(resMark$Gene.Symbol),]
resMark<-resMark[resMark$Mechanism!="U",]
colnames(resMark)[colnames(resMark)=="Gene.Symbol"]<-"Symbol"

normals <- gtex_tpm[rownames(gtex_tpm) %in% resMark$Symbol ,]
normals<-normals[order(rownames(normals)),]
exprs(normals) <- log2(exprs(normals) + 1)

th<-apply(exprs(normals), 1, quantile, probs = 0.80)
detection<-apply(exprs(normals), 2, function(x)x>th)
detection<-apply(detection, 1, function(x)tapply(x, normals$study, sum)/table(normals$study)*100)
df.det<-data.frame(Tissue = rownames(detection), detection)
df.det = reshape::melt(df.det, id = c("Tissue"))
colnames(df.det)[3]<-"Perc_samples"
df.det$ID<-paste(df.det$Tissue, df.det$variable, sep="_")

fc<-data.frame(Tissue=rep(unique(normals$study), each=nrow(normals)),
               Symbol=rep(rownames(normals), length(unique(normals$study))), logFC=0, p.value=0, Bonferroni=0)
for(i in 1:nrow(normals)){
  print(i)
  for(j in unique(normals$study)){
    statTest <- t.test(exprs(normals)[i,normals$study==j], exprs(normals)[i,normals$study!=j], alternative = "greater")
    fc$logFC[fc$Symbol==rownames(normals)[i] & fc$Tissue==j] <- -diff(statTest$estimate)
    fc$p.value[fc$Symbol==rownames(normals)[i] & fc$Tissue==j] <- statTest$p.value
  }
}
fc$Bonferroni<-p.adjust(fc$p.value, method="bonferroni")
fc$logBonferroni <- -log10(fc$Bonferroni)
fc$ID<-paste(fc$Tissue, fc$Symbol, sep="_")

fcMat<-reshape(data = fc[,1:3], idvar="Tissue", timevar="Symbol", direction="wide")
rownames(fcMat)<-fcMat$Tissue
colnames(fcMat)<-gsub("logFC\\.", "", colnames(fcMat))
fcMat<-as.matrix(fcMat[,-1])
pMat<-reshape(data = fc[,c(1,2,5)], idvar="Tissue", timevar="Symbol", direction="wide")
rownames(pMat)<-pMat$Tissue
colnames(pMat)<-gsub("Bonferroni\\.", "", colnames(pMat))
pMat<-as.matrix(pMat[,-1])
fcMat[fcMat<log2(2) | pMat >= 0.05]<-0

hcrow<-hclust(dist(fcMat, method="euclidean"), method="ward.D2")
hccol<-hclust(dist(t(fcMat), method="euclidean"), method="ward.D2")
fcMat<-fcMat[hcrow$order, hccol$order]
df<-data.frame(Tissue = rownames(fcMat), fcMat)
df = reshape::melt(df, id = c("Tissue"))
colnames(df)[2]<-"Symbol"
colnames(df)[3]<-"logFC"
df$ID<-paste(df$Tissue, df$Symbol, sep="_")
df$Index<-1:nrow(df)
df<-merge(df, fc[,-c(1:3)],by="ID",all=T)
df<-merge(df, df.det[,-c(1:2)],by="ID",all=T)
nTissueSamples<-as.data.frame(table(normals$study))
colnames(nTissueSamples)<-c("Tissue","N")
df<-merge(df, nTissueSamples, by="Tissue",all=T)
df<-merge(df, resMark[,2:6], by="Symbol")
df$`Non-cleavable.linker.specific`<-factor(df$`Non-cleavable.linker.specific`, levels=c("No","Yes","UK"))
df$Mechanism<-ifelse(df$Mechanism=="R","Resistance","Sensitivity")
df<-df[order(df$Index),]
df$TissueN<-paste0(df$Tissue, " (N = ", df$N, ")")
df$TissueN<-factor(df$TissueN, levels=unique(df$TissueN))
df$Symbol<-factor(df$Symbol, levels=unique(df$Symbol))
df$Bonferroni_cat<-1
df$Bonferroni_cat[df$Bonferroni<0.05]<-2
df$Bonferroni_cat[df$logFC<0]<-1
df$Bonferroni_cat[df$Bonferroni<0.01]<-3
df$Bonferroni_cat[df$Bonferroni<0.001]<-4
df$All_higher<-factor(ifelse(df$Perc_samples==100,"100%","< 100%"))
df$All_higher_num<-ifelse(df$Perc_samples==100,1,0.25)

toxic<-read.xlsx("~/ADC/data/ADC_toxicities.xlsx", sheet = 1)
toxic<-toxic[,c(3,7:ncol(toxic))]
ngenes<-unique(c(unlist(strsplit(toxic$Predictor[1],",")), unlist(strsplit(toxic$Predictor[2],","))))
myMat<-matrix(0, nrow=length(ngenes), ncol=(ncol(toxic)-1))
rownames(myMat)<-ngenes
colnames(myMat)<-colnames(toxic)[-1]
for(i in 1:nrow(myMat)){
  myMat[i,]<-ifelse(colSums(toxic[grep(rownames(myMat)[i],toxic$Predictor),-1])>0,1,0)
}
myMat<-data.frame(Symbol=rownames(myMat), myMat)
myMat<-reshape2::melt(myMat, id="Symbol")
myMat$ID<-paste(myMat$variable, myMat$Symbol, sep="_")
colnames(myMat)[3]<-"Adverse"
myMat$Adverse<-ifelse(myMat$Adverse==1,"Frequency of adverse\nreaction > 9%", NA)
df<-merge(df, myMat[,3:4], by="ID", all.x=T)

pl1<-ggplot(data = df, aes(x = Symbol, y = TissueN)) +
  geom_tile(color="black", aes(fill=Adverse)) +
  scale_fill_manual(name="", values = c("olivedrab2"), na.translate = F) +
  #scale_fill_manual(name="Clinical status", breaks=c("Phase 2","Phase 3","Approved"), values=c(carto_pal(7,"Emrld")[c(2,4,6)]),na.value="white") +
  theme_pubr() +
  facet_grid(. ~ Mechanism, scales = "free", space='free') +
  theme(legend.key=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), 
        #axis.text.x = element_text(angle=45, hjust=1), 
        axis.text.y = element_text(colour = "black", face = "plain", size = 11), 
        legend.text = element_text(size = 10, face ="plain", colour ="black"), 
        legend.title = element_text(size = 12, face = "plain"), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "top",
        plot.margin=unit(c(0.5, 0.25, 0, 0.15), "cm"),
        strip.text.x = element_text(size=14, face="bold")) +
  xlab("") +
  ylab("") +
  new_scale_fill() +
  geom_point(aes(fill = logFC,  size = Perc_samples, color=All_higher, stroke=All_higher_num), shape = 21) +
  scale_size_continuous(name="% samples", range = c(-1,4), breaks=seq(20,100,20), limits = c(0,100)) +
  scale_color_manual(name="", values=c("black","red"), guide="none") +
  binned_scale(aesthetics = "fill",
               scale_name = "stepsn",
               name=expression(log['2']*FC),
               #palette = function(x)inferno(n=5),
               palette = function(x)c("grey90",carto_pal(7,"ag_Sunset")[5:1]),
               breaks = c(0,seq(1, ceiling(max(df$logFC)), 2)),
               limits = c(0, ceiling(max(df$logFC))),
               labels=c("< 1", seq(1, ceiling(max(df$logFC)), 2)),
               guide = guide_colorsteps(frame.colour = "black")
  ) +
  guides(size = guide_legend(override.aes = list(stroke = c(0.25,0.25,0.25,0.25,1), color=c("black","black","black","black","red"))))


tmp1<-resMark[,c(2,4)]
colnames(tmp1)[2]<-"value"
tmp1$variable<-"Potential bystander effect"
tmp2<-resMark[,c(2,5)]
colnames(tmp2)[2]<-"value"
tmp2$variable<-"Mechanism"
tmp3<-resMark[,c(2,6)]
colnames(tmp3)[2]<-"value"
tmp3$variable<-"Non cleavable linker specific"
tmp<-rbind(tmp1,tmp2, tmp3)
tmp$Symbol<-factor(tmp$Symbol, levels=levels(df$Symbol))
tmp$value[tmp$value=="UK"]<-"Undetermined"
tmp$value<-factor(tmp$value, levels=c("Yes","No","Undetermined",sort(unique(df$Process))))
tmp<-merge(tmp, df[,c("Symbol","Mechanism")], by="Symbol")


pl2<-ggplot(data = tmp, aes(x = Symbol, y = variable)) +
  geom_tile(aes(fill = value),color="black") +
  theme_pubr() +
  xlab("") +
  ylab("") +
  scale_y_discrete(expand=c(0,0)) +
  scale_fill_manual(name="", values=c("grey30","grey80", "white", carto_pal(12, "Prism")[1:8])) +
  facet_grid(. ~ Mechanism, scales = "free", space='free') +
  theme(legend.key=element_blank(),
        axis.text.x = element_text(colour = "black", size = 10, face = "plain", angle = 45, vjust = 1, hjust = 1),
        axis.ticks = element_blank(),
        #axis.text.y = element_blank(), 
        legend.text = element_text(size = 10, face ="plain", colour ="black"), 
        legend.title = element_text(size = 12, face = "plain"), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "bottom",
        plot.margin=unit(c(-0.5, 0.25, 0.2, 0.15), "cm"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  guides(fill=guide_legend(ncol=4))

pdf("~/ADC/results/supp_figure_8.pdf", width=13, height=9.6)
plot_grid(pl1, pl2,ncol=1, rel_heights = c(1,0.32))
dev.off()

