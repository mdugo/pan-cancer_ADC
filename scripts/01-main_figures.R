
library(Biobase)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(ggnewscale)
library(ggridges)
library(geomtextpath)
library(rcartocolor)
library(RColorBrewer)
library(cowplot)
library(openxlsx)
library(viridis)

dir.create("~/ADC/results", recursive=T)

#************************************
# Figure 1
#************************************

load("/Users/dugo.matteo/Desktop/hsr_bioinfo/databases/TCGA/pancancer/data/rnaseq/eSet_TCGA_recount3_tpm_GDCmatched_noDup.RData")
load("/Users/dugo.matteo/Desktop/hsr_bioinfo/databases/GTEx/data/eSet_GTEx_recount3_tpm.RData")

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

th<-apply(exprs(normals), 1, quantile, probs = 0.80)
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
    statTest <- t.test(exprs(tumors)[i,tumors$study==j], exprs(normals)[i,], alternative = "greater")
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
               palette = function(x)c("grey90",carto_pal(7,"ag_Sunset")[5:1]),
               breaks = c(0,seq(1, ceiling(max(df$logFC)), 2)),
               limits = c(0, ceiling(max(df$logFC))),
               labels=c("< 1", seq(1, ceiling(max(df$logFC)), 2)),
               guide = guide_colorsteps(frame.colour = "black", label.theme = element_text(size=10))
  ) +
  guides(size = guide_legend(override.aes = list(stroke = c(0.25,0.25,0.25,0.25,1), color=c("black","black","black","black","red"))))


diseaseCol<-structure(c(carto_pal(12,"Safe")[-12],carto_pal(12,"Vivid")[-12],carto_pal(12,"Bold")[-c(10:12)],"grey80"), names=c("ACC","BLCA","LGG","BRCA","CESC", "CHOL", "COAD", "ESCA", "GBM",
                                                                                                                                "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "MESO",
                                                                                                                                "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
                                                                                                                                "TGCT", "THYM", "THCA", "UCS", "UCEC", "UVM", "Normals"))

# density ERBB2

df2<-data.frame(Expression=c(exprs(tumors)["ERBB2",tumors$study%in%c("BRCA","UCEC","STAD","BLCA")],exprs(normals)["ERBB2",]),
                Tissue=c(tumors$study[tumors$study%in%c("BRCA","UCEC","STAD","BLCA")],rep("Normals",ncol(normals))))
mediane<-sort(tapply(df2$Expression, df2$Tissue, median), decreasing=F)
df2$Tissue<-factor(df2$Tissue, levels=names(mediane))
pl2<-ggplot(data=df2, aes(x=Expression, y=Tissue, fill=Tissue)) +
  geom_density_ridges(alpha=0.5) +
  scale_fill_manual(values=diseaseCol[levels(df2$Tissue)], guide="none") +
  theme_pubr() +
  xlab(expression(log['2']*"(TPM+1)")) +
  ylab("") +
  geom_vline(xintercept = quantile(df2$Expression[df2$Tissue=="Normals"], probs=0.8), linetype="dashed") +
  labs(title = "ERBB2") +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=10))

# density FOLR1

df2<-data.frame(Expression=c(exprs(tumors)["FOLR1",tumors$study%in%c("OV","LUAD","THCA","KIRC")],exprs(normals)["FOLR1",]),
                Tissue=c(tumors$study[tumors$study%in%c("OV","LUAD","THCA","KIRC")],rep("Normals",ncol(normals))))
mediane<-sort(tapply(df2$Expression, df2$Tissue, median), decreasing=F)
df2$Tissue<-factor(df2$Tissue, levels=names(mediane))
pl3<-ggplot(data=df2, aes(x=Expression, y=Tissue, fill=Tissue)) +
  geom_density_ridges(alpha=0.5) +
  scale_fill_manual(values=diseaseCol[levels(df2$Tissue)], guide="none") +
  theme_pubr() +
  xlab(expression(log['2']*"(TPM+1)")) +
  ylab("") +
  geom_vline(xintercept = quantile(df2$Expression[df2$Tissue=="Normals"], probs=0.8), linetype="dashed") +
  labs(title = "FOLR1") +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=10))

# density NECTIN4

df2<-data.frame(Expression=c(exprs(tumors)["NECTIN4",tumors$study%in%c("BLCA","BRCA","HNSC","CESC")],exprs(normals)["NECTIN4",]),
                Tissue=c(tumors$study[tumors$study%in%c("BLCA","BRCA","HNSC","CESC")],rep("Normals",ncol(normals))))
mediane<-sort(tapply(df2$Expression, df2$Tissue, median), decreasing=F)
df2$Tissue<-factor(df2$Tissue, levels=names(mediane))
pl4<-ggplot(data=df2, aes(x=Expression, y=Tissue, fill=Tissue)) +
  geom_density_ridges(alpha=0.5) +
  scale_fill_manual(values=diseaseCol[levels(df2$Tissue)], guide="none") +
  theme_pubr() +
  xlab(expression(log['2']*"(TPM+1)")) +
  ylab("") +
  geom_vline(xintercept = quantile(df2$Expression[df2$Tissue=="Normals"], probs=0.8), linetype="dashed") +
  labs(title = "NECTIN4") +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=10))

# density TACSTD2

df2<-data.frame(Expression=c(exprs(tumors)["TACSTD2",tumors$study%in%c("BRCA","PAAD","HNSC","BLCA")],exprs(normals)["TACSTD2",]),
                Tissue=c(tumors$study[tumors$study%in%c("BRCA","PAAD","HNSC","BLCA")],rep("Normals",ncol(normals))))
mediane<-sort(tapply(df2$Expression, df2$Tissue, median), decreasing=F)
df2$Tissue<-factor(df2$Tissue, levels=names(mediane))
pl5<-ggplot(data=df2, aes(x=Expression, y=Tissue, fill=Tissue)) +
  geom_density_ridges(alpha=0.5) +
  scale_fill_manual(values=diseaseCol[levels(df2$Tissue)], guide="none") +
  theme_pubr() +
  xlab(expression(log['2']*"(TPM+1)")) +
  ylab("") +
  geom_vline(xintercept = quantile(df2$Expression[df2$Tissue=="Normals"], probs=0.8), linetype="dashed") +
  labs(title = "TACSTD2") +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=10))

plDensity<-plot_grid(pl2, pl3, pl4, pl5, ncol=4)


pdf("~/ADC/results/figure_1.pdf", width=13, height=11)
plot_grid(pl1, plDensity, labels = c("A","B"), label_size = 18, ncol=1, rel_heights = c(1, 0.55))
dev.off()

#******************************************************************
# Figure 2
#******************************************************************

adc<-read.xlsx("~/ADC/data/ADC_target_list.xlsx", sheet = 1)
adc<-adc[adc$Tested.in.hematological.malignancies.only==0,]

fc<-read.table("~/ADC/data/PrimaryVsMet_ADC_target_logFCscaled.txt", header=T, sep="\t", as.is=T, row.names=1)
fdr<-read.table("~/ADC/data/PrimaryVsMet_ADC_target_FDR.txt", header=T, sep="\t", as.is=T)
colnames(fc)<-c("Breast","Kidney","Colon","Liver","Pancreas","Ovarian","Esophagus","Skin")
colnames(fdr)<-c("Symbol","Breast","Kidney","Colon","Liver","Pancreas","Ovarian","Esophagus","Skin")
fc<-fc[rownames(fc)%in%adc$Gene_Symbol,]
fdr<-fdr[fdr$Symbol%in%adc$Gene_Symbol,]
colSums(fdr[,-1]<0.05)/nrow(fdr)*100
rowSums(fdr[,-1]<0.05)/ncol(fdr)*100
hcrow<-hclust(dist(fc, method="euclidean"), method="ward.D2")
hccol<-hclust(dist(t(fc), method="euclidean"), method="ward.D2")
fc<-fc[hcrow$order, hccol$order]
fc<-data.frame(Symbol=rownames(fc),fc)

df.fc<-reshape2::melt(fc, id="Symbol")
df.fc$Index<-1:nrow(df.fc)
colnames(df.fc)[3]<-"logFC"
df.fc$ID<-paste(df.fc$Symbol, df.fc$variable,sep="_")
df.fdr<-reshape2::melt(fdr, id="Symbol")
colnames(df.fdr)[3]<-"FDR"
df.fdr$ID<-paste(df.fdr$Symbol, df.fdr$variable,sep="_")
df<-merge(df.fc, df.fdr[,c(3,4)], by="ID")
df$Significant<-ifelse(df$FDR<0.05,"FDR < 0.05","n.s.")
df$Significant_num<-ifelse(df$FDR<0.05,1.4,0.3)
df<-df[order(df$Index),]
df$Symbol<-factor(df$Symbol, levels=unique(df$Symbol))
df$variable<-factor(df$variable, levels=unique(df$variable))

pl1<-ggplot(data = df, aes(x = Symbol, y = variable)) +
  geom_tile(color="black", fill="white") +
  theme_pubr() +
  theme(legend.key=element_blank(),
        axis.text.x = element_text(colour = "black", size = 11, face = "plain", angle = 45, vjust = 1, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "plain", size = 13), 
        legend.text = element_text(size = 11, face ="plain", colour ="black", angle=0), 
        legend.title = element_text(size = 12, face = "plain"), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "bottom",
        plot.margin=unit(c(6, 6, 6, 6), "pt")) +
  xlab("") +
  ylab("") +
  new_scale_fill() +
  geom_point(aes(fill = logFC,  size = abs(logFC), stroke=Significant_num), shape = 21) +
  scale_size_continuous(name="Standardised\n|log(FC)|", range = c(0,4), breaks=seq(1,4,1), limits = c(0,4)) +
  scale_color_manual(name="", values=c("black","red")) +
  binned_scale(aesthetics = "fill",
               scale_name = "stepsn",
               name="Standardised\nlog(FC)",
               palette = function(x)c(brewer.pal(9,"Blues")[8:2],brewer.pal(9,"RdPu")[2:8]),
               breaks = seq(-3.5, 3.5, 0.5),
               labels = c("","-3","","-2","","-1","","0","","1","","2","","3",""),
               limits = c(-3.5, 3.5),
               #show.limits = TRUE,
               guide = guide_colorsteps(frame.colour = "black", label.theme = element_text(angle=90, size=11, hjust=0.5, vjust=0.5)))

pdf("~/ADC/results/figure_2.pdf", width=12, height=3.6)
pl1
dev.off()



#*********************************************************************************
# Figure 3
#*********************************************************************************

load("/Users/dugo.matteo/Desktop/hsr_bioinfo/databases/TCGA/pancancer/data/rnaseq/eSet_TCGA_recount3_tpm_GDCmatched_noDup.RData")
load("/Users/dugo.matteo/Desktop/hsr_bioinfo/databases/GTEx/data/eSet_GTEx_recount3_tpm.RData")

resMark<-read.xlsx("~/ADC/data/1-s2.0-S0959804923006810-mmc2.xlsx", sheet = 8, startRow = 2)
resMark<-resMark[!is.na(resMark$Gene.Symbol),]
resMark<-resMark[resMark$Mechanism!="U",]
colnames(resMark)[colnames(resMark)=="Gene.Symbol"]<-"Symbol"
tumors <- pancan_tpm[rownames(pancan_tpm) %in% resMark$Symbol , pancan_tpm$tcga.gdc_cases.samples.sample_type %in%c("Primary Tumor")]
tumors<-tumors[,!tumors$tcga.cgc_sample_sample_type%in%""]
tumors<-tumors[,!is.na(tumors$tcga.cgc_sample_sample_type)]
tumors<-tumors[order(rownames(tumors)),]
dim(tumors)
exprs(tumors) <- log2(exprs(tumors) + 1)
nTumorSamples<-as.data.frame(table(tumors$study))
colnames(nTumorSamples)<-c("Cancer","N")

normals <- gtex_tpm[rownames(gtex_tpm) %in% resMark$Symbol ,]
normals<-normals[order(rownames(normals)),]
dim(normals)
identical(rownames(tumors), rownames(normals))
exprs(normals) <- log2(exprs(normals) + 1)

th<-apply(exprs(normals), 1, quantile, probs = 0.8)
detection<-apply(exprs(tumors), 2, function(x)x>th)
detection<-apply(detection, 1, function(x)tapply(x, tumors$study, sum)/table(tumors$study)*100)
df.det<-data.frame(Cancer = rownames(detection), detection)
df.det = reshape::melt(df.det, id = c("Cancer"))
colnames(df.det)[3]<-"Perc_samples"
df.det$ID<-paste(df.det$Cancer, df.det$variable, sep="_")

fc<-data.frame(Cancer=rep(unique(tumors$study), each=nrow(tumors)),
               Symbol=rep(rownames(tumors), length(unique(tumors$study))), logFC=0, p.value=0, Bonferroni=0)
for(i in 1:nrow(tumors)){
  print(i)
  for(j in unique(tumors$study)){
    statTest <- t.test(exprs(tumors)[i,tumors$study==j], exprs(normals)[i,], alternative = "greater")
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
df$All_higher<-factor(ifelse(df$Perc_samples==100,"100%","< 100%"))
df$All_higher_num<-ifelse(df$Perc_samples==100,1,0.25)
df<-merge(df, resMark[,c(2:6)], by="Symbol")
df$`Non-cleavable.linker.specific`[df$`Non-cleavable.linker.specific`=="UK"]<-"Undetermined"
df$`Non-cleavable.linker.specific`<-factor(df$`Non-cleavable.linker.specific`, levels=c("No","Yes","Undetermined"))
df$Mechanism<-ifelse(df$Mechanism=="R","Resistance","Sensitivity")

pl1<-ggplot(data = df, aes(x = Symbol, y = CancerN)) +
  geom_tile(fill = "white",color="black") +
  scale_fill_manual(name="Clinical status", breaks=c("Phase 2","Phase 3","Approved"), values=c(carto_pal(7,"Emrld")[c(2,4,6)]),na.value="white") +
  theme_pubr() +
  facet_grid(. ~ Mechanism, scales = "free", space='free') +
  theme(legend.key=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", face = "plain", size = 11), 
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 10, face ="plain", colour ="black"), 
        legend.title = element_text(size = 12, face = "plain"), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right",
        plot.margin=unit(c(1, 0.235, 0, 2.47), "cm"),
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
               labels = c("< 1", seq(1, ceiling(max(df$logFC)), 2)),
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
        axis.text.x = element_text(colour = "black", size = 9, face = "plain", angle = 45, vjust = 1, hjust = 1),
        axis.ticks = element_blank(),
        #axis.text.y = element_blank(), 
        legend.text = element_text(size = 10, face ="plain", colour ="black"), 
        legend.title = element_text(size = 12, face = "plain"), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "bottom",
        plot.margin=unit(c(-0.5, 3.13, 0.2, 0.2), "cm"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  guides(fill=guide_legend(ncol=4))


pdf("~/ADC/results/figure_3.pdf", width=13, height=7.9)
plot_grid(pl1, pl2,ncol=1, rel_heights = c(1,0.39))
dev.off()


#**************************************************************************************
# Figure 4
#**************************************************************************************

resMark<-read.xlsx("~/ADC/data/1-s2.0-S0959804923006810-mmc2.xlsx", sheet = 8, startRow = 2)
resMark<-resMark[!is.na(resMark$Gene.Symbol),]
resMark<-resMark[resMark$Mechanism!="U",]
colnames(resMark)[colnames(resMark)=="Gene.Symbol"]<-"Symbol"

fc<-read.table("~/ADC/data/PrimaryVsMet_ADC_predictors_logFCscaled.txt", header=T, sep="\t", as.is=T, row.names=1)
fdr<-read.table("~/ADC/data/PrimaryVsMet_ADC_predictors_FDR.txt", header=T, sep="\t", as.is=T)
colnames(fc)<-c("Breast","Kidney","Colon","Liver","Pancreas","Ovarian","Esophagus","Skin")
colnames(fdr)<-c("Symbol","Breast","Kidney","Colon","Liver","Pancreas","Ovarian","Esophagus","Skin")
cg<-intersect(rownames(fc), resMark$Symbol)
resMark<-resMark[match(cg,resMark$Symbol),]
fc<-fc[cg,]
fdr<-fdr[match(cg,fdr$Symbol),]
identical(resMark$Symbol, fdr$Symbol)
identical(resMark$Symbol, rownames(fc))

hcrow<-hclust(dist(fc, method="euclidean"), method="ward.D2")
hccol<-hclust(dist(t(fc), method="euclidean"), method="ward.D2")
fc<-fc[hcrow$order, hccol$order]
fc<-data.frame(Symbol=rownames(fc),fc)

df.fc<-reshape2::melt(fc, id="Symbol")
df.fc$Index<-1:nrow(df.fc)
colnames(df.fc)[3]<-"logFC"
df.fc$ID<-paste(df.fc$Symbol, df.fc$variable,sep="_")
df.fdr<-reshape2::melt(fdr, id="Symbol")
colnames(df.fdr)[3]<-"FDR"
df.fdr$ID<-paste(df.fdr$Symbol, df.fdr$variable,sep="_")
df<-merge(df.fc, df.fdr[,c(3,4)], by="ID")
df$Significant<-ifelse(df$FDR<0.05,"FDR < 0.05","n.s.")
df$Significant_num<-ifelse(df$FDR<0.05,1.4,0.3)
df<-df[order(df$Index),]
df$Symbol<-factor(df$Symbol, levels=unique(df$Symbol))
df$variable<-factor(df$variable, levels=unique(df$variable))
df<-merge(df, resMark[,2:6], by="Symbol")
df$`Non-cleavable.linker.specific`[df$`Non-cleavable.linker.specific`=="UK"]<-"Undetermined"
df$`Non-cleavable.linker.specific`<-factor(df$`Non-cleavable.linker.specific`, levels=c("No","Yes","Undetermined"))
df$Mechanism<-ifelse(df$Mechanism=="R","Resistance","Sensitivity")

pl1<-ggplot(data = df, aes(x = Symbol, y = variable)) +
  geom_tile(color="black", fill="white") +
  theme_pubr() +
  facet_grid(. ~ Mechanism, scales = "free", space='free') +
  theme(legend.key=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(colour = "black", face = "plain", size = 13), 
        legend.text = element_text(size = 11, face ="plain", colour ="black", angle=0), 
        legend.title = element_text(size = 12, face = "plain"), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "top",
        plot.margin=unit(c(1, 0.25, 0, 3.16), "cm"),
        strip.text.x = element_text(size=14, face="bold")) +
  xlab("") +
  ylab("") +
  new_scale_fill() +
  geom_point(aes(fill = logFC,  size = abs(logFC), stroke=Significant_num), shape = 21) +
  scale_size_continuous(name="Standardised\n|log(FC)|", range = c(0,6), breaks=seq(1,6,1), limits = c(0,6)) +
  scale_color_manual(name="", values=c("black","red")) +
  binned_scale(aesthetics = "fill",
               scale_name = "stepsn",
               name="Standardised\nlog(FC)",
               palette = function(x)c(brewer.pal(9,"Blues")[8:3],brewer.pal(9,"RdPu")[2:8]),
               breaks = seq(-6, 6, 1),
               labels = c("-6","","-4","","-2","","0","","2","","4","","6"),
               limits = c(-6, 6),
               #show.limits = TRUE,
               guide = guide_colorsteps(frame.colour = "black", label.theme = element_text(angle=0, size=11, hjust=0.5, vjust=0.5)))


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

pdf("~/ADC/results/figure_4.pdf", width=13, height=6.2)
plot_grid(pl1, pl2,ncol=1, rel_heights = c(1,0.66))
dev.off()


#********************************************************
# Figure 5
#********************************************************

#### A - Bubble plot of ADC targets in TCGA-BRCA vs GTEx

diseaseCol<-structure(c(carto_pal(12,"Safe")[-12],carto_pal(12,"Vivid")[-12],carto_pal(12,"Bold")[-c(10:12)],"grey80"), names=c("ACC","BLCA","LGG","BRCA","CESC", "CHOL", "COAD", "ESCA", "GBM",
                                                                                                                                "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "MESO",
                                                                                                                                "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
                                                                                                                                "TGCT", "THYM", "THCA", "UCS", "UCEC", "UVM", "Normals"))

load("/Users/dugo.matteo/Desktop/hsr_bioinfo/databases/TCGA/pancancer/data/rnaseq/eSet_TCGA_recount3_tpm_GDCmatched_noDup.RData")
load("/Users/dugo.matteo/Desktop/hsr_bioinfo/databases/GTEx/data/eSet_GTEx_recount3_tpm.RData")

adc<-read.xlsx("~/ADC/data/ADC_target_list.xlsx", sheet = 1)
adc<-adc[adc$Tested.in.hematological.malignancies.only==0,]
adc$phase2_count<-sapply(adc$condition.phase2, function(x)length(unlist(strsplit(x, ","))))
adc$phase2_count[is.na(adc$condition.phase2)]<-0
adc$phase3_count<-sapply(adc$condition.phase3, function(x)length(unlist(strsplit(x, ","))))
adc$phase3_count[is.na(adc$condition.phase3)]<-0
adc$approved_count<-sapply(adc$approved_indications, function(x)length(unlist(strsplit(x, ","))))
adc$approved_count[is.na(adc$approved_indications)]<-0
tumors <- pancan_tpm[rownames(pancan_tpm) %in% adc$Gene_Symbol , pancan_tpm$tcga.gdc_cases.samples.sample_type %in%c("Primary Tumor") & pancan_tpm$study%in%"BRCA"]
tumors<-tumors[,!tumors$tcga.cgc_sample_sample_type%in%""]
tumors<-tumors[,!is.na(tumors$tcga.cgc_sample_sample_type)]
tumors<-tumors[order(rownames(tumors)),]
subty<-read.xlsx("~/ADC/data/1-s2.0-S0959804923006810-mmc2.xlsx", sheet = 9, startRow = 2)
colnames(subty)<-c("tcga.xml_bcr_patient_barcode","Subtype")
pdata<-pData(tumors)
pdata<-merge(pdata,subty,by="tcga.xml_bcr_patient_barcode", all.x=T)
rownames(pdata)<-pdata$external_id
pdata<-pdata[colnames(tumors),]
identical(rownames(pdata), colnames(tumors))
pData(tumors)<-pdata
tumors<-tumors[,tumors$Subtype!="Undetermined"]
dim(tumors)
exprs(tumors) <- log2(exprs(tumors) + 1)
nTumorSamples<-as.data.frame(table(tumors$Subtype))
colnames(nTumorSamples)<-c("Subtype","N")

normals <- gtex_tpm[rownames(gtex_tpm) %in% adc$Gene_Symbol ,]
normals<-normals[order(rownames(normals)),]
dim(normals)
exprs(normals) <- log2(exprs(normals) + 1)
identical(rownames(tumors), rownames(normals))
adc<-adc[match(rownames(tumors), adc$Gene_Symbol),]
identical(rownames(tumors), adc$Gene_Symbol)

th<-apply(exprs(normals), 1, quantile, probs = 0.80)
detection<-apply(exprs(tumors), 2, function(x)x>th)
detection<-apply(detection, 1, function(x)tapply(x, tumors$Subtype, sum)/table(tumors$Subtype)*100)
df.det<-data.frame(Subtype = rownames(detection), detection)
df.det = reshape::melt(df.det, id = c("Subtype"))
colnames(df.det)[3]<-"Perc_samples"
df.det$ID<-paste(df.det$Subtype, df.det$variable, sep="_")

fc<-data.frame(Subtype=rep(unique(tumors$Subtype), each=nrow(tumors)),
               Symbol=rep(rownames(tumors), length(unique(tumors$Subtype))), logFC=0, p.value=0, Bonferroni=0)
for(i in 1:nrow(tumors)){
  print(i)
  for(j in unique(tumors$Subtype)){
    statTest <- t.test(exprs(tumors)[i,tumors$Subtype==j], exprs(normals)[i,], alternative = "greater")
    fc$logFC[fc$Symbol==rownames(tumors)[i] & fc$Subtype==j] <- -diff(statTest$estimate)
    fc$p.value[fc$Symbol==rownames(tumors)[i] & fc$Subtype==j] <- statTest$p.value
  }
}
fc$Bonferroni<-p.adjust(fc$p.value, method="bonferroni")
fc$logBonferroni <- -log10(fc$Bonferroni)
fc$ID<-paste(fc$Subtype, fc$Symbol, sep="_")

fcMat<-reshape(data = fc[,1:3], idvar="Subtype", timevar="Symbol", direction="wide")
rownames(fcMat)<-fcMat$Subtype
colnames(fcMat)<-gsub("logFC\\.", "", colnames(fcMat))
fcMat<-as.matrix(fcMat[,-1])
pMat<-reshape(data = fc[,c(1,2,5)], idvar="Subtype", timevar="Symbol", direction="wide")
rownames(pMat)<-pMat$Subtype
colnames(pMat)<-gsub("Bonferroni\\.", "", colnames(pMat))
pMat<-as.matrix(pMat[,-1])
fcMat[fcMat<log2(2) | pMat >= 0.05]<-0

hcrow<-hclust(dist(fcMat, method="euclidean"), method="ward.D2")
hccol<-hclust(dist(t(fcMat), method="euclidean"), method="ward.D2")
fcMat<-fcMat[hcrow$order, hccol$order]
df<-data.frame(Subtype = rownames(fcMat), fcMat)
df = reshape::melt(df, id = c("Subtype"))
colnames(df)[2]<-"Symbol"
colnames(df)[3]<-"logFC"
df$ID<-paste(df$Subtype, df$Symbol, sep="_")
df$Index<-1:nrow(df)
df<-merge(df, fc[,-c(1:3)],by="ID",all=T)
df<-merge(df, df.det[,-c(1:2)],by="ID",all=T)
df<-merge(df, nTumorSamples, by="Subtype",all=T)
df<-df[order(df$Index),]
df$SubtypeN<-paste0(df$Subtype, " (N = ", df$N, ")")
df$SubtypeN<-factor(df$SubtypeN, levels=unique(df$SubtypeN))
df$Symbol<-factor(df$Symbol, levels=unique(df$Symbol))
df$Bonferroni_cat<-1
df$Bonferroni_cat[df$Bonferroni<0.05]<-2
df$Bonferroni_cat[df$logFC<0]<-1
df$Bonferroni_cat[df$Bonferroni<0.01]<-3
df$Bonferroni_cat[df$Bonferroni<0.001]<-4
df$All_higher<-factor(ifelse(df$Perc_samples==100,"100%","< 100%"))
df$All_higher_num<-ifelse(df$Perc_samples==100,1,0.25)

pl1<-ggplot(data = df, aes(x = Symbol, y = SubtypeN)) +
  geom_tile(fill="white",color="black") +
  theme_pubr() +
  theme(legend.key=element_blank(),
        axis.text.x = element_text(colour = "black", size = 9, face = "plain", angle = 45, vjust = 1, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "plain", size = 11), 
        legend.text = element_text(size = 10, face ="plain", colour ="black"), 
        legend.title = element_text(size = 12, face = "plain"), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "top",
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
               palette = function(x)c("grey90",carto_pal(7,"ag_Sunset")[5:1]),
               breaks = c(0,seq(1, ceiling(max(df$logFC)), 2)),
               limits = c(0, ceiling(max(df$logFC))),
               labels=c("< 1", seq(1, ceiling(max(df$logFC)), 2)),
               guide = guide_colorsteps(frame.colour = "black", label.theme = element_text(size=10))
  ) +
  guides(size = guide_legend(override.aes = list(stroke = c(0.25,0.25,0.25,0.25,1), color=c("black","black","black","black","red"))))

#### B - Density plot of selected targets in TCGA-BRCA vs GTEx

# density ERBB2

df2<-data.frame(Expression=c(exprs(tumors)["ERBB2",],exprs(normals)["ERBB2",]),
                Tissue=c(tumors$Subtype,rep("Normals",ncol(normals))))
mediane<-sort(tapply(df2$Expression, df2$Tissue, median), decreasing=F)
df2$Tissue<-factor(df2$Tissue, levels=names(mediane))
pl2<-ggplot(data=df2, aes(x=Expression, y=Tissue, fill=Tissue)) +
  geom_density_ridges(alpha=0.5) +
  scale_fill_manual(values=c("grey80",carto_pal(12, "Antique")[c(8:11)]), guide="none") +
  theme_pubr() +
  xlab(expression(log['2']*"(TPM+1)")) +
  ylab("") +
  geom_vline(xintercept = quantile(df2$Expression[df2$Tissue=="Normals"], probs=0.8), linetype="dashed") +
  labs(title = "ERBB2") +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, size=10))

# density TACSTD2

df2<-data.frame(Expression=c(exprs(tumors)["TACSTD2",],exprs(normals)["TACSTD2",]),
                Tissue=c(tumors$Subtype,rep("Normals",ncol(normals))))
mediane<-sort(tapply(df2$Expression, df2$Tissue, median), decreasing=F)
df2$Tissue<-factor(df2$Tissue, levels=names(mediane))
pl3<-ggplot(data=df2, aes(x=Expression, y=Tissue, fill=Tissue)) +
  geom_density_ridges(alpha=0.5) +
  scale_fill_manual(values=c("grey80",carto_pal(12, "Antique")[c(8:11)]), guide="none") +
  theme_pubr() +
  xlab(expression(log['2']*"(TPM+1)")) +
  ylab("") +
  geom_vline(xintercept = quantile(df2$Expression[df2$Tissue=="Normals"], probs=0.8), linetype="dashed") +
  labs(title = "TACSTD2") +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, size=10))

# density NECTIN4

df2<-data.frame(Expression=c(exprs(tumors)["NECTIN4",],exprs(normals)["NECTIN4",]),
                Tissue=c(tumors$Subtype,rep("Normals",ncol(normals))))
mediane<-sort(tapply(df2$Expression, df2$Tissue, median), decreasing=F)
df2$Tissue<-factor(df2$Tissue, levels=names(mediane))
pl4<-ggplot(data=df2, aes(x=Expression, y=Tissue, fill=Tissue)) +
  geom_density_ridges(alpha=0.5) +
  scale_fill_manual(values=c("grey80",carto_pal(12, "Antique")[c(8:11)]), guide="none") +
  theme_pubr() +
  xlab(expression(log['2']*"(TPM+1)")) +
  ylab("") +
  geom_vline(xintercept = quantile(df2$Expression[df2$Tissue=="Normals"], probs=0.8), linetype="dashed") +
  labs(title = "NECTIN4") +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, size=10))

plDensity<-plot_grid(pl2, pl3, pl4, ncol=3)

#### C - bubble plot of predictors in TCGA-BRCA vs GTEx

resMark<-read.xlsx("~/ADC/data/1-s2.0-S0959804923006810-mmc2.xlsx", sheet = 8, startRow = 2)
resMark<-resMark[!is.na(resMark$Gene.Symbol),]
resMark<-resMark[resMark$Mechanism!="U",]
colnames(resMark)[colnames(resMark)=="Gene.Symbol"]<-"Symbol"
tumors <- pancan_tpm[rownames(pancan_tpm) %in% resMark$Symbol , pancan_tpm$tcga.gdc_cases.samples.sample_type %in%c("Primary Tumor") & pancan_tpm$study%in%"BRCA"]
tumors<-tumors[,!tumors$tcga.cgc_sample_sample_type%in%""]
tumors<-tumors[,!is.na(tumors$tcga.cgc_sample_sample_type)]
tumors<-tumors[order(rownames(tumors)),]
subty<-read.xlsx("~/ADC/data/1-s2.0-S0959804923006810-mmc2.xlsx", sheet = 9, startRow = 2)
colnames(subty)<-c("tcga.xml_bcr_patient_barcode","Subtype")
pdata<-pData(tumors)
pdata<-merge(pdata,subty,by="tcga.xml_bcr_patient_barcode", all.x=T)
rownames(pdata)<-pdata$external_id
pdata<-pdata[colnames(tumors),]
identical(rownames(pdata), colnames(tumors))
pData(tumors)<-pdata
tumors<-tumors[,tumors$Subtype!="Undetermined"]
dim(tumors)
exprs(tumors) <- log2(exprs(tumors) + 1)
nTumorSamples<-as.data.frame(table(tumors$Subtype))
colnames(nTumorSamples)<-c("Subtype","N")

normals <- gtex_tpm[rownames(gtex_tpm) %in% resMark$Symbol ,]
normals<-normals[order(rownames(normals)),]
dim(normals)
identical(rownames(tumors), rownames(normals))
exprs(normals) <- log2(exprs(normals) + 1)

th<-apply(exprs(normals), 1, quantile, probs = 0.8)
detection<-apply(exprs(tumors), 2, function(x)x>th)
detection<-apply(detection, 1, function(x)tapply(x, tumors$Subtype, sum)/table(tumors$Subtype)*100)
df.det<-data.frame(Subtype = rownames(detection), detection)
df.det = reshape::melt(df.det, id = c("Subtype"))
colnames(df.det)[3]<-"Perc_samples"
df.det$ID<-paste(df.det$Subtype, df.det$variable, sep="_")

fc<-data.frame(Subtype=rep(unique(tumors$Subtype), each=nrow(tumors)),
               Symbol=rep(rownames(tumors), length(unique(tumors$Subtype))), logFC=0, p.value=0, Bonferroni=0)
for(i in 1:nrow(tumors)){
  print(i)
  for(j in unique(tumors$Subtype)){
    statTest <- t.test(exprs(tumors)[i,tumors$Subtype==j], exprs(normals)[i,], alternative = "greater")
    fc$logFC[fc$Symbol==rownames(tumors)[i] & fc$Subtype==j] <- -diff(statTest$estimate)
    fc$p.value[fc$Symbol==rownames(tumors)[i] & fc$Subtype==j] <- statTest$p.value
  }
}
fc$Bonferroni<-p.adjust(fc$p.value, method="bonferroni")
fc$logBonferroni <- -log10(fc$Bonferroni)
fc$ID<-paste(fc$Subtype, fc$Symbol, sep="_")

fcMat<-reshape(data = fc[,1:3], idvar="Subtype", timevar="Symbol", direction="wide")
rownames(fcMat)<-fcMat$Subtype
colnames(fcMat)<-gsub("logFC\\.", "", colnames(fcMat))
fcMat<-as.matrix(fcMat[,-1])
pMat<-reshape(data = fc[,c(1,2,5)], idvar="Subtype", timevar="Symbol", direction="wide")
rownames(pMat)<-pMat$Subtype
colnames(pMat)<-gsub("Bonferroni\\.", "", colnames(pMat))
pMat<-as.matrix(pMat[,-1])
fcMat[fcMat<log2(2) | pMat >= 0.05]<-0

hcrow<-hclust(dist(fcMat, method="euclidean"), method="ward.D2")
hccol<-hclust(dist(t(fcMat), method="euclidean"), method="ward.D2")
fcMat<-fcMat[hcrow$order, hccol$order]
df<-data.frame(Subtype = rownames(fcMat), fcMat)
df = reshape::melt(df, id = c("Subtype"))
colnames(df)[2]<-"Symbol"
colnames(df)[3]<-"logFC"
df$ID<-paste(df$Subtype, df$Symbol, sep="_")
df$Index<-1:nrow(df)
df<-merge(df, fc[,-c(1:3)],by="ID",all=T)
df<-merge(df, df.det[,-c(1:2)],by="ID",all=T)
df<-merge(df, nTumorSamples, by="Subtype",all=T)
df<-df[order(df$Index),]
df$SubtypeN<-paste0(df$Subtype, " (N = ", df$N, ")")
df$SubtypeN<-factor(df$SubtypeN, levels=unique(df$SubtypeN))
df$Symbol<-factor(df$Symbol, levels=unique(df$Symbol))
df$Bonferroni_cat<-1
df$Bonferroni_cat[df$Bonferroni<0.05]<-2
df$Bonferroni_cat[df$logFC<0]<-1
df$Bonferroni_cat[df$Bonferroni<0.01]<-3
df$Bonferroni_cat[df$Bonferroni<0.001]<-4
df$All_higher<-factor(ifelse(df$Perc_samples==100,"100%","< 100%"))
df$All_higher_num<-ifelse(df$Perc_samples==100,1,0.25)
df<-merge(df, resMark[,c(2:6)], by="Symbol")
df$`Non-cleavable.linker.specific`[df$`Non-cleavable.linker.specific`=="UK"]<-"Undetermined"
df$`Non-cleavable.linker.specific`<-factor(df$`Non-cleavable.linker.specific`, levels=c("No","Yes","Undetermined"))
df$Mechanism<-ifelse(df$Mechanism=="R","Resistance","Sensitivity")

pl5<-ggplot(data = df, aes(x = Symbol, y = SubtypeN)) +
  geom_tile(fill = "white",color="black") +
  scale_fill_manual(name="Clinical status", breaks=c("Phase 2","Phase 3","Approved"), values=c(carto_pal(7,"Emrld")[c(2,4,6)]),na.value="white") +
  theme_pubr() +
  facet_grid(. ~ Mechanism, scales = "free", space='free') +
  theme(legend.key=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", face = "plain", size = 11), 
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 10, face ="plain", colour ="black"), 
        legend.title = element_text(size = 12, face = "plain"), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "top",
        plot.margin=unit(c(1, 0.25, 0, 1.6), "cm"),
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
               labels = c("< 1", seq(1, ceiling(max(df$logFC)), 2)),
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

pl6<-ggplot(data = tmp, aes(x = Symbol, y = variable)) +
  geom_tile(aes(fill = value),color="black") +
  theme_pubr() +
  xlab("") +
  ylab("") +
  scale_y_discrete(expand=c(0,0)) +
  scale_fill_manual(name="", values=c("grey30","grey80", "white", carto_pal(12, "Prism")[1:8])) +
  facet_grid(. ~ Mechanism, scales = "free", space='free') +
  theme(legend.key=element_blank(),
        axis.text.x = element_text(colour = "black", size = 9, face = "plain", angle = 45, vjust = 1, hjust = 1),
        axis.ticks = element_blank(),
        #axis.text.y = element_blank(), 
        legend.text = element_text(size = 10, face ="plain", colour ="black"), 
        legend.title = element_text(size = 12, face = "plain"), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "bottom",
        plot.margin=unit(c(-0.5, 0.25, 0.2, 0.07), "cm"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  guides(fill=guide_legend(ncol=4))


pPred<-plot_grid(pl5, pl6,ncol=1, rel_heights = c(1,0.87))

pdf("~/ADC/results/figure_5.pdf", width=13, height=10.5)
plot_grid(pl1, plDensity, pPred, labels = c("A","B","C"), label_size = 18, ncol=1, rel_heights = c(0.55, 0.65,1))
dev.off()


#********************************************************
# Figure 6B
#********************************************************

selTar<-data.frame(TCGA_code=c("BRCA","LUAD","BLCA"), Gene_Symbol=c("TACSTD2,ERBB2,ERBB3,CD276,NECTIN4,VTCN1,MUC1,ITGB6", "ERBB2,MET,CEACAM5,CEACAM6,FOLR1,CD276,TACSTD2,MUC1", "NECTIN4,ERBB2,TACSTD2,VTCN1,CD276,MUC1,PTK7,ITGB6"))
resMark<-read.xlsx("~/ADC/data/1-s2.0-S0959804923006810-mmc2.xlsx", sheet = 8, startRow = 2)
resMark<-resMark[!is.na(resMark$Gene.Symbol),]
resMark<-resMark[resMark$Mechanism!="U",]
colnames(resMark)[colnames(resMark)=="Gene.Symbol"]<-"Symbol"

load("/Users/dugo.matteo/Desktop/hsr_bioinfo/databases/TCGA/pancancer/data/rnaseq/eSet_TCGA_recount3_tpm_GDCmatched_noDup.RData")
exprs(pancan_tpm)<-log2(exprs(pancan_tpm)+1)

tumors <- pancan_tpm[rownames(pancan_tpm) %in% unlist(strsplit(selTar$Gene_Symbol[1],",")) , pancan_tpm$tcga.gdc_cases.samples.sample_type %in%c("Primary Tumor")]
tumors<-tumors[,!tumors$tcga.cgc_sample_sample_type%in%""]
tumors<-tumors[,!is.na(tumors$tcga.cgc_sample_sample_type)]
predSens<-pancan_tpm[rownames(pancan_tpm) %in% resMark$Symbol[resMark$Mechanism=="S"] , colnames(tumors)]
predRes<-pancan_tpm[rownames(pancan_tpm) %in% resMark$Symbol[resMark$Mechanism=="R"] , colnames(tumors)]
medianSens<-apply(exprs(predSens), 2, median)
medianRes<-apply(exprs(predRes), 2, median)
brca <- tumors[,tumors$study%in%"BRCA"]
duplic<-names(which(table(brca$tcga.xml_bcr_patient_barcode)>1))
brca<-brca[,!brca$tcga.xml_bcr_patient_barcode%in%duplic]
percMat<-exprs(brca)
for(i in 1:nrow(brca)){
  percentile <- ecdf(exprs(tumors)[i,])
  percMat[i,]<-percentile(exprs(brca)[i,])*100
}
percSens <- ecdf(medianSens)
percMatSens<-percSens(medianSens[colnames(brca)])*100
percRes <- ecdf(medianRes)
percMatRes<-percRes(medianRes[colnames(brca)])*100
colnames(brca)<-brca$tcga.xml_bcr_patient_barcode
colnames(percMat)<-brca$tcga.xml_bcr_patient_barcode
names(percMatSens)<-brca$tcga.xml_bcr_patient_barcode
names(percMatRes)<-brca$tcga.xml_bcr_patient_barcode
dfPred2<-data.frame(Var2=names(percMatSens),Sens=percMatSens,Res=percMatRes)
df2<-reshape2::melt(percMat)
df2$groups<-cut(df2$value, breaks = seq(0,100,25), include.lowest = T)
df2$study<-"BRCA"
set.seed(123)
df2<-df2[df2$Var2%in%sample(colnames(brca),3),]
dfPred2<-dfPred2[dfPred2$Var2%in%df2$Var2,]

tumors <- pancan_tpm[rownames(pancan_tpm) %in% unlist(strsplit(selTar$Gene_Symbol[2],",")) , pancan_tpm$tcga.gdc_cases.samples.sample_type %in%c("Primary Tumor")]
tumors<-tumors[,!tumors$tcga.cgc_sample_sample_type%in%""]
tumors<-tumors[,!is.na(tumors$tcga.cgc_sample_sample_type)]
predSens<-pancan_tpm[rownames(pancan_tpm) %in% resMark$Symbol[resMark$Mechanism=="S"] , colnames(tumors)]
predRes<-pancan_tpm[rownames(pancan_tpm) %in% resMark$Symbol[resMark$Mechanism=="R"] , colnames(tumors)]
medianSens<-apply(exprs(predSens), 2, median)
medianRes<-apply(exprs(predRes), 2, median)
luad <- tumors[,tumors$study%in%"LUAD"]
duplic<-names(which(table(luad$tcga.xml_bcr_patient_barcode)>1))
luad<-luad[,!luad$tcga.xml_bcr_patient_barcode%in%duplic]
percMat<-exprs(luad)
for(i in 1:nrow(luad)){
  percentile <- ecdf(exprs(tumors)[i,])
  percMat[i,]<-percentile(exprs(luad)[i,])*100
}
percSens <- ecdf(medianSens)
percMatSens<-percSens(medianSens[colnames(luad)])*100
percRes <- ecdf(medianRes)
percMatRes<-percRes(medianRes[colnames(luad)])*100
colnames(luad)<-luad$tcga.xml_bcr_patient_barcode
colnames(percMat)<-luad$tcga.xml_bcr_patient_barcode
names(percMatSens)<-luad$tcga.xml_bcr_patient_barcode
names(percMatRes)<-luad$tcga.xml_bcr_patient_barcode
dfPred3<-data.frame(Var2=names(percMatSens),Sens=percMatSens,Res=percMatRes)
df3<-reshape2::melt(percMat)
df3$groups<-cut(df3$value, breaks = seq(0,100,25), include.lowest = T)
df3$study<-"LUAD"
set.seed(123)
df3<-df3[df3$Var2%in%sample(colnames(luad),3),]
dfPred3<-dfPred3[dfPred3$Var2%in%df3$Var2,]

tumors <- pancan_tpm[rownames(pancan_tpm) %in% unlist(strsplit(selTar$Gene_Symbol[3],",")) , pancan_tpm$tcga.gdc_cases.samples.sample_type %in%c("Primary Tumor")]
tumors<-tumors[,!tumors$tcga.cgc_sample_sample_type%in%""]
tumors<-tumors[,!is.na(tumors$tcga.cgc_sample_sample_type)]
predSens<-pancan_tpm[rownames(pancan_tpm) %in% resMark$Symbol[resMark$Mechanism=="S"] , colnames(tumors)]
predRes<-pancan_tpm[rownames(pancan_tpm) %in% resMark$Symbol[resMark$Mechanism=="R"] , colnames(tumors)]
medianSens<-apply(exprs(predSens), 2, median)
medianRes<-apply(exprs(predRes), 2, median)
blca <- tumors[,tumors$study%in%"BLCA"]
duplic<-names(which(table(blca$tcga.xml_bcr_patient_barcode)>1))
blca<-blca[,!blca$tcga.xml_bcr_patient_barcode%in%duplic]
percMat<-exprs(blca)
for(i in 1:nrow(blca)){
  percentile <- ecdf(exprs(tumors)[i,])
  percMat[i,]<-percentile(exprs(blca)[i,])*100
}
percSens <- ecdf(medianSens)
percMatSens<-percSens(medianSens[colnames(blca)])*100
percRes <- ecdf(medianRes)
percMatRes<-percRes(medianRes[colnames(blca)])*100
colnames(blca)<-blca$tcga.xml_bcr_patient_barcode
colnames(percMat)<-blca$tcga.xml_bcr_patient_barcode
names(percMatSens)<-blca$tcga.xml_bcr_patient_barcode
names(percMatRes)<-blca$tcga.xml_bcr_patient_barcode
dfPred4<-data.frame(Var2=names(percMatSens),Sens=percMatSens,Res=percMatRes)
df4<-reshape2::melt(percMat)
df4$groups<-cut(df4$value, breaks = seq(0,100,25), include.lowest = T)
df4$study<-"BLCA"
set.seed(123)
df4<-df4[df4$Var2%in%sample(colnames(blca),3),]
dfPred4<-dfPred4[dfPred4$Var2%in%df4$Var2,]


g1<-ggplot(df2) +
  geom_hline(yintercept = seq(0,100,25), color = "grey60") +
  geom_col(aes(x = Var1, y = value, fill=groups),
           position = "dodge2",
           show.legend = TRUE,
           alpha = 0.7,
           color="black"
  ) +
  geom_vline(xintercept = 1:8, color = "black", linetype="dotted") +
  scale_fill_manual(name="Quartiles", breaks = levels(df2$groups), labels=c("Q1","Q2","Q3","Q4"), values=carto_pal(7,"Sunset")[c(1,3,5,7)], drop=F) +
  facet_wrap(~Var2, ncol=1) +
  coord_curvedpolar(clip = "off") +
  theme_pubr(border=T) +
  xlab("") +
  ylab("") +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size=10, face="bold"),
        plot.title = element_text(size=15, face="bold"),
        strip.text.x = element_text(size=13, face="bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=13)) +
  geom_hline(data = dfPred2, aes(yintercept = Sens, color="Sensitivity markers\nexpression"), linetype="solid", linewidth=1) +
  geom_hline(data = dfPred2, aes(yintercept = Res, color="Resistance markers\nexpression"), linetype="solid", linewidth=1) +
  scale_color_manual(name="",values=c("coral3", "darkgreen")) +
  guides(color = guide_legend(override.aes = list(fill = NA)))

g2<-ggplot(df3) +
  geom_hline(yintercept = seq(0,100,25), color = "grey60") +
  geom_col(aes(x = Var1, y = value, fill=groups),
           position = "dodge2",
           show.legend = TRUE,
           alpha = 0.7,
           color="black"
  ) +
  geom_vline(xintercept = 1:8, color = "black", linetype="dotted") +
  scale_fill_manual(name="Quartiles", breaks = levels(df3$groups), labels=c("Q1","Q2","Q3","Q4"), values=carto_pal(7,"Sunset")[c(1,3,5,7)], drop=F) +
  facet_wrap(~Var2, ncol=1) +
  coord_curvedpolar(clip = "off") +
  theme_pubr(border=T) +
  xlab("") +
  ylab("") +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size=10, face="bold"),
        plot.title = element_text(size=15, face="bold"),
        strip.text.x = element_text(size=13, face="bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=13)) +
  geom_hline(data = dfPred3, aes(yintercept = Sens, color="Sensitivity markers\nexpression"), linetype="solid", linewidth=1) +
  geom_hline(data = dfPred3, aes(yintercept = Res, color="Resistance markers\nexpression"), linetype="solid", linewidth=1) +
  scale_color_manual(name="",values=c("coral3", "darkgreen")) +
  guides(color = guide_legend(override.aes = list(fill = NA)))

g3<-ggplot(df4) +
  geom_hline(yintercept = seq(0,100,25), color = "grey60") +
  geom_col(aes(x = Var1, y = value, fill=groups),
           position = "dodge2",
           show.legend = TRUE,
           alpha = 0.7,
           color="black"
  ) +
  geom_vline(xintercept = 1:8, color = "black", linetype="dotted") +
  scale_fill_manual(name="Quartiles", breaks = levels(df4$groups), labels=c("Q1","Q2","Q3","Q4"), values=carto_pal(7,"Sunset")[c(1,3,5,7)], drop=F) +
  facet_wrap(~Var2, ncol=1) +
  coord_curvedpolar(clip = "off") +
  theme_pubr(border=T) +
  xlab("") +
  ylab("") +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size=10, face="bold"),
        plot.title = element_text(size=15, face="bold"),
        strip.text.x = element_text(size=13, face="bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=13)) +
  geom_hline(data = dfPred4, aes(yintercept = Sens, color="Sensitivity markers\nexpression"), linetype="solid", linewidth=1) +
  geom_hline(data = dfPred4, aes(yintercept = Res, color="Resistance markers\nexpression"), linetype="solid", linewidth=1) +
  scale_color_manual(name="",values=c("coral3", "darkgreen")) +
  guides(color = guide_legend(override.aes = list(fill = NA))) 

pdf("~/ADC/results/figure_6B.pdf", width=11, height=10)
ggarrange(g3, g1, g2, ncol=3, nrow=1, common.legend = T)
dev.off()

