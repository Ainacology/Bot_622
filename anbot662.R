#install.packages("BiocManager")
#install.packages("phyloseq")
#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("bipartite")
#install.packages("permute")
#install.packages("lattice")
#install.packages("statnet.common")
#install.packages("network")
#install.packages("sna")
#install.packages("dplyr")
#install.packages("plyr")
#install.packages("reshape2")
#install.packages("glmmTMB")
#install.packages("lme4")
#install.packages("TMB")
#install.packages("stats")
#install.packages("bbmle")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("phyloseq")

#Load Libraries
library(BiocManager)
library(phyloseq)
library(tidyverse)
library(vegan)
library(bipartite)
library(permute)
library(lattice)
library(statnet.common)
library(network)
library(sna)
library(dplyr)
library(plyr)              
library(reshape2)
library(glmmTMB)
library(lme4)
library(TMB)
library(stats)
library(Matrix)
#library(quanteda)
library(ggplot2)
library(dplyr)
library(plyr)
library(reshape2)
library(tidyr)
library(glmmTMB)
library(bbmle)
library(ggeffects)
#library(epitools)


###############           Edit on basis of control samples     ############################
setwd("~/Dropbox/ANBOT662/16S")
dir()
abun16S<-read.csv("abundance_table_99.csv", header=TRUE)

#get the correction vector
negs<-subset(abun16S, Type==c("N","P"))
rownames(negs)<-negs$Sample_ID
negs = negs[,-c(1:2)]
correct2<-apply(negs, 2, function(x){mean(x)+(2*sd(x))})
#apply vector
rownames(abun16S)<-abun16S$Sample_ID
abun16S2<-abun16S[,-c(1:2)]
abun16S2dm<-data.matrix(abun16S2)
corrector<-as.vector(correct2)
corrector<-round(corrector)
data16S<-sweep(abun16S2dm, 2, corrector, "-")
#transform all negative values to zero
abundance_table16S <-ifelse(data16S < 0, 0, data16S)

#############            Transpose data          #################################
abun16Stotal <- as.data.frame(t(as.matrix(abundance_table16S)))
abun16S<-abun16Stotal[,c(-1:-13,-26)]

#############            Create Phyloseq object           #################################

tax16S<-read.csv("annotations_99.csv", header=TRUE)
metadata<-read.csv("metadata.csv")
row.names(tax16S) <- tax16S$OTU
tax16S <- tax16S %>% dplyr::select(-OTU)
row.names(metadata) <- metadata$Sample_ID
metadata <- metadata %>% dplyr::select (-Sample_ID)
OTU_mat<-as.matrix(abun16S)
tax_mat<-as.matrix(tax16S)
OTU = otu_table(OTU_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
TAX = subset_taxa(TAX, !Genus %in% "Chloroplast_ge") #filter out chloroplast from tax table
TAX = subset_taxa(TAX, !Genus %in% "Mitochondria_ge") #filter out mitochondria from tax table
TAX = subset_taxa(TAX, !Genus %in% "Wolbachia") #filter out wolbachia from tax table
TAX = subset_taxa(TAX, !Kingdom %in% "unknown") #filter out unknown kingdom from tax table
SAMPLES = sample_data(metadata)
physeq16S <- phyloseq(OTU, TAX, SAMPLES) #create phyloseq object

#############       Filter all samples          #################################
physeq16S_prune <- prune_samples(sample_sums(physeq16S) >= 1, physeq16S)
physeq16S_trim = filter_taxa(physeq16S_prune, function(x) sum(x > 1) > (0.05*length(x)), TRUE) 
ae16Sdf=as.data.frame(t(otu_table(physeq16S_trim)))
#ae16Sdf.na <- na.omit(t(ae16Sdf))
#ae16Sdf.mat<-as.matrix(ae16Sdf.na)
#ae16Sdf.na <- na.omit(ae16Sdf.mat)
#ae16Sdf.OTU = otu_table(ae16Sdf.na, taxa_are_rows = TRUE)
#physeq16S <- phyloseq(OTU, TAX, SAMPLES) #create phyloseq object


#############       Ordinate         #####################
GP.ord <- ordinate(physeq16S_trim, "NMDS", "horn")
p1 = plot_ordination(physeq16S_trim, GP.ord, type="samples", color="Treatment")
print(p1)

#############     Visualizations    ######################

p<-plot_bar(physeq16S_trim, fill="Family", x="Ntreatment")+facet_wrap(~Wtreatment, ncol=1)
p + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")

#############       Build dataset for GLMM          #################################
ae16Sdf$sum<-rowSums(ae16Sdf)
ae16Sdf<-subset(ae16Sdf, sum>5)
ae16Sdf$Sample_ID<-row.names(ae16Sdf)

metadata2<-read.csv("metadata.csv")
ma16Sdata2<-join(ae16Sdf,metadata2, by="Sample_ID", type="left")

ma16Sdata<-subset(ma16Sdata2, Wtreatment!="Natural")

colnames(ma16Sdata)
wide16S<-ma16Sdata[,c(1:65,67,69,72:74)] #pulls out only relevant metadata
colnames(wide16S)
longdata16S1<-melt(wide16S, id.vars=c("sum","Sample_ID","Sex","Treatment","Mesocosm","Wtreatment","Ntreatment"))
longdata16S1$rowid<-1:nrow(longdata16S1) 

Px16S<-read.csv("16SPx.csv", header=TRUE)
Px16S$Pnf<-Px16S$Pn+Px16S$Pf
longdata16S<-join(longdata16S1, Px16S, by="variable", type="left")
longdata16S$Wtreatment<-relevel(as.factor(longdata16S$Wtreatment), ref="Lab")
longdata16S$Ntreatment<-relevel(as.factor(longdata16S$Ntreatment), ref="Sterile")

#Danya's model

dweb2<-glmmTMB(value~offset(log(sum))+Sex+Wtreatment+Ntreatment
               +(1|variable)+(1|variable:Sex)+(1|variable:Wtreatment)
               +(1|variable:Ntreatment),  ziformula=~(1|variable), 
               family=nbinom2, data=longdata16S)

dweb2nn<-glmmTMB(value~offset(log(sum))+Sex+Wtreatment+Ntreatment
                 +(1|variable)+(1|variable:Sex)+(1|variable:Wtreatment)
                 ,  ziformula=~(1|variable), 
                 family=nbinom2, data=longdata16S)

dweb2nw<-glmmTMB(value~offset(log(sum))+Sex+Wtreatment+Ntreatment
                 +(1|variable)+(1|variable:Sex)+(1|variable:Ntreatment)
                 ,  ziformula=~(1|variable), 
                 family=nbinom2, data=longdata16S)


dweb2ns<-glmmTMB(value~offset(log(sum))+Sex+Wtreatment+Ntreatment
               +(1|variable)+(1|variable:Wtreatment)
               +(1|variable:Ntreatment),  ziformula=~(1|variable), 
               family=nbinom2, data=longdata16S)

### Log likelihood ratio of nested models
anova(dweb2,dweb2nn) #p=2.2e-16
anova(dweb2,dweb2nw) #p=2.2e-16
anova(dweb2,dweb2ns) #p=1.585e-08

summary(dweb2)

### Extract and view conditional modes
bestfitw<-ranef(dweb2)
Ev<-bestfitw[["cond"]][["variable"]]
wT<-bestfitw[["cond"]][["variable:Wtreatment"]]
nT<-bestfitw[["cond"]][["variable:Ntreatment"]]

nT$mode<-row.names(nT)
x<-str_split_fixed(nT$mode, ":", 2)
nT<-cbind(nT,x)
colnames(nT)<- c("value","mode","asv","type")
nTdcasted<-dcast(nT, asv ~ type)

Ev$asv<-row.names(Ev)

modes<-join(nTdcasted,Ev, by="asv", type="left")
modes$NSresponse<-modes$Nonsterile-modes$Sterile

tax16Scb2<-tax16S
tax16Scb2$asv<-row.names(tax16Scb2)

fullmodes<-join(modes,tax16Scb2, by="asv", type="left")

