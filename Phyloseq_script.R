
###################### Mike Henson
###################### Started 7-13-2015 by A. Webber
###################### This is how we import our data and analyze it with the R package "phyloseq"
###################### This is an example code for what we should do
###################### For first time, a lot of packages need to be installed and required
###################### This needs a Mothur shared.file to be converted to a .csv file with seperate Tax and OTU tables.
cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                "#CC79A7", "#F0E442")

############ Extra just google "install phyloseq"
# May need new version of R
# Need to get new Rtools package as well


rm(list=ls())
require(ggplot2)
require(grid)
require(plyr)
require(lubridate)
require(phyloseq)
require(devtools)
require(knitr)
require(reshape2)
require(vegan)
require(data.table)


#Set working directory
setwd("path/to/dirctory")


# Gets in the Nutrient Data (csv file format)
NUT<- read.csv("path/to/nutrient/and/sample/data", colClasses = "character",row.names = 1)



###############################################################################################
# IMPORTANT!!!
# Sample data and OTU table need to have matching sample names in the same order!!!!!!!!!!!!!!!
###############################################################################################

# Read in just the OTU for each sample (csv file with OTUs in rows and samples in columns)
OTU<- read.csv("path/to/OTUtable")


# Create matrix from data frame
USC<- data.matrix(OTU, rownames.force = NA)


# Name the Columns based on OTU numer
rownames(USC)<- paste0("OTU", 1:nrow(USC))


# Read in Taxonomy Table (csv with columns Kingdom, Phylum...)
TAX<- read.csv("path", colClasses = "character")


# Create another matrix from data frame and label first column as OTU number
LSU<- as.matrix(TAX, rownames.force = NA)
rownames(LSU)<- paste0("OTU", 1:nrow(LSU))


#Tell phyloseq what is what (otu and tax tables need to be matrixes)
OTU = otu_table(VOL, taxa_are_rows = TRUE)
TAX = tax_table(LSU)


# Read the data into phyloseq
physeq = phyloseq(OTU, TAX)
#Check that it looks correct
physeq


# Make nutrient data the sample data
sampledata = sample_data(NUT)


# Add row names thats correspond to other phyloseq data
rownames(sampledata) = sample_names(physeq)


# MAKE SURE THE SAMPLE NAMES MATCH EACH OTHER
sampledata


# Nutrient data can be in data.frame
SAM = sample_data(sampledata)

# Read the nutrient data into phyloseq also
ALL = phyloseq(OTU, TAX, SAM)
ALL


#Rarefy dataset to even depth
set.seed(200)
OTUSTER_rarefy<- rarefy_even_depth(ALL, rngseed=200)
OTUSTER_rarefy


# Filter taxa that arent seen at least twice in 11% of the data
OTUSTER_rarefy<-filter_taxa(ALL, function(x) sum(x >= 2) >= (0.11*length(x)), TRUE)
OTUSTER_rarefy
#write out table to have data
write.csv(cbind(tax_table(OTUSTER_rarefy), otu_table(OTUSTER_rarefy)), "OTU_rarefy.csv")

# Run NMDS
iDist<- phyloseq::distance(OTUSTER_rarefy, method= 'bray')
MDS2<- ordinate(OTUSTER_rarefy, 'NMDS', distance = iDist) 
MDS2
p5<- plot_ordination(OTUSTER_rarefy, MDS2)
p5 + geom_point(size=3) + theme_bw()


# envfit of NMDS data with NUT data
NUT_ster_full<-as.data.frame((sample_data(OTUSTER_rarefy)))
# making sure you only have nutrient dataset and not metadata
NUT_ster<- as.data.frame(NUT[,c(6,7,8,9,11,12,13,14)])
# Steps to get data readable for envfit function in vegan
NUT_ster[]<-lapply(NUT_ster, as.numeric)
nmds.envfit<-envfit(MDS2, NUT_ster, permu=999, na.rm = TRUE)
nmds.envfit
#Pulling out data from envfit 
vec.sp.df<-as.data.frame(cbind((nmds.envfit$vectors$arrows),pvals=nmds.envfit$vectors$pvals, rsqr=nmds.envfit$vectors$r))
#Select only those variables that are significant
env.scores.nmds=as.data.frame(vec.sp.df[vec.sp.df$pvals<0.05,])
#check what variables remain
env.scores.nmds
#write out dataset
write.csv(env.scores.nmds, "env.scores.nmds.csv")



#cultivar RA plots
cult_RAplot<-read.csv("path/to/RA/dataset")
ggplot(cult_RAplot, aes(Salinity, RA, color=OTU_type)) + geom_point(aes(shape=Was.it.cultivated)) + facet_wrap(~Group, scales = "free_y") +
  theme_bw() + scale_color_manual(values=cbbPalette) +geom_smooth(se=FALSE)


# Top50 RA Plots
cbbPalette <- c("#9ad0f3", "#e79f00","#009E73", "#D55E00", "brown")
otu_top50<-read.csv("path/to/TOP50RA/dataset", header=T)
ra_melted<- melt(otu_top50, id=c("OTU", "Phylum", "LSUCC"))
write.csv(ra_melted, "otu_Overalltop50_03182020")
#cleaned up dataset outside R to rename varaibles
melted_OTUtop50Overall<-read.csv("path/to/cleanedup/data")
melted_OTUtop50Overall$OTU<-factor(melted_OTUtop50Overall$OTU, levels=unique(melted_OTUtop50Overall$OTU))
ggplot(melted_OTUtop50Overall, aes(OTU, RA)) +geom_boxplot(aes(fill=Phylum), outlier.shape=NA) + geom_point(aes(color=Salinity)) + theme_bw() + theme(axis.text.x = element_text(angle = 90), legend.position="none") + scale_fill_manual(values=cbbPalette) + scale_color_continuous(low="red", high="blue")


# Top50 Fresh RA plots
cbbPalette <- c("#000000", "#9ad0f3", "#e79f00","#009E73", "#808080","#D55E00",
                "#F0E442", "brown")
otu_freshtop50<-read.csv("path/to/TO50FreshRA/data", header=T)
ra_melted<- melt(otu_freshtop50, id=c("OTU", "Phylum", "LSUCC"))
write.csv(ra_melted, "FreshOTURank_melted_03182020.csv")
#cleaned up dataset outside R to rename varaibles
melted_OTUtop50Fresh<-read.csv("path/to/TO50FreshRA/cleaneddata")
melted_OTUtop50Fresh$OTU<-factor(melted_OTUtop50Fresh$OTU, levels=unique(melted_OTUtop50Fresh$OTU))
ggplot(melted_OTUtop50Fresh, aes(OTU, RA)) +geom_boxplot(aes(fill=Phylum), outlier.shape=NA) + geom_point(aes(color=Salinity)) + theme_bw() + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + scale_fill_manual(values=cbbPalette) + scale_color_continuous(limits=c(0,36), low="red", high="blue")


#Top 50 Salt RA plots
cbbPalette <- c("#9ad0f3", "#e79f00", "#009E73", "#0072B2", "#D55E00", 
                "#CC79A7", "#F0E442")
otu_Salttop50<-read.csv("path/to/TO50SaltRA/data", header=T)
ra_melted<- melt(otu_Salttop50, id=c("OTU", "Phylum", "LSUCC"))
write.csv(ra_melted, "SaltOTURank_melted_03182020.csv")
melted_OTUtop50Salt<-read.csv("path/to/TO50SaltRA/cleaneddata")
melted_OTUtop50Salt$OTU<-factor(melted_OTUtop50Salt$OTU, levels=unique(melted_OTUtop50Salt$OTU))
ggplot(melted_OTUtop50Salt, aes(OTU, RA)) +geom_boxplot(aes(fill=Phylum), outlier.shape=NA) + geom_point(aes(color=Salinity)) + theme_bw() + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + scale_fill_manual(values=cbbPalette) + scale_color_continuous(limits=c(0,36), low="red", high="blue")


