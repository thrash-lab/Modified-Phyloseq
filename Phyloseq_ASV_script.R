
###################### Micbael W. Henson
###################### Started 7-13-2015 by A. Webber
###################### This is how we import our data and analyze it with the R package "phyloseq"
###################### This is an example code for what we should do for ASVs
###################### For first time, a lot of packages need to be installed and required
###################### This needs a MED file to be converted to a .csv file with seperate Tax and ASV tables.


############ Extra just google "install phyloseq"
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
setwd("path/to/ASV/directory") 


# Gets in the Nutrient Data
NUT<- read.csv("path/to/nutrient/data", colClasses = "character")



###############################################################################################

# IMPORTANT!!!
# Sample data and OTU table need to have matching sample names in the same order!!!!!!!!!!!!!!!
###############################################################################################

# Read in just the ASV sample file
ASV<- read.csv("path/to/asV file")


# Create matrix from data frame because phyloseq is prissy
USC<- data.matrix(ASV, rownames.force = NA)


# Name the Columns based on ASV numer
rownames(USC)<- paste0("OTU", 1:nrow(VOL))


# Read in Taxonomy Table
TAX<- read.csv("path/to/taxonomy/file", colClasses = "character")


# Create another matrix from data frame and label first column as OTU number
LSU<- as.matrix(TAX, rownames.force = NA)
rownames(LSU)<- paste0("OTU", 1:nrow(LSU))


# Load phyloseq
library("phyloseq")

#Tell phyloseq what is what otu and tax tables need to be matrix
OTU = otu_table(USC, taxa_are_rows = TRUE)
TAX = tax_table(LSU)



# Read the data into phyloseq
physeq = phyloseq(OTU, TAX)
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

set.seed(200)
ALL_rarefy<- rarefy_even_depth(ALL, rngseed = 200)


# Filter taxa that arent seen at least twice in 20% of the data. Per PHyloseq
ASVSTER_rarefy<-filter_taxa(ALL, function(x) sum(x >= 2) >= (0.11*length(x)), TRUE)
ASVSTER_rarefy
# Read out ASV file
write.csv(cbind(tax_table(OTUSTER_rarefy), otu_table(OTUSTER_rarefy)), "ASV_phyloseq_Rarefy-Trim.csv")


# Run our NMDS now with only our 16S Sterivex filter data
iDist<- phyloseq::distance(ASVSTER_rarefy, method= 'bray')
MDS2<- ordinate(ASVSTER_rarefy, 'NMDS', distance = iDist) 
MDS2
p5<- plot_ordination(OTUSTER_rarefy, MDS2)
p5 + geom_point(size=3) + theme_bw()


# envfit of NMDS data with NUT data
NUT_ster_full<-as.data.frame((sample_data(OTUSTER_rarefy)))
# making sure you only have nutrient dataset and not metadata
NUT_ster<- as.data.frame(NUT[,c(7,8,9,11,12,13,14,15)])
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


# Top50 RA Plots
cbbPalette <- c("#9ad0f3", "#e79f00", "#009E73", "#D55E00", 
                "limegreen")
asv_top50<-read.csv("path/to/Top50ASV_RA/data", header=T)
ra_melted<- melt(asv_top50, id=c("ASV", "Phylum", "LSUCC"))
write.csv(ra_melted, "asv_Overalltop50_03052020.csv")
#cleaned up dataset outside R to rename varaibles
melted_ASVtop50Overall<-read.csv("path/to/Top50fixedASV_RA/data")
melted_ASVtop50Overall$ASV<-factor(melted_ASVtop50Overall$ASV, levels=unique(melted_ASVtop50Overall$ASV))
ggplot(melted_ASVtop50Overall, aes(ASV, RA)) +geom_boxplot(aes(fill=Phylum), outlier.shape=NA) + geom_point(aes(color=Salinity)) + theme_bw() + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + scale_fill_manual(values=cbbPalette) + scale_color_continuous(low="red", high="blue")

# Top50 Fresh RA plots
cbbPalette <- c("#000000","#9ad0f3" , "#e79f00", "#009E73", "#D55E00", 
                 "#F0E442", "brown")
asv_freshtop50<-read.csv("path/to/Top50FreshASV_RA/data", header=T)
ra_melted<- melt(asv_freshtop50, id=c("ASV","Phylum", "LSUCC"))
write.csv(ra_melted, "asv_Freshtop50_03132020.csv")
#cleaned up dataset outside R to rename varaibles
melted_ASVtop50Fresh<-read.csv("path/to/Top50fixedFreshASV_RA/data")
melted_ASVtop50Fresh$ASV<-factor(melted_ASVtop50Fresh$ASV, levels=unique(melted_ASVtop50Fresh$ASV))
ggplot(melted_ASVtop50Fresh, aes(ASV, RA)) +geom_boxplot(aes(fill=Phylum), outlier.shape=NA) + geom_point(aes(color=Salinity)) + theme_bw() + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + scale_fill_manual(values=cbbPalette) + scale_color_continuous(limits=c(0,36), low="red", high="blue")

#Top 50 Salt RA plots

cbbPalette <- c("#9ad0f3", "#e79f00", "#009E73", "#0072B2", "#D55E00")
asv_Salttop50<-read.csv("path/to/Top50SaltASV_RA/data", header=T)
ra_melted<- melt(asv_Salttop50, id=c("ASV","Phylum", "LSUCC"), value.names=c("Salinity", "RA"))
write.csv(ra_melted, "asv_SaltTop50_03132020.csv")
#cleaned up dataset outside R to rename varaibles
melted_ASVtop50Salt<-read.csv("path/to/Top50fixedSaltASV_RA/data")
melted_ASVtop50Salt$ASV<-factor(melted_ASVtop50Salt$ASV, levels=unique(melted_ASVtop50Salt$ASV))
ggplot(melted_ASVtop50Salt, aes(ASV, RA)) +geom_boxplot(aes(fill=Phylum), outlier.shape=NA) + geom_point(aes(color=Salinity)) + theme_bw() + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + scale_fill_manual(values=cbbPalette) + scale_color_continuous(limits=c(0,36), low="red", high="blue")

#cultivar RA plots
cbbPalette <- c("#808080", "#ff7f50", "#F0E442", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#333366")
cult_RAplot<-read.csv("path/to/ASV-Cultivar/RAdata")
ggplot(cult_RAplot, aes(Salinity, RA, color=Type)) + geom_point(aes(shape=Was.it.cultivated)) + facet_wrap(~Group, scales = "free_y") +
  theme_bw() + scale_color_manual(values=cbbPalette) +geom_smooth(se=FALSE)

#cultivar Selected RA plots
cbbPalette <- c("#808080", "#ff7f50", "#F0E442", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#333366")
cult_RAplot<-read.csv("path/to/seleced_ASV-Cultivar/RAdata")
ggplot(cult_RAplot, aes(Salinity, RA, color=Type)) + geom_point(aes(shape=Was.it.cultivated)) + facet_wrap(~Group, scales = "free_y") +
  theme_bw() + scale_color_manual(values=cbbPalette) +geom_smooth(se=FALSE)







