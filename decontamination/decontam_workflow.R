library(ape)
library(decontam)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(tidyverse)

### After exporting data from QIIME2.
# Metadata
metadata=read.csv("metadata.txt",sep = "\t")
metadata$sample.id=as.character(metadata$sample.id)
metadata$DNAconc=as.numeric(metadata$DNAconc)
metadata$SampleOrControl=as.factor(metadata$SampleOrControl)

metadata=metadata %>%
unique()
metadata=as.data.frame(metadata)
row.names(metadata)=metadata[,1]
metadata[,1]=NULL

# ASV-table
## Beware the "#" included in the header
asvtable=read.csv("asv-table.txt", sep="\t", check.names = F)
numsamples=length(asvtable)
otumat=as.matrix(asvtable[,2:numsamples])
rownames(otumat) <- as.data.frame(asvtable)[,1]

# IMPORT DATA TO PHYLOSEQ
OTU = otu_table(otumat, taxa_are_rows = TRUE)
samplesdata <- sample_data(metadata)
phy=phyloseq(OTU, samplesdata)
ps=phy

# INSPECT LIBRARY SIZES
df <- as.data.frame(sample_data(ps))
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=SampleOrControl)) + geom_point()

# IDENTIFY CONTAMINANTS - PREVALENCE
## In the prevalence test there is a special value worth knowing, threshold=0.5, that will identify as contaminants all sequences thare are more prevalent in negative controls than in positive samples.
contamdf.prev <- isContaminant(phy, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

# LIST CONTAMINANTS
contaminants.prev=subset(contamdf.prev, contamdf.prev$contaminant == TRUE)
#write.table(contaminants.prev, "contaminants-prev05-Caecum.txt", col.names = NA, quote = F, row.names = T, sep = "\t")
#write.table(contamdf.prev, "decontam-p05-NMQ.txt", col.names = NA, quote = F, row.names = T, sep = "\t")
#contamdf.prev=read.table("contaminants-prev05-Blood-Plate1.txt",sep = "\t", row.names = 1, dec = ".", header = T)


# PLOT CONTAMINANTS
ps.pa <- transform_sample_counts(phy, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$SampleOrControl == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$SampleOrControl == "Sample", ps.pa)
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
#ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
 # xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
#And/Or
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_jitter() + 
  ggtitle("Contamination - TRUE vs. FALSE") +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
