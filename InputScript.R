library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)

#Reading RNAseq TPM file
seq <- read.csv("C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/rnaseq_donor9861/RNAseqTPM.csv")
#Reading Sample Annotation file
names <- read.csv("C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/rnaseq_donor9861/SampleAnnot.csv")

#Taking the column names to add them back later, as they are currently the first row of the TPM file
rowRecovery <- c(colnames(seq))
namesVec <- c(names$ontology_structure_acronym)

#Changing column names to what we need
colnames(seq) <- namesVec

#Final step
seq <- rbind(rowRecovery, seq)

#Put path after coma; R doesn't like using backslashes (\), so replace them with forward slashes(/), and 
#then direct
write.csv(seq, )

