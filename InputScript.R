library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)

#Reading RNAseq TPM file
seq <- read.csv("C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/rnaseq_donor9861/RNAseqTPM.csv")
#Reading Sample Annotation file
names <- read.csv("C:/Users/Shixt/OneDrive/Desktop/INTERNFILES2/rnaseq_donor9861/SampleAnnot.csv")

seq <- rbind(names$ontology_structure_acronym, seq)

#Put path after coma; R doesn't like using backslashes (\), so replace them with forward slashes(/), and 
#then direct
write.csv(seq, )

