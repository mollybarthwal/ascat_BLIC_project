## LogR and BAF Files for E34 and E43 

library(ASCAT)
library(tidyverse)
library(dplyr)

## Need to ask how to integrate cosmetic change and confirm which files to use

## E34 ASCAT Files

# Read E34 files 
E34_B.logR <- read.table("E34_baseline/E34_B.LogR", sep = "\t", header = TRUE)
E34_B.BAF <- read.table("E34_baseline/E34_B.BAF", sep = "\t", header = TRUE)
E34_CL.logR <- read.table("E34_baseline/E34_CL.LogR", sep = "\t", header = TRUE)
E34_CL.BAF <- read.table("E34_baseline/E34_CL.BAF", sep = "\t", header = TRUE)

##delete the "chr" prefix to run the latest version of ASCAT
E34_B.logR %>%mutate(Chr = sub("chr", "", Chr)) %>% write.table(., "E34_B.logR", sep = "\t", quote = F, row.names = F)
E34_B.BAF %>% mutate(Chr = sub("chr", "", Chr)) %>% write.table(., "E34_B.BAF", sep = "\t", quote = F, row.names = F)
E34_CL.logR  %>% mutate(Chr = sub("chr", "", Chr)) %>% write.table(., "E34_CL.logR", sep = "\t", quote = F, row.names = F)
E34_CL.BAF  %>% mutate(Chr = sub("chr", "", Chr)) %>% write.table(., "E34_CL.BAF", sep = "\t", quote = F,  row.names = F)

# Create files after
ascat.bc = ascat.loadData(Tumor_LogR_file = "E34_CL.logR", Tumor_BAF_file = "E34_CL.BAF", Germline_LogR_file = "E34_B.logR", Germline_BAF_file = "E34_B.BAF", gender = rep('XX',100))
ascat.plotRawData(ascat.bc, img.prefix = "E34_Baseline_new")

ascat.bc = ascat.loadData(Tumor_LogR_file = "E34_CL.logR", Tumor_BAF_file = "E34_CL.BAF.BAF", Germline_LogR_file = "E34_B.logR", Germline_BAF_file = "E34_B.BAF", gender = rep('XX',100), genomeVersion = "hg38")


###########

## E43 ASCAT Files

# Read E34 files 
E43_B.logR <- read.table("E43_CNA/E43_B.LogR", sep = "\t", header = TRUE)
E43_B.BAF <- read.table("E43_CNA/E43_B.BAF", sep = "\t", header = TRUE)
E43_CL.logR <- read.table("E43_CNA/E43_CL.LogR", sep = "\t", header = TRUE)
E43_CL.BAF <- read.table("E43_CNA/E43_CL.BAF", sep = "\t", header = TRUE)

##delete the "chr" prefix to run the latest version of ASCAT
E43_B.logR %>%mutate(Chr = sub("chr", "", Chr)) %>% write.table(., "E43_B.logR", sep = "\t", quote = F, row.names = F)
E43_B.BAF %>% mutate(Chr = sub("chr", "", Chr)) %>% write.table(., "E43_B.BAF", sep = "\t", quote = F, row.names = F)
E43_CL.logR  %>% mutate(Chr = sub("chr", "", Chr)) %>% write.table(., "E43_CL.logR", sep = "\t", quote = F, row.names = F)
E43_CL.BAF  %>% mutate(Chr = sub("chr", "", Chr)) %>% write.table(., "E43_CL.BAF", sep = "\t", quote = F,  row.names = F)

# Create files
ascat.bc = ascat.loadData(Tumor_LogR_file = "E43_CL.logR", Tumor_BAF_file = "E43_CL.BAF.BAF", Germline_LogR_file = "E43_baseline/E34_B.LogR", Germline_BAF_file = "E43_B.BAF", gender = rep('XX',100))
ascat.plotRawData(ascat.bc, img.prefix = "E43_Baseline_ASCAT")
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc, img.prefix = "E43_ASCAT_Plot")# penalty=25 for targeted sequencing data
ascat.output = ascat.runAscat(ascat.bc, gamma = 1) # gamma=1 for HTS data


### 4. Running ASCAT for BAF and LogR plots
read.table %>% mutate(Chr = sub("chr", "", Chr)) %>%
  ascat.bc = ascat.loadData(Tumor_LogR_file = "E34_baseline/E34_CL.LogR",
                            Tumor_BAF_file = "E34_baseline/E34_CL.BAF", 
                            Germline_LogR_file = "E34_baseline/E34_B.LogR", 
                            Germline_BAF_file = "E34_baseline/E34_B.BAF", 
                            gender = rep('XX',100))
ascat.plotRawData(ascat.bc, img.prefix = "E34_ASCAT_Plot")
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc, img.prefix = "E34_ASCAT_Plot")# penalty=25 for targeted sequencing data
ascat.output = ascat.runAscat(ascat.bc, gamma = 1) # gamma=1 for HTS data

