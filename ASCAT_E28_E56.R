library(ASCAT)
library(dplyr)
setwd("~/Downloads/ascat_intro/")

### 1. Running ASCAT with E28 only
ascat.bc = ascat.loadData(Tumor_LogR_file = "E28_chr10_deletion/E28_CL_new.LogR",
                          Tumor_BAF_file = "E28_chr10_deletion//E28_CL_new.BAF", 
                          Germline_LogR_file = "E28_chr10_deletion/E28_B_new.LogR", 
                          Germline_BAF_file = "E28_chr10_deletion/E28_B_new.BAF", 
                          gender = rep('XX',100))
ascat.plotRawData(ascat.bc, img.prefix = "Non_mixed_28")
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc, img.prefix = "Non_mixed_28")# penalty=25 for targeted sequencing data
ascat.output = ascat.runAscat(ascat.bc, gamma = 1) # gamma=1 for HTS data

### 2. Running ASCAT with E56 only
ascat.bc = ascat.loadData(Tumor_LogR_file = "E56_chr10_neutral/E56_CL_new.LogR",
                          Tumor_BAF_file = "E56_chr10_neutral/E56_CL_new.BAF", 
                          Germline_LogR_file = "E56_chr10_neutral/E56_B_new.LogR", 
                          Germline_BAF_file = "E56_chr10_neutral/E56_B_new.BAF", 
                          gender = rep('XX',100))
ascat.plotRawData(ascat.bc, img.prefix = "Non_mixed_56_checking")
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc, img.prefix = "Non_mixed_56")# penalty=25 for targeted sequencing data
ascat.output = ascat.runAscat(ascat.bc, gamma = 1) # gamma=1 for HTS data

### 3. Mixing data -- currently we are trying to add portions of chr10 from E28_CL into E56_CL
### Loading data 
##germline normal chr10 (i.e. E56_B)
gl_n10.logR <- read.table("./E56_chr10_neutral/E56_B_new.LogR", sep = "\t", header = TRUE)
gl_n10.BAF <- read.table("./E56_chr10_neutral/E56_B_new.BAF", sep = "\t", header = TRUE)
##cell line normal chr10 (i.e. E56_CL)
cl_n10.logR <- read.table("./E56_chr10_neutral/E56_CL_new.LogR", sep = "\t", header = TRUE)
cl_n10.BAF <- read.table("./E56_chr10_neutral/E56_CL_new.BAF", sep = "\t", header = TRUE)
##germline deleted chr10 (i.e. E28_B)
gl_d10.logR <- read.table("./E28_chr10_deletion/E28_B_new.LogR", sep = "\t", header = TRUE)
gl_d10.BAF <- read.table("./E28_chr10_deletion/E28_B_new.BAF", sep = "\t", header = TRUE)
##cell line deleted chr10 (i.e. E28_CL)
cl_d10.logR <- read.table("./E28_chr10_deletion/E28_CL_new.LogR", sep = "\t", header = TRUE)
cl_d10.BAF <- read.table("./E28_chr10_deletion/E28_CL_new.BAF", sep = "\t", header = TRUE)

min(gl_n10.logR[gl_n10.logR$Chr == "10", ]$Position)
max(gl_n10.logR[gl_n10.logR$Chr == "10", ]$Position)
## rows refers to the position numbers to be replaced, regardless of the sample. 
rows <- which(gl_n10.logR$Chr == "10" & gl_n10.logR$Position >= 26760000 & gl_n10.logR$Position <= 53520000)
# subset <- gl_n10.logR[gl_n10.logR$Chr == "chr10",]
# which(subset$Chr == "chr10" & subset$Position == 18828)

gl_mix.logR <- gl_n10.logR
gl_mix.BAF <- gl_n10.BAF
cl_mix.logR <- cl_n10.logR
cl_mix.BAF <- cl_n10.BAF

gl_mix.logR[rows, ] <- gl_d10.logR[rows, ]
gl_mix.BAF[rows, ] <- gl_d10.BAF[rows, ]
cl_mix.logR[rows, ] <- cl_d10.logR[rows, ]
cl_mix.BAF[rows, ] <- cl_d10.BAF[rows, ]

gl_mix.logR %>%  write.table(.,"gl_mix.logR", sep = "\t", quote = F, row.names = F)
gl_mix.BAF %>% write.table(., "gl_mix.BAF", sep = "\t", quote = F, row.names = F)
cl_mix.logR %>% write.table(., "cl_mix.logR", sep = "\t", quote = F, row.names = F)
cl_mix.BAF %>% write.table(., "cl_mix.BAF", sep = "\t", quote = F,  row.names = F)

ascat.bc = ascat.loadData(Tumor_LogR_file = "cl_mix.logR", 
                          Tumor_BAF_file = "cl_mix.BAF", 
                          Germline_LogR_file = "gl_mix.logR", 
                          Germline_BAF_file = "gl_mix.BAF", 
                          gender = rep('XX',100))

ascat.plotRawData(ascat.bc, img.prefix = "Mixed_twofifthE28_inE56")
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc, img.prefix = "Mixed_twofifthE28_inE56")# penalty=25 for targeted sequencing data
ascat.output = ascat.runAscat(ascat.bc, gamma = 1) 

sessionInfo()
