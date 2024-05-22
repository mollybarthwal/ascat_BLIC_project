library(ASCAT)
setwd("~/Downloads/ascat_intro/")

##germline normal chr10 (i.e. E56_B)
gl_n10.logR <- read.table("./E56_chr10_neutral/E56_B.LogR", sep = "\t", header = TRUE)
gl_n10.BAF <- read.table("./E56_chr10_neutral/E56_B.BAF", sep = "\t", header = TRUE)
##cell line normal chr10 (i.e. E56_CL)
cl_n10.logR <- read.table("./E56_chr10_neutral/E56_CL.LogR", sep = "\t", header = TRUE)
cl_n10.BAF <- read.table("./E56_chr10_neutral/E56_CL.BAF", sep = "\t", header = TRUE)
##germline deleted chr10 (i.e. E28_B)
gl_d10.logR <- read.table("./E28_chr10_deletion/E28_B.LogR", sep = "\t", header = TRUE)
gl_d10.BAF <- read.table("./E28_chr10_deletion/E28_B.BAF", sep = "\t", header = TRUE)
##cell line deleted chr10 (i.e. E28_CL)
cl_d10.logR <- read.table("./E28_chr10_deletion/E28_CL.LogR", sep = "\t", header = TRUE)
cl_d10.BAF <- read.table("./E28_chr10_deletion/E28_CL.BAF", sep = "\t", header = TRUE)


##germline normal chr13 (i.e. E28_B)
gl_n10.logR <- read.table("./E56_chr10_neutral/E56_B.LogR", sep = "\t", header = TRUE)
gl_n10.BAF <- read.table("./E56_chr10_neutral/E56_B.BAF", sep = "\t", header = TRUE)
##cell line normal chr13 (i.e. E43_CL)
cl_n10.logR <- read.table("./E56_chr10_neutral/E56_CL.LogR", sep = "\t", header = TRUE)
cl_n10.BAF <- read.table("./E56_chr10_neutral/E56_CL.BAF", sep = "\t", header = TRUE)
##germline deleted chr10 (i.e. E28_B)
gl_d10.logR <- read.table("./E28_chr10_deletion/E28_B.LogR", sep = "\t", header = TRUE)
gl_d10.BAF <- read.table("./E28_chr10_deletion/E28_B.BAF", sep = "\t", header = TRUE)
##cell line deleted chr10 (i.e. E28_CL)
cl_d10.logR <- read.table("./E28_chr10_deletion/E28_CL.LogR", sep = "\t", header = TRUE)
cl_d10.BAF <- read.table("./E28_chr10_deletion/E28_CL.BAF", sep = "\t", header = TRUE)

## Checking the length of chr10  
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



##sliding window attempt


#1. initialise a sliding window at 1Mbps long
sw_start <- 1
sw_end <- 1000000


# ## we propose i to be the number of sliding window.
# ## deciding how many windows do we need (134). and we know that the length od chr10 is 135M.
n <- round(max(gl_n10.logR[gl_n10.logR$Chr == "10", ]$Position)/1000000)

gl_mix.logR <- list()
gl_mix.BAF <- list()
cl_mix.logR <- list()
cl_mix.BAF <- list()

for (i in 1:5) {
  rows <- which(gl_n10.logR$Chr == "chr10" & gl_n10.logR$Position >=  sw_start & gl_n10.logR$Position <= sw_end)
  gl_mix.logR[[i]] <- gl_n10.logR
  gl_mix.BAF[[i]] <- gl_n10.BAF
  cl_mix.logR[[i]] <- cl_n10.logR
  cl_mix.BAF[[i]] <- cl_n10.BAF

  gl_mix.logR[[i]][rows, ] <- gl_d10.logR[rows, ]
  gl_mix.BAF[[i]][rows, ] <- gl_d10.BAF[rows, ]
  cl_mix.logR[[i]][rows, ] <- cl_d10.logR[rows, ]
  cl_mix.BAF[[i]][rows, ] <- cl_d10.BAF[rows, ]

  sw_start <- sw_end + 1
  sw_end <-  sw_end + 1000000
  i <- i + 1
}







##delete the "chr" prefix to run the latest version of ASCAT
library(dplyr)
gl_mix.logR %>% mutate(Chr = sub("chr", "", Chr)) %>% write.table(., "gl_mix.logR", sep = "\t", quote = F, row.names = F)
gl_mix.BAF %>% mutate(Chr = sub("chr", "", Chr)) %>% write.table(., "gl_mix.BAF", sep = "\t", quote = F, row.names = F)
cl_mix.logR  %>% mutate(Chr = sub("chr", "", Chr)) %>% write.table(., "cl_mix.logR", sep = "\t", quote = F, row.names = F)
cl_mix.BAF  %>% mutate(Chr = sub("chr", "", Chr)) %>% write.table(., "cl_mix.BAF", sep = "\t", quote = F,  row.names = F)

library(ASCAT)
ascat.bc = ascat.loadData(Tumor_LogR_file = "cl_mix.logR", Tumor_BAF_file = "cl_mix.BAF", Germline_LogR_file = "gl_mix.logR", Germline_BAF_file = "gl_mix.BAF", gender = rep('XX',100), genomeVersion = "hg38")
ascat.plotRawData(ascat.bc, img.prefix = "Mixed_twofifthE28_inE56")

## Skipping GC correction at the moment as we don't have the correct ref file 
# ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "GC_example.txt", replictimingfile = "RT_example.txt")
# ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")

ascat.bc = ascat.aspcf(ascat.bc, out.dir = NA) # penalty=25 for targeted sequencing data


ascat.output = ascat.runAscat(ascat.bc) # gamma=1 for HTS data
save(ascat.bc,ascat.output,file='ASCAT_objects_mixed_twofifthE28_inE56.Rdata')



# 
ascat.bc = ascat.loadData(Tumor_LogR_file = "E28_chr10_deletion/E28_CL.LogR", Tumor_BAF_file = "E56_chr10_neutral/E56_CL.BAF", Germline_LogR_file = "gl_mix.logR", Germline_BAF_file = "gl_mix.BAF", gender = rep('XX',100))
ascat.plotRawData(ascat.bc, img.prefix = "Non_mixed_56")
# 
# ## Skipping GC correction at the moment as we don't have the correct ref file 
# # ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "GC_example.txt", replictimingfile = "RT_example.txt")
# # ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
# 
# ascat.bc = ascat.aspcf(ascat.bc, out.dir = NA) # penalty=25 for targeted sequencing data
# 
# 
# ascat.output = ascat.runAscat(ascat.bc) # gamma=1 for HTS data


