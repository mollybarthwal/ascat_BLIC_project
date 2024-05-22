library(ASCAT)
library(tidyverse)

args <- commandArgs(trailingOnly = T)
window_filename <- "sliding_window.chr10_1Mb.csv"
window_num <- 4
folder_name <- "test"
window_filename <- args[1]
window_num <- args[2]
folder_name <- args[3]
normal_sample <- args[4]
cna_sample <- args[5]
corrections <- args[6]

if (corrections == "GC") {
  do_all_corrections <- F
  do_gc_correction <- T
} else if (corrections == "all") {
  do_all_corrections <- T
  do_gc_correction <- F # ignored if do_all_corrections is TRUE
} else {
  do_all_corrections <- F
  do_gc_correction <- F
}

cat("Reading", window_filename, "\n")
cat("Starting run k =", window_num, "\n")
cat("Saving to", folder_name, "\n")

sliding_window <- read_csv(window_filename, col_types = "ccc")

molly_ascat <- function(chromosome = "10", window_start, window_end, folder) {
  
  ### 3. Mixing data -- currently we are trying to add portions of chr10 from E28_CL into E56_CL
  
  ### Loading data 
  data_folder <- file.path("..", "rawdata", "gcgr_samples")
  normal_germline_logR_file <- file.path(data_folder, paste0(normal_sample, "_B_new.LogR"))
  normal_germline_BAF_file <- file.path(data_folder, paste0(normal_sample, "_B_new.BAF"))
  normal_tumour_logR_file <- file.path(data_folder, paste0(normal_sample, "_CL_new.LogR"))
  normal_tumour_BAF_file <- file.path(data_folder, paste0(normal_sample, "_CL_new.BAF"))
  cna_germline_logR_file <- file.path(data_folder, paste0(cna_sample, "_B_new.LogR"))
  cna_germline_BAF_file <- file.path(data_folder, paste0(cna_sample, "_B_new.BAF"))
  cna_tumour_logR_file <- file.path(data_folder, paste0(cna_sample, "_CL_new.LogR"))
  cna_tumour_BAF_file <- file.path(data_folder, paste0(cna_sample, "_CL_new.BAF"))
  
  normal_germline_logR <- read.table(normal_germline_logR_file, sep = "\t", header = TRUE)
  normal_germline_BAF <- read.table(normal_germline_BAF_file, sep = "\t", header = TRUE)
  normal_tumour_logR <- read.table(normal_tumour_logR_file, sep = "\t", header = TRUE)
  normal_tumour_BAF <- read.table(normal_tumour_BAF_file, sep = "\t", header = TRUE)
  cna_germline_logR <- read.table(cna_germline_logR_file, sep = "\t", header = TRUE)
  cna_germline_BAF <- read.table(cna_germline_BAF_file, sep = "\t", header = TRUE)
  cna_tumour_logR <- read.table(cna_tumour_logR_file, sep = "\t", header = TRUE)
  cna_tumour_BAF <- read.table(cna_tumour_BAF_file, sep = "\t", header = TRUE)
  
  ## rows refers to the position numbers to be replaced, regardless of the sample. 
  rows <- which(normal_germline_logR$Chr == chromosome &
                  normal_germline_logR$Position >= as.numeric(window_start) &
                  normal_germline_logR$Position <= as.numeric(window_end))
  
  mix_germline_logR <- normal_germline_logR
  mix_germline_BAF <- normal_germline_BAF
  mix_tumour_logR <- normal_tumour_logR
  mix_tumour_BAF <- normal_tumour_BAF
  
  mix_germline_logR[rows, ] <- cna_germline_logR[rows, ]
  mix_germline_BAF[rows, ] <- cna_germline_BAF[rows, ]
  mix_tumour_logR[rows, ] <- cna_tumour_logR[rows, ]
  mix_tumour_BAF[rows, ] <- cna_tumour_BAF[rows, ]
  
  
  # Define the sub-folder where to store the results
  output_folder = file.path(folder,
                            paste(window_start, window_end, sep = "_"))
  
  # Create the sub-folder where to store the results
  dir.create(output_folder,
             recursive = T, # Create the parent folders if needs by
             showWarnings = F) # Avoid to compain if the folder already exists
  
  mix_germline_logR_filename <- file.path(output_folder, "mix_germline.logR.gz")
  mix_germline_BAF_filename <- file.path(output_folder, "mix_germline.BAF.gz")
  mix_tumour_logR_filename <- file.path(output_folder, "mix_tumour.logR.gz")
  mix_tumour_BAF_filename <- file.path(output_folder, "mix_tumour.BAF.gz")
  
  write_tsv(mix_germline_logR, file = mix_germline_logR_filename)
  write_tsv(mix_germline_BAF, file = mix_germline_BAF_filename)
  write_tsv(mix_tumour_logR, file = mix_tumour_logR_filename)
  write_tsv(mix_tumour_BAF, file = mix_tumour_BAF_filename)
  
  cat("Loading data", fill = T)
  ascat.bc = ascat.loadData(Tumor_LogR_file = mix_tumour_logR_filename, 
                            Tumor_BAF_file = mix_tumour_BAF_filename, 
                            Germline_LogR_file = mix_germline_logR_filename, 
                            Germline_BAF_file = mix_germline_BAF_filename)
  
  ascat.plotRawData(ascat.bc,
                    img.dir = output_folder)
  
  if (do_all_corrections) {
    ascat.bc = ascat.correctLogR(ascat.bc,
                                 GCcontentfile = "../rawdata/LogRcorrection/GCcontent_SNPloci.txt.gz",
                                 replictimingfile = "../rawdata/LogRcorrection/ReplicationTiming_SNPloci.txt.gz")
    
    ascat.plotRawData(ascat.bc,
                      img.dir = output_folder,
                      img.prefix = "post_ALLcorrection-")
    
  } else if (do_gc_correction) {
    ascat.bc = ascat.correctLogR(ascat.bc,
                                 GCcontentfile = "../rawdata/LogRcorrection/GCcontent_SNPloci.txt.gz")
    
    ascat.plotRawData(ascat.bc,
                      img.dir = output_folder,
                      img.prefix = "post_GCcorrection-")
    
  }
  
  cat("Segmenting data", fill = T)
  ascat.bc = ascat.aspcf(ascat.bc, out.dir = NA)
  
  ascat.plotSegmentedData(ascat.bc,
                          img.dir = output_folder)
  
  cat("Running ASCAT (fitting)", fill = T)
  ascat.output = ascat.runAscat(ascat.bc,
                                img.dir = output_folder,
                                gamma = 1) 
  
  write_tsv(ascat.output$segments_raw,
            file = file.path(output_folder, "segments_raw.txt.gz"))
  write(paste0("Ploidy:\t", ascat.output$ploidy),
        file = file.path(output_folder, "stats.txt"))
  write(paste0("Aberrant cell fraction:\t", ascat.output$aberrantcellfraction),
        file = file.path(output_folder, "stats.txt"),
        append = T)
  write(paste0("Goodness of fit:\t", ascat.output$goodnessOfFit),
        file = file.path(output_folder, "stats.txt"),
        append = T)
  
  return(invisible(ascat.output))
}


molly_ascat(chromosome = sliding_window[window_num, ] %>% pull(chr),
            window_start = sliding_window[window_num, ] %>% pull(window_start),
            window_end = sliding_window[window_num, ] %>% pull(window_end),
            folder = folder_name)

