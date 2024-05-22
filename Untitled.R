## Proof of concept

## sliding window of 1/5 (of the length of chr10)
chr10_logR <- cl_mix.logR[cl_mix.logR$Chr == "chr10", ]
chr10_original <- cl_mix.logR[rows, ] %>% filter(Chr == "chr10")
chr10_remaining <- cl_mix.logR[-rows, ] %>% filter(Chr == "chr10")

chr10_BAF <- cl_mix.BAF[cl_mix.logR$Chr == "chr10", ]
chr10_BAF_original <- cl_mix.BAF[rows, ] %>% filter(Chr == "chr10")
chr10_BAF_remaining <- cl_mix.BAF[-rows, ] %>% filter(Chr == "chr10")


## Full chromosome
mean(chr10_logR$E56_CL)
sd(chr10_logR$E56_CL)

## Deletion chromosome
mean(chr10_original$E56_CL)
sd(chr10_original$E56_CL)

## Remaining normal 
mean(chr10_remaining$E56_CL)
sd(chr10_remaining$E56_CL)



ggplot(chr10_BAF, aes(x=Position, y = E56_CL)) + geom_point()






## after running ascat, comparing those of segmented logRs and BAFs. 


