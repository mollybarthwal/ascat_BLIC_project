# Explaination of this very messy project
## Acknowledgements 
I am incredibly grateful for my supervisor, Dr. Javier Herrero, for putting up with me through this project. Managing a BSc student with no coding experience is a challenge I would struggle to take. Big thank you to Chuling Ding, PhD student at the BLIC who showed me endless photos of her cat, explained my own scripts to me, taught me some brilliant shortcuts, helped me deal with UCL's HPC and most importantly, showed me how to change my R studio theme. Thank you to Dan Champion for teaching me how to make this README file and put this project all on GitHub. 

## Project Summary 
I took a healthy patient's chromosome as my base line. I then took a cancer patient's chromosome, cut into into different lengths to test, how sensitive the cancer signal was to getting "lost in the noise". By generating a value for this, we can provide a numeric evaluation of the impact of noise, rather than just saying, "noise can impact our findings". I've taken excerpts from my final dissertation to explain what is going on in all these files. Bare with me as this is my first time with GitHub, have taken note to organise better for future projects :)

## Research Question 
What is the error range on the raw copy-number estimates generated by ASCAT?

## Abstract
Copy number aberrations (CNAs) play a crucial role in understanding cancer. The Allele- Specific Copy number Analysis of Tumours (ASCAT) tool is widely used to infer allele-specific copy number profiles from whole-genome sequencing data. However, the impact of noisy data on the accuracy of ASCAT results has not been quantified in previous literature.
In our study, we aimed to measure the error range in the raw copy-number (CN) estimates generated by ASCAT and evaluate the effect of GC content and replication timing bias corrections on the data.
Our findings revealed a noise margin ranging from 7.65% to 16.85%, determined by sliding simulated windows of CN events along assumed normal chromosomes. Notably, the application of GC and replication timing corrections significantly improved the detection of windows, particularly for smaller focal events, with an increase of up to 7.28%. However, we observed that bias corrections did not always result in decreased variability across all data sets, indicating the complexities of these correction methods and the potential for inadvertently highlighting or introducing additional noise.
Our study highlights the importance of minimising biases and quantifying noise for accurate CNA detection. It provides valuable insights for the refinement of ASCAT and similar computational tools. Further research is needed to validate and expand upon our findings by incorporating different sequencing platforms and studying various species. By doing so, we can enhance our understanding of the impact of noise and biases on genomic data analysis, ultimately improving the accuracy of CNA detection in cancer research.

## Process
It is well established in cancer genomics research that the size of CNAs can vary widely, with alterations ranging from a few kilobases affecting a single gene to several megabases encompassing multiple genes or even whole chromosome arms. Small-scale CNAs (typically a few kilobases to a few hundred kilobases) can involve key genes driving oncogenesis and cancer progression, while large-scale CNAs (megabases to tens of megabases) can drastically alter the genetic landscape of cancer cells and contribute to genomic instability (Beroukhim et al., 2010; Zack et al., 2013). Therefore, we chose our focal lengths to capture this wide range of potential CNAs, from small, gene-specific alterations to larger, chromosomal-scale events.

We employed a sliding window technique, a common method in bioinformatics, to analyse focal regions within larger datasets. This strategy allowed for the systematic replacement of sections of the normal dataset with the simulated copy-number aberrations, thus enabling an examination of the implications of these aberrations at different genomic scales.

This process was executed by developing a custom function in R named 'mixed_ascat'. This managed the data loading, data mixing, substituting the normal (E34) data with aberrations (E43), and the segmentation and analysis of the mixed data. During the mixing phase, segments of the normal data (i.e., E43 blood and cell lines) were replaced with corresponding sections from the aberrant data (i.e., E34 blood and cell lines), based on the predefined chromosome and window positions. The resulting mixed data were subsequently loaded into the ASCAT model for analysis. Initially, the raw data were plotted, followed by data segmentation using the aspcf function from the ASCAT package. The segmented data were also visualized for inspection. This iterative approach was applied to each window in the sliding window dataset, enabling a comprehensive and systematic analysis of the impact of various copy-number aberrations.



