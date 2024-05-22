
setwd("/Users/mollybarthwal/downloads/ascat_intro")


library(tidyverse)
library(ggbeeswarm)
library(EnvStats)
library(Gviz)
library(ggplot2)

  
  add_rosner_outlier_info <- function(data) {
    # Add a flag to say the data are not an outlier (to be changed for outliers later on)
    data$is_outlier <- FALSE
    
    for (this_window_size in levels(data$window_size)) {
      print(this_window_size)
      all_nraws <- data %>%
        filter(window_size == this_window_size) %>%
        filter(!is.na(nraw)) %>%
        pull(nraw)
      
      if (length(all_nraws) < 20) {
        next # Skip the test if not enough values
      }
      outliers <- rosnerTest(all_nraws, k = 15) %>%
        .$all.stats %>%
        filter(Outlier) %>%
        select(Value, Outlier)
      print(outliers)
      data <- data %>%
        mutate(is_outlier = case_when(
          window_size == this_window_size & nraw %in% outliers$Value ~ TRUE,
          window_size == this_window_size ~ FALSE,
          TRUE ~ is_outlier
        ))
    }
    
    return(data)
  }
  
  
  ##Specifying amplification location with filtered outliers 
  
  all_windows1GCRT %>%
    filter(!is_outlier) %>%
    filter(window_start > 75622137 & window_end < 179902577)
  
  #chr1 GCRT corrected 
  all_windows1GCRT <- read_tsv('./E43_Chr1/chr1_GCRTcorrected-all_windows_results.txt.gz',
                                col_types = "cfdddddd")
  
  all_windows1GCRT <- add_rosner_outlier_info(all_windows1GCRT)
  
  
  # geom_segment plot of window start and nraw
  all_windows1GCRT %>% 
    filter(window_start > 75622137 & window_end < 179902577) %>%
    # filter(window_size == "500000") %>%
    ggplot(aes(x = window_start, y = nraw, xend = window_end, yend = nraw)) +
    # geom_rect(data = . %>% filter(n_segments == 0),
    #           aes(xmin = window_start, xmax = window_start + 1000000, ymin = -Inf, ymax = Inf),
    #           fill = "lightgrey") +
    # geom_rect(data = . %>% filter(n_segments == 0) %>%
    #             mutate(window_end = ifelse(window_end < window_start + 200000, window_start + 00000, window_start)),
    #           aes(xmin = window_start, xmax = window_end, ymin = -Inf, ymax = Inf),
    #           fill = "lightgrey") +
    geom_rect(data = . %>% filter(n_segments == 0),
              aes(xmin = window_start, xmax = window_end, ymin = -Inf, ymax = Inf,
                  fill = "undetected"),
              col = "lightgrey") +
    geom_segment(data = . %>% filter(is_outlier),
                 col = "#777777", linewidth = 1, alpha = 1,
                 arrow = arrow(angle = 90, length = unit(1, "mm"), ends = "both"),
                 show.legend = F) +
    geom_segment(data = . %>% filter(!is_outlier),
                 aes(col = window_size), linewidth = 1, alpha = 1,
                 arrow = arrow(angle = 90, length = unit(1, "mm"), ends = "both"),
                 show.legend = F) +
    scale_fill_manual(name = NULL, values = c("undetected" = "lightgrey")) +
    labs(x = "Chr 1 (Mbs)", y = "Average raw CN value") +
    facet_grid(rows = vars(window_size)) +
    scale_x_continuous(breaks = as.integer(c(75622137, 1:12 * 20000000, 179902577)),
                       labels = c(1, paste0(1:12 * 20, "Mb"), "179 Mb")) +
    theme_bw() +
    theme(legend.position = "top")
  
  ##Summary for Chr1GCRT
  summary_all_windows1GCRT <- all_windows1GCRT %>%
    filter(!is_outlier) %>%
    filter(window_start > 75622137 & window_end < 179902577) %>%
    group_by(window_size) %>%
    summarise(q05 = quantile(nraw, 0.05, na.rm = T),
              q95 = quantile(nraw, 0.95, na.rm = T)) %>%
    mutate(iqr5.95 = q95 - q05)
  
  
  #improved violin plot 
  all_windows1GCRT %>%
    filter(window_start > 75622137 & window_end < 179902577) %>%
    filter(!is_outlier) %>%
    ggplot(aes(x = factor(window_size), y = nraw)) +
    # geom_boxplot(fill = "white", width = 0.2, outlier.color = NA, notch = TRUE, ) +
    geom_violin(trim = FALSE, scale = "count", size = 0.2) +
    geom_segment(data = summary_all_windows1GCRT %>% filter(window_size != "10000"),
                 aes(x = as.numeric(window_size) - 1.3,
                     xend = as.numeric(window_size) - 1.3, y = q05, yend = q95),
                 arrow = arrow(angle = 90, ends = "both", length = unit(0.01, "native")),
                 size = 1
    ) +
    geom_text(data = summary_all_windows1GCRT %>% filter(window_size != "10000"),
              aes(x = as.numeric(window_size) - 1.4,
                  y = (q05 + q95) / 2,
                  label = round(iqr5.95, digits = 4)),
              angle = 90
    ) +
    # geom_segment(data = all_windows13GCRT_summary %>% filter(window_size != "10000"),
    #              aes(x = as.numeric(window_size) - 1.5, xend = as.numeric(window_size) - 0.5, y = q05, yend = q05),
    #              size = 1.25
    # ) +
    # geom_segment(data = all_windows13GCRT_summary %>% filter(window_size != "10000"),
    #              aes(x = as.numeric(window_size) - 1.5, xend = as.numeric(window_size) - 0.5, y = q95, yend = q95),
    #              size = 1.25
    # ) +
    # geom_tile(data = all_windows13GCRT_summary %>% filter(window_size != "10000"),
    #           aes(x = window_size, y = (q05 + q95) / 2, width = 0.9, height = q95 - q05,
    #               fill = window_size), col = "black", size = 0.5
  # ) +
  geom_beeswarm(aes(fill = window_size), size = 2, alpha = 0.7, color = "black", shape = 21) +
    scale_fill_discrete(limits = c("10000", "50000", "100000", "200000", "300000", "400000", "500000", "1000000", "5000000")) +
    labs(x = "Window size (Mbs)", y = "Raw Copy-Number value") +
    scale_x_discrete(limits = c("50000", "100000", "200000", "300000", "400000", "500000", "1000000", "5000000")) +
    #scale_y_continuous(limits = c(0.75, 1.25)) +
    theme_bw() +
    theme(legend.position = "none")
  
  
  