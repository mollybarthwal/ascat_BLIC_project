library(tidyverse)
library(ggbeeswarm)
library(EnvStats)
library(Gviz)
library(ggplot2)
library(wesanderson)

setwd("/Users/mollybarthwal/downloads/ascat_intro")

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



#chr10 non corrected 
all_segments10 <- read_tsv('./plots_set1/all_segments.txt',
                         col_types = "cfdddddd")
all_windows10 <- read_tsv('./plots_set1/all_windows_results.txt',
                        col_types = "cfdddddd")

#chr13 
all_segments13 <- read_tsv('./plots_set3/chr13-all_segments.txt.gz',
                             col_types = "cfdddddd")
all_windows13 <- read_tsv('./plots_set3/chr13-all_windows_results.txt.gz',
                            col_types = "cfdddddd")

#chr13 GC corrected 
all_segments13GC <- read_tsv('./plots_set3/chr13_GCcorrected-all_segments.txt.gz',
                         col_types = "cfdddddd")
all_windows13GC <- read_tsv('./plots_set3/chr13_GCcorrected-all_windows_results.txt.gz',
                        col_types = "cfdddddd")
#chr13 GCRT corrected 
all_segments13GCRT <- read_tsv('./plots_set3/chr13_GCRTcorrected-all_segments.txt.gz',
                             col_types = "cfdddddd")
all_windows13GCRT <- read_tsv('./plots_set3/chr13_GCRTcorrected-all_windows_results.txt.gz',
                            col_types = "cfdddddd")

all_windows13GCRT <- add_rosner_outlier_info(all_windows13GCRT)

#chr13-vs-E28 GCRT corrected 
all_segments13_E28GCRT <- read_tsv('./plots_set2/chr13-vs-E28_GCRTcorrected-all_segments.txt.gz',
                               col_types = "cfdddddd")
all_windows13_E28GCRT <- read_tsv('./plots_set2/chr13-vs-E28_GCRTcorrected-all_windows_results.txt.gz',
                              col_types = "cfdddddd")


plot(coverage, nraw)

nraw_wsize <- ggplot(windows, aes(x = window_size, y = nraw))
    nraw_wsize +
    +     geom_violin()
    
ggplot(windows,aes(y = nraw,group=window_size))+
      geom_boxplot()

# column plot of nraw by window size
all_segments13GCRT %>%
  ggplot(aes(x = window_size, y = nraw)) +
  geom_col(fill = "blue") +
  labs(x = "Window size", y = "Total number of segments")

#boxplot avg number of segments by window size
all_segments13GCRT %>%
  group_by(window_size, window_start) %>%
  summarise(num_segments = n()) %>%
  filter(window_size != "5000000") %>%
  ggplot(aes(x = window_size, y = num_segments)) +
  geom_boxplot(aes(fill = window_size), show.legend = FALSE) +
  scale_fill_discrete(limits = c("10000", "50000", "100000", "200000", "300000", "400000", "500000", "1000000", "5000000")) +
  labs(x = "Window size", y = "Average number of segments") +
  theme_bw()

#boxplot avg number of segments by window size including 5000000
all_segments13GCRT %>%
  group_by(window_size, window_start) %>%
  summarise(num_segments = n()) %>%
  ggplot(aes(x = window_size, y = num_segments)) +
  geom_boxplot(aes(fill = window_size), show.legend = FALSE) +
  scale_fill_discrete(limits = c("10000", "50000", "100000", "200000", "300000", "400000", "500000", "1000000", "5000000")) +
  labs(x = "Window size", y = "Average number of segments") +
  theme_bw()

#violin plot of nraw by window size (final)
all_segments13GCRT %>% 
  ggplot(aes(x = factor(window_size), y = nraw)) +
  geom_violin(aes(fill = window_size), trim = FALSE, show.legend = F) +
  geom_boxplot(fill = "white", width = 0.2, outlier.color = NA) +
  geom_quasirandom(size = 1, alpha = 0.5) +
  scale_fill_discrete(limits = c("10000", "50000", "100000", "200000", "300000", "400000", "500000", "1000000", "5000000")) +
  labs(x = "Window size", y = "raw CN value") +
  theme_bw()

##SUMMARY
all_windows13GCRT_summary <- all_windows13GCRT %>%
  filter(!is_outlier) %>%
  group_by(window_size) %>%
  summarise(q05 = quantile(nraw, 0.05, na.rm = T),
            q95 = quantile(nraw, 0.95, na.rm = T)) %>%
  mutate(iqr5.95 = q95 - q05)
#improved violin plot 
all_windows13GCRT %>%
  filter(!is_outlier) %>%
  ggplot(aes(x = factor(window_size), y = nraw)) +
  # geom_boxplot(fill = "white", width = 0.2, outlier.color = NA, notch = TRUE, ) +
  geom_violin(trim = FALSE, scale = "count", size = 0.2) +
  geom_segment(data = all_windows13GCRT_summary %>% filter(window_size != "10000"),
               aes(x = as.numeric(window_size) - 1.3,
                   xend = as.numeric(window_size) - 1.3, y = q05, yend = q95),
               arrow = arrow(angle = 90, ends = "both", length = unit(0.01, "native")),
               size = 1
  ) +
  geom_text(data = all_windows13GCRT_summary %>% filter(window_size != "10000"),
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
  scale_y_continuous(limits = c(0.75, 1.25)) +
  theme_bw() +
  theme(legend.position = "none")

## wes anderson color pallet 
palette <- wesanderson::wes_palette(name = "Zissou1", n = 8, type = "continuous")
all_windows13GCRT %>%
  +     filter(!is_outlier) %>%
  +     ggplot(aes(x = factor(window_size), y = nraw)) +
  +     geom_violin(aes(fill = window_size), trim = FALSE, scale = "count") +
  +     geom_boxplot(fill = "white", width = 0.2, outlier.color = NA, notch = TRUE) +
  +     geom_beeswarm(size = 2, alpha = 0.5, color = "black") +
  +     scale_fill_manual(values = palette) +
  +     labs(x = "Window size (Mbs)", y = "Raw CN value") +
  +     scale_x_discrete(limits = c("50000", "100000", "200000", "300000", "400000", "500000", "1000000", "5000000")) +
  +     scale_y_continuous(breaks = seq(0.5, 2, by = 0.5), limits = c(0.75, 1.25)) +
  +     theme_bw() +
  +     theme(legend.position = "none")

## Summary Table
all_windows13GCRT_summary <- all_windows13GCRT %>%
  filter(!is_outlier) %>%
  group_by(window_size) %>%
  summarize(
    Q1 = quantile(nraw, 0.25, na.rm = TRUE),
    Median = quantile(nraw, 0.5, na.rm = TRUE),
    Q3 = quantile(nraw, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    LowerRange = min(nraw, na.rm = TRUE),
    UpperRange = max(nraw, na.rm = TRUE)
  )

# Print the summary
all_windows13GCRT_summary

# violin plot of nraw by window size, limited to 50000, 300000 and 5000000 for visibility
all_windows13GCRT %>%
  filter(!is_outlier) %>%
  ggplot(aes(x = factor(window_size), y = nraw)) +
  geom_violin(aes(fill = window_size), trim = FALSE, show.legend = F) +
  geom_boxplot(width = 0.2, fill = "white", outlier.color = NA) +
  scale_fill_discrete(limits = c("10000", "50000", "100000", "200000", "300000", "400000", "500000", "1000000", "5000000")) +
  labs(x = "Window size", y = "nraw") +
  scale_x_discrete(limits = c("50000", "300000", "5000000")) +
  theme_bw()

# geom_segment plot of window start and nraw
all_windows13GCRT %>% 
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
  labs(x = "Chr 13 GCRT Corrected", y = "Average raw CN value") +
  facet_grid(rows = vars(window_size)) +
  scale_x_continuous(breaks = as.integer(c(1, 1:6 * 20000000, 133797422)),
                     labels = c(1, paste0(1:6 * 20, "Mb"), "133.8Mb")) +
  theme_bw() +
  theme(legend.position = "top")

# geom_segment plot of window start and nraw
all_windows %>% 
  filter(n_segments == 0) %>%
  ggplot(aes(x = window_start, y = window_size, xend = window_end, yend = window_size)) +
  geom_segment(aes(col = window_size), linewidth = 1, alpha = 1,
               arrow = arrow(angle = 90, length = unit(1, "mm"), ends = "both")) +
  # labs(x = "Window start", y = "nraw") +
  # facet_grid(rows = vars(window_size)) +
  theme_bw()

#why this x axis?
all_windows %>% 
  filter(n_segments == 0) %>%
  ggplot() +
  geom_histogram(aes(x = (window_start + window_end) / 2, fill = window_size), bins = 134) +
  facet_grid(rows = vars(window_size)) +
  scale_y_continuous(breaks = c(0,1))

all_windows %>% 
  filter(n_segments == 0) %>%
  ggplot() +
  geom_density(aes(x = (window_start + window_end) / 2, fill = window_size))
  facet_wrap(vars(window_size), ncol = 1, strip.position = "right", scales = "free_y")

all_windows13GCRT %>% 
  group_by(window_size) %>%
  summarise(number_missed_windows = sum((n_segments == 0) == FALSE))
all_windows13GCRT %>%
  mutate(failed = n_segments == 0) %>%
  ggplot() +
  geom_bar(aes(x = window_size, fill = failed), position = "fill") +
  scale_fill_manual(name = NULL,
                    values = c("FALSE" = "darkgrey", "TRUE" = scales::muted("red")),
                    labels = c("detected", "undetected")) +
  labs(y = "Number of windows") +
  scale_y_continuous(labels = scales::percent, expand = expansion()) +
  theme_bw()

### rosnerTest Script 

# Add a flag to say the data are not an outlier (to be changed for outliers later on)
all_windows13GCRT$is_outlier <- FALSE

for (this_window_size in levels(all_segments13GCRT$window_size)) {
  print(this_window_size)
  all_nraws <- all_windows13GCRT%>%
    filter(window_size == this_window_size) %>%
    pull(nraw)
  
  if (length(all_nraws) < 20) {
    next # Skip the test if not enough values
  }
  outliers <- rosnerTest(all_nraws, k = 10) %>%
    .$all.stats %>%
    filter(Outlier) %>%
    select(Value, Outlier)
  print(outliers)
  all_windows13GCRT <- all_windows13GCRT %>%
    mutate(is_outlier = case_when(
      window_size == this_window_size & nraw %in% outliers$Value ~ TRUE,
      window_size == this_window_size ~ FALSE,
      TRUE ~ is_outlier
     ))
}

all_windows13GCRT%>% filter(!is_outlier)

#need to find a way to group by window without creating assignments for each group (i.e. ws_50000)
all_segments13GCRT %>% group_by(window_size) %>% rosnerTest(nraw)

all_segments13GCRT %>%  group_by(window_size) %>% mutate(nraw_numeric = as.numeric(nraw)) %>% rosnerTest(nraw_numeric)

all_segments13GCRT %>%  
  group_by(window_size) %>% 
  rosnerTest(all_segments13GCRT$nraw)

for (this_window_size in levels(all_windows_summary$window_size)) { all_nraws <- all_windows_summary %>% 
  filter(window_size == this_window_size) %>% filter(!is.na(nraw)) %>% 
  pull(nraw) if (length(all_nraws) < 20) { next } rosnerTest(all_nraws, k = 5) }

###plot showing differences in noise amplitude 

set.seed(123)

# Generate a hypothetical dataset with true copy number states
n_points <- 100
true_copy_number <- rep(c(2, 3, 1, 4), each = n_points / 4)
Mbs <- 1:n_points

simulated_data <- data.frame(Mbs, true_copy_number)

# Add different levels of noise
noise_levels <- c(0.1, 0.5, 1, 2)
simulated_data_with_noise <- lapply(noise_levels, function(noise_amplitude) {
  simulated_data %>%
    mutate(
      observed_copy_number = true_copy_number + runif(n(), min = -noise_amplitude, max = noise_amplitude),
      noise_level = paste0("Noise amplitude = ", noise_amplitude)
    )
}) %>% bind_rows()

# Plot the data
ggplot(simulated_data_with_noise, aes(x = Mbs, y = observed_copy_number, color = as.factor(true_copy_number))) +
  geom_point() +
  facet_wrap(~ noise_level) +
  labs(x = "Mbs", y = "Observed Copy Number", color = "True Copy Number") +
  theme_bw() +
  theme(legend.position = "bottom")

# Calculate the mean absolute deviation for each window size
deviation_data <- all_windows13GCRT %>%
  group_by(window_size) %>%
  summarize(mean_absolute_deviation = mean(abs(nraw - 1), na.rm = TRUE))

# Plot the mean absolute deviation for each window size
ggplot(deviation_data, aes(x = window_size, y = mean_absolute_deviation)) +
  geom_point(stat = "identity", fill = "steelblue") +
  labs(x = "Window Size", y = "Mean Absolute Deviation from 1") +
  
  theme_bw()

#chromosome13ideogram
itrack <- IdeogramTrack(genome = "hg38", chromosome = 13)
plotTracks(itrack, from = 1.5e+07, to = 71500000)

#Detected windows table for Chr13 Loss

summary_table <- all_windows13GCRT %>%
  group_by(window_size) %>%
  summarize(detected_windows = sum(!is.na(nraw)))

formatted_table <- summary_table %>%
  gt() %>%
  cols_label(
    window_size = "Window Size (Mbs)",
    detected_windows = "Detected Chr13 Loss by Window Size"
  ) %>%
  tab_header(title = "Detected Chr13 Loss by Window Size") %>%
  data_color(
    columns = c(detected_windows),
    fn = scales::col_numeric(
      palette = RColorBrewer::brewer.pal(9, "Blues"),
      domain = NULL
    )
  )

formatted_table

#relationship between IQR value and window size 
ggplot(all_windows13GCRT_summary, aes(x = WINDOW_SIZE, y = IQR5.95)) +
  geom_point() +
  scale_y_log10() +
  labs(x = "Window Size", y = "IQR Value") +
  theme_minimal()

##creating tribble
all_windows13GCRT_summarylog <- tribble(
  ~window_size, ~q05,   ~q95,   ~iqr5.95, ~log_iqr,
  10000,        0.949,  1.16,   0.209,    -1.57,
  50000,        0.924,  1.09,   0.168,    -1.78,
  100000,       0.932,  1.09,   0.162,    -1.82,
  200000,       0.944,  1.09,   0.147,    -1.92,
  300000,       0.949,  1.08,   0.134,    -2.01,
  400000,       0.941,  1.08,   0.136,    -2.00,
  500000,       0.946,  1.07,   0.128,    -2.05,
  1000000,      0.955,  1.56,   0.601,    -0.510,
  5000000,      0.980,  1.54,   0.557,    -0.586
)

## weird but close exponential line of best fit for IQR vs window size
ggplot(all_windows13GCRT_summarylog, aes(x = window_size, y = iqr5.95)) +
  geom_point() +
  geom_smooth(method = "lm", formula = -y ~ log(x), se = FALSE, linetype = "solid", color = "red") +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Exponential Line of Best Fit for IQR vs Window Size",
       x = "Window Size",
       y = "IQR") +
  theme_minimal()
