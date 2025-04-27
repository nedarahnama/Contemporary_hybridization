library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(gridExtra)
library(patchwork)
library(cowplot)

seg.dis.df.cor_map4 <- read.csv("./seg.dis.df.cor_map4.csv")
str(seg.dis.df.cor_map4)

### Segregation distortion----
#### Plot----
seg.dis <- 
  ggplot(seg.dis.df.cor_map4, aes(x = marker.pos, y = marker.neglog10P, color = marker.crit.val)) +
  geom_point(size=1, alpha=0.7) +
  labs(
    x = "Position",
    y = expression(-log[10](P))
  ) +
  scale_color_manual(
    name = "", 
    values = c("FALSE" = "#FF5555", "TRUE" = "gray"),
    labels = c("TRUE" = "Not-distorted", "FALSE" = "Distorted")
  ) +
  theme_bw() +
  facet_grid(. ~ marker.chr, scales = "free_x") +
  theme(
  legend.position = "bottom",
  legend.text = element_text(size = 12),  # Increase legend text size
  legend.title = element_text(size = 14), # Increase legend title size
  axis.title.x = element_text(size = 12, margin = margin(t = 10)),
  axis.title.y = element_text(size = 12, margin = margin(r = 10)),
  axis.text.x = element_blank(),          # Remove x-axis text
  axis.ticks.x = element_blank(),         # Remove x-axis ticks
  axis.text.y = element_text(size = 10),
  strip.text.x = element_text(size = 12)  # Increase facet label text size
) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 0.95)))  # Increase legend dot size




### Selection coefficient----
p or N = 0.5
q or S = 0.5 

selection coefficient = s
delta p = s * p (1-p) = s * (pq)
s = delta p / pq
pq = genetic variation 
0.5 * 0.5 = 0.25

# Set the constants
exp_freq <- 0.5
N <- 0.5

# Calculate the observed frequency
obs_freq <- ((2 * (seg.dis.df.cor_map4$marker.NN * 742)) + (seg.dis.df.cor_map4$marker.NS * 742)) / (2 * 742)

# Calculate the selection coefficient for each row
seg.dis.df.cor_map4$selection_coefficient <- (obs_freq - exp_freq) / (N * (1 - N))

# Assuming 'data' is your dataframe containing the selection coefficients
seg.dis.df.cor_map4$selection_sign <- ifelse(seg.dis.df.cor_map4$selection_coefficient > 0, "Positive", "Negative")

#### Significant test - Z score----
##### ALL----
# Mean and standard deviation under the null hypothesis (assuming the null mean is 0)
null_mean <- 0
null_sd <- sd(seg.dis.df.cor_map4$selection_coefficient) # or specify a known standard deviation

# Z-score calculation
seg.dis.df.cor_map4$z_score <- (seg.dis.df.cor_map4$selection_coefficient - null_mean) / null_sd

# P-value calculation (two-tailed test)
seg.dis.df.cor_map4$p_value <- 2 * pnorm(-abs(seg.dis.df.cor_map4$z_score))

# Determine significance (e.g., p < 0.05)
seg.dis.df.cor_map4$significant <- seg.dis.df.cor_map4$p_value < 0.05


##### Per chr----
# Calculate standard deviation per chromosome
sd_per_chromosome <- seg.dis.df.cor_map4 %>%
  group_by(marker.chr) %>%
  summarise(sd_selection_coefficient = sd(selection_coefficient, na.rm = TRUE))

# Merge the standard deviation per chromosome back into the original dataframe
seg.dis.df.cor_map4 <- seg.dis.df.cor_map4 %>%
  left_join(sd_per_chromosome, by = "marker.chr")

# Z-score calculation using chromosome-specific standard deviation
seg.dis.df.cor_map4$z_score <- (seg.dis.df.cor_map4$selection_coefficient - null_mean) / seg.dis.df.cor_map4$sd_selection_coefficient

# P-value calculation (two-tailed test)
seg.dis.df.cor_map4$p_value <- 2 * pnorm(-abs(seg.dis.df.cor_map4$z_score))

# Determine significance (e.g., p < 0.05)
seg.dis.df.cor_map4$significant <- seg.dis.df.cor_map4$p_value < 0.05


#### Plot----
selection.co <- 
  ggplot(seg.dis.df.cor_map4, aes(x = marker.pos, y = selection_coefficient, color = significant)) +
  geom_point(size = 1, alpha = 0.7) +
  labs(
    #title = "Selection Coefficient",
    x = "Position",
    y = "Selection Coefficient (s)"
  ) +
  scale_color_manual(
    name = "", 
    values = c("FALSE" = "gray", "TRUE" = "#FF5555"),
    labels = c("FALSE" = "Not-significant s", "TRUE" = "Significant s"),
    breaks = c("TRUE", "FALSE")  # Explicitly set the order of the legend
  ) +
  theme_bw() +
  facet_grid(. ~ marker.chr, scales = "free_x") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14), # Increase legend title size
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.x = element_blank(),          # Remove x-axis text
    axis.ticks.x = element_blank(),         # Remove x-axis ticks
    axis.text.y = element_text(size = 10),
    strip.text.x = element_text(size = 12)  # Increase facet label text size
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 0.95)))  # Increase legend dot size


# Combined
combined_grid <- (seg.dis) /
  (selection.co)
combined_grid

# Combined label 
combined_grid <- plot_grid(seg.dis, selection.co,
                           labels = "AUTO", ncol = 1, nrow = 2)
combined_grid



### Save data----
save.image("./selection_coeffiecient.RData")