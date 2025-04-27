# < Genetic architecture of phenotypic differences between endangered hybridizing Arabis floodplain species >
# < PhD Thesis >
# < Neda Rahnamae> 

# < Chapter 2 - Genetic architecture of phenotypic differences between endangered hybridizing Arabis floodplain species >
# < Parental Phenotypic Differences + Cross Direction >
# March 2024
# 742 individuals - 2082 markers - 22 Phenotypes


# Load R libraries----
library(agricolae) # 1.3-7
library(cowplot) # 1.1.3
library(DescTools) # 0.99.57
library(dplyr) # 1.1.4
library(GGally) # 2.2.1
library(ggcorrplot) # 0.1.4.1
library(ggplot2) # 3.5.1
library(Hmisc) # 5.1-3
library(hrbrthemes) # 0.8.7
library(lme4) # 1.1-35.5
library(magrittr) # 2.0.3
library(patchwork) # 1.3.0
library(ppcor) 
library(psych) # 2.4.6.26
library(qgraph) # 1.9.8
library(tidyr) # 1.3.1
library(tidyverse) # 2.0.0
library(viridis) # 0.6.5



# Data manipulation----
## Parental phenotypes
# parents_pheno_clean = parents_pheno_clean.csv

## F2 phenotypes
# pheno_f2_clean = pheno_f2_clean.csv

## Plot----
# Reshape the data to a long format for easier plotting
# This assumes that "Gen" and all trait columns are in the dataframe "parents_pheno_clean"
long_data <- melt(parents_pheno_clean, id.vars = "Gen", 
                  measure.vars = names(parents_pheno_clean)[1:22], 
                  variable.name = "Trait", value.name = "Value")

#### Boxplots per trait-----
# Plot with customized boxplots for each species
ggplot(long_data, aes(x = Gen, y = Value, color = Gen, fill = Gen)) +
  geom_jitter(width = 0.2, alpha = 0.7) + # Jitter for individual points
  geom_boxplot(alpha = 0.6, width = 0.5, position = position_dodge(width = 0.75)) + # Boxplots with dodge to align
  scale_fill_manual(values = c("A.nem" = "#D98D96", "A.sag" = "#91B3D6"), 
                    labels = c(expression(italic("Arabis nemorensis")), expression(italic("Arabis sagittata")))) + # Italicized legend labels for fill
  scale_color_manual(values = c("A.nem" = "#d1495b", "A.sag" = "#00798c"), 
                     labels = c(expression(italic("Arabis nemorensis")), expression(italic("Arabis sagittata")))) + # Italicized legend labels for color
  facet_wrap(~ Trait, scales = "free") +
  labs(x = "Species", y = "Trait Value") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),  
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 11),
    legend.title = element_blank(), # Remove legend title for a cleaner look
    legend.text = element_text(size = 11)
  )


#### Histogram & Boxplots per trait-----
# Loop over each trait in columns_of_interest and create a plot
for (trait in columns_of_interest) {
  combined_plot_all <- ggplot(pheno_f2_clean, aes_string(x = trait)) +
    geom_histogram(color = "#808080", fill = "#E6E6E6", alpha = 0.9, bins = 30) +
    geom_boxplot(data = filter(parents_pheno_clean, Gen == "A.nem"), 
                 aes_string(x = trait), fill = "#D98D96", color = "#d1495b", alpha = 0.6, 
                 width = 5, position = position_nudge(y = 70)) +
    geom_boxplot(data = filter(parents_pheno_clean, Gen == "A.sag"), 
                 aes_string(x = trait), fill = "#91B3D6", color = "#00798c", alpha = 0.6, 
                 width = 5, position = position_nudge(y = 55)) +
    labs(x = trait, y = "Count") +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 12, margin = margin(t = 10)),  
      axis.title.y = element_text(size = 12, margin = margin(r = 10)),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11)  
    )
  
  # Print or save the plot for each trait
  print(combined_plot_all)
}



# Initialize a list to store plots
plot_list <- list()

# Loop over each trait and generate individual plots
for (trait in columns_of_interest) {
  # Create the plot for each trait
  combined_plot_all <- ggplot(pheno_f2_clean, aes_string(x = trait)) +
    geom_histogram(color = "#808080", fill = "#E6E6E6", alpha = 0.9, bins = 30) +
    geom_boxplot(data = filter(parents_pheno_clean, Gen == "A.nem"), 
                 aes_string(x = trait), fill = "#D98D96", color = "#d1495b", alpha = 0.6, 
                 width = 5, position = position_nudge(y = 70)) +
    geom_boxplot(data = filter(parents_pheno_clean, Gen == "A.sag"), 
                 aes_string(x = trait), fill = "#91B3D6", color = "#00798c", alpha = 0.6, 
                 width = 5, position = position_nudge(y = 55)) +
    labs(x = trait, y = "Count") +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 10.5, margin = margin(t = 5)),  
      #axis.title.y = element_text(size = 12, margin = margin(r = 10)),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9)  
    )
  
  # Add the plot to the list
  plot_list[[trait]] <- combined_plot_all
}

# Combine all plots into a grid
combined_grid_traits <- plot_grid(plotlist = plot_list, ncol = 5, nrow = 5)
combined_grid_traits


## Residuals----
### Extract residuals----
# Select columns representing traits
traits <- names(parents_pheno_clean)[1:22]  # Adjust based on trait column range

# Initialize an empty dataframe to store residuals, with NA to handle missing values
# Adding an 'Individual' column to store row names
residuals_df <- data.frame(Id = rownames(parents_pheno_clean), 
                           matrix(NA, ncol = length(traits), nrow = nrow(parents_pheno_clean)))
colnames(residuals_df)[-1] <- traits  # Set column names for traits

# Create the Gen column based on pattern in Individual names and convert it to factor
residuals_df$Gen <- factor(ifelse(grepl("^nem_", residuals_df$Id), "A.nem", 
                                  ifelse(grepl("^sag_", residuals_df$Id), "A.sag", NA)))

# Loop over each trait to fit the model and extract residuals
for (trait in traits) {
  # Define formula dynamically
  formula <- as.formula(paste(trait, "~ Tray"))
  
  # Fit the quasi-Poisson model, handling missing data with na.omit
  model <- glm(formula, data = parents_pheno_clean, family = quasipoisson, na.action = na.omit)
  
  # Extract residuals
  residual_values <- residuals(model)
  
  # Find non-missing indices for the current trait
  non_missing_indices <- which(!is.na(parents_pheno_clean[[trait]]))
  
  # Store residuals only for non-missing rows in residuals_df
  residuals_df[non_missing_indices, trait] <- residual_values
}



### t-Test----
# Select columns representing traits, excluding "Id"
traits <- setdiff(names(residuals_df), c("Id", "Gen"))

# Initialize a dataframe to store t-test results
t_test_results <- data.frame(Trait = character(), 
                             p_value = numeric(), 
                             t_statistic = numeric(), 
                             stringsAsFactors = FALSE)

# Loop over each trait to perform t-test
for (trait in traits) {
  # Subset residuals based on Gen groups
  group_nem <- residuals_df[[trait]][residuals_df$Gen == "A.nem"]
  group_sag <- residuals_df[[trait]][residuals_df$Gen == "A.sag"]
  
  # Ensure both groups contain non-missing numeric values and more than one unique value
  if (all(is.na(group_nem)) | all(is.na(group_sag)) |
      length(unique(na.omit(group_nem))) <= 1 | length(unique(na.omit(group_sag))) <= 1) {
    # Skip this trait if one of the groups is entirely NA or constant
    next
  }
  
  # Perform t-test if both groups contain numeric data and more than one unique value
  t_test <- t.test(group_nem, group_sag, na.rm = TRUE)
  
  # Store the results
  t_test_results <- rbind(t_test_results, 
                          data.frame(Trait = trait, 
                                     p_value = t_test$p.value, 
                                     t_statistic = t_test$statistic))
}



# View the t-test results
t_test_results

# Add significance codes based on p-values
t_test_results$Significance <- with(t_test_results, ifelse(p_value < 0.001, "***",
                                                           ifelse(p_value < 0.01, "**",
                                                           ifelse(p_value < 0.05, "*",
                                                           ifelse(p_value < 0.1, ".", "")))))


write.csv(t_test_results, "./parental_pheno_res_ttest.csv", row.names = F)