# < Genetic architecture of phenotypic differences between endangered hybridizing Arabis floodplain species >
# < PhD Thesis >
# < Neda Rahnamae> 

# < Chapter 2 - Genetic architecture of phenotypic differences between endangered hybridizing Arabis floodplain species >
# < F2 Phenotypic Differences >
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



# path_network = map6j_742_2082_path_network.csv
# First column that does not have a name is ID column

# Mixed Model----
## GLM----
# Mixed model for cross direction
# Select only numeric columns
numeric_vars <- path_network[, c(1, 6, 8:14, 16:20, 22:26, 31:34)]

# Step 1: Initialize list to store residuals for each variable
variables <- names(numeric_vars)  # Adjust to include only the relevant numeric variables in `numeric_vars`
residuals_list <- list()

# Step 2: Create an initial data frame with all sample IDs
residuals_df <- data.frame(row = rownames(path_network))

# Step 3: Fit models for each variable, controlling for Cross/TrayBlock
for (var in variables) {
  # Define the model formula
  formula <- as.formula(paste(var, "~ Cross/TrayBlock"))
  
  # Fit a quasi-Poisson model on each variable using the defined formula
  model <- glm(formula, family = quasipoisson, data = path_network, na.action = na.exclude)
  
  # Extract residuals and create a data frame with row names
  temp_df <- data.frame(row = rownames(model$model), residual = model$residuals)
  colnames(temp_df) <- c("row", var)
  
  # Merge the residuals into the main residuals_df by "row"
  residuals_df <- merge(residuals_df, temp_df, by = "row", all.x = TRUE)
}

# Step 4: Reorder rows to match the original order in path_network
residuals_df <- residuals_df[match(rownames(path_network), residuals_df$row), ]

# Set the "row" column as row names and then remove it from the columns
rownames(residuals_df) <- residuals_df$row
residuals_df <- residuals_df[-1]  # Remove the "row" column now that it is set as row names
write.csv(residuals_df, "./output/cross_control/map6j_742_2082_cross_control_res.csv")
# First column that does not have a name is ID column

# Step 5: Initialize matrices to store pairwise correlations and p-values
cor_matrix <- cor(residuals_df, method = "spearman", use = "pairwise.complete.obs")
write.csv(cor_matrix, "./output/cross_control/cor_matrix_res.csv")



## Summary stats----
# Directory to save the model summaries
output_dir <- "./output/cross_control/model_summaries"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Step 1: Fit models and save results to text files
for (var in variables) {
  # Define the model formula
  formula <- as.formula(paste(var, "~ Cross/TrayBlock"))
  
  # Fit a quasi-Poisson model on each variable using the defined formula
  model <- glm(formula, family = quasipoisson, data = path_network, na.action = na.exclude)
  
  # Create a file path for saving the summary
  file_path <- file.path(output_dir, paste0(var, "_model_summary.txt"))
  
  # Write the model summary to a text file
  summary_text <- capture.output(summary(model))
  writeLines(summary_text, file_path)
}


# Phenotype Network----
# Mixed-model residuals
# Plot the network with significant correlations using qgraph
spearman.res_net <- qgraph(
  matrix_data,
  graph = "cor",
  sampleSize = nrow(path_network),
  alpha = 0.05, # significance 
  bonf = T,
  layout = "spring",
  minimum = "sig",
  maximum = 1,
  cut = 0.5,
  theme = "colorblind",
  posCol = "#018571",
  negCol = "#a6611a",
  vsize = 4.5,
  shape = "circle",
  border.width = 1.5,
  graph = "association",
  details = TRUE,
  label.scale = TRUE,
  label.cex = 1,
  label.font = 1.2,        
  label.color = "#000000", 
  color = "#FFFFFF", 
  curveAll = TRUE,
  fade = F
)


# Clustered Heatmap----
heat_map <- 
  ggcorrplot(matrix_data, 
           colors = color_p,
           outline.col = "white",
           hc.order = T,
           lab = T,
           lab_col = "#313131",
           lab_size = 3,
           sig.level = 0.05,
           tl.cex = 9,
           tl.col = "#313131",
           tl.srt = 45,
           show.legend = T)


