# < Genetic architecture of phenotypic differences between endangered hybridizing Arabis floodplain species >
# < PhD Thesis >
# < Neda Rahnamae> 

# < Chapter 2 - Genetic architecture of phenotypic differences between endangered hybridizing Arabis floodplain species >
# < QTL mapping - Stepwise QTL Mapping on Residuals >
# November 2024
# 742 individuals - 2082 markers - 56 QTLs for 19 traits


########################## 
# 1. Load R libraries----
library(ASMap) # 1.0-8
library(dplyr) # 1.1.4
library(ggplot2) # 3.5.1
library(qtl) # 1.70
library(qtlTools) # 1.2.0
library(reshape) # 0.8.9
library(reshape2) # 1.4.4
library(snow) # 0.4-4
library(tidyr) # 1.3.1
library(tidyverse) # 2.0.0

# '-': 0 = #DBDBD7 = Gray
# 'NN': 1 = #FC8D62 = Orange
# 'NS': 2 = #8DA0CB = Blue
# 'SS': 3 = #66C2A5 = Green

########################## 



########################## 
# 2. QTL Mapping----
map6 <- subset(map5_7, chr=1:8)
target_id <- map6$pheno$ID
summary(map6)
# BC(0)F(2) cross
# No. individuals:    742 
# No. phenotypes:     52 
# Percent phenotyped: 96 
# No. chromosomes:    8 
# Autosomes:          1 2 3 4 5 6 7 8 
# Total markers:      2082 
# No. markers:        304 232 244 369 202 160 220 351 
# Percent genotyped:  99.7 
# Genotypes (%):      NN:19.5  NS:48.7  SS:31.8  not SS:0.0  not NN:0.0 

## Add new structure checked pheno data
new_pheno <- read.table("./map6j_742_2082_cross_control_res.csv", sep = ",", header = T, na.strings = "NA", check.names = F, quote = "" )
new_pheno <- new_pheno[match(target_id, new_pheno$ID),]
map6$pheno <- merge(map6$pheno, new_pheno,by = 'ID', all.y = TRUE)
map6$pheno <- subset(map6$pheno, select = -c(2:52))
# Remove ".y" suffix from column names
# names(map6$pheno) <- gsub("\\.y$", "", names(map6$pheno))
map6$pheno <- map6$pheno[match(target_id, map6$pheno$ID),]


######################################################
## 0. Recombination fractions----
# Estimate pairwise recombination fractions
mapthis <- est.rf(map6)
# Pull out recombination fractions/LOD scores from cross
rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
#par(mar=c(4.1,4.1,0.6,0.6), las=1, cex=0.8)
#plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")
# Plot recombination fractions
plotRF(mapthis, what = "rf")

# Jitter markers
map6_jitter <- jittermap(map6)

write.cross(
  map6_jitter,
  format = "csvsr",
  filestem = "./output/map6j_742_2082_cross_control_res",
  chr = c("1", "2", "3", "4", "5", "6", "7", "8"),
  digits = NULL
)


######################################################
## 1. QTL analysis using stepwiseqtl ----
# Let’s proceed to QTL mapping via a Multiple-QTL model.
# Stepwise-QTL----
# This loop goes through all traits and runs qtl mapping pipeline for each of them. 

# List of trait names from map7$pheno, skipping the first two columns (ID and DNA) and filtering for your test
trait_names <- colnames(map7$pheno)[-c(1, 2)]  # Skip the first two columns (ID and DNA)

# For testing, if you want to use only specific traits (e.g., columns 3 and 4), specify them:
trait_names <- trait_names[c(9:22)]  # Adjust this to select your desired traits

# Create an output directory if it doesn't exist
output_dir <- "QTL_Results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Initialize the result list
result_list <- list()

# Loop through the filtered traits
for (i in seq_along(trait_names)) {
  trait_name <- trait_names[i]
  trait_col <- which(colnames(map7$pheno) == trait_name)  # Get the column index for the trait
  
  # 1. QTL genotype probabilities
  #map7 <- calc.genoprob(map6_jitter, map.function = "kosambi", step = 1)
  plotPheno(map7, pheno.col = trait_col, main = paste("Trait", trait_name))
  
  # Save the phenotype plot as SVG
  svg(filename = file.path(output_dir, paste0("phenotype_plot_", trait_name, ".svg")), width = 827 / 72, height = 583 / 72)
  plotPheno(map7, pheno.col = trait_col, main = paste("Trait", trait_name))
  dev.off()
  
  # 2. Perform a 2D-QTL genome scan
  penalties <- scantwo(map7,
                       method = "hk",
                       n.perm = 1000,
                       pheno.col = trait_col,
                       batchsize = 100,
                       n.cluster = 10)
  assign(paste0("penalties_", trait_name), penalties)
  print(summary(penalties))
   
  # Extract and format penalties summary
  penalties_summary_values <- summary(penalties)
  penalties_summary_df <- data.frame(
    FDR = c("5%", "10%"),
    full = penalties_summary_values$full[, 1],
    fv1 = penalties_summary_values$fv1[, 1],
    int = penalties_summary_values$int[, 1],
    add = penalties_summary_values$add[, 1],
    av1 = penalties_summary_values$av1[, 1],
    one = penalties_summary_values$one[, 1]
  )

  # Save the penalties summary as a CSV file
  write.csv(penalties_summary_df, file = file.path(output_dir, paste0("penalties_2D_summary_", trait_name, ".csv")), row.names = FALSE)

  # 3. Calculate penalties
  pen <- calc.penalties(penalties)
  assign(paste0("pen_", trait_name), pen)
  print(pen)

  # Save penalties as a table
  pen_df <- as.data.frame(t(pen))
  write.csv(pen_df, file = file.path(output_dir, paste0("penalties_", trait_name, ".csv")), row.names = TRUE)

  # 4. Perform step-wise QTL analyses
  stepout <- stepwiseqtl(
    map7,
    pheno.col = trait_col,
    max.qtl = 5,
    method = "hk",
    incl.markers = TRUE,
    refine.locations = TRUE,
    additive.only = FALSE,
    penalties = pen,
    keeplodprofile = TRUE,
    keeptrace = TRUE
  )
  
  # Check if stepout is valid and contains QTLs
  if (!is.list(stepout) || is.null(stepout$pos) || length(stepout$pos) == 0) {
    warning(paste("No QTLs found for", trait_name))
    next  # Skip to the next trait if no QTLs are found
  }
  
  assign(paste0("stepout_", trait_name), stepout)
  print(stepout)

  # Save stepwise QTL results as a table
  stepout_summary <- summary(stepout)
  write.csv(stepout_summary, file = file.path(output_dir, paste0("stepwise_qtl_", trait_name, ".csv")), row.names = FALSE)

  # Plot QTL model and LOD profiles
  svg(filename = file.path(output_dir, paste0("qtl_model_plot_", trait_name, ".svg")), width = 827 / 72, height = 583 / 72)
  plot(stepout, alpha = 0.05, chr = c(1, 2, 3, 4, 5, 6, 7, 8), pvalues = TRUE, gap = 50, bandcol = "gray90", main = paste("Trait", trait_name))
  dev.off()

  svg(filename = file.path(output_dir, paste0("lod_profile_plot_", trait_name, ".svg")), width = 827 / 72, height = 583 / 72)
  plotLodProfile(stepout)
  dev.off()

  # 5. Fit QTL positions
  stepout_fit <- fitqtl(
    map7,
    pheno.col = trait_col,
    qtl = stepout,
    formula = as.formula(paste("y ~", paste(paste0("Q", 1:length(stepout$pos)), collapse = " + "))),
    method = "hk",
    dropone = TRUE,
    get.ests = TRUE
  )
  assign(paste0("stepout_fit_", trait_name), stepout_fit)

  # Save the summary.fitqtl object as a text file
  fit_summary <- capture.output(summary(stepout_fit))
  writeLines(fit_summary, file.path(output_dir, paste0("fit_summary_", trait_name, ".txt")))

  # 6. Confidence interval calculations using lodint
  if (length(stepout$pos) > 0) {  # Check if there are QTLs identified
    for (qtl_index in 1:length(stepout$pos)) {
      lod_interval <- lodint(stepout, qtl.index = qtl_index, expandtomarkers = TRUE)

      # Convert the marker name factor to character and include as a new column
      lod_interval_df <- data.frame(
        Marker = rownames(lod_interval),
        chr = as.character(lod_interval$chr),
        pos = lod_interval$pos,
        lod = lod_interval$lod
      )

      # Save the lod interval with marker names as a CSV file
      write.csv(lod_interval_df, file = file.path(output_dir, paste0("lod_interval_", trait_name, "_QTL", qtl_index, ".csv")), row.names = FALSE)
    }
  } else {
    warning(paste("No QTLs found for", trait_name))
  }
  
  # Store results in the result list
  result_list[[trait_name]] <- list(
    penalties = penalties,
    pen = pen,
    stepout = stepout,
    stepout_fit = stepout_fit
  )
}

# Save the entire result list if needed
save(result_list, file = file.path(output_dir, "all_traits_QTL_results.RData"))


## 2. QTL peaks positions ----
lodint(stepout_B.T, qtl.index = 1, expandtomarkers = TRUE)
marker_B.T_Q1 <- find.marker(map7, chr = 3, pos = 163.00) # chr3_3_20.040893 - 
lodint(stepout_B.T, qtl.index = 3, expandtomarkers = TRUE)
marker_B.T_Q3 <- find.marker(map7, chr = 8, pos = 12.23969) # chr8_8_2.061023 - 
# chr8_8_0.656384 - chr8_8_3.118527
# 0.000001 - 17.4757188188002

marker_F.T_Q1 <- find.marker(map7, chr = 3, pos = 190.00) # chr3_3_26.88284 - 
marker_F.T_Q2 <- find.marker(map7, chr = 5, pos = 111.00) # chr5_5_23.418119 - 
marker_F.T_Q4 <- find.marker(map7, chr = 8, pos = 13.00) # chr8_8_2.061022 - 
marker_F.T_Q5 <- find.marker(map7, chr = 8, pos = 133.00) # chr8_8_18.832553 - 


lodint(stepout_L.L, qtl.index = 2, expandtomarkers = TRUE)
marker_L.L_Q2 <- find.marker(map7, chr = 3, pos = 136.2193) # chr3_3_14.255522
# chr3_3_12.540387 - chr3_3_19.199451
# 122.827699700737 - 161.028047128389
lodint(stepout_L.L, qtl.index = 3, expandtomarkers = TRUE)
marker_L.L_Q3 <- find.marker(map7, chr = 3, pos = 201.4139) # chr3_3_29.031007
# chr3_3_26.88284 - chr3_3_29.703667
# 193.967573263987 - 207.661200044656
marker_L.L_Q5 <- find.marker(map7, chr = 8, pos = 35.00) # chr8_8_5.624776

marker_L.W_Q1 <- find.marker(map7, chr = 1, pos = 116.00) # chr1_1_21.722429
marker_L.W_Q2 <- find.marker(map7, chr = 2, pos = 87.00) # chr2_2_19.435694
marker_L.W_Q3 <- find.marker(map7, chr = 2, pos = 129.00) # chr2_2_25.576894
marker_L.W_Q4 <- find.marker(map7, chr = 8, pos = 32.00) # chr8_8_4.356973

marker_Lam_Q1 <- find.marker(map7, chr = 8, pos = 36.00) # chr8_8_5.481058

marker_LLW_Q1 <- find.marker(map7, chr = 3, pos = 61.00) # chr3_3_6.744303
marker_LLW_Q2 <- find.marker(map7, chr = 4, pos = 163.00) # chr4_4_38.897287

marker_N.L_Q2 <- find.marker(map7, chr = 4, pos = 18.00) # chr4_4_7.558305
marker_N.L_Q4 <- find.marker(map7, chr = 8, pos = 7.00) # chr8_8_1.180103

marker_P.H_Q2 <- find.marker(map7, chr = 3, pos = 27.00) # chr3_3_2.664384
marker_P.H_Q3 <- find.marker(map7, chr = 8, pos = 6.00) # chr8_8_1.180103
marker_P.H_Q4 <- find.marker(map7, chr = 8, pos = 101.00) # chr8_8_10.128705

marker_Pet_Q1 <- find.marker(map7, chr = 1, pos = 51.00) # chr1_1_9.074627

marker_RD1_Q1 <- find.marker(map7, chr = 1, pos = 95.00) # chr1_1_17.879753
marker_RD1_Q3 <- find.marker(map7, chr = 8, pos = 28.00) # chr8_8_4.362085

marker_RD2_Q1 <- find.marker(map7, chr = 1, pos = 115.00) # chr1_1_21.954996
marker_RD2_Q3 <- find.marker(map7, chr = 8, pos = 31.00) # chr8_8_4.356973

marker_RD3_Q2 <- find.marker(map7, chr = 2, pos = 118.00) # chr2_2_23.54193
marker_RD3_Q3 <- find.marker(map7, chr = 8, pos = 18.00) # chr8_8_2.718077

marker_RD4_Q2 <- find.marker(map7, chr = 8, pos = 32.00) # chr8_8_4.356973

marker_S.H_Q1 <- find.marker(map7, chr = 1, pos = 162.00) # chr1_1_29.847516
marker_S.H_Q2 <- find.marker(map7, chr = 3, pos = 36.00) # chr3_3_4.035375
marker_S.H_Q3 <- find.marker(map7, chr = 4, pos = 160.00) # chr4_4_38.897287
marker_S.H_Q5 <- find.marker(map7, chr = 8, pos = 102.00) # chr8_8_10.506892

marker_SLD_Q1 <- find.marker(map7, chr = 3, pos = 29.00) # chr3_3_3.491752

marker_SLL_Q2 <- find.marker(map7, chr = 8, pos = 15.00) # chr8_8_3.118527

marker_W.S_Q1 <- find.marker(map7, chr = 3, pos = 63.00) # chr3_3_8.36208
marker_W.S_Q2 <- find.marker(map7, chr = 6, pos = 10.00) # chr6_6_1.052484
marker_W.S_Q3 <- find.marker(map7, chr = 7, pos = 31.00) # chr7_7_3.192272
marker_W.S_Q4 <- find.marker(map7, chr = 8, pos = 17.00) # chr8_8_3.118527


########################## 
# 3. Visualization of QTLs----
# Plot QTL regions on the linkage map
phe_names <- c("Days to Bolting", "Days to Flowering", "Fertility Score", "Inflorescence Height",
               "Lamina Length", "Lamina L/W", "Leaf Length", "Leaf Width",
               "Number of Leaves", "Petal Length", "Petiole Length",
               "Rosette Diameter 1", "Rosette Diameter 2", "Rosette Diameter 3", "Rosette Diameter 4",
               "Side Shoots", "Stem Height", "Stem Leaf Density",
               "Stem Leaf Length", "Stem Leaf Width")

phe_qtlobjects <- list(B.T=stepout_B.T, F.T=stepout_F.T, W.S=stepout_W.S, P.H=stepout_P.H,
                       Lam=stepout_Lam, LLW=stepout_LLW, L.L=stepout_L.L, L.W=stepout_L.W,
                       N.L=stepout_N.L, Pet=stepout_Pet, Pti=stepout_Pti,
                       RD1=stepout_RD1, RD2=stepout_RD2, RD3=stepout_RD3, RD4=stepout_RD4,
                       Ssh=stepout_Ssh, S.H=stepout_S.H, SLD=stepout_SLD,
                       SLL=stepout_SLL, SLW=stepout_SLW)

names(phe_qtlobjects) <- phe_names

phe_colorp <- c("#b00058", "#ff0000", "#ffc000", "#7f7f7f",
                "#0070c0", "#004272", "#002060", "#87afff",
                "#266f8b", "#ff89c4", "#00b0f0",
                "#502273", "#7030a0", "#b17ed8", "#dac2ec",
                "#7fbd98", "#007b32", "#b9e08c",
                "#533b1e", "#cfa879")

map_vis_out <- phe_qtlobjects

cis_map <- lapply(map_vis_out, function(x) {
  calcCis(map7,mod=x, qtlnames=x$altname, returnChr=TRUE, returnMaxLod=TRUE, expandtomarkers=T)
})

cis_map.df <- plyr::ldply(cis_map, data.frame)
colnames(cis_map.df)[1] <- "phe"
for(i in c("lowposition","highposition")) cis_map.df[,i] <- as.numeric(cis_map.df[,i])

# Add marker names:
cis_map.df[3, "lowmarker"] <- "chr8_8_0.656384"
cis_map.df[3, "highmarker"] <- "chr8_8_3.118527"
cis_map.df[3, "lowposition"] <- "0.000001"
cis_map.df[3, "highposition"] <- "17.4757188188002"

cis_map.df[22, "lowmarker"] <- "chr3_3_12.540387"
cis_map.df[22, "highmarker"] <- "chr3_3_19.199451"
cis_map.df[22, "lowposition"] <- "122.827699700737"
cis_map.df[22, "highposition"] <- "161.028047128389"

cis_map.df[23, "lowmarker"] <- "chr3_3_26.88284"
cis_map.df[23, "highmarker"] <- "chr3_3_29.703667"
cis_map.df[23, "lowposition"] <- "193.967573263987"
cis_map.df[23, "highposition"] <- "207.661200044656"

str(cis_map.df)
cis_map.df$lowposition <- as.numeric(cis_map.df$lowposition)
cis_map.df$highposition <- as.numeric(cis_map.df$highposition)

# QTL regions on genetic map
segmentsOnMap(
  cross = map7,
  phe = cis_map.df$phe,
  chr = cis_map.df$chr,
  l = cis_map.df$lowposition,
  h = cis_map.df$highposition,
  peaklod = cis_map.df$maxLod,
  peakcM = cis_map.df$pos,
  showPeaks = TRUE,
  lwd = 2.5,
  leg.lwd	= 4,
  orderBy = cis_map.df$phe,
  segSpread = 1.5,
  legendCex = 1,
  chrBuffer = c(0.2, 0),
  leg.inset = 0.0,
  legendPosition = NULL,
  col = phe_colorp
)

cis_map.df.pos <- cis_map.df %>%
  mutate(qtl_start = as.numeric(sub('.*_(\\d+\\.?\\d*)$', '\\1', lowmarker)) * 1000000,
         qtl_end = as.numeric(sub('.*_(\\d+\\.?\\d*)$', '\\1', highmarker)) * 1000000)
write.csv(cis_map.df.pos, "./QTL_Results/seg_plot/map7_all_qtls_cross_control.csv", row.names = F)



map7_all_qtls_cross_control_peaks <- read.csv("./QTL_Results/seg_plot/map7_all_qtls_cross_control_peaks.csv", check.names = F)
cis_map.df.pos.peaks <- map7_all_qtls_cross_control_peaks %>%
  mutate(qtl_peak_pos = as.numeric(sub('.*_(\\d+\\.?\\d*)$', '\\1', qtl_peak)) * 1000000)
write.csv(cis_map.df.pos.peaks, "./QTL_Results/seg_plot/map7_all_qtls_cross_control_peaks_pos.csv", row.names = F)



# # Plot legend
# # Get unique phenotypes and their corresponding colors
# unique_phe <- unique(cis_map.df$phe)
# # Define colors
# phe_colors <- phe_colorp
# # Adjust transparency
# alpha <- 0.85  # Set the transparency level (0 = fully transparent, 1 = fully opaque)
# adjusted_colors <- adjustcolor(phe_colors, alpha = alpha)
# # Determine the number of rows and columns for the legend
# num_phe <- length(unique_phe)
# num_cols <- 5  # Number of columns for the legend
# num_rows <- ceiling(num_phe / num_cols)  # Calculate the number of rows
# # Create a plot with no axes or labels and no visible grid
# plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 0.3), axes = F, xaxs = "i", yaxs = "i")
# # Calculate box dimensions and text position
# box_width <- 0.18
# box_height <- 0.02
# text_offset <- 0.01
# # Loop through unique phenotypes
# for (i in 1:num_phe) {
#   # Calculate row and column indices
#   row_index <- num_rows - ((i - 1) %/% num_cols)
#   col_index <- (i - 1) %% num_cols + 1
#   # Calculate box position
#   box_x <- (col_index - 1) * (box_width + text_offset)
#   box_y <- (row_index - 1) * (box_height + text_offset)
#   # Draw colored box
#   rect(box_x, box_y, box_x + box_width, box_y + box_height, col = adjusted_colors[i])
#   # Add phenotype name inside the box
#   text((box_x ), (box_y + box_height / 2), unique_phe[i], pos = 4, cex = 0.8)
# }


# Select rows of ftime
ftime_cis_map.df <- subset(cis_map.df, phe == "Days to Flowering")
# QTL regions on genetic map
segmentsOnMap(
  cross = map7,
  phe = ftime_cis_map.df$phe,
  chr = ftime_cis_map.df$chr,
  l = ftime_cis_map.df$lowposition,
  h = ftime_cis_map.df$highposition,
  peaklod = ftime_cis_map.df$maxLod,
  peakcM = ftime_cis_map.df$pos,
  showPeaks = TRUE,
  #lwd = 8,
  leg.lwd	= 5,
  orderBy = ftime_cis_map.df$phe,
  segSpread = 5,
  legendCex = 1,
  chrBuffer = c(0.2, 0),
  leg.inset = 0.0,
  legendPosition = NULL,
  col = c("#ff0000", "#ff0000", "#ff0000", "#ff0000", "#ff0000")
)

# Select rows of fertility
fertility_cis_map.df <- subset(cis_map.df, phe == "Fertility Score")
# QTL regions on genetic map
segmentsOnMap(
  cross = map7,
  phe = fertility_cis_map.df$phe,
  chr = fertility_cis_map.df$chr,
  l = fertility_cis_map.df$lowposition,
  h = fertility_cis_map.df$highposition,
  peaklod = fertility_cis_map.df$maxLod,
  peakcM = fertility_cis_map.df$pos,
  showPeaks = TRUE,
  #lwd = 8,
  leg.lwd	= 5,
  orderBy = fertility_cis_map.df$phe,
  segSpread = 5,
  legendCex = 1,
  chrBuffer = c(0.2, 0),
  leg.inset = 0.0,
  legendPosition = NULL,
  col = c("#ffc000", "#ffc000", "#ffc000", "#ffc000")
)

########################## 
# 4. Visualization of QTL effect sizes----
# Plotting maxLod values
qtls_effect <- cis_map.df.pos %>%
  ggplot(aes(x = maxLod, fill = phe)) +
  geom_histogram(color = "#e9ecef", alpha = 0.9, position = 'stack', binwidth = 1) +
  scale_fill_manual(values = c("#b00058", "#ff0000", "#ffc000", "#7f7f7f",
                               "#0070c0", "#004272", "#002060", "#87afff",
                               "#266f8b", "#ff89c4", "#00b0f0",
                               "#502273", "#7030a0", "#b17ed8", "#dac2ec",
                               "#7fbd98", "#007b32", "#b9e08c",
                               "#533b1e", "#cfa879")) +
  labs(fill = "Phenotypes", #title = "QTLs and their maxLod values", 
       y = "Number of QTLs", x = "LOD Score") + 
  theme_bw() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10.5),         
    legend.title = element_text(size = 11, face = "bold"),        
    axis.title.x = element_text(size = 12,  margin = margin(t = 10)),        
    axis.title.y = element_text(size = 12, margin = margin(r = 10))
  )

# Display the plot
print(qtls_effect)




# Plot "effectplot" nicer!
# Neda

fertility_m_1 <- find.marker(map7, chr = 3, pos = 63.0)
fertility_m_2 <- find.marker(map7, chr = 5, pos = 0.0)
fertility_m_3 <- find.marker(map7, chr = 6, pos = 11.0)
fertility_m_4 <- find.marker(map7, chr = 7, pos = 50.0)
fertility_m_5 <- find.marker(map7, chr = 8, pos = 27.0)
par(mfrow=c(2,5))
effectplot(map7, pheno.col = 8, mname1 = "3@63.0")
effectplot(map7, pheno.col = 8, mname1 = "5@0.0")
effectplot(map7, pheno.col = 8, mname1 = "6@11.0")
effectplot(map7, pheno.col = 8, mname1 = "7@50.0")
effectplot(map7, pheno.col = 8, mname1 = "8@27.0")
plotPXG(map7, marker = c(fertility_m_1), pheno.col = 8, jitter = 1)
plotPXG(map7, marker = c(fertility_m_2), pheno.col = 8, jitter = 1)
plotPXG(map7, marker = c(fertility_m_3), pheno.col = 8, jitter = 1)
plotPXG(map7, marker = c(fertility_m_4), pheno.col = 8, jitter = 1)
plotPXG(map7, marker = c(fertility_m_5), pheno.col = 8, jitter = 1)


# plot effectplot
efp_f <- effectplot(map7, 
                    pheno.col = 3, 
                    mname1="3@63.0", mname2="7@31.0",
                    main = "Fertility Score QTLs (3,7)",
                    xlab = "Marker on Chr7",
                    ylab = "W.S = Total weight of seeds [g]",
                    col = c("#d1495b","#984ea3","#00798c"),
                    legend.lab = "Marker on Chr3")

# convert matrix data to data frame
fer_means <- efp_f$Means
fer_means <- as.data.frame(fer_means)
fer_se <- efp_f$SEs
fer_se <- as.data.frame(fer_se)

# melt data
dataframe_fer_means=as.data.frame(fer_means)
print(dataframe_fer_means)
print("Original data frame:\n") 
print(dataframe_fer_means) 
melt_data <- melt(dataframe_fer_means, id = dataframe_fer_means[,1]) 
colnames(melt_data) <- c("chr7","W.S")
print("Reshaped data frame:\n") 
print(melt_data) 
chr3 = c("chr3_3_8.36208.NN","chr3_3_8.36208.NS","chr3_3_8.36208.SS")
melt_data1 <- cbind(melt_data, chr3)
melt_data1

# plot
melt_data1 %>%
  ggplot( aes(x=chr7, y=W.S)) +
  geom_line(aes(group=chr3, color=as.factor(chr3))) +
  geom_point(aes( color=as.factor(chr3)))

# melt sd data
dataframe_fer_se=as.data.frame(fer_se)
print("Original data frame:\n") 
print(dataframe_fer_se) 
melt_data_sd <- melt(dataframe_fer_se, id = dataframe_fer_se[,1]) 
colnames(melt_data_sd) <- c("chr7","sd")
print("Reshaped data frame:\n") 
print(melt_data_sd) 
colsd <- melt_data_sd$sd 
melt_data2 <- cbind(melt_data1, colsd)

# extract Hexadecimal color
display.brewer.pal(3, name = "Dark2")
brewer.pal(3, name = "Dark2")
#  "#1B9E77" "#D95F02" "#7570B3"
# or
# "#edae49","#d1495b","#00798c" # I used these

# plot +sd
finalefp <-
  melt_data2 %>%
  ggplot(aes(x=chr7, y=W.S, color=as.factor(chr3))) +
  geom_line(aes(group=chr3), size = 1) +
  geom_point() +
  geom_errorbar(aes(ymin=W.S-colsd, ymax=W.S+colsd), width=0.05, alpha=0.8, size=0.6, linetype = "twodash") +
  #scale_colour_brewer(type="qual", palette="Dark2") +
  scale_x_discrete("\nGenotypes of Marker on Chr7", labels = c("chr7_7_8.505683.SS" = "SS","chr7_7_8.505683.NS" = "NS", "chr7_7_8.505683.NN" = "NN")) + 
  labs(title = "Interaction Plot for Fertility Score Markers on Chr3 and Chr7\n", y = "W.S = Total weight of seeds [g]\n", color = "Genotypes of Marker on Chr3\n") +
  scale_color_manual(labels = c("chr3_3_8.36208.NN"="NN","chr3_3_8.36208.NS"="NS","chr3_3_8.36208.SS"="SS"), values = c("#d1495b","#984ea3","#00798c")) +
  theme_light()
finalefp


# Only chromosome 3
# plot effectplot
efp_f_chr3 <- effectplot(map7,
                         pheno.col = 3, 
                         mname1="3@63.0",
                         #main = "Fertility Score QTLs (3)",
                         xlab = "Marker Genotype on Chromosome 3",
                         ylab = "Seed Production",
                         col = c("#ffc000"))

# convert matrix data to data frame
fer_means_chr3 <- efp_f_chr3$Means
fer_means_chr3 <- as.data.frame(fer_means_chr3)
fer_se_chr3 <- efp_f_chr3$SEs
fer_se_chr3 <- as.data.frame(fer_se_chr3)

# Merge the data frames by their row names
merged_df_chr3 <- merge(fer_means_chr3, fer_se_chr3, by = "row.names")

# Rename the row names back
rownames(merged_df_chr3) <- merged_df_chr3$Row.names

# Remove the Row.names column if you don't need it
merged_df_chr3 <- merged_df_chr3[, -1]

# First, prepare the data frame for plotting
merged_df_chr3$Genotypechr3 <- rownames(merged_df_chr3)
merged_df_chr3$Group <- c("NN", "NS", "SS")

# Create the line plot with error bars
ggplot(merged_df_chr3, aes(x = Group, y = fer_means_chr3, group = 1)) + 
  geom_line(color = "#ffc000", alpha = 0.9, size = 1.5) + 
  geom_point(color = "#ffc000", size = 3) +
  geom_errorbar(aes(ymin = fer_means_chr3 - fer_se_chr3, ymax = fer_means_chr3 + fer_se_chr3), 
                color = "#CC9A00", width=0.05, alpha=0.7, size=1, linetype = "longdash") +
  labs(x = "Marker Genotype on Chromosome 3", y = "Mean Seed Production") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),  
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.x = element_text(size = 11, margin = margin(t = 5),
                               color = c("#d1495b","#984ea3","#00798c"), face = "bold"),
    axis.text.y = element_text(size = 10)  
    ) +
  annotate("segment", x = 0.70, xend = 0.70, y = 0.12, yend = -0.08, 
           arrow = arrow(length = unit(0.2, "inches")), color = "#FF5555", size = 2) +
  annotate("text", x = 0.78, y = 0.011, label = "30%", color = "#FF5555", size = 5, hjust = 0)




########################## 
# 5. Save workspace----
save.image("single_qtl_mapping_residuals.RData")