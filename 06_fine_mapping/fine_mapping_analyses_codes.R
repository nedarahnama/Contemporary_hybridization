# < Genetic architecture of phenotypic differences between endangered hybridizing Arabis floodplain species >
# < PhD Thesis >
# < Neda Rahnamae> 

# < Chapter 3 - Flowering time QTL fine-mapping >
# < Flowering time fine-mapping experiment on F3 hybrids >
# October 2024
# 422 individuals - 15 F3 Families


########################## 
# 1. Load R libraries----
library(cowplot) # 1.1.3
library(dplyr) # 1.1.4
library(ggplot2) # 3.5.1
library(gridExtra) # 2.3
library(hrbrthemes) # 0.8.7
library(lme4) # 1.1-35.5
library(lmerTest) # 3.1-3
library(performance) # 0.12.3
library(reshape2) # 1.4.4
library(see) # 0.9.0
library(smplot2) # 0.2.4
library(tidyverse) # 2.0.0
library(viridis) # 0.6.5


# map6j_742_2082_gen_fine_mapping_families-173.csv = genetic data of F3 families 

########################## 
# 2. Data manipulation----
## Load data----
new.ftime_pcr_f3_nona <- read.csv("./new.ftime_pcr_f3_nona.csv", check.names = F)

## Infer genotypes----
# Define the list of columns you want to modify
columns_to_infer <- c("up_Fp", "int", "02int1.2", "DW2p", "p8p")
# Use mutate and across to apply the logic to those columns based on F3_genotype
new.ftime_pcr_f3_nona <- new.ftime_pcr_f3_nona %>%
  mutate(across(all_of(columns_to_infer), ~ ifelse(F3_genotype == "S", "S", 
                                                   ifelse(F3_genotype == "N", "N", 
                                                          ifelse(F3_genotype == "HH", "H", .)))))

new.ftime_pcr_f3_nona <- new.ftime_pcr_f3_nona %>%
  mutate(int = ifelse(is.na(int) & up_Fp == DW2p, up_Fp, int))

new.ftime_pcr_f3_nona <- new.ftime_pcr_f3_nona %>%
  mutate(`02int1.2` = ifelse(is.na(`02int1.2`) & int == DW2p, int, `02int1.2`))

new.ftime_pcr_f3_nona <- new.ftime_pcr_f3_nona %>%
  mutate(`02int1.2` = ifelse(is.na(`02int1.2`) & int == p8p, int, `02int1.2`))

new.ftime_pcr_f3_nona <- new.ftime_pcr_f3_nona %>%
  mutate(DW2p = ifelse(is.na(DW2p) & int == p8p, int, DW2p))

str(new.ftime_pcr_f3_nona)


## Count recombinations----
# Define the columns in which you want to count the changes
columns_to_check <- c("up_Fp", "int", "02int1.2", "DW2p", "p8p")

# Function to count changes in a row across specified columns
count_changes <- function(row) {
  # Remove NA values from the row
  row <- na.omit(row)
  
  # Count the number of transitions (changes) between consecutive elements in the row
  sum(row[-length(row)] != row[-1])
}

# Apply the function to each row of the dataframe and store in a new column 'nr_recom'
new.ftime_pcr_f3_nona$nr_recom <- apply(new.ftime_pcr_f3_nona[, columns_to_check], 1, count_changes)
str(new.ftime_pcr_f3_nona)
new.ftime_pcr_f3_nona$nr_recom <-as.factor(new.ftime_pcr_f3_nona$nr_recom)
summary(new.ftime_pcr_f3_nona$nr_recom)
# 0   1   2   3 
# 272 138   7   5 

# Filter the dataframe to keep only rows where nr_recom is less than or equal to 1
# Convert nr_recom to numeric if it's a factor (optional, depends on the situation)
new.ftime_pcr_f3_nona$nr_recom <- as.numeric(as.character(new.ftime_pcr_f3_nona$nr_recom))

# Filter the dataframe to keep only rows where nr_recom is not equal to 2 or 3
new.ftime_pcr_f3_nona.1recom <- new.ftime_pcr_f3_nona[!(new.ftime_pcr_f3_nona$nr_recom %in% c(2, 3)), ]
new.ftime_pcr_f3_nona.1recom$nr_recom <- as.factor(new.ftime_pcr_f3_nona.1recom$nr_recom)
summary(new.ftime_pcr_f3_nona.1recom$nr_recom)
# 0   1 
# 272 138 
str(new.ftime_pcr_f3_nona.1recom)

## Creat recom subset----
# Create a subset where p8p and up_Fp are different
new.ftime_pcr_f3_nona.1recom.H <- new.ftime_pcr_f3_nona.1recom[new.ftime_pcr_f3_nona.1recom$p8p != new.ftime_pcr_f3_nona.1recom$up_Fp, ]
str(new.ftime_pcr_f3_nona.1recom.H)


## Split experiments----
exp1_nona.1recom <- new.ftime_pcr_f3_nona.1recom %>%
  filter(exp == '1')
summary(exp1_nona.1recom$flowering_time)

exp2_nona.1recom <- new.ftime_pcr_f3_nona.1recom %>%
  filter(exp == '2')
summary(exp2_nona.1recom$flowering_time)

## Mean Ftime----
# Calculate the mean flowering time for each genotype group
# exp1
exp1_mean_ftime <- exp1_nona.1recom %>%
  group_by(genotype) %>%
  summarise(exp1_mean_ftime = mean(flowering_time, na.rm = TRUE))
# genotype exp1_mean_ftime
#   <fct>              <dbl>
# 1 170                 226.
# 2 173                 235.
# 3 263                 233.
# 4 460                 225.
# 5 523                 227.
# 6 643                 224.
# 7 719                 233.
# 8 741                 226.
# 9 825                 234.
# 10 885                231.
# 11 1047               228.
# 12 1094               222.
# 13 1101               225.
# 14 1110               233.
# 15 1117               224.
exp1_nona.1recom.1094 <- exp1_nona.1recom %>%
  filter(genotype == '1094')
mean(exp1_nona.1recom.1094$flowering_time)
sd(exp1_nona.1recom.1094$flowering_time)
exp1_nona.1recom.825 <- exp1_nona.1recom %>%
  filter(genotype == '825')
mean(exp1_nona.1recom.825$flowering_time)
sd(exp1_nona.1recom.825$flowering_time)

# exp2
exp2_mean_ftime <- exp2_nona.1recom %>%
  group_by(genotype) %>%
  summarise(exp2_mean_ftime = mean(flowering_time, na.rm = TRUE))
# genotype exp2_mean_ftime
# <fct>              <dbl>
# 1 173                 200.
# 2 263                 202.
# 3 523                 204.
# 4 719                 202.
# 5 741                 198.
# 6 885                 204.
# 7 1094                198.
# 8 1101                201.
# 9 1110                201.
# 10 1117               201.
exp2_nona.1recom.1094 <- exp2_nona.1recom %>%
  filter(genotype == '1094')
mean(exp2_nona.1recom.1094$flowering_time)
sd(exp2_nona.1recom.1094$flowering_time)
exp2_nona.1recom.885 <- exp2_nona.1recom %>%
  filter(genotype == '885')
str(exp2_nona.1recom.885)
mean(exp2_nona.1recom.885$flowering_time, na.rm = T)
sd(exp2_nona.1recom.885$flowering_time, na.rm = T)
exp2_nona.1recom.173 <- exp2_nona.1recom %>%
  filter(genotype == '173')
mean(exp2_nona.1recom.173$flowering_time)
sd(exp2_nona.1recom.173$flowering_time)
exp2_nona.1recom.263 <- exp2_nona.1recom %>%
  filter(genotype == '263')
mean(exp2_nona.1recom.263$flowering_time)
sd(exp2_nona.1recom.263$flowering_time)


########################## 
# 3. Residuals----
# Define the conversion function
convert_values <- function(value) {
  if (is.na(value)) {
    return(NA)
  } else if (value == "N") {
    return(0)
  } else if (value == "S") {
    return(2)
  } else if (value %in% c("H", "NH", "SH", "HN", "HS", "NS", "SN")) {
    return(1)
  } else {
    return(NA)  # Return NA for any other unexpected values
  }
}

## ftime----
new.ftime_pcr_f3_nona.1recomu <- new.ftime_pcr_f3_nona.1recom[!is.na(new.ftime_pcr_f3_nona.1recom$DW2p), ]
new.ftime_pcr_f3_nona.1recomu <- new.ftime_pcr_f3_nona.1recom[!is.na(new.ftime_pcr_f3_nona.1recom$`02int1.2`), ]
new.ftime_pcr_f3_nona.1recomu <- new.ftime_pcr_f3_nona.1recom[!is.na(new.ftime_pcr_f3_nona.1recom$flowering_time), ]

modres_ftime <- glm(flowering_time ~ genotype * (exp/tray_garden), 
                    family = "quasipoisson", data = new.ftime_pcr_f3_nona.1recomu)
summary(modres_ftime)
res_ftime <- modres_ftime$residuals

hist(res_ftime, 
     main = "Histogram of Residuals (res_ftime) with Normal Curve", 
     xlab = "Residuals", 
     ylab = "Frequency", 
     col = "lightblue", 
     border = "black", 
     probability = TRUE)
curve(dnorm(x, mean=mean(res_ftime), sd=sd(res_ftime)), 
      col = "red", lwd = 2, add = TRUE)

# Perform a Shapiro-Wilk test for normality on res_ftime
shapiro_test <- shapiro.test(res_ftime) # not normal

### QTL Region----
# Create the new columns based on the conditions
# Primers: "up_Fp", "int", "02int1.2", "DW2p", "p8p"
str(new.ftime_pcr_f3_nona.1recomu)

new.ftime_pcr_f3_nona.1recomu.i <- new.ftime_pcr_f3_nona.1recomu %>%
  mutate(
    interval1 = case_when(
      is.na(up_Fp) ~ int,
      is.na(int) ~ up_Fp,
      up_Fp == int ~ up_Fp,
      TRUE ~ paste(up_Fp, int, sep = "")
    ),
    interval2 = case_when(
      is.na(int) ~ `02int1.2`,
      is.na(`02int1.2`) ~ int,
      int == `02int1.2` ~ int,
      TRUE ~ paste(int, `02int1.2`, sep = "")
    ),
    interval3 = case_when(
      is.na(`02int1.2`) ~ DW2p,
      is.na(DW2p) ~ `02int1.2`,
      `02int1.2` == DW2p ~ `02int1.2`,
      TRUE ~ paste(`02int1.2`, DW2p, sep = "")
    ),
    interval4 = case_when(
      is.na(DW2p) ~ p8p,
      is.na(p8p) ~ DW2p,
      DW2p == p8p ~ DW2p,
      TRUE ~ paste(DW2p, p8p, sep = "")
    )
  )


# Convert the values in interval columns to numbers
new.ftime_pcr_f3_nona.1recomu.i <- new.ftime_pcr_f3_nona.1recomu.i %>%
  mutate_at(vars(interval1, interval2, interval3, interval4), ~sapply(., convert_values))

# Run the linear model
model.resftime.int <- 
  glm(res_ftime ~ interval1 + interval2 + interval3 + interval4, 
      family = "gaussian", data = new.ftime_pcr_f3_nona.1recomu.i)
summary(model.resftime.int)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  0.0027195  0.0011501   2.364   0.0186 *
#   interval1   -0.0016936  0.0018858  -0.898   0.3697  
#   interval2   -0.0132025  0.0052176  -2.530   0.0118 *
#   interval3    0.0120917  0.0056111   2.155   0.0318 *
#   interval4    0.0003302  0.0027229   0.121   0.9035  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 0.0001210603)
# 
# Null deviance: 0.047290  on 379  degrees of freedom
# Residual deviance: 0.045398  on 375  degrees of freedom
# (27 observations deleted due to missingness)
# AIC: -2341.9
# 
# Number of Fisher Scoring iterations: 2


# Define the intervals and corresponding p-values
intervals <- c("Interval1", "Interval2", "Interval3", "Interval4")
p_values <- c(0.3697, 0.0118, 0.0318, 0.9035)

# Calculate -log10 of p-values
neg_log_p_values <- -log10(p_values)

# Create a data frame to use in ggplot
data.ftime.int <- data.frame(
  intervals = factor(intervals, levels = intervals),  # Convert to factor for ggplot
  neg_log_p_values = neg_log_p_values
)

# Create the ggplot
ggplot(data.ftime.int, aes(x = intervals, y = neg_log_p_values)) +
  geom_point(size = 4, color = "#404040") +
  geom_line(group = 1, color = "#404040") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#FF5555") +  # Significance line
  theme_bw() +
  theme(
    #plot.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 12), # Increase legend title size
    axis.title.x = element_text(size = 12, margin = margin(t = 2)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  labs(x = "Intervals", y = expression(-log[10](P)))





## ph----
modres_ph <- glm(plant_height ~ genotype * (exp/tray_garden), 
                 family = "quasipoisson", data = new.ftime_pcr_f3_nona.1recomu)
summary(modres_ph)
res_ph <- modres_ph$residuals
new.ftime_pcr_f3_nona.1recomu.ph <- new.ftime_pcr_f3_nona.1recom[!is.na(new.ftime_pcr_f3_nona.1recom$plant_height), ]


## sh----
new.ftime_pcr_f3_nona.1recomu.sh <- new.ftime_pcr_f3_nona.1recom[!is.na(new.ftime_pcr_f3_nona.1recom$stem_height), ]
modres_sh <- glm(stem_height ~ genotype * (exp/tray_garden), 
                 family = "quasipoisson", data = new.ftime_pcr_f3_nona.1recomu.sh)
summary(modres_sh)
res_sh <- modres_sh$residuals

### QTL Region----
# Create the new columns based on the conditions
# Primers: "up_Fp", "int", "02int1.2", "DW2p", "p8p"

new.ftime_pcr_f3_nona.1recomu.sh.i <- new.ftime_pcr_f3_nona.1recomu.sh %>%
  mutate(
    interval1 = case_when(
      is.na(up_Fp) ~ int,
      is.na(int) ~ up_Fp,
      up_Fp == int ~ up_Fp,
      TRUE ~ paste(up_Fp, int, sep = "")
    ),
    interval2 = case_when(
      is.na(int) ~ `02int1.2`,
      is.na(`02int1.2`) ~ int,
      int == `02int1.2` ~ int,
      TRUE ~ paste(int, `02int1.2`, sep = "")
    ),
    interval3 = case_when(
      is.na(`02int1.2`) ~ DW2p,
      is.na(DW2p) ~ `02int1.2`,
      `02int1.2` == DW2p ~ `02int1.2`,
      TRUE ~ paste(`02int1.2`, DW2p, sep = "")
    ),
    interval4 = case_when(
      is.na(DW2p) ~ p8p,
      is.na(p8p) ~ DW2p,
      DW2p == p8p ~ DW2p,
      TRUE ~ paste(DW2p, p8p, sep = "")
    )
  )


# Convert the values in interval columns to numbers
new.ftime_pcr_f3_nona.1recomu.sh.i <- new.ftime_pcr_f3_nona.1recomu.sh.i %>%
  mutate_at(vars(interval1, interval2, interval3, interval4), ~sapply(., convert_values))

# Run the linear model
model.ressh.int <- 
  glm(res_sh ~ interval1 + interval2 + interval3 + interval4, 
      family = "gaussian", data = new.ftime_pcr_f3_nona.1recomu.sh.i)
summary(model.ressh.int)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  0.032621   0.015562   2.096   0.0367 *
# interval1    0.012197   0.025549   0.477   0.6334  
# interval2    0.017870   0.070687   0.253   0.8006  
# interval3   -0.050986   0.076543  -0.666   0.5058  
# interval4   -0.004991   0.037939  -0.132   0.8954  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 0.02221938)
# 
# Null deviance: 8.5302  on 380  degrees of freedom
# Residual deviance: 8.3545  on 376  degrees of freedom
# (27 observations deleted due to missingness)
# AIC: -362.19
# 
# Number of Fisher Scoring iterations: 2



## inh----
new.ftime_pcr_f3_nona.1recomu.inh <- new.ftime_pcr_f3_nona.1recom[!is.na(new.ftime_pcr_f3_nona.1recom$inflorescence_height), ]
modres_inh <- glm(inflorescence_height ~ genotype * (exp/tray_garden), 
                  family = "gaussian", data = new.ftime_pcr_f3_nona.1recomu.inh)
summary(modres_inh)
res_inh <- modres_inh$residuals

### QTL Region----
# Create the new columns based on the conditions
# Primers: "up_Fp", "int", "02int1.2", "DW2p", "p8p"

new.ftime_pcr_f3_nona.1recomu.inh.i <- new.ftime_pcr_f3_nona.1recomu.inh %>%
  mutate(
    interval1 = case_when(
      is.na(up_Fp) ~ int,
      is.na(int) ~ up_Fp,
      up_Fp == int ~ up_Fp,
      TRUE ~ paste(up_Fp, int, sep = "")
    ),
    interval2 = case_when(
      is.na(int) ~ `02int1.2`,
      is.na(`02int1.2`) ~ int,
      int == `02int1.2` ~ int,
      TRUE ~ paste(int, `02int1.2`, sep = "")
    ),
    interval3 = case_when(
      is.na(`02int1.2`) ~ DW2p,
      is.na(DW2p) ~ `02int1.2`,
      `02int1.2` == DW2p ~ `02int1.2`,
      TRUE ~ paste(`02int1.2`, DW2p, sep = "")
    ),
    interval4 = case_when(
      is.na(DW2p) ~ p8p,
      is.na(p8p) ~ DW2p,
      DW2p == p8p ~ DW2p,
      TRUE ~ paste(DW2p, p8p, sep = "")
    )
  )


# Convert the values in interval columns to numbers
new.ftime_pcr_f3_nona.1recomu.inh.i <- new.ftime_pcr_f3_nona.1recomu.inh.i %>%
  mutate_at(vars(interval1, interval2, interval3, interval4), ~sapply(., convert_values))

# Run the linear model
model.resinh.int <- 
  glm(res_inh ~ interval1 + interval2 + interval3 + interval4, 
      family = "gaussian", data = new.ftime_pcr_f3_nona.1recomu.inh.i)
summary(model.resinh.int)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)  -0.9908     0.6584  -1.505    0.133
# interval1    -1.1274     1.0811  -1.043    0.298
# interval2     0.2495     2.9909   0.083    0.934
# interval3     0.7909     3.2387   0.244    0.807
# interval4     0.9603     1.6053   0.598    0.550
# 
# (Dispersion parameter for gaussian family taken to be 39.78041)
# 
# Null deviance: 15246  on 380  degrees of freedom
# Residual deviance: 14957  on 376  degrees of freedom
# (27 observations deleted due to missingness)
# AIC: 2491.6
# 
# Number of Fisher Scoring iterations: 2



## rd----
new.ftime_pcr_f3_nona.1recomu.rd <- new.ftime_pcr_f3_nona.1recom[!is.na(new.ftime_pcr_f3_nona.1recom$rd_cm), ]
modres_rd <- glm(rd_cm ~ genotype * (exp/tray_garden), 
                 family = "quasipoisson", data = new.ftime_pcr_f3_nona.1recomu.rd)
summary(modres_rd)
res_rd <- modres_rd$residuals

### QTL Region----
# Create the new columns based on the conditions
# Primers: "up_Fp", "int", "02int1.2", "DW2p", "p8p"

new.ftime_pcr_f3_nona.1recomu.rd.i <- new.ftime_pcr_f3_nona.1recomu.rd %>%
  mutate(
    interval1 = case_when(
      is.na(up_Fp) ~ int,
      is.na(int) ~ up_Fp,
      up_Fp == int ~ up_Fp,
      TRUE ~ paste(up_Fp, int, sep = "")
    ),
    interval2 = case_when(
      is.na(int) ~ `02int1.2`,
      is.na(`02int1.2`) ~ int,
      int == `02int1.2` ~ int,
      TRUE ~ paste(int, `02int1.2`, sep = "")
    ),
    interval3 = case_when(
      is.na(`02int1.2`) ~ DW2p,
      is.na(DW2p) ~ `02int1.2`,
      `02int1.2` == DW2p ~ `02int1.2`,
      TRUE ~ paste(`02int1.2`, DW2p, sep = "")
    ),
    interval4 = case_when(
      is.na(DW2p) ~ p8p,
      is.na(p8p) ~ DW2p,
      DW2p == p8p ~ DW2p,
      TRUE ~ paste(DW2p, p8p, sep = "")
    )
  )


# Convert the values in interval columns to numbers
new.ftime_pcr_f3_nona.1recomu.rd.i <- new.ftime_pcr_f3_nona.1recomu.rd.i %>%
  mutate_at(vars(interval1, interval2, interval3, interval4), ~sapply(., convert_values))

# Run the linear model
model.resrd.int <- 
  glm(res_rd ~ interval1 + interval2 + interval3 + interval4, 
      family = "gaussian", data = new.ftime_pcr_f3_nona.1recomu.rd.i)
summary(model.resrd.int)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.004927   0.011046   0.446    0.656
# interval1    0.015881   0.018330   0.866    0.387
# interval2   -0.043532   0.047631  -0.914    0.361
# interval3    0.016161   0.052112   0.310    0.757
# interval4    0.009660   0.027276   0.354    0.723
# 
# (Dispersion parameter for gaussian family taken to be 0.009939647)
# 
# Null deviance: 3.2569  on 330  degrees of freedom
# Residual deviance: 3.2403  on 326  degrees of freedom
# (22 observations deleted due to missingness)
# AIC: -580.02
# 
# Number of Fisher Scoring iterations: 2


## plot----
# Define the intervals and corresponding p-values for each model
intervals <- c("Interval1", "Interval2", "Interval3", "Interval4")
interval_length_kb <- c(1150000, 294003, 234479, 751879) / 1000

# P-values for the four models
p_values_ftime <- c(0.3697, 0.0118, 0.0318, 0.9035)
p_values_sh <- c(0.6334, 0.8006, 0.5058, 0.8954)
p_values_inh <- c(0.298, 0.934, 0.807, 0.550)
p_values_rd <- c(0.387, 0.361, 0.757, 0.723)

# Calculate -log10 of p-values
neg_log_p_values_ftime <- -log10(p_values_ftime)
neg_log_p_values_sh <- -log10(p_values_sh)
neg_log_p_values_inh <- -log10(p_values_inh)
neg_log_p_values_rd <- -log10(p_values_rd)

# Combine data into a single data frame for ggplot
data <- data.frame(
  intervals = factor(rep(intervals, 4), levels = intervals),
  neg_log_p_values = c(neg_log_p_values_ftime, neg_log_p_values_sh, neg_log_p_values_inh, neg_log_p_values_rd),
  model = factor(rep(c("Flowering Time", "Stem Height", "Inflorescence Height", "Rosette Diameter"), each = 4))
)


# Create the ggplot with specified colors for each model
ggplot(data, aes(x = intervals, y = neg_log_p_values, color = model, group = model)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_line(size = 1, alpha = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "twodash", color = "#808080", size = 0.5) +  # Significance line
  theme_bw() +
  theme(
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 12), # Increase legend title size
    axis.title.x = element_text(size = 12, margin = margin(t = 2)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  labs(x = "Intervals", y = expression(-log[10](P)), color = "Trait") +
  scale_color_manual(values = c(
    "Flowering Time" = "#ff0000",        
    "Inflorescence Height" = "#7f7f7f", 
    "Rosette Diameter" = "#7030a0",  
    "Stem Height"  = "#007b32"       
  ))



## Segment plot----

# Define the intervals and corresponding lengths (in kb)
intervals <- c("Interval1", "Interval2", "Interval3", "Interval4")
interval_length_kb <- c(1165727, 293759, 233614, 748207) / 1000

# P-values for the four models
p_values_ftime <- c(0.3697, 0.0118, 0.0318, 0.9035)
p_values_sh <- c(0.6334, 0.8006, 0.5058, 0.8954)
p_values_inh <- c(0.298, 0.934, 0.807, 0.550)
p_values_rd <- c(0.387, 0.361, 0.757, 0.723)

# Calculate -log10 of p-values
neg_log_p_values_ftime <- -log10(p_values_ftime)
neg_log_p_values_sh <- -log10(p_values_sh)
neg_log_p_values_inh <- -log10(p_values_inh)
neg_log_p_values_rd <- -log10(p_values_rd)

# Define the starting position
start_position_bp <- 665597
# Calculate cumulative positions for each interval in bp
cumulative_lengths_bp <- cumsum(c(0, head(interval_length_kb * 1000, -1)))
start_positions_bp <- start_position_bp + cumulative_lengths_bp
end_positions_bp <- start_positions_bp + (interval_length_kb * 1000)

# Combine data into a single data frame for ggplot
new.data <- data.frame(
  intervals = intervals,
  start_position = start_positions_bp,
  end_position = end_positions_bp,
  neg_log_p_values_ftime = neg_log_p_values_ftime,
  neg_log_p_values_sh = neg_log_p_values_sh,
  neg_log_p_values_inh = neg_log_p_values_inh,
  neg_log_p_values_rd = neg_log_p_values_rd
)

# Reshape the data for ggplot
data_long <- pivot_longer(
  new.data,
  cols = starts_with("neg_log_p_values"),
  names_to = "model",
  values_to = "neg_log_p_values"
)

# Rename model values for better readability
data_long$model <- factor(data_long$model, levels = c(
  "neg_log_p_values_ftime",
  "neg_log_p_values_sh",
  "neg_log_p_values_inh",
  "neg_log_p_values_rd"
), labels = c("Flowering Time", "Stem Height", "Inflorescence Height", "Rosette Diameter"))


# Create the ggplot with connected line segments and intervals as x-axis labels
ggplot(data_long, aes(x = start_position, y = neg_log_p_values, color = model, group = model)) +
  geom_segment(aes(xend = end_position, yend = neg_log_p_values), size = 1.2, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#808080", size = 0.5) +  # Significance line
  #annotate("text", x = 2850000, y = -log10(0.05) + 0.1, label = "p = 0.05", 
          # size = 4.5, color = "#808080", hjust = 0) +  # Annotation for the significance line
  theme_bw() +
  theme(
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 12), # Increase legend title size
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.x = element_text(size = 12, angle = 30, hjust = 1),  # Rotate x-axis text
    axis.text.y = element_text(size = 11),
    axis.ticks.x = element_blank()  # Remove x-axis ticks
  ) +
  labs(x = "Intervals", y = expression(-log[10](P)), color = "Trait") +
  scale_color_manual(values = c(
    "Flowering Time" = "#ff0000",        
    "Inflorescence Height" = "#7f7f7f", 
    "Rosette Diameter" = "#7030a0",  
    "Stem Height" = "#007b32"
  )) +
  scale_x_continuous(
    breaks = (start_positions_bp + end_positions_bp) / 2,  # Center label at midpoint of each segment
    labels = intervals
  )


# 4. Ftime Data visualization----
# ftime
barplot_ftime_pcr_f3 <- 
  new.ftime_pcr_f3_nona.1recomu.i %>%
  ggplot(aes(x=flowering_time, fill = exp, color = exp)) +
  geom_bar(aes(color = exp), alpha = 0.3) +
  scale_fill_manual(values = c("1" = "#e7d4e8", "2" = "#a6dba0")) +  
  scale_color_manual(values = c("1" = "#762a83", "2" = "#1b7837")) +  
  theme_bw() +
  theme(
    #plot.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 12), # Increase legend title size
    axis.title.x = element_text(size = 12, margin = margin(t = 2)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  #ggtitle(expression(bold("Flowering time distribution in"~italic("Arabis")~"F3 hybrids"))) 
  labs(y = "Count", x = "Flowering Time", fill = "Trial", color = "Trial")
barplot_ftime_pcr_f3


boxplot_ftime_pcr_f3 <-
  new.ftime_pcr_f3_nona.1recomu.i %>%
  ggplot(aes(x = reorder(genotype, flowering_time, na.rm = TRUE), 
             y = flowering_time, fill = as.factor(exp), color = as.factor(exp))) +  
  geom_boxplot(aes(color = exp), alpha = 0.1, outlier.size = -1, position = position_dodge(width = 0.10)) +  
  geom_jitter(aes(fill = as.factor(ifelse(new.ftime_pcr_f3_nona.1recomu.i$F3_genotype == "H", exp, NA))),
              height = 0.5, width = 0.2, shape = 21, size = 1.5, stroke = 0.5, alpha = ifelse(new.ftime_pcr_f3_nona.1recomu.i$F3_genotype == "H", 1, 0.7)) +  
  scale_fill_manual(values = c("1" = "#762a83", "2" = "#1b7837"), na.value = "white", name = "Recombinant") +  
  scale_color_manual(values = c("1" = "#762a83", "2" = "#1b7837"), name = "Trial") +  
  theme_bw() +
  theme(
    #plot.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 12), # Increase legend title size
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  #ggtitle(expression(bold("Flowering time distribution in"~italic("Arabis")~"F3 hybrids"))) 
  labs(y = "Flowering Time", x = "F3 Families") 
boxplot_ftime_pcr_f3

mixplot_ftime_pcr_f3 <- grid.arrange(barplot_ftime_pcr_f3, boxplot_ftime_pcr_f3, nrow=2)

# Combined label 
combined_grid <- plot_grid(barplot_ftime_pcr_f3, boxplot_ftime_pcr_f3,
                           labels = "AUTO", ncol = 1, nrow = 2)
combined_grid



# 5. Phenotypes Data analysis----
model_ftime_ph <- lmerTest::lmer(flowering_time ~ plant_height * genotype + (1|exp) , data = new.ftime_pcr_f3_nona.1recomu.i)
summary(model_ftime_ph)
# the relationship between plant height and flowering time is not significantly different across genotypes
anova(model_ftime_ph)
# Type III Analysis of Variance Table with Satterthwaite's method
#                       Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# plant_height           14.61  14.614     1 418.17  0.5988 0.4394667    
# genotype              959.89  68.564    14 418.01  2.8095 0.0004973 ***
# plant_height:genotype 737.34  52.667    14 418.01  2.1581 0.0086970 ** -643
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
pf(2.16, df1 = 14, df2 = 418, lower.tail = FALSE)
# 0.008628692

# Get the rows used in the model (non-missing values)
used_rows <- rownames(ftime_pcr_f3)[!is.na(predict(model_ftime_ph, ftime_pcr_f3))]
# Subset the data
all_f3_data_used <- ftime_pcr_f3[used_rows, ]
# ake the predictions
all_f3_data_used$pred <- predict(model_ftime_ph, all_f3_data_used)

mm_plot <- ggplot(all_f3_data_used, aes(x = plant_height, y = flowering_time, colour = exp)) +
  facet_wrap(~genotype, nrow = 2) +   # a panel for each genotype
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_line(aes(y = pred)) +  # use pred directly from all_f3_data_used
  theme(legend.position = "none",
        panel.spacing = unit(2, "lines"))  # adding space between panels


model_ftime_sh <- lmerTest::lmer(flowering_time ~ stem_height * genotype + (1|exp) , data = new.ftime_pcr_f3_nona.1recomu.i)
summary(model_ftime_sh)
# the relationship between stem height and flowering time is not significantly different across genotypes
anova(model_ftime_sh)
# Type III Analysis of Variance Table with Satterthwaite's method
#                      Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# stem_height          124.78 124.779     1 418.09  5.0987 0.024459 * 
# genotype             832.92  59.494    14 418.01  2.4310 0.002709 **
# stem_height:genotype 642.91  45.922    14 418.01  1.8765 0.027120 * 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
pf(1.89, df1 = 14, df2 = 418, lower.tail = FALSE)
# 0.02572451


model_ftime_inh <- lmerTest::lmer(flowering_time ~ inflorescence_height * genotype + (1|exp) , data = new.ftime_pcr_f3_nona.1recomu.i)
summary(model_ftime_inh)
# the relationship between inflorescence height and flowering time is not significantly different across genotypes
anova(model_ftime_inh)
# Type III Analysis of Variance Table with Satterthwaite's method
#                               Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# inflorescence_height          258.78 258.784     1 418.13 10.5117 0.001281 **
# genotype                      845.75  60.411    14 418.01  2.4539 0.002452 **
# inflorescence_height:genotype 316.62  22.616    14 418.01  0.9186 0.538535   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
pf(0.9186, df1 = 14, df2 = 418, lower.tail = FALSE)
# 0.5385798


model_ftime_shn<- lmerTest::lmer(flowering_time ~ shoot_number * genotype + (1|exp) , data = ftime_pcr_f3)
summary(model_ftime_shn)
# the relationship between shoot number and flowering time is not significantly different across genotypes
anova(model_ftime_shn)
# Type III Analysis of Variance Table with Satterthwaite's method
#                       Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# shoot_number          365.13  365.13     1 417.03 14.9066 0.0001309 ***
# genotype              727.62   51.97    14 417.00  2.1218 0.0101173 *  
# shoot_number:genotype 358.77   25.63    14 417.00  1.0462 0.4059386    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
pf(1.0462, df1 = 14, df2 = 417, lower.tail = FALSE)
# 0.4059386


model_ftime_ld <- lmerTest::lmer(flowering_time ~ leaves_distance_cm * genotype + (1|exp) , data = ftime_pcr_f3)
summary(model_ftime_ld)
# the relationship between leaves distance and flowering time is not significantly different across genotypes
anova(model_ftime_ld)
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# leaves_distance_cm            4.84   4.840     1 409.01  0.1967 0.657624   
# genotype                    855.96  61.140    14 409.00  2.4847 0.002155 **
# leaves_distance_cm:genotype 398.53  28.467    14 409.00  1.1569 0.306393   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
pf(1.1569, df1 = 14, df2 = 409, lower.tail = FALSE)
# 0.3063633


model_ftime_sld <- lmerTest::lmer(flowering_time ~ stem_leaf_density * genotype + (1|exp) , data = ftime_pcr_f3)
summary(model_ftime_sld)
# the relationship between stem leaf density and flowering time is not significantly different across genotypes
anova(model_ftime_sld)
# Type III Analysis of Variance Table with Satterthwaite's method
#                            Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# stem_leaf_density            0.18   0.181     1 408.01  0.0074 0.93150  
# genotype                   637.45  45.532    14 408.01  1.8589 0.02911 *
# stem_leaf_density:genotype 414.57  29.612    14 408.00  1.2090 0.26556  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
pf(1.2090, df1 = 14, df2 = 408, lower.tail = FALSE)
# 0.2655361


model_ftime_rd <- glm(flowering_time ~ rd_cm * genotype * exp, family = "quasipoisson" , data = new.ftime_pcr_f3_nona.1recomu.i)
summary(model_ftime_rd)
# the relationship between rd and flowering time is not significantly different across genotypes
anova(model_ftime_rd)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# rd_cm          130.75 130.751     1 338.24  5.5777 0.018757 * 
# genotype       727.49  51.964    14 338.01  2.2167 0.007135 **
# rd_cm:genotype 380.92  27.208    14 338.01  1.1607 0.304265   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
pf(1.1607, df1 = 14, df2 = 338, lower.tail = FALSE)
# 0.3042513





# 6. Correlation plots----
# all
# Reshape data to long format
f3_long <- new.ftime_pcr_f3_nona.1recomu.i %>%
  pivot_longer(cols = c(plant_height, stem_height, inflorescence_height, shoot_number, leaves_distance_cm, stem_leaf_density, rd_cm),
               names_to = "variable",
               values_to = "value") %>%
  mutate(variable = factor(variable, levels = c("plant_height", "stem_height", "inflorescence_height", "shoot_number", "leaves_distance_cm", "stem_leaf_density", "rd_cm")))

# Create combined plot
ggplot(data = f3_long, mapping = aes(x = value, y = flowering_time, color = exp)) +
  geom_point(aes(color = exp), fill = "white", shape = 21, size = 2, alpha = 0.5) +
  sm_statCorr(text_size = 3) +
  scale_color_manual(values = c("1" = "#762a83", "2" = "#1b7837")) +  
  facet_wrap(~variable, scales = "free_x", nrow = 2, 
             labeller = labeller(variable = c(plant_height = "Plant Height", 
                                              stem_height = "Stem Height",
                                              inflorescence_height = "Inflorescence Height",
                                              shoot_number = "Number of Shoots",
                                              leaves_distance_cm = "Internode Length",
                                              stem_leaf_density = "Stem Leaf Density", 
                                              rd_cm = "Rosette Diameter"))) + 
  labs(y = "Flowering Time", x = "Phenotypes", color = "Trial" ) +
  theme_bw()


# Create combined plot
ggplot(data = f3_long, mapping = aes(x = value, y = flowering_time, color = as.factor(exp))) +
  geom_point(aes(fill = as.factor(ifelse(f3_long$F3_genotype == "H", exp, NA))),
             shape = 21, size = 2, alpha = ifelse(f3_long$F3_genotype == "H", 1, 0.5), stroke = 0.5) +
  sm_statCorr(text_size = 4) +
  scale_fill_manual(values = c("1" = "#762a83", "2" = "#1b7837"), na.value = "white", name = "Recombinant") +
  scale_color_manual(values = c("1" = "#762a83", "2" = "#1b7837"), name = "Trial") +  
  facet_wrap(~variable, scales = "free_x", nrow = 2, 
             labeller = labeller(variable = c(plant_height = "Plant Height", 
                                              stem_height = "Stem Height",
                                              inflorescence_height = "Inflorescence Height",
                                              shoot_number = "Number of Shoots",
                                              leaves_distance_cm = "Internode Length",
                                              stem_leaf_density = "Stem Leaf Density", 
                                              rd_cm = "Rosette Diameter"))) + 
  labs(y = "Flowering Time", x = "Phenotypes") +
  theme_bw() + 
  theme(
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 12), # Increase legend title size
    axis.title.x = element_text(size = 12, margin = margin(t = 5)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.y = element_text(size = 10),
    strip.text.x = element_text(size = 12)  # Increase facet label text size
  )


# ftime - plant height
# ggplot(data = new.ftime_pcr_f3_nona.1recomu.i, mapping = aes(x = plant_height, y = flowering_time, colour = exp)) +
#   geom_point(aes(color = exp), fill = "white", shape = 21, size = 2, alpha = 0.5) +
#   sm_statCorr(text_size = 3) +
#   scale_color_manual(values = c("1" = "#762a83", "2" = "#1b7837")) +  
#   facet_wrap(~genotype, nrow = 3) +
#   labs(y = "Flowering Time", x = "Plant Height", color = "Trial") +
#   theme_bw()
# 
# 
# 
# ggplot(data = new.ftime_pcr_f3_nona.1recomu.i, mapping = aes(x = plant_height, y = flowering_time, colour = exp)) +
#   geom_point(aes(fill = as.factor(ifelse(new.ftime_pcr_f3_nona.1recomu.i$F3_genotype == "H", exp, NA))),
#              shape = 21, size = 2, alpha = ifelse(new.ftime_pcr_f3_nona.1recomu.i$F3_genotype == "H", 1, 0.5), stroke = 0.5) +
#   sm_statCorr(text_size = 3) +
#   scale_fill_manual(values = c("1" = "#762a83", "2" = "#1b7837"), na.value = "white") +
#   scale_color_manual(values = c("1" = "#762a83", "2" = "#1b7837")) +
#   facet_wrap(~genotype, nrow = 3) +
#   labs(y = "Flowering Time", x = "Plant Height", color = "Trial", fill = "Recombinant") +
#   theme_bw()


# ftime - stem height
ggplot(data = new.ftime_pcr_f3_nona.1recomu.i, mapping = aes(x = stem_height, y = flowering_time, colour = exp)) +
  geom_point(aes(color = exp), fill = "white", shape = 21, size = 2, alpha = 0.5) +
  sm_statCorr(text_size = 3) +
  scale_color_manual(values = c("1" = "#762a83", "2" = "#1b7837")) +  
  facet_wrap(~genotype, nrow = 3) +
  labs(y = "Flowering Time", x = "Stem Height", color = "Trial") +
  theme_bw()

ggplot(data = new.ftime_pcr_f3_nona.1recomu.i, mapping = aes(x = stem_height, y = flowering_time, colour = exp)) +
  geom_point(aes(fill = as.factor(ifelse(new.ftime_pcr_f3_nona.1recomu.i$F3_genotype == "H", exp, NA))),
             shape = 21, size = 2, alpha = ifelse(new.ftime_pcr_f3_nona.1recomu.i$F3_genotype == "H", 1, 0.5), stroke = 0.5) +
  sm_statCorr(text_size = 3) +
  scale_fill_manual(values = c("1" = "#762a83", "2" = "#1b7837"), na.value = "white") +
  scale_color_manual(values = c("1" = "#762a83", "2" = "#1b7837")) +
  facet_wrap(~genotype, nrow = 3) +
  labs(y = "Flowering Time", x = "Stem Height", color = "Trial", fill = "Recombinant") +
  theme_bw()


# ftime - inflorescence height
ggplot(data = new.ftime_pcr_f3_nona.1recomu.i, mapping = aes(x = inflorescence_height, y = flowering_time, colour = exp)) +
  geom_point(aes(color = exp), fill = "white", shape = 21, size = 2, alpha = 0.5) +
  sm_statCorr(text_size = 3) +
  scale_color_manual(values = c("1" = "#762a83", "2" = "#1b7837")) +  
  facet_wrap(~genotype, nrow = 3) +
  labs(y = "Flowering Time", x = "Inflorescence Height", color = "Trial") +
  theme_bw()

ggplot(data = new.ftime_pcr_f3_nona.1recomu.i, mapping = aes(x = inflorescence_height, y = flowering_time, colour = exp)) +
  geom_point(aes(fill = as.factor(ifelse(new.ftime_pcr_f3_nona.1recomu.i$F3_genotype == "H", exp, NA))),
             shape = 21, size = 2, alpha = ifelse(new.ftime_pcr_f3_nona.1recomu.i$F3_genotype == "H", 1, 0.5), stroke = 0.5) +
  sm_statCorr(text_size = 3) +
  scale_fill_manual(values = c("1" = "#762a83", "2" = "#1b7837"), na.value = "white") +
  scale_color_manual(values = c("1" = "#762a83", "2" = "#1b7837")) +
  facet_wrap(~genotype, nrow = 3) +
  labs(y = "Flowering Time", x = "Inflorescence Height", color = "Trial", fill = "Recombinant") +
  theme_bw()


# ftime - shoot number 
ggplot(data = new.ftime_pcr_f3_nona.1recomu.i, mapping = aes(x = shoot_number, y = flowering_time, colour = exp)) +
  geom_point(aes(color = exp), fill = "white", shape = 21, size = 2, alpha = 0.5) +
  sm_statCorr(text_size = 3) +
  scale_color_manual(values = c("1" = "#762a83", "2" = "#1b7837")) +  
  facet_wrap(~genotype, nrow = 3) +
  labs(y = "Flowering Time", x = "Number of Shoots", color = "Trial") +
  theme_bw()

ggplot(data = new.ftime_pcr_f3_nona.1recomu.i, mapping = aes(x = shoot_number, y = flowering_time, colour = exp)) +
  geom_point(aes(fill = as.factor(ifelse(new.ftime_pcr_f3_nona.1recomu.i$F3_genotype == "H", exp, NA))),
             shape = 21, size = 2, alpha = ifelse(new.ftime_pcr_f3_nona.1recomu.i$F3_genotype == "H", 1, 0.5), stroke = 0.5) +
  sm_statCorr(text_size = 3) +
  scale_fill_manual(values = c("1" = "#762a83", "2" = "#1b7837"), na.value = "white") +
  scale_color_manual(values = c("1" = "#762a83", "2" = "#1b7837")) +
  facet_wrap(~genotype, nrow = 3) +
  labs(y = "Flowering Time", x = "Number of Shoots", color = "Trial", fill = "Recombinant") +
  theme_bw()


# ftime - internode length
# ggplot(data = new.ftime_pcr_f3_nona.1recomu.i, mapping = aes(x = leaves_distance_cm, y = flowering_time, colour = exp)) +
#   geom_point(aes(color = exp), fill = "white", shape = 21, size = 2, alpha = 0.5) +
#   sm_statCorr(text_size = 3) +
#   scale_color_manual(values = c("1" = "#762a83", "2" = "#1b7837")) +  
#   facet_wrap(~genotype, nrow = 3) +
#   labs(y = "Flowering Time", x = "Internode Length", color = "Trial") +
#   theme_bw()
# 
# ggplot(data = new.ftime_pcr_f3_nona.1recomu.i, mapping = aes(x = leaves_distance_cm, y = flowering_time, colour = exp)) +
#   geom_point(aes(fill = as.factor(ifelse(new.ftime_pcr_f3_nona.1recomu.i$F3_genotype == "H", exp, NA))),
#              shape = 21, size = 2, alpha = ifelse(new.ftime_pcr_f3_nona.1recomu.i$F3_genotype == "H", 1, 0.5), stroke = 0.5) +
#   sm_statCorr(text_size = 3) +
#   scale_fill_manual(values = c("1" = "#762a83", "2" = "#1b7837"), na.value = "white") +
#   scale_color_manual(values = c("1" = "#762a83", "2" = "#1b7837")) +
#   facet_wrap(~genotype, nrow = 3) +
#   labs(y = "Flowering Time", x = "Internode Length", color = "Trial", fill = "Recombinant") +
#   theme_bw()


# ftime - stem leaf density 
# ggplot(data = new.ftime_pcr_f3_nona.1recomu.i, mapping = aes(x = stem_leaf_density, y = flowering_time, colour = exp)) +
#   geom_point(aes(color = exp), fill = "white", shape = 21, size = 2, alpha = 0.5) +
#   sm_statCorr(text_size = 3) +
#   scale_color_manual(values = c("1" = "#762a83", "2" = "#1b7837")) +  
#   facet_wrap(~genotype, nrow = 3) +
#   labs(y = "Flowering Time", x = "Stem Leaf Density", color = "Trial") +
#   theme_bw()
# 
# ggplot(data = new.ftime_pcr_f3_nona.1recomu.i, mapping = aes(x = stem_leaf_density, y = flowering_time, colour = exp)) +
#   geom_point(aes(fill = as.factor(ifelse(new.ftime_pcr_f3_nona.1recomu.i$F3_genotype == "H", exp, NA))),
#              shape = 21, size = 2, alpha = ifelse(new.ftime_pcr_f3_nona.1recomu.i$F3_genotype == "H", 1, 0.5), stroke = 0.5) +
#   sm_statCorr(text_size = 3) +
#   scale_fill_manual(values = c("1" = "#762a83", "2" = "#1b7837"), na.value = "white") +
#   scale_color_manual(values = c("1" = "#762a83", "2" = "#1b7837")) +
#   facet_wrap(~genotype, nrow = 3) +
#   labs(y = "Flowering Time", x = "Stem Leaf Density", color = "Trial", fill = "Recombinant") +
#   theme_bw()


# ftime - rd
ggplot(data = new.ftime_pcr_f3_nona.1recomu.i, mapping = aes(x = rd_cm, y = flowering_time, colour = exp)) +
  geom_point(aes(color = exp), fill = "white", shape = 21, size = 2, alpha = 0.5) +
  sm_statCorr(text_size = 3) +
  scale_color_manual(values = c("1" = "#762a83", "2" = "#1b7837")) +  
  facet_wrap(~genotype, nrow = 3) +
  labs(y = "Flowering Time", x = "Rosette Diameter", color = "Trial") +
  theme_bw()

ggplot(data = new.ftime_pcr_f3_nona.1recomu.i, mapping = aes(x = rd_cm, y = flowering_time, colour = exp)) +
  geom_point(aes(fill = as.factor(ifelse(new.ftime_pcr_f3_nona.1recomu.i$F3_genotype == "H", exp, NA))),
             shape = 21, size = 2, alpha = ifelse(new.ftime_pcr_f3_nona.1recomu.i$F3_genotype == "H", 1, 0.5), stroke = 0.5) +
  sm_statCorr(text_size = 3) +
  scale_fill_manual(values = c("1" = "#762a83", "2" = "#1b7837"), na.value = "white") +
  scale_color_manual(values = c("1" = "#762a83", "2" = "#1b7837")) +
  facet_wrap(~genotype, nrow = 3) +
  labs(y = "Flowering Time", x = "Rosette Diameter", color = "Trial", fill = "Recombinant") +
  theme_bw()


# 7. Save data---- 
save.image("./new_analysis.RData")