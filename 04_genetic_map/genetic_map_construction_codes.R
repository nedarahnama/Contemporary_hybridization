# < Genetic architecture of phenotypic differences between endangered hybridizing Arabis floodplain species >
# < PhD Thesis >
# < Neda Rahnamae> 

# < Chapter 2 - Genetic architecture of phenotypic differences between endangered hybridizing Arabis floodplain species >
# < Genetic Map Construction >
# Feb 2024
# 742 individuals - 5630 markers



########################## 
# 1. Load R libraries----
library(ASMap) # 1.0-8
library(dplyr) # 1.1.4
library(ggplot2) # 3.5.1
library(gridExtra) # 2.3
library(qtl) # 1.70
library(readr) # 2.1.5
library(tibble) # 3.2.1
library(tidyr) # 1.3.1
library(tidyverse) # 2.0.0


########################## 
# 2. Functions----
## check marker distribution on physical maps after linkage group separation
## help assign LG to chromosomes
# marker name should be like this: IWB65373_1A_3.8 = marker name, physical chromosome and Mb, separated by "_"
check_chrom = function(maps) {
  aa = strsplit(rownames(maps), "_")
  cc = sapply(aa, function(x)
    x[2]) # chrom names
  maps$chrom = cc
  # xtabs(~chr+chrom,maps)
  write.table(xtabs( ~ chr + chrom, maps), "chromosome_check.txt", sep =
                "\t")
}

## check genetic orientation compared to physical order
# marker name should be like this: IWB65373_1A_3.8 = marker name, physical chromosome and Mb, separated by "_"
check_orientation = function(maps) {
  aa = strsplit(rownames(maps), "_")
  bb = sapply(aa, function(x)
    as.numeric(x[3]))
  cc = sapply(aa, function(x)
    x[2]) # chrom names
  maps$Mb = bb
  print(xyplot(Mb ~ pos |
                 chr, data = maps[maps$chr == cc, ], as.table = TRUE)) # check orientation and make sure to use markers on the same chrom
  write.table(maps,
              "orientation_check.txt",
              col.names = NA,
              quote = F)
}


## function to reorder a map manually
reorder.marker = function(cross, chr, new.marker.order) {# new.marker.order can be either numbers or marker names
  genodata = cross$geno[[chr]]$data
  mapdata = cross$geno[[chr]]$map
  new.geno = genodata[,new.marker.order]
  new.map = mapdata[new.marker.order]
  cross$geno[[chr]]$data = new.geno
  cross$geno[[chr]]$map = new.map
  cross2 = quickEst(cross, chr)
  print(pull.map(cross2, chr, as.table = T))
  #new.pos = est.map(cross, chr, map.function="kosambi")[[chr]]
  #print(new.pos)
  #ll = length(new.pos) # n marker
  #cross$geno[[chr]]$map[1:ll] = new.pos[1:ll]
  return(cross2)
}




########################## 
# 3. Load data as cross object----
cross781 <- read.cross(format = "csvsr",
                       dir = "../input_data",
                       genfile = "13_geno781.csv",
                       phefile = "13_pheno781.csv",
                       na.strings = c("./.","NA"), #./. is NA in genotype file, NA is NA in phenotype file
                       genotypes = c("0/0","0/1","1/1"),  
                       alleles = c("N","S"), 
                       BC.gen = 0, F.gen = 2,
                       estimate.map = F,
                       sep = ",")


########################## 
# 4. Analyze the cross object----
## 1. Calculate observed het across markers
ch1_8 <-
  geno.table(cross781,
             chr = c(1, 2, 3, 4, 5, 6, 7, 8),
             scanone.output = F)   #for chr 1 to 8
write.csv(ch1_8,
          "output/cross781/cross781_gen_distribution.csv")

## 2. Look at summary statistics
summary(cross781)
# BC(0)F(2) cross
# No. individuals:    781 
# No. phenotypes:     52 
# Percent phenotyped: 56.3 
# No. chromosomes:    8 
# Autosomes:          1 2 3 4 5 6 7 8 
# Total markers:      5360 
# No. markers:        732 548 697 859 599 499 588 838 
# Percent genotyped:  73.6 
# Genotypes (%):      NN:22.0  NS:46.9  SS:31.1  not SS:0.0  not NN:0.0 

## 3. Plot preliminary results
plotMissing(cross781)
geno.image(
  cross781,
  reorder = TRUE,
  col = c("#DBDBD7", "#FC8D62", "#8DA0CB", "#66C2A5")) #missing, NN, NS, SS


########################## 
# 5. Pre-construction filtering for markers----
## 1. Missing values: genotypes----
# identify the genotypes with a certain number of missing values
sg_miss <- statGen(cross781,
              bychr = FALSE,
              stat.type = "miss",
              id = "ID")
write.csv(sg_miss,
          "output/cross781/cross781_missing_genotypes_per_indv.csv")

summary(sg_miss$miss)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 91     694    1203    1413    1830    5098 
hist(sg_miss$miss)
hist(sg_miss$miss[sg_miss$miss < 3500])
sg_miss$miss[sg_miss$miss > 3500] 
# 142  240  965   79  294 1073  693  987  771  168  299  790  773  173  248  403  943  438 1082   61  966  237  711   68  842  275 
# 3739 3587 4225 4225 5030 3743 3739 4337 4585 3505 3519 4216 4860 4045 4815 3815 4139 4752 3791 3503 3800 3945 4121 3793 3678 3829 
# 277  969 1087  298  714  341  357  695  308  311 1095 
# 4258 3615 3819 3976 3663 4777 3527 3513 5098 4026 3734  
sum(sg_miss$miss > 3500) # 37
sum(sg_miss$miss < 3500) # 744
# only keep individuals with missing data < 3500 -> 781 - 37 = 744
map1 <- subsetCross(cross781, ind = sg_miss$miss < 3500)
summary(map1)
# BC(0)F(2) cross
# No. individuals:    744 
# No. phenotypes:     52 
# Percent phenotyped: 56.4 
# No. chromosomes:    8 
# Autosomes:          1 2 3 4 5 6 7 8 
# Total markers:      5360 
# No. markers:        732 548 697 859 599 499 588 838 
# Percent genotyped:  76.1 
# Genotypes (%):      NN:22.0  NS:46.9  SS:31.0  not SS:0.0  not NN:0.0  

plotMissing(map1)
geno.image(
  map1,
  reorder = TRUE,
  col = c("#DBDBD7", "#FC8D62", "#8DA0CB", "#66C2A5")) #missing, NN, NS, SS

## 2. Genetic clones----
# check highly related individuals
gc <- genClones(map1, tol = 0.95, id = "ID") 
# There are no genotype pairs with matching allele proportions greater than 0.95.

## 3. Check null and duplicates markers----
map1 <- drop.nullmarkers(map1)
totmar(map1) # 5360
dupmar <- findDupMarkers(map1, exact.only=FALSE, adjacent.only=FALSE) # 46
dupmar.adjonly <- findDupMarkers(map1, adjacent.only = TRUE) # 22
dupmar.exconly <- findDupMarkers(map1, exact.only = TRUE) # 23
dupmar.nexact <-
  findDupMarkers(map1, exact.only = FALSE, adjacent.only = TRUE) # 38
# one might consider dropping the extra markers
totmar(map1) # 5360 markers
map2 <- drop.markers(map1, unlist(dupmar, dupmar.exconly, dupmar.adjonly)) # -47
summary(map2)
# BC(0)F(2) cross
# No. individuals:    744 
# No. phenotypes:     52 
# Percent phenotyped: 56.4 
# No. chromosomes:    8 
# Autosomes:          1 2 3 4 5 6 7 8 
# Total markers:      5313 
# No. markers:        731 544 688 848 593 493 584 832 
# Percent genotyped:  76.1 
# Genotypes (%):      NN:22.0  NS:47.0  SS:31.0  not SS:0.0  not NN:0.0 
plotMissing(map2)
geno.image(
  map2,
  reorder = TRUE,
  col = c("#DBDBD7", "#FC8D62", "#8DA0CB", "#66C2A5")) #missing, NN, NS, SS

## 4. Marker profiles----
profileMark(
  map2,
  stat.type = c("seg.dist", "lod", "erf", "miss"),
  #segregation distortion, allelic proportion, missing value proportion of markers
  crit.val = "bonf",
  #chr = "8",
  layout = c(1, 4),
  type = "l",
  #"l" gives you info as line,"p" as data points
  cex = 0.5
)

sg_dis_map2 <-
  profileMark(
    map2,
    stat.type = c("seg.dist"),
    crit.val = "bonf",
    #chr = "8",
    layout = c(1, 1),
    type = "l",
    #"l" gives you info as line,"p" as data points
    cex = 0.5
  )
df.sg_dis_map2 <- as.data.frame(sg_dis_map2)
write.csv(df.sg_dis_map2, "./output/df.sg_dis_map2.csv", row.names = T)

## 5. Check individuals' genotype frequencies----
g <- pull.geno(map2)
gfreq <- apply(g, 1, function(a)
  table(factor(a, levels = 1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))

par(mfrow = c(1, 3), las = 1)
for (i in 1:3)
  plot(
    gfreq[i,],
    ylab = "Genotype frequency",
    main = c("NN", "SN", "SS")[i],
    ylim = c(0, 1)
  )
par(mfrow = c(1, 1))

## 6. Check segregation distortion----
# we can use bonferroni criteria = chi squre pvalue / number of markers < 0.05 to remove distorted markers or chr specific thresholds
# profileMark(map2, stat.type = c("seg.dist", "lod", "erf", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l", cex = 0.5)
totmar(map2) # 5313
gt <- geno.table(map2)
nrow(gt[gt$P.value < 0.5 / totmar(map2),]) # number of distorted markers with bonferroni criteria 3762
sgd <- (gt[gt$P.value < 0.05 / totmar(map2),])

## 7. Pull out 3 types of markers before making genetic maps----
### 1. colocated markers: markers with identical genotying data----
map3 <- pullCross(map2, type = "co.located")
#No markers were found with these characteristics.
totmar(map3) # 5313
nmar(map3)
# 1   2   3   4   5   6   7   8 
# 731 544 688 848 593 493 584 832

### 2. markers with too many missing data----
#I only keep markers with missing data <10%, we can push back markers with more missing data to the map after making the genetic maps
map3 <-
  pullCross(map3, type = "missing", pars = list(miss.thresh = 0.235)) # using high stringency here
totmar(map3) # 2164
nmar(map3)
summary(map3)
# Total markers:      2164 
# No. markers:        310 241 247 391 206 165 240 364 
# Percent genotyped:  78.6 
# Genotypes (%):      NN:20.9  NS:47.8  SS:31.4  not SS:0.0  not NN:0.0 
plotMissing(map3)
geno.image(
  map3,
  reorder = T ,
  col = c("#DBDBD7", "#FC8D62", "#8DA0CB", "#66C2A5")
)
profileMark(
  map3,
  stat.type = c("seg.dist", "miss"),
  #segregation distortion, allelic proportion, missing value proportion of markers
  crit.val = "bonf",
  layout = c(1, 2),
  type = "l",
  #"l" gives you info as line,"p" as data points
  cex = 0.5
)

profileMark(
  map3,
  stat.type = c("seg.dist", "prop"),
  #segregation distortion, allelic proportion, missing value proportion of markers
  crit.val = "bonf",
  layout = c(1, 6),
  type = "l",
  #"l" gives you info as line,"p" as data points
  cex = 0.5
)

### 3. Filter out markers that do not fall into some ranges----
#(not good enough proportion of alleles)
totmar(map3) # 2164
nmar(map3)
# 1   2   3   4   5   6   7   8 
# 310 241 247 391 206 165 240 364 

## Filter markers according to their segregation distortion
# map4 <- pullCross(map3, type = "seg.distortion", 
#                   pars = list(seg.ratio = "10:80:0.3"))  
# totmar(map4) # 2274
# nmar(map3)
# nmar(map4)

# Plot ration along the genome 
# profileMark(
#   map4,
#   stat.type = c("seg.dist", "prop"),
#   #segregation distortion, allelic proportion, missing value proportion of markers
#   crit.val = "bonf",
#   layout = c(1, 6),
#   #chr = "1",
#   type = "l",
#   #"l" gives you info as line,"p" as data points
#   cex = 0.5
# )
# distortion <- map4$seg.distortion$table

## Filter markers according to their allele proportion
nmar(map3)
# 1   2   3   4   5   6   7   8 
# 310 241 247 391 206 165 240 364 

mm_ns <- statMark(map3, stat.type = "marker")$marker$NS
map4 <- drop.markers(map3, c(markernames(map3)[mm_ns > 0.75],   
                               markernames(map3)[mm_ns < 0.25])) 
totmar(map4) # 2096
summary(map4)
# Total markers:      2096 
# No. markers:        304 232 244 369 205 161 227 354 
# Percent genotyped:  78.6 
# Genotypes (%):      NN:20.9  NS:47.1  SS:32.0  not SS:0.0  not NN:0.0 

profileMark(
  map4,
  stat.type = c("seg.dist", "prop"),
  #segregation distortion, allelic proportion, missing value proportion of markers
  crit.val = "bonf",
  layout = c(1, 6),
  type = "l",
  #"l" gives you info as line,"p" as data points
  cex = 0.5
)
plotMissing(map4)
geno.image(
  map4,
  reorder = T ,
  col = c("#DBDBD7", "#FC8D62", "#8DA0CB", "#66C2A5")
)

sg_dis_map4 <-
  profileMark(
    map4,
    stat.type = c("seg.dist"),
    crit.val = "bonf",
    #chr = "8",
    layout = c(1, 1),
    type = "l",
    #"l" gives you info as line,"p" as data points
    cex = 0.5,
    display.markers = F
  )
df.sg_dis_map4 <- as.data.frame(sg_dis_map4)
write.csv(df.sg_dis_map4, "./output/df.sg_dis_map4.csv", row.names = T)


## 8. How many markers have been removed?----
names(map4)
sum(
  ncol(map4$missing$data)
)                  # 3149 markers have been removed due to missing
totmar(map4) # 2096
nind(map4) # 744
## 9. Write cross (to rename marker name column)----
write.cross(
  map4,
  format = "csvsr",
  filestem = "cross_data/map4/map4_744_2096",
  chr = c("1", "2", "3", "4", "5", "6", "7", "8"),
  digits = NULL
)


########################## 
# 6. Map construction----
## 1. Threshold marker clustering----
totmar(map4) # 2096
nind(map4)   # 744
pValue(
  dist = seq(5, 40, by = 5),
  pop.size = 1:744,
  map.function = "kosambi",
  LOD = F
)
summary(map4)
# Total markers:      2096 
# No. markers:        304 232 244 369 205 161 227 354 
# Percent genotyped:  78.6 
# Genotypes (%):      NN:20.9  NS:47.1  SS:32.0  not SS:0.0  not NN:0.0 

## 2. Load corrected data as cross object----
# '-': 0, 'NN': 1, 'NS': 2, 'SS': 3
cor_map4 <- read.cross(format = "csvsr",
                       dir = "./cross_data",
                       genfile = "map4_744_2096_cor_genotypes_2mb_0.5mb.csv",
                       phefile = "map4_744_2096_phe.csv",
                       na.strings = c("0","NA"), #- is NA in genotype file, NA is NA in phenotype file
                       genotypes = c("1","2","3"),  
                       alleles = c("N","S"), 
                       BC.gen = 0, F.gen = 2,
                       estimate.map = F,
                       sep = ",")

summary(cor_map4)
# BC(0)F(2) cross
# No. individuals:    744 
# No. phenotypes:     52 
# Percent phenotyped: 96 
# No. chromosomes:    8 
# Autosomes:          1 2 3 4 5 6 7 8 
# Total markers:      2096 
# No. markers:        304 232 244 369 205 161 227 354 
# Percent genotyped:  99.6 
# Genotypes (%):      NN:19.6  NS:48.8  SS:31.6  not SS:0.0  not NN:0.0 

## 3. Plot preliminary results----
plotMissing(cor_map4)
geno.image(
  cor_map4,
  reorder = TRUE,
  col = c("#DBDBD7", "#FC8D62", "#8DA0CB", "#66C2A5")) #missing, NN, NS, SS
plotMap(cor_map4)

sg_dis_cor_map4 <-
  profileMark(
    cor_map4,
    stat.type = c("seg.dist"),
    crit.val = "bonf",
    #chr = "8",
    layout = c(1, 1),
    type = "l",
    #"l" gives you info as line,"p" as data points
    cex = 0.5
  )
df.sg_dis_cor_map4 <- as.data.frame(sg_dis_cor_map4)
write.csv(df.sg_dis_cor_map4, "./output/df.sg_dis_cor_map4.csv", row.names = T)

## 4. First map construction----
map5 <-
  mstmap(
    cor_map4,
    id = "ID",
    bychr = T,
    trace = TRUE,
    suffix = "numeric",
    dist.fun = "kosambi",
    objective.fun = "COUNT",
    p.value = 1e-50,
    noMap.dist = 30,
    noMap.size = 30,
    miss.thresh = 1,
    mvest.bc = T,
    detectBadData	= T,
    return.imputed	= T
  )

nchr(map5)
totmar(map5)
nind(map5)
nmar(cor_map4)
nmar(map5)
# 1   2   3   4   5 6.1 6.2 7.1 7.2 7.3   8 
# 304 232 244 369 205 160   1 221   1   5 354
plotMap(map5)
maps = pull.map(map5, as.table=T)
check_orientation(maps)


### 1. Check orientation----
physical_distance <- read.csv("./cross_data/map4/map4_744_2096_gen_phy_dis.csv",header=T, dec = ".",na.strings = "NA")

#### Chr 1
genetic_distance1 <- map5$geno$'1'$map
gd1 <- data.frame(genetic_distance1)
gd1 <- rownames_to_column(gd1, "ID")
names(gd1)[names(gd1)=="genetic_distance1"] <- "genetic_distance"
marker_info1 <- inner_join(physical_distance, gd1, by = "ID")
or_c1_map5 <-
  ggplot(marker_info1, aes(POS,genetic_distance))+
  geom_point()+
  geom_smooth(se = F)+
  labs(title = "Chromosome 1",
       x = "Physical distance [bp]",
       y = "Genetic distance [cM]")

#### Chr 2
genetic_distance2 <- map5$geno$'2'$map
gd2 <- data.frame(genetic_distance2)
gd2 <- rownames_to_column(gd2, "ID")
names(gd2)[names(gd2)=="genetic_distance2"] <- "genetic_distance"
marker_info2 <- inner_join(physical_distance, gd2, by = "ID")
or_c2_map5 <-
  ggplot(marker_info2, aes(POS,genetic_distance))+
  geom_point()+
  geom_smooth(se = F)+
  labs(title = "Chromosome 2",
       x = "Physical distance [bp]",
       y = "Genetic distance [cM]")

#### Chr 3
genetic_distance3 <- map5$geno$'3'$map
gd3 <- data.frame(genetic_distance3)
gd3 <- rownames_to_column(gd3, "ID")
names(gd3)[names(gd3)=="genetic_distance3"] <- "genetic_distance"
marker_info3 <- inner_join(physical_distance, gd3, by = "ID")
or_c3_map5 <-
  ggplot(marker_info3, aes(POS,genetic_distance))+
  geom_point()+
  geom_smooth(se = F)+
  labs(title = "Chromosome 3",
       x = "Physical distance [bp]",
       y = "Genetic distance [cM]")

#### Chr 4
genetic_distance4 <- map5$geno$'4'$map
gd4 <- data.frame(genetic_distance4)
gd4 <- rownames_to_column(gd4, "ID")
names(gd4)[names(gd4)=="genetic_distance4"] <- "genetic_distance"
marker_info4 <- inner_join(physical_distance, gd4, by = "ID")
or_c4_map5 <-
  ggplot(marker_info4, aes(POS,genetic_distance))+
  geom_point()+
  geom_smooth(se = F)+
  labs(title = "Chromosome 4",
       x = "Physical distance [bp]",
       y = "Genetic distance [cM]")

#### Chr 5
genetic_distance5 <- map5$geno$'5'$map
gd5 <- data.frame(genetic_distance5)
gd5 <- rownames_to_column(gd5, "ID")
names(gd5)[names(gd5)=="genetic_distance5"] <- "genetic_distance"
marker_info5 <- inner_join(physical_distance, gd5, by = "ID")
or_c5_map5 <-
  ggplot(marker_info5, aes(POS,genetic_distance))+
  geom_point()+
  geom_smooth(se = F)+
  labs(title = "Chromosome 5",
       x = "Physical distance [bp]",
       y = "Genetic distance [cM]")

#### Chr 6
genetic_distance6 <- map5$geno$'6.1'$map
gd6 <- data.frame(genetic_distance6)
gd6 <- rownames_to_column(gd6, "ID")
names(gd6)[names(gd6)=="genetic_distance6"] <- "genetic_distance"
marker_info6 <- inner_join(physical_distance, gd6, by = "ID")
or_c6_map5 <-
  ggplot(marker_info6, aes(POS,genetic_distance))+
  geom_point()+
  geom_smooth(se = F)+
  labs(title = "Chromosome 6",
       x = "Physical distance [bp]",
       y = "Genetic distance [cM]")

#### Chr 7
genetic_distance7 <- map5$geno$'7.1'$map
gd7 <- data.frame(genetic_distance7)
gd7 <- rownames_to_column(gd7, "ID")
names(gd7)[names(gd7)=="genetic_distance7"] <- "genetic_distance"
marker_info7 <- inner_join(physical_distance, gd7, by = "ID")
or_c7_map5 <-
  ggplot(marker_info7, aes(POS,genetic_distance))+
  geom_point()+
  geom_smooth(se = F)+
  labs(title = "Chromosome 7",
       x = "Physical distance [bp]",
       y = "Genetic distance [cM]")

#### Chr 8
genetic_distance8 <- map5$geno$'8'$map
gd8 <- data.frame(genetic_distance8)
gd8 <- rownames_to_column(gd8, "ID")
names(gd8)[names(gd8)=="genetic_distance8"] <- "genetic_distance"
marker_info8 <- inner_join(physical_distance, gd8, by = "ID")
or_c8_map5 <-
  ggplot(marker_info8, aes(POS,genetic_distance))+
  geom_point()+
  geom_smooth(se = F)+
  labs(title = "Chromosome 8",
       x = "Physical distance [bp]",
       y = "Genetic distance [cM]")

#### Genome
grid.arrange(or_c1_map5, or_c2_map5, or_c3_map5, or_c4_map5, or_c5_map5, or_c6_map5, or_c7_map5, or_c8_map5, ncol=4)
par(mfrow = c(1, 1))

### 2. Merge, rename and flip chr 1, 2, 3, 4, 5, 6, 7 & 8 and check orientation again----
# merge
nmar(map5)
map5_1 <- mergeCross(map5, merge = list (
  "unlink7" = c("7.2","7.3")
))
nmar(map5_1)
nind(map5_1)

# rename chr
names(map5_1$geno) = c("1","2","3","4","5","6","unlink6","7","8","unlink7")
nmar(map5_1)
maps = pull.map(map5_1, chr=(nmar(map5_1)>6), as.table=T)
write.table(maps,"maps_map5_1_for_checking.txt",sep="\t", col.names=NA,quote=F)
check_chrom(maps) # check output file "chromosome_check.txt" to see whether all chromosomes have been sepated and which LG is which chromosome, so we can regroup them together.
check_orientation(maps)

# flip chr order
map5_2 = flip.order(map5_1, c("1","2","3","4","5","6","7","8"))
nmar(map5_2)
nind(map5_2)
maps = pull.map(map5_2, chr=(nmar(map5_2)>6), as.table=T)
write.table(maps,"maps_map5_2_for_checking.txt",sep="\t", col.names=NA,quote=F)
check_chrom(maps) # check output file "chromosome_check.txt" to see whether all chromosomes have been sepated and which LG is which chromosome, so we can regroup them together.
check_orientation(maps)

# For checking order of markers
totmar(map5_2) # 2096
nmar(map5_2)
write.cross(
  map5_2,
  format = "csvsr",
  filestem = "cross_data/map5_2/map5_2_744_2096",
  chr = c("1", "2", "3", "4", "5", "6", "7", "8"),
  digits = NULL
)


### 3. Re-order the inverted regions of chr 3 and 6----
# reorder markers on chr 3
nmar(map5_2)
map5_3 = reorder.marker(map5_2, "3", c(1:78, 244:79))
maps = pull.map(map5_3, chr=(nmar(map5_3)>6), as.table=T)
write.table(maps,"maps_map5_3_for_checking.txt",sep="\t", col.names=NA,quote=F)
check_chrom(maps) # check output file "chromosome_check.txt" to see whether all chromosomes have been sepated and which LG is which chromosome, so we can regroup them together.
check_orientation(maps)
nmar(map5_3)

# reorder markers on chr 4
map5_4 = reorder.marker(map5_3, "4", c(1:72, 139:73, 140:369))
nmar(map5_4)
maps = pull.map(map5_4, chr=(nmar(map5_4)>6), as.table=T)
write.table(maps,"maps_map5_4_for_checking.txt",sep="\t", col.names=NA,quote=F)
check_chrom(maps) # check output file "chromosome_check.txt" to see whether all chromosomes have been sepated and which LG is which chromosome, so we can regroup them together.
check_orientation(maps)
nmar(map5_4)

# reorder markers on chr 6
map5_5 = reorder.marker(map5_4, "6", c(1:50, 160:51))
nmar(map5_5)
maps = pull.map(map5_5, chr=(nmar(map5_5)>6), as.table=T)
write.table(maps,"maps_map5_5_for_checking.txt",sep="\t", col.names=NA,quote=F)
check_chrom(maps) # check output file "chromosome_check.txt" to see whether all chromosomes have been sepated and which LG is which chromosome, so we can regroup them together.
check_orientation(maps)

# reorder markers on chr 7
map5_6 = reorder.marker(map5_5, "7", c(1:79, 153:80, 154:221))
nmar(map5_6)
maps = pull.map(map5_6, chr=(nmar(map5_6)>6), as.table=T)
write.table(maps,"maps_map5_6_for_checking.txt",sep="\t", col.names=NA,quote=F)
check_chrom(maps) # check output file "chromosome_check.txt" to see whether all chromosomes have been sepated and which LG is which chromosome, so we can regroup them together.
check_orientation(maps)


### 4. Drop markers----
# markers with overlap with repeat regions:
#                     ID sequence position
# 1294 chr5_5_16.352767     chr5 16352767
# 1295 chr5_5_16.352789     chr5 16352789
# 1296 chr5_5_16.352808     chr5 16352808
map5_7 <- drop.markers(map5_6, c("chr5_5_16.352767","chr5_5_16.352789","chr5_5_16.352808",
                                 "chr7_7_9.43433","chr7_7_9.43433",
                                 "chr8_8_38.085842","chr8_8_38.705893","chr8_8_39.77297"))
maps = pull.map(map5_7, chr=(nmar(map5_7)>6), as.table=T)
check_orientation(maps)
plotMap(map5_7, chr = (nmar(map5_7)>6))
nmar(map5_7)
totmar(map5_7) 
maps = pull.map(map5_7, chr=(nmar(map5_7)>6), as.table=T)
write.table(maps,"maps_map5_7_for_checking.txt",sep="\t", col.names=NA,quote=F)
check_chrom(maps) # check output file "chromosome_check.txt" to see whether all chromosomes have been sepated and which LG is which chromosome, so we can regroup them together.
check_orientation(maps)

sg_dis_map5_7 <-
  profileMark(
    map5_7,
    stat.type = c("seg.dist"),
    crit.val = "bonf",
    chr = c("1","2","3","4","5","6","7","8"),
    layout = c(1, 1),
    type = "l",
    #"l" gives you info as line,"p" as data points
    cex = 0.5
  )
df.sg_dis_map5_7 <- as.data.frame(sg_dis_map5_7)
write.csv(df.sg_dis_map5_7, "./output/df.sg_dis_map5_7.csv", row.names = T)

### 4. Write cross----
totmar(map5_7) # 2089 - 7 = 2082
nmar(map5_7)
# 1       2       3       4       5       6 unlink6       7       8 unlink7 
# 304     232     244     369     202     160       1     220     351       6 
nind(map5_7) # 744
# Remove individuals with wrong information
ind2drop <- c("10101010","69696969")
map5_7 <- subset(map5_7, ind= -match(ind2drop, map5_7$pheno$ID))
nind(map5_7) # 742

write.cross(
  map5_7,
  format = "csvsr",
  filestem = "cross_data/map5_7/map5_7_742_2082",
  chr = c("1", "2", "3", "4", "5", "6", "7", "8"),
  digits = NULL
)


geno.image(
  cor_map4,
  reorder = TRUE,
  col = c("#DBDBD7", "#FC8D62", "#8DA0CB", "#66C2A5")) #missing, NN, NS, SS


########################## 
save.image("./genetic_map_final.RData")