####READ ME####
#FILE NAME: Community_data.R     
#DATE: last update 17/03/2024
#TITLE: Oxygen and pH fluxes in shallow bay habitats: Evaluating the effectiveness of a macroalgal forest restoration
#AUTHORS: Cristina Galobart, Cèlia Sitjà, Sònia de Caralt, Jorge Santamaría, Alba Vergés, Jordi Boada, Emma Cebrian
#SCRIPT: C Galobart (cgalobart@ceab.csic.es; cgalobart@gmail.com)
#JOURNAL: Journal of Phycology - under review

# DISCLAMER: This script has been developed by an ecologist, not a programmer, 
# so please take into account that the code may have room to be optimised. 
# Positive feedback will always be more than welcome.

# Script Content
# 1. Load data and data arrangements
# 2. Macroalgal biomass comparison
# 3. Macroinvertebrate biomass comparison
# 4. Macroinvertebrate number of individuals comparison
# 5. Macroalgal vs. macroinvertebrate biomass proportion (Figure 3A)
# 6. nMDS macroalgae (Figure 3B)
# 7. nMDS macroinvertebrates (Figure 3C)


#Load libraries
library(tidyverse)
library(vegan)
library(ggplot2)
library(car)
library(ggplot2)
library(ggrepel)


####1. Load data and data arrangements####
#Set working directory 
setwd("C:/Community_data")
options(scipen=999)

#Import metadata
metadata <- read.csv('metadata.csv', sep = ";", dec = ".")
#Reorder levels
metadata$assemblage <- factor(metadata$assemblage, levels = c("degraded", "forest", "restored_forest"))

#Import summary data of macroalgae biomass, macroinvertebrate biomass, macroivertebrate number of individuals, etc
data <- read.csv('summary_data.csv', sep = ";", dec = ".", stringsAsFactors = T) 
str(data)

####2. Macroalgal biomass comparison####
#Mean values
biomass_sum <- data %>%
  group_by(assemblage) %>%
    summarise(biomass = mean(algae_biomass_m2, na.rm = TRUE),
              sd = sd(algae_biomass_m2, na.rm = TRUE))
biomass_sum

#One-way ANOVA
mod <- aov(algae_biomass_m2 ~ assemblage, data = data)
summary(mod)

#Anova assumptions
#Normality and Homocedasticity
qqPlot(residuals(mod))
plot(fitted(mod), residuals(mod), xlab="Fitted Values", ylab="Residuals")
abline(h=0, lty=2)
lines(lowess(fitted(mod), residuals(mod)))

shapiro.test(mod$residuals) # formal normality statistical test
leveneTest(algae_biomass_m2 ~ assemblage, data = data) #formal statistical test for homogeinity of variances

#Tukey test for pair-wise comparison
TukeyHSD(mod)

####3. Macroinvertebrate biomass comparison####
#Remove the incubations with no data (macroinvertebrate were sorted in three incubations per assemblage, as indicated in the manuscript)
data_invert <- na.omit(data)
rownames(data_invert) <- NULL

#Mean values
biomass_sum <- data_invert %>%
  group_by(assemblage) %>%
  summarise(biomass = mean(invert_biomass_m2, na.rm = TRUE),
    sd = sd(invert_biomass_m2, na.rm = TRUE))
biomass_sum

#One-way ANOVA
mod <- aov(invert_biomass_m2 ~ assemblage, data = data_invert)
summary(mod)

#Anova assumptions
#Normality and Homocedasticity
qqPlot(residuals(mod))
plot(fitted(mod), residuals(mod), xlab="Fitted Values", ylab="Residuals")
abline(h=0, lty=2)
lines(lowess(fitted(mod), residuals(mod)))

shapiro.test(mod$residuals) # formal normality statistical test
leveneTest(invert_biomass_m2 ~ assemblage, data = data_invert) #formal statistical test for homogeinity of variance

#Tukey test for pair-wise comparison
TukeyHSD(mod)

####4. Macroinvertebrate number of individuals comparison####
#Mean values
indiv_sum <- data_invert %>%
  group_by(assemblage) %>%
  summarise(ind = mean(invert_ind_m2_log, na.rm = TRUE),
            sd = sd(invert_ind_m2_log, na.rm = TRUE))
indiv_sum

#One-way ANOVA
mod <- aov(invert_ind_m2_log ~ assemblage, data = data_invert)
summary(mod)

#Anova assumptions
#Normality and Homocedasticity
qqPlot(residuals(mod))
plot(fitted(mod), residuals(mod), xlab="Fitted Values", ylab="Residuals")
abline(h=0, lty=2)
lines(lowess(fitted(mod), residuals(mod)))

shapiro.test(mod$residuals) #formal normality statistical test
leveneTest(invert_ind_m2_log ~ assemblage, data = data_invert) #formal statistical test for homogeinity of variance

#Tukey test for pair-wise comparison
TukeyHSD(mod)

####5. Macroalgal vs. macroinvertebrate biomass####
#Mean values
proportion_sum <- data_invert %>%
  group_by(assemblage) %>%
  summarise(b_alg = mean(algae_biomass_m2, na.rm = TRUE),
            sd_alg = sd(algae_biomass_m2, na.rm = TRUE),
            b_invert = mean(invert_biomass_m2, na.rm = TRUE),
            sd_invert = sd(invert_biomass_m2, na.rm = TRUE))
proportion_sum

#plot
ggplot(proportion_sum, aes(x = assemblage)) +
  geom_bar(aes(y = b_alg, fill = "alga"), stat = "identity", position = "stack", width = 0.5) +
  geom_bar(aes(y = b_invert, fill = "invert"), stat = "identity", position = "stack", width = 0.5) +
  geom_errorbar(aes(ymin = b_alg - sd_alg, ymax = b_alg + sd_alg), width = 0.1, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = b_invert - sd_invert, ymax = b_invert + sd_invert), width = 0.1, position = position_dodge(width = 0.5)) +
  theme_classic() +
  scale_fill_manual(values = c(alga = "#727A54", invert = "#E15F30"))  # Specify fill colors

ggsave("Figure_3A.pdf", width=12, height=15, units= "cm") 


#Remove all objects in the Environment tab
####7. NMDs macroalgae####

#Import data and data arrangements
options(scipen=999)
metadata <- read.csv('metadata.csv', sep = ";", dec = ".")

#Reorder the levels
metadata$assemblage <- factor(metadata$assemblage, levels = c("degraded", "forest", "restored_forest"))

#Macroalgae data
abund_raw <- read.csv('biomass_algae.csv', sep = ";", dec = ".")

abund <- abund_raw %>% 
  remove_rownames %>% 
  column_to_rownames(var="id")

#Perform NMDs
mds_data <- metaMDS(abund, trymax = 300, k=2) 

#Save information
mds_points <- as.data.frame(mds_data$points)
mds_points <- tibble::rownames_to_column(mds_points, "id")
mds_data$stress
stress <- mds_data$stress

nMDS_info <- merge(metadata, mds_points, by = "id")

#Calculate the ellipses for the plot
plot.new()
ord <- ordiellipse(mds_data, as.factor(nMDS_info$assemblage), display="sites", kind="sd", conf=0.95, label=T)

#Create dataframe with the data for the ellipses
veganCovEllipse<-function(cov, center=c(0,0), scale=1, npoints=100)
{
  theta<-(0:npoints)*2*pi/npoints
  Circle<-cbind(cos(theta), sin(theta))
  t(center+scale*t(Circle%*%chol(cov)))
  
}
#Generate the ellipse points
ellipse <- data.frame()
for (g in levels(nMDS_info$assemblage)){
  if(g!="" && (g %in% names(ord))){
    ellipse<-rbind(ellipse, cbind(as.data.frame(with(nMDS_info[nMDS_info$assemblage==g,],
                                                     veganCovEllipse(ord[[g]]$cov, ord[[g]]$center, ord[[g]]$scale))), assemblage=g))
  }
}

#nMDS
ggplot(nMDS_info,aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color = assemblage), size = 3) + 
  geom_label(label = 0.08, x = 0.9, y = 1.1) +
  scale_color_manual(values = c("#EAA54A", "#6B896E", "#A3C580")) +
  theme_light() +
  coord_fixed() +
  geom_path(data = ellipse, aes(x=NMDS1, y=NMDS2, colour=assemblage), size=0.8, linetype=1) +
  xlim (-2.5, 1.5) + 
  ylim (-1, 1.20)

ggsave("Figure_3B.pdf", width=22, height=15, units= "cm")


#Remove all objects in the Environment tab
####8. nMDS macroinvertebrates####

#Import data and data arrangements
metadata <- read.csv('metadata.csv', sep = ";", dec = ".")

#Reorder levels
metadata$assemblage <- factor(metadata$assemblage, levels = c("degraded", "forest", "restored_forest"))

#Macroinvertebrate data (already ln (x+1) transformed)
abund <- read.csv('indiv_macroinvertebrates.csv', sep = ";", dec = ".", stringsAsFactors = T) 
str(abund)

abund <- abund %>% 
  remove_rownames %>% 
  column_to_rownames(var="id")

#Perform NMDs
mds_data <- metaMDS(abund, trymax = 300, k=2) 

#Save information
mds_points <- as.data.frame(mds_data$points)
mds_points <- tibble::rownames_to_column(mds_points, "id")
mds_data$stress
stress <- mds_data$stress

nMDS_info <- merge(metadata, mds_points, by = "id")

#Calculate the ellipses for the plot
plot.new()
ord <- ordiellipse(mds_data, as.factor(nMDS_info$assemblage), display="sites", kind="sd", conf=0.95, label=T)

#Create dataframe with the data for the ellipses
veganCovEllipse<-function(cov, center=c(0,0), scale=1, npoints=100)
{
  theta<-(0:npoints)*2*pi/npoints
  Circle<-cbind(cos(theta), sin(theta))
  t(center+scale*t(Circle%*%chol(cov)))
  
}

#Generate the ellipse points
ellipse <- data.frame()
for (g in levels(nMDS_info$assemblage)){
  if(g!="" && (g %in% names(ord))){
    ellipse<-rbind(ellipse, cbind(as.data.frame(with(nMDS_info[nMDS_info$assemblage==g,],
                    veganCovEllipse(ord[[g]]$cov, ord[[g]]$center, ord[[g]]$scale))), assemblage=g))
    }
}

#nMDS plot
ggplot(nMDS_info,aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color = assemblage), size = 3) + 
  geom_label(label = 0.01, x = 0.5, y = 0.7) +
  scale_color_manual(values = c("#EAA54A", "#6B896E", "#A3C580")) +
  theme_light() +
  coord_fixed() +
  geom_path(data = ellipse, aes(x=NMDS1, y=NMDS2, colour=assemblage), size=0.8, linetype=1) +
  xlim (-1.3, 1.3) +
  ylim (-0.85, 0.85)

ggsave("Figure_3C.pdf", width=22, height=15, units= "cm")

