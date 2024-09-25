####READ ME####
#FILE NAME: Metabolic_fluxes.R     
#DATE: last update 17/03/2024
#TITLE: Oxygen and pH fluxes in shallow bay habitats: Evaluating the effectiveness of a macroalgal forest restoration
#AUTHORS: Cristina Galobart, Cèlia Sitjà, Sònia de Caralt, Jorge Santamaría, Alba Vergés, Jordi Boada, Emma Cebrian
#SCRIPT: C Galobart (cgalobart@ceab.csic.es; cgalobart@gmail.com)
#JOURNAL: Journal of Phycology  - under review

# DISCLAMER: This script has been developed by an ecologist, not a programmer, 
# so please take into account that the code may have room to be optimised. 
# Positive feedback will always be more than welcome.

# Script Content
# 1. Load data and data arrangements
# 2. Oxygen flux (Figure 4)
# 3. pH variation at the start and end of incubations (Figure 5)
 

#Load libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(tidyr)
require(gridExtra)
library(car)


####1. Load data and data arrangements####
#Set working directory 
setwd("C:/Metabolic_fluxes")
options(scipen=999)

data <- read.csv2("data_O2.csv", sep = ";", dec = ".", stringsAsFactors = T)
data$replicate <- as.factor(data$replicate)
str(data)

#Reorder the levels
data$assemblage <- factor(data$assemblage, levels = c("Degraded", "Forest", "Restored_forest"))

#Transformation of oxygen values to oxygen fluxes
data$do_diff <- (data$do_t1 - data$do_t0) #difference between t1 and t0 (mg O2/L)
data$do_diff_ch <- (data$do_diff * 17.5) #difference t1 and t0 in the chamber (17.5L volume) (mg O2)
data$do_diff_ch_A <- (data$do_diff_ch / 0.138544) #considering the area of the chamber (0.138544) (mg O2/m2)
data$do_flux <- (data$do_diff_ch_A * 60 / data$inc_time) #considering the minutes that have been incubating (mg O2/m2 h)
data$do_flux_mmol <- (data$do_flux / 32) #oxygen to moles (mmol O2/m2 h)

####2. Oxygen flow (Figure 4)####
#Mean values
df.summary <- data %>%
  group_by(assemblage, treatment) %>%
  summarise(do_flux_mmol_m = mean(do_flux_mmol),
            sd = sd(do_flux_mmol, na.rm = TRUE))
df.summary

data %>%
  group_by(treatment) %>%
  summarise(do_flux_mmol_m = mean(do_flux_mmol),
    sd = sd(do_flux_mmol, na.rm = TRUE))

#Figure 4A - Community Net Primary Production (NPP) and Community Respiration (R)
do <- ggplot(df.summary, aes(fill=treatment, y=do_flux_mmol_m, x=assemblage)) + 
  geom_col(position="stack", width = 0.6) +
  geom_hline(yintercept = 0) +
  theme(legend.position = "top") +
  geom_errorbar(aes(ymin = do_flux_mmol_m-sd, ymax = do_flux_mmol_m+sd), width = 0.1) +
  geom_point(data = data, aes(y = do_flux_mmol, x=assemblage, colour = treatment)) +
  theme_classic() + 
  ylim (-6 , 8) + 
  ylab("DO flux (mmol O2/m2 h)") +
  scale_x_discrete(limits = c("Degraded", "Forest", "Restored_forest")) +
  scale_fill_manual(values=c("#ac9d93", '#ffd737'))+
  scale_colour_manual(values = c("#58504b","#86711d")) 

print(do)

ggsave("Figure_4A.pdf", do, width=15, height=20, units= "cm")

#Test differences between NPP among assemblages
#Subset of light incubations
data_test_light <- subset(data, treatment=="light")

#One-way ANOVA
mod <- aov(do_flux_mmol ~ assemblage, data = data_test_light)
summary(mod)

#Anova assumptions
#Normality and Homocedasticity
qqPlot(residuals(mod))
plot(fitted(mod), residuals(mod), xlab="Fitted Values", ylab="Residuals")
abline(h=0, lty=2)
lines(lowess(fitted(mod), residuals(mod)))

shapiro.test(mod$residuals) #formal normality statistical test 
leveneTest(do_flux_mmol ~ assemblage, data = data_test_light) #formal statistical test for homogeneity of variance

#Tukey test for pair-wise comparison
TukeyHSD(mod)

#Test differences between R_O2 among assemblages
#Subset of dark incubations
data_test_dark <- subset(data, treatment=="dark")

#One way ANOVA
mod <- aov(do_flux_mmol ~ assemblage, data = data_test_dark)
summary(mod)

#Anova assumptions
#Normality and Homocedasticity
qqPlot(residuals(mod))
plot(fitted(mod), residuals(mod), xlab="Fitted Values", ylab="Residuals")
abline(h=0, lty=2)
lines(lowess(fitted(mod), residuals(mod)))

shapiro.test(mod$residuals) #formal normality statistical test 
leveneTest(do_flux_mmol ~ assemblage, data = data_test_dark) #formal statistical test for homogeneity of variance

#Figure 4B - Community Gross Primary Production (GPP)
#This is NPP + |R| 
data_light_GPP <-  data %>% 
  subset(treatment == "light") %>% 
  select(id, assemblage, replicate, treatment, do_flux_mmol)

data_dark_GPP <- data %>% 
  subset(treatment == "dark") %>% 
  select(id, assemblage, replicate, treatment, do_flux_mmol)

#Merge both databases by "id" and rearrange the database
GPP_data <- merge(data_light_GPP, data_dark_GPP, by = "id")

GPP_data <- GPP_data %>% 
  select(id, assemblage.x, replicate.x, do_flux_mmol.x, do_flux_mmol.y)

colnames(GPP_data)[c(1, 2, 3, 4, 5)] <- c("id", "assemblage", "replicate", "do_flux_mmol_light", "do_flux_mmol_dark")

#Compute GPP_O2
GPP_data$GPP <- abs(GPP_data$do_flux_mmol_light) + abs(GPP_data$do_flux_mmol_dark)

#Mean values
GPP_data %>%
  group_by(assemblage) %>%
  summarise(GPP_m = mean(GPP),
    sd = sd(GPP, na.rm = TRUE))

#Figure 3B - Community Gross Primary Production (GPP_O2)
do_gpp <- ggplot(GPP_data, aes(x = assemblage, y = GPP)) +
  geom_jitter(width = 0, pch = 17, cex = 8, color = "#0B445F") +
  theme_classic() +
  stat_summary(fun = mean, pch = 17, size = 3, color = "#1791C9") +
  scale_x_discrete(limits = c("Degraded", "Forest", "Restored_forest")) +
  ylim (0 , 12)

print(do_gpp)
#Warning message, 3 rows removed are because they overlap with the mean triangles

ggsave("Figure_4B.pdf", do_gpp, width=10, height=20, units= "cm")

#One-way ANOVA
mod <- aov(GPP ~ assemblage, data = GPP_data)
summary(mod)

#Anova assumptions
#Normality and Homocedasticity
qqPlot(residuals(mod))
plot(fitted(mod), residuals(mod), xlab="Fitted Values", ylab="Residuals")
abline(h=0, lty=2)
lines(lowess(fitted(mod), residuals(mod)))

shapiro.test(mod$residuals) #formal normality statistical test 
leveneTest(GPP ~ assemblage, data = GPP_data) #formal statistical test for homogeneity of variance

#Tukey test for pair-wise comparison
TukeyHSD(mod)

#Remove all objects in the Environment tab
#### 3. pH variation at the start and end of incubations (Figure 5) ####

data_pH <- read.csv2("data_pH.csv", sep = ";", dec = ".", stringsAsFactors = T)
data_pH$replicate <- as.factor(data_pH$replicate)
data_pH$time <- as.factor(data_pH$time)
str(data_pH)

data_pH$assemblage <- factor(data_pH$assemblage, levels= c("Degraded", "Restored_forest", "Forest"))
data_pH$treatment <- factor(data_pH$treatment, levels= c("light", "dark"))
data_pH$time <- factor(data_pH$time, levels= c("0", "1"))


#Figure 4 - Change in pH at the start and end of the incubations
ph_point <- ggplot(data_pH, aes(x = time, y = pH)) + 
  geom_point(colour = "#486b69", size = 3)+ 
  geom_path(aes(group = id), colour = "#486b69", size = 1.2) +
  facet_wrap(~ treatment + assemblage) +
  ylim (8.1 , 8.5) +
  theme_bw()+
  theme(axis.text.y = element_text(size = 15)) 

print(ph_point)
ggsave("Figure_5.pdf", ph_point, width=27, height=18, units= "cm")

#Divide the dataset for statistical analyses
#Repeated measures ANOVA to compate pH between the start and end of incubations

#Light & degraded
L_D <- data_pH %>% 
  subset(treatment == "light") %>% 
  subset(assemblage == "Degraded")

#Repeated measures ANOVA
mod <- aov(pH ~ time + Error(id/time), data = L_D)
summary(mod)
#Normality test
shapiro.test(mod$`id:time`$residuals)
leveneTest(pH ~ time, data = L_D) 

#Light & forest
L_F <- data_pH %>% 
  subset(treatment == "light") %>% 
  subset(assemblage == "Forest")

#Repeated measures ANOVA
mod <- aov (pH~time+Error(id/time), data = L_F)
summary(mod)
#Normality test
shapiro.test(mod$`id:time`$residuals)
leveneTest(pH ~ time, data = L_F) 

#Light & restored forest
L_RF <- data_pH %>% 
  subset(treatment == "light") %>% 
  subset(assemblage == "Restored_forest")

#Repeated measures ANOVA
mod <- aov (pH~time+Error(id/time), data = L_RF)
summary(mod)
#Normality test
shapiro.test(mod$`id:time`$residuals)
leveneTest(pH ~ time, data = L_RF) 

#Dark & degraded
D_D <- data_pH %>% 
  subset(treatment == "dark") %>% 
  subset(assemblage == "Degraded")

#Repeated measures ANOVA
mod <- aov (pH~time+Error(id/time), data = D_D)
summary(mod)
#Normality test
shapiro.test(mod$`id:time`$residuals)
leveneTest(pH ~ time, data = D_D) 

#Dark & forest
D_F <- data_pH %>% 
  subset(treatment == "dark") %>% 
  subset(assemblage == "Forest")

#Repeated measures ANOVA
mod <- aov (pH~time+Error(id/time), data = D_F)
summary(mod)
#Normality test
shapiro.test(mod$`id:time`$residuals)
leveneTest(pH ~ time, data = D_F) 

#Dark & restored forest
D_RF <- data_pH %>% 
  subset(treatment == "dark") %>% 
  subset(assemblage == "Restored_forest")

#Repeated measures ANOVA
mod <- aov (pH~time+Error(id/time), data = D_RF)
summary(mod)
#Normality test
shapiro.test(mod$`id:time`$residuals)
leveneTest(pH ~ time, data = D_RF) 

