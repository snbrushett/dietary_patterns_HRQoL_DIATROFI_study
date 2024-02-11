###     SCRIPT: DIATROFI FFQ DIETARY PATTERNS VIA PCA VARIMAX ROTATION
###     AUTHOR(S): SIOBHAN BRUSHETT 
###     DESCRIPTION: VARIMAC PCA AND ROBUST PCA ON FFQ ITEMS FROM THE BASELINE 2015-2016 AND 2017-2018 DIATROFI DATA
###     PROJECT: MERGED DIATROFI 2015-2018
###     NOTE(S): PCA: https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome/blob/master/Projects/IBD_onset_Diet/05_PCA.R
###              PAPERS: https://doi.org/10.1093/ecco-jcc/jjab219 (Genetics Department, Groningen); https://doi.org/10.1093/eurpub/ckaa178 (Prolepsis)
###              TESTS: https://www.statology.org/bartletts-test-of-sphericity/; https://www.statisticshowto.com/kaiser-meyer-olkin/


#libraries
library(tidyverse)
library(DataExplorer)
library(psych)
library(corrplot)#version 0.92
library(GPArotation)
library(rrcov) #robust PCA
library(pracma) #robust PCA
library(factoextra) #PCA HC
library(FactoMineR) #PCA HC
library(matrixStats) #cluster analysis
library(dendextend) #cluster analysis
library(circlize) #cluster analysis
library(RColorBrewer)
library(vegan) #version 2.6-2; Bray Curtis distance matrix
 
##Contents of this file

## 0. LOAD DATA
## 1. CLEAN FFQ ITEMS BY EXCLUSION CRITERIA
## 2. DERIVE FFQ WEIGHTED FREQUENCIES
## 3. TESTS TO VALIDATE PC FACTOR ANALYSIS
## 4. PCA USING VARIMAX ROTATION

        #### =============== 0. LOAD DATA   ================ ####

dia_merged <- read.delim("diatrofi_merged_2015_2018_baseline.txt", sep = "\t", header = T, stringsAsFactors = T)
#6583  117
row.names(dia_merged) <- dia_merged$pseudo_ids

        #### =============== 1. CLEAN FFQ ITEMS BY EXCLUSION CRITERIA ================ ####

#subset ffq items
dia_merged_ffq <- dia_merged[,c("pseudo_ids",
                                "child_ffq_bread_white",
                                "child_ffq_bread_wholegrain",
                                "child_ffq_bread_bagel_burger",
                                "child_ffq_bread_crumbs",
                                "child_ffq_cereal_products",
                                "child_ffq_rice_white",
                                "child_ffq_rice_wholegrain",
                                "child_ffq_pasta",
                                "child_ffq_pasta_wholegrain",
                                "child_ffq_potatoes",
                                "child_ffq_french_fries",
                                "child_ffq_tomato_cucumber_carrot_pepper",
                                "child_ffq_lettuce_cabbage_spinach_rocket",
                                "child_ffq_broccoli_pumpkin_cauliflower",
                                "child_ffq_leeks_spinach_celery",
                                "child_ffq_fruit_orange",
                                "child_ffq_fruit_apple_pear",
                                "child_ffq_fruit_winter_other",
                                "child_ffq_fruit_banana",
                                "child_ffq_fruit_summer_other",
                                "child_ffq_fruit_dried",
                                "child_ffq_lentils_beans_chickpeas",
                                "child_ffq_olives",
                                "child_ffq_nuts",
                                "child_ffq_milk_yoghurt_full_fat",
                                "child_ffq_milk_yoghurt_low_fat",
                                "child_ffq_chocolate_milk_drink",
                                "child_ffq_yellow_or_cream_cheese",
                                "child_ffq_feta_anthotyros_cheese",
                                "child_ffq_cottage_light_cheese_low_fat",
                                "child_ffq_egg",
                                "child_ffq_beef",
                                "child_ffq_burger_meatballs_minced_meat",
                                "child_ffq_chicken_turkey",
                                "child_ffq_pork",
                                "child_ffq_pork_souvlaki_pita",
                                "child_ffq_chicken_souvlaki_pita",
                                "child_ffq_hamburger",
                                "child_ffq_hotdog",
                                "child_ffq_lamb_goat_ribs",
                                "child_ffq_meat_cold_cuts",
                                "child_ffq_meat_cold_cuts_low_no_fat",
                                "child_ffq_sausage_bacon",
                                "child_ffq_fish_small",
                                "child_ffq_fish_big",
                                "child_ffq_seafood",
                                "child_ffq_veg_risotto_stuff_veg",
                                "child_ffq_mousaka",
                                "child_ffq_grean_beans_peas_okra",
                                "child_ffq_homemade_pies",
                                "child_ffq_standard_pies",
                                "child_ffq_toast_sandwich",
                                "child_ffq_pizza",
                                "child_ffq_pastries",
                                "child_ffq_dessert",
                                "child_ffq_jelly_composte",
                                "child_ffq_cakes",
                                "child_ffq_chocolate_milk",
                                "child_ffq_chocolate_dark",
                                "child_ffq_honey_marmalade_sugar",
                                "child_ffq_ice_cream_milkshake",
                                "child_ffq_crisps_popcorn",
                                "child_ffq_coffee",
                                "child_ffq_fruit_juice_fresh",
                                "child_ffq_fruit_juice_standard",
                                "child_ffq_herbal_teas",
                                "child_ffq_soft_drinks",
                                "child_ffq_soft_drinks_light",
                                "child_ffq_energy_drinks",
                                "child_ffq_mayo",
                                "child_ffq_mayo_light",
                                "child_ffq_olive_oil",
                                "child_ffq_seed_oil",
                                "child_ffq_margarine",
                                "child_ffq_butter")]
#6583 76
rownames(dia_merged_ffq) <- dia_merged_ffq$pseudo_ids
dia_merged_ffq$pseudo_ids <- NULL

#based on following study: https://doi.org/10.1186/s12966-016-0353-2 
dia_merged_ffq_cri <- dia_merged_ffq

#criteria 1: participants with >10% of the FFQ data missing are considered invalid and participants were omitted from analysis
dia_merged_ffq_cri$count_na <- rowSums(is.na(dia_merged_ffq_cri))
sum(dia_merged_ffq_cri$count_na==0)
#3251 (49.4%) participants answered every item of the FFQ
range(dia_merged_ffq_cri$count_na) #0-75
sum(dia_merged_ffq_cri$count_na>7.5) #participants with more than 10% of the items missing
#499 (7.6%) of participants did not complete >10% of the questionnaire

dia_merged_ffq_cri$count_na[dia_merged_ffq_cri$count_na>7.5] <- NA

for(i in 1:ncol(dia_merged_ffq_cri)) {
        dia_merged_ffq_cri[,i][is.na(dia_merged_ffq_cri$count_na)] <- NA
}

#criteria 2: all other missings (<10% missing) are considered as not consumed
dia_merged_ffq_cri <- dia_merged_ffq_cri %>%
        filter(!is.na(count_na))
#6084 76 (as expected)

dia_merged_ffq_cri$count_na <- NULL

summary(dia_merged_ffq_cri)
for(i in 1:ncol(dia_merged_ffq_cri)) {
        print(table(dia_merged_ffq_cri[,i], useNA = "ifany"))
        dia_merged_ffq_cri[,i][is.na(dia_merged_ffq_cri[,i])] <- "0_never_less_than_once_per_month"
        print(table(dia_merged_ffq_cri[,i], useNA = "ifany"))
}

summary(dia_merged_ffq_cri)

#included a 3rd criteria, where if >90% of the items indicated as '0_never_less_than_once_per_month', then omit these participants
#given that all NA's are accounted for, make all answers of "0_never_less_than_once_per_month" equivalent to NA, count the NAs and remove those >90%
dia_merged_ffq_cri_add <- dia_merged_ffq_cri
colnames(dia_merged_ffq_cri_add)
head(rownames(dia_merged_ffq_cri_add))
summary(dia_merged_ffq_cri_add) #QC to confirm that there are no NAs in the dataset

for (i in 1:ncol(dia_merged_ffq_cri_add)) {
        dia_merged_ffq_cri_add[,i][dia_merged_ffq_cri_add[,i]=="0_never_less_than_once_per_month"] <- NA
}

dia_merged_ffq_cri_add$count_na <- rowSums(is.na(dia_merged_ffq_cri_add))
range(dia_merged_ffq_cri_add$count_na) #0-75
sum(dia_merged_ffq_cri_add$count_na>67.5) #participants with more than 90% of the items indicated recorded as '0_never_less_than_once_per_month'
#9 (0.15%) of participants with more than 90% of the items indicated as a frequency of 0.01

dia_merged_ffq_cri_add <- dia_merged_ffq_cri_add %>% filter (
        count_na<67.5
)
#6075  76

summary(dia_merged_ffq_cri_add)

for(i in 1:ncol(dia_merged_ffq_cri_add)) {
        dia_merged_ffq_cri_add[,i][is.na(dia_merged_ffq_cri_add[i])] <- '0_never_less_than_once_per_month'
}
summary(dia_merged_ffq_cri_add)

dia_merged_ffq_cri_add$count_na <- NULL

dia_merged_ffq <- dia_merged_ffq_cri_add

#4th criteria:
#retrospectively added a 4th criteria where if all participants indicate consuming every item in the FFQ '7_times_per_day_more_than_2',
#these participants were omitted

colnames(dia_merged_ffq)
head(rownames(dia_merged_ffq))
summary(dia_merged_ffq) #QC to confirm that there are no NAs in the dataset

for (i in 1:ncol(dia_merged_ffq)) {
        dia_merged_ffq[,i][dia_merged_ffq[,i]=='7_times_per_day_more_than_2'] <- NA
}

dia_merged_ffq$count_na <- rowSums(is.na(dia_merged_ffq))
range(dia_merged_ffq$count_na) #0-69
sum(dia_merged_ffq$count_na>67.5) #participants with more than 90% of the items indicated recorded as '0_never_less_than_once_per_month'
#1 participants with more than 90% of the items indicated with a frequency of 2

dia_merged_ffq <- dia_merged_ffq %>% filter (
        count_na<67.5
)
#6074  76

summary(dia_merged_ffq)

for(i in 1:ncol(dia_merged_ffq)) {
        dia_merged_ffq[,i][is.na(dia_merged_ffq[i])] <- "7_times_per_day_more_than_2"
        
}
summary(dia_merged_ffq)

dia_merged_ffq$count_na <- NULL

#clean RStudio
rm(dia_merged_ffq_cri_add, dia_merged_ffq_cri)

#for steps below:
#dia_merged_ffq_bray_curtis <- dia_merged_ffq

        #### =============== 2. DERIVE FFQ WEIGHTED FREQUENCIES   ================ ####

for (i in 1:ncol(dia_merged_ffq)){
        print(table(dia_merged_ffq[,i], useNA = "ifany"))
        dia_merged_ffq[,i] <- recode_factor(dia_merged_ffq[,i], 
                                            "0_never_less_than_once_per_month"="0.01",
                                            "1_times_per_month_1_3"="0.07",
                                            "2_times_per_week_1"="0.14",
                                            "3_times_per_week_2"="0.29",
                                            "4_times_per_week_3_4"="0.5",
                                            "5_times_per_week_5_6"="0.79",
                                            "6_times_per_day_1"="1",
                                            "7_times_per_day_more_than_2"="2")
        dia_merged_ffq[,i] <- as.numeric(as.character(dia_merged_ffq[,i]))
        print(table(dia_merged_ffq[,i], useNA = "ifany"))
}

setwd("C:/Users/Siobhan Brushett/OneDrive - UMCG/Prolepsis_internship/merged_2015_2018/dietary_patterns_ffq/")
pdf(file = "2023_02_20_ffq_items_criteria1_4.pdf", useDingbats = F, onefile = T, width = 20, height=12)
dia_merged_ffq %>% plot_histogram(nrow = 4L, ncol=4L)
dev.off()

        #### =============== 3. TESTS TO VALIDATE PC FACTOR ANALYSIS   ================ ####

cor_matrix = cor(dia_merged_ffq)
print(round(dia_merged_ffq, 3))

#Bartlett’s Test of Sphericity
t_bart_1 <- psych::cortest.bartlett(cor_matrix, n=nrow(dia_merged_ffq))
#Barlett results:
t_bart_1
#$chisq
#[1] 152546.7

#$p.value
#[1] 0

#$df
#[1] 2775
#Conclusion: the test statistic is high and the p-value is == 0 which is less than 0.05; thus the data is suitable for PCA/factor analysis

t_kmo <- psych::KMO(cor_matrix)
t_kmo
summary(t_kmo$MSAi)
#the minimum value is 0.8399 indicating that the variables are adequate for PCA/factor analysis

        #### =============== 4. PCA USING VARIMAX ROTATION   ================ ####
setwd("C:/Users/Siobhan Brushett/OneDrive - UMCG/Prolepsis_internship/merged_2015_2018/dietary_patterns_ffq/pca_varimax_robust/")
#https://pablobernabeu.github.io/2018/naive-principal-component-analysis-in-r/

#pca visualizations
res.pca <- PCA(dia_merged_ffq,  graph = FALSE)
get_eig(res.pca)
eigen_values <- as.data.frame(get_eig(res.pca))
eigen_values$eigenvalue <- round(eigen_values$eigenvalue, digits = 1)
eigen_values
pdf(file = "2023_02_20_pca_scree_plot.pdf", useDingbats = F, onefile = T, width = 20, height=12)
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 20))
#1st component: 18.8%
#2nd component: 5.9%
#3rd component: 3.8%
#4th component: 2.9%
#5th component: 2.6%
#>> 33.9% variance explained
dev.off()


        #4.1. pca varimax rotation (default: weighted sums)
pca_varimax <- psych::principal(dia_merged_ffq, rotate="varimax", nfactors=5, scores=TRUE) #first 5 components
summary(pca_varimax)

pdf(file = "2023_02_16_diatrofi_ffq_pca_varimax_correlation.pdf", useDingbats = F, onefile = T, width = 20, height=12)
corrplot(pca_varimax$loadings,is.corr=F)
#fa.diagram(pca_varimax)
dev.off()

#pdf(file = "2023_02_16_diatrofi_ffq_pca_varimax_biplots.pdf", useDingbats = F, onefile = T, width = 20, height=12)
#biplot(princomp(dia_merged_ffq))
#psych::biplot.psych(pca_varimax)
#dev.off()

pc_loadings <- pca_varimax$loadings

ffq_pc <- as.matrix(dia_merged_ffq)
ffq_pc <- ffq_pc %*% pc_loadings

write.table(ffq_pc, "2023_02_15_ffq_pca_varimax_loadings.txt", sep="\t", row.names=T, quote = F)

#defining key items of loadings
varimax_item_loadings <- as.data.frame(loadings(pca_varimax)[])
summary(varimax_item_loadings)

varimax_item_loadings$items <- rownames(varimax_item_loadings)
RC1_loadings <- varimax_item_loadings %>%
        filter(RC1 > 0.45) %>%
        select (c(items, RC1)) %>%
        na.omit(RC1_loadings) %>%
        arrange(desc(RC1))
#19 items

RC2_loadings <- varimax_item_loadings %>% 
        filter(RC2 > 0.45) %>%
        select (c(items, RC2)) %>%
        na.omit(RC2_loadings) %>%
        arrange(desc(RC2))
#11 items

RC3_loadings <- varimax_item_loadings %>% 
        filter(RC3 > 0.45) %>%
        select (c(items, RC3)) %>%
        na.omit(RC3_loadings) %>%
        arrange(desc(RC3))
#4 items

RC4_loadings <- varimax_item_loadings %>% 
        filter(RC4 > 0.45) %>%
        select (c(items, RC4)) %>%
        na.omit(RC4_loadings) %>%
        arrange(desc(RC4))
#6 items

RC5_loadings <- varimax_item_loadings %>% 
        filter(RC5 > 0.45) %>%
        select (c(items, RC5)) %>%
        na.omit(RC5_loadings) %>%
        arrange(desc(RC5))
#4 items

#clean RStudio
rm(cor_matrix, pca_varimax, ffq_pc, varimax_item_loadings, RC1_loadings, RC2_loadings, RC3_loadings, RC4_loadings, RC5_loadings)

#generate loadings according to methodology of Prolepsis article: 

        #4.2. varimax rotation with regression method
pca_varimax_regression <- factanal(dia_merged_ffq, 5, scores= c("regression"), rotation= "varimax")

pdf(file = "2023_06_16_diatrofi_ffq_pca_varimax_reg_correlation.pdf", useDingbats = F, onefile = T, width = 20, height=12)
corrplot(pca_varimax_regression$loadings,is.corr=F, tl.col = "black")
#fa.diagram(pca_varimax_regression)
dev.off()

pc_loadings_regression <- pca_varimax_regression$loadings

ffq_pc_regression <- as.matrix(dia_merged_ffq)
ffq_pc_regression <- ffq_pc_regression %*% pc_loadings_regression

write.table(ffq_pc_regression, "2023_02_20_ffq_pca_varimax_regression_loadings.txt", sep="\t", row.names=T, quote = F)

#defining key items of loadings
regression_item_loadings <- as.data.frame(loadings(pca_varimax_regression)[])
summary(regression_item_loadings)

regression_item_loadings$items <- rownames(regression_item_loadings)
F1_loadings <- regression_item_loadings %>%
        filter(Factor1 > 0.4) %>%
        select (c(items, Factor1)) %>%
        na.omit(F1_loadings) %>%
        arrange(desc(Factor1))
#15 items

F2_loadings <- regression_item_loadings %>% 
        filter(Factor2 > 0.4) %>%
        select (c(items, Factor2)) %>%
        na.omit(F2_loadings) %>%
        arrange(desc(Factor2))
#7 items

F3_loadings <- regression_item_loadings %>% 
        filter(Factor3 > 0.4) %>%
        select (c(items, Factor3)) %>%
        na.omit(F3_loadings) %>%
        arrange(desc(Factor3))
#8 items

F4_loadings <- regression_item_loadings %>% 
        filter(Factor4 > 0.4) %>%
        select (c(items, Factor4)) %>%
        na.omit(F4_loadings) %>%
        arrange(desc(Factor4))
#6 items

F5_loadings <- regression_item_loadings %>% 
        filter(Factor5 > 0.4) %>%
        select (c(items, Factor5)) %>%
        na.omit(F5_loadings) %>%
        arrange(desc(Factor5))
#5 items

rm(pca_varimax_regression, ffq_pc_regression, regression_item_loadings, F1_loadings, F2_loadings, F3_loadings, F4_loadings, F5_loadings)
