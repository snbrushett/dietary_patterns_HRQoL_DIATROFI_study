###     SCRIPT: TABLE 2: QoL IN RELATION TO DIETARY FEATURES
###     AUTHOR(S): SIOBHAN BRUSHETT 
###     DESCRIPTION: TOTAL HEALTH RELATED QUALITY OF LIFE SCORE AND SCALE SCORES IN RELATION TO DIETARY INDEXES AND DIETARY PATTERNS
###     PROJECT: DIATROFI_BASELINE_2015_2018

###     LAST UPDATED: 16 MAY 2023

#libraries
library(tidyverse)

## 0. LOAD DATA
## 1. TABLE 2: QoL AND DIET

          #### =============== 0. LOAD DATA ================ ####
#MAC-book location
setwd("~/OneDrive - UMCG/Prolepsis_internship/merged_2015_2018/masterfile/2_tidy_masterfile_incl_marg_add_hrqol_percentiles/")
masterfile <- read.delim("2023_05_09_diatrofi_2015_2018_masterfile_updated.txt", sep = "\t", header = T, stringsAsFactors = T)
#6583 66
str(masterfile)
masterfile$child_demo_school_code <- as.factor(masterfile$child_demo_school_code)

        #### =============== 1. TABLE 1: DESCRIPTIVE STATISTICS   ================ ####
#NOTE: THIS SCRIPT IS EXACTLY THE SAME AS THE SCRIPT 'derive_and_tidy_masterfile.R' FOR THE RELEVANT SECTION, WITH THE 
#EXCEPTION THAT THE PSYCHO-SOCIAL FUNCTIONING SCALE WAS REMOVED AS IT IS HIGHLY CORRELATED WITH ALL OTHER QOL OUTCOMES
        ####1.1 Total summary statistics####
mean(masterfile$child_ffq_hpdi_score, na.rm=T)
sd(masterfile$child_ffq_hpdi_score, na.rm=T)

mean(masterfile$child_ffq_animal_score, na.rm=T)
sd(masterfile$child_ffq_animal_score, na.rm=T)

mean(masterfile$child_ffq_diet_quality_score, na.rm=T)
sd(masterfile$child_ffq_diet_quality_score, na.rm=T)

median(masterfile$child_ffq_pattern_1_meat_seafood_prepd_meals, na.rm=T)
IQR(masterfile$child_ffq_pattern_1_meat_seafood_prepd_meals, na.rm=T)

median(masterfile$child_ffq_pattern_2_cooked_veg_grains_legumes, na.rm=T)
IQR(masterfile$child_ffq_pattern_2_cooked_veg_grains_legumes, na.rm=T)

median(masterfile$child_ffq_pattern_3_fruits_raw_veg_cheese, na.rm=T)
IQR(masterfile$child_ffq_pattern_3_fruits_raw_veg_cheese, na.rm=T)

median(masterfile$child_ffq_pattern_4_confectioneries_pizza, na.rm=T)
IQR(masterfile$child_ffq_pattern_4_confectioneries_pizza, na.rm=T)

median(masterfile$child_ffq_pattern_5_starchy_foods_sweetened_bev, na.rm=T)
IQR(masterfile$child_ffq_pattern_5_starchy_foods_sweetened_bev, na.rm=T)

        ####1.2 Q15: HRQoL####
hrqol_numeric <- masterfile %>%
  filter(is.na(child_pedsql_health_related_qol_score_total_catg_Q15) == F) %>%
  group_by(child_pedsql_health_related_qol_score_total_catg_Q15) %>%
  select (child_ffq_hpdi_score, 
          child_ffq_animal_score, 
          child_ffq_diet_quality_score,
          child_ffq_pattern_1_meat_seafood_prepd_meals,
          child_ffq_pattern_2_cooked_veg_grains_legumes,
          child_ffq_pattern_3_fruits_raw_veg_cheese,
          child_ffq_pattern_4_confectioneries_pizza,
          child_ffq_pattern_5_starchy_foods_sweetened_bev) %>%
  summarise_all(list(mean = ~ mean(., na.rm=T), sd = ~ sd (., na.rm=T),
                     median = ~ median (., na.rm=T), 
                     #min = ~ min (., na.rm=T), max = ~ max (., na.rm=T),
                     IQR = ~ IQR (., na.rm = T)))

#HRQoL: Q15 p-values
x <- t.test(child_ffq_hpdi_score ~ child_pedsql_health_related_qol_score_total_catg_Q15, data=masterfile)
x$p.value

x <- t.test(child_ffq_animal_score ~ child_pedsql_health_related_qol_score_total_catg_Q15, data=masterfile)
x$p.value

x <- t.test(child_ffq_diet_quality_score ~ child_pedsql_health_related_qol_score_total_catg_Q15, data=masterfile)
x$p.value

wilcox.test(child_ffq_pattern_1_meat_seafood_prepd_meals ~ child_pedsql_health_related_qol_score_total_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_2_cooked_veg_grains_legumes ~ child_pedsql_health_related_qol_score_total_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_3_fruits_raw_veg_cheese ~ child_pedsql_health_related_qol_score_total_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_4_confectioneries_pizza ~ child_pedsql_health_related_qol_score_total_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_5_starchy_foods_sweetened_bev ~ child_pedsql_health_related_qol_score_total_catg_Q15, data=masterfile)

        ####1.3 Q15: Physical Functioning####
phy_numeric <- masterfile %>%
  filter(is.na(child_pedsql_phy_score_catg_Q15) == F) %>%
  group_by(child_pedsql_phy_score_catg_Q15) %>%
  select (child_ffq_hpdi_score, 
          child_ffq_animal_score, 
          child_ffq_diet_quality_score,
          child_ffq_pattern_1_meat_seafood_prepd_meals,
          child_ffq_pattern_2_cooked_veg_grains_legumes,
          child_ffq_pattern_3_fruits_raw_veg_cheese,
          child_ffq_pattern_4_confectioneries_pizza,
          child_ffq_pattern_5_starchy_foods_sweetened_bev) %>%
  summarise_all(list(mean = ~ mean(., na.rm=T), sd = ~ sd (., na.rm=T),
                     median = ~ median (., na.rm=T), 
                     #min = ~ min (., na.rm=T), max = ~ max (., na.rm=T),
                     IQR = ~ IQR (., na.rm = T)))

#phy: Q15 p-values
t.test(child_ffq_hpdi_score ~ child_pedsql_phy_score_catg_Q15, data=masterfile)
t.test(child_ffq_animal_score ~ child_pedsql_phy_score_catg_Q15, data=masterfile)
t.test(child_ffq_diet_quality_score ~ child_pedsql_phy_score_catg_Q15, data=masterfile)

wilcox.test(child_ffq_pattern_1_meat_seafood_prepd_meals ~ child_pedsql_phy_score_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_2_cooked_veg_grains_legumes ~ child_pedsql_phy_score_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_3_fruits_raw_veg_cheese ~ child_pedsql_phy_score_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_4_confectioneries_pizza ~ child_pedsql_phy_score_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_5_starchy_foods_sweetened_bev ~ child_pedsql_phy_score_catg_Q15, data=masterfile)

        ####1.4 Q15: Emotional Functioning####
#emo: Q15 summary stats
emo_numeric <- masterfile %>%
  filter(is.na(child_pedsql_emo_score_catg_Q15) == F) %>%
  group_by(child_pedsql_emo_score_catg_Q15) %>%
  select (child_ffq_hpdi_score, 
          child_ffq_animal_score, 
          child_ffq_diet_quality_score,
          child_ffq_pattern_1_meat_seafood_prepd_meals,
          child_ffq_pattern_2_cooked_veg_grains_legumes,
          child_ffq_pattern_3_fruits_raw_veg_cheese,
          child_ffq_pattern_4_confectioneries_pizza,
          child_ffq_pattern_5_starchy_foods_sweetened_bev) %>%
  summarise_all(list(mean = ~ mean(., na.rm=T), sd = ~ sd (., na.rm=T),
                     median = ~ median (., na.rm=T), 
                     #min = ~ min (., na.rm=T), max = ~ max (., na.rm=T),
                     IQR = ~ IQR (., na.rm = T)))

#emo: Q15 p-values
t.test(child_ffq_hpdi_score ~ child_pedsql_emo_score_catg_Q15, data=masterfile)
t.test(child_ffq_animal_score ~ child_pedsql_emo_score_catg_Q15, data=masterfile)
t.test(child_ffq_diet_quality_score ~ child_pedsql_emo_score_catg_Q15, data=masterfile)

wilcox.test(child_ffq_pattern_1_meat_seafood_prepd_meals ~ child_pedsql_emo_score_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_2_cooked_veg_grains_legumes ~ child_pedsql_emo_score_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_3_fruits_raw_veg_cheese ~ child_pedsql_emo_score_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_4_confectioneries_pizza ~ child_pedsql_emo_score_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_5_starchy_foods_sweetened_bev ~ child_pedsql_emo_score_catg_Q15, data=masterfile)

        ####1.5 Q15: Social Functioning####
#soc: Q15 summary stats
soc_numeric <- masterfile %>%
  filter(is.na(child_pedsql_soc_score_catg_Q15) == F) %>%
  group_by(child_pedsql_soc_score_catg_Q15) %>%
  select (child_ffq_hpdi_score, 
          child_ffq_animal_score, 
          child_ffq_diet_quality_score,
          child_ffq_pattern_1_meat_seafood_prepd_meals,
          child_ffq_pattern_2_cooked_veg_grains_legumes,
          child_ffq_pattern_3_fruits_raw_veg_cheese,
          child_ffq_pattern_4_confectioneries_pizza,
          child_ffq_pattern_5_starchy_foods_sweetened_bev) %>%
  summarise_all(list(mean = ~ mean(., na.rm=T), sd = ~ sd (., na.rm=T),
                     median = ~ median (., na.rm=T), 
                     #min = ~ min (., na.rm=T), max = ~ max (., na.rm=T),
                     IQR = ~ IQR (., na.rm = T)))

#soc: Q15 p-values
t.test(child_ffq_hpdi_score ~ child_pedsql_soc_score_catg_Q15, data=masterfile)
t.test(child_ffq_animal_score ~ child_pedsql_soc_score_catg_Q15, data=masterfile)
t.test(child_ffq_diet_quality_score ~ child_pedsql_soc_score_catg_Q15, data=masterfile)

wilcox.test(child_ffq_pattern_1_meat_seafood_prepd_meals ~ child_pedsql_soc_score_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_2_cooked_veg_grains_legumes ~ child_pedsql_soc_score_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_3_fruits_raw_veg_cheese ~ child_pedsql_soc_score_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_4_confectioneries_pizza ~ child_pedsql_soc_score_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_5_starchy_foods_sweetened_bev ~ child_pedsql_soc_score_catg_Q15, data=masterfile)

        ####1.6 Q15: School Functioning####
#sch: Q15 summary stats
sch_numeric <- masterfile %>%
  filter(is.na(child_pedsql_sch_score_catg_Q15) == F) %>%
  group_by(child_pedsql_sch_score_catg_Q15) %>%
  select (child_ffq_hpdi_score, 
          child_ffq_animal_score, 
          child_ffq_diet_quality_score,
          child_ffq_pattern_1_meat_seafood_prepd_meals,
          child_ffq_pattern_2_cooked_veg_grains_legumes,
          child_ffq_pattern_3_fruits_raw_veg_cheese,
          child_ffq_pattern_4_confectioneries_pizza,
          child_ffq_pattern_5_starchy_foods_sweetened_bev) %>%
  summarise_all(list(mean = ~ mean(., na.rm=T), sd = ~ sd (., na.rm=T),
                     median = ~ median (., na.rm=T), 
                     #min = ~ min (., na.rm=T), max = ~ max (., na.rm=T),
                     IQR = ~ IQR (., na.rm = T)))

#sch: Q15 p-values
t.test(child_ffq_hpdi_score ~ child_pedsql_sch_score_catg_Q15, data=masterfile)
t.test(child_ffq_animal_score ~ child_pedsql_sch_score_catg_Q15, data=masterfile)
t.test(child_ffq_diet_quality_score ~ child_pedsql_sch_score_catg_Q15, data=masterfile)

wilcox.test(child_ffq_pattern_1_meat_seafood_prepd_meals ~ child_pedsql_sch_score_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_2_cooked_veg_grains_legumes ~ child_pedsql_sch_score_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_3_fruits_raw_veg_cheese ~ child_pedsql_sch_score_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_4_confectioneries_pizza ~ child_pedsql_sch_score_catg_Q15, data=masterfile)
wilcox.test(child_ffq_pattern_5_starchy_foods_sweetened_bev ~ child_pedsql_sch_score_catg_Q15, data=masterfile)