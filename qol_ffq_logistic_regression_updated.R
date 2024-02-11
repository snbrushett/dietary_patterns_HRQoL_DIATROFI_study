###     SCRIPT: LOGISTIC REGRESSION OF QOL (CATEGORIZED) AND DIETARY FEATURES 
###     AUTHOR(S): SIOBHAN BRUSHETT 
###     DESCRIPTION: FITTING LOGISTIC REGRESSION OF Q15 HRQOL OUTCOMES AND TOTAL DIETARY INDEX; ENTIRE DATASET USED (NO DOWNSAMPLING)
###     PROJECT: DIATROFI_BASELINE_2015_2018
###     NOTE(S): https://juliasilge.com/blog/student-debt/
###              https://bookdown.org/pdr_higgins/rmrwr/logistic-regression-and-broom-for-tidying-models.html#evaluating-your-model-assumptions
###              https://juliasilge.com/blog/bird-baths/
###              https://jtools.jacob-long.com/articles/summ.html

#libraries
library(tidyverse)
library(foreach)
library(basecamb) #to add confidence intervals to summary table output

## 0. LOAD DATA
## 1. FIT LOGISTIC MODELS
## 2. MODELS FOR MANUSCRIPT - FDR CORRECTION

                #### =============== 0. LOAD DATA ================ ####
#MAC-book location
setwd("~/OneDrive - UMCG/Prolepsis_internship/merged_2015_2018/masterfile/2_tidy_masterfile_incl_marg_add_hrqol_percentiles/")
setwd("C:/Users/Siobhan Brushett/OneDrive - UMCG/Prolepsis_internship/merged_2015_2018/masterfile/2_tidy_masterfile_incl_marg_add_hrqol_percentiles/")
masterfile <- read.delim("2023_05_09_diatrofi_2015_2018_masterfile_updated.txt", sep = "\t", header = T, stringsAsFactors = T)
#6583 66
str(masterfile)
masterfile$child_demo_school_code <- as.factor(masterfile$child_demo_school_code)

#retrospectively updated dietary indexes to reflect a 10-unit increase with regards to OR (7th July 2023)
table(masterfile$child_ffq_hpdi_score, useNA="ifany")
masterfile$child_ffq_hpdi_score <- masterfile$child_ffq_hpdi_score / 10

table(masterfile$child_ffq_animal_score, useNA="ifany")
masterfile$child_ffq_animal_score <- masterfile$child_ffq_animal_score / 10

                #### =============== 1. FIT LOGISTIC MODELS   ================ ####
#NOTE: for dietary indexes AND NOT DIETARY PATTERNS, models need to be corrected for margarine intake as it was NOT included in the indexes (as in previous studies)

                ####1.1 prepare vectors for linear models ####
#prepare different HRQoL measurements
Q15_measurements <- masterfile [, c("pseudo_ids", 
                                    "child_pedsql_health_related_qol_score_total_catg_Q15",
                                    #"child_pedsql_psysoc_score_catg_Q15",
                                    "child_pedsql_phy_score_catg_Q15",
                                    "child_pedsql_emo_score_catg_Q15",
                                    "child_pedsql_soc_score_catg_Q15",
                                    "child_pedsql_sch_score_catg_Q15")]
rownames(Q15_measurements) <- masterfile$pseudo_ids
head(Q15_measurements)
Q15_measurements$pseudo_ids <- NULL

                ####1.2 model 1: crude model of QoL and diet ####
#dietary indexes
model_1_indexes <- foreach(j = 1:ncol(Q15_measurements),.combine = rbind) %do% {
                glm1 = glm(Q15_measurements[,j] ~ 
                                   masterfile$child_ffq_margarine + 
                                   masterfile$child_ffq_hpdi_score +
                                   masterfile$child_ffq_animal_score, family='binomial')
                sum1 = summary(glm1)
                conf = or_model_summary(glm1) #package: basecamb
                #child_ffq_hpdi_score
                summary_result_hpdi = sum1$coef[nrow(sum1$coef)-1,]
                conf_result_hpdi = conf[nrow(conf)-1,]
                #child_ffq_animal_score
                summary_result_animal = sum1$coef[nrow(sum1$coef),]
                conf_result_animal = conf[nrow(conf),]
                data.frame(
                        #child_ffq_hpdi_score
                        QoL_outcomes_hpdi = colnames(Q15_measurements)[j],
                        Dietary_features_hpdi = "child_ffq_hpdi_score",
                        Odds_ratio_hpdi = conf_result_hpdi[1],
                        Standard_error_hpdi = summary_result_hpdi[2],
                        Statistic_hpdi = summary_result_hpdi[3],
                        CI_lower_hpdi = conf_result_hpdi[2],
                        CI_upper_hpdi = conf_result_hpdi[3],
                        P_value_hpdi = summary_result_hpdi[4],
                        #child_ffq_animal_score
                        QoL_outcomes_animal = colnames(Q15_measurements)[j],
                        Dietary_features_animal = "child_ffq_animal_score",
                        Odds_ratio_animal = conf_result_animal[1],
                        Standard_error_animal = summary_result_animal[2],
                        Statistic_animal = summary_result_animal[3],
                        CI_lower_animal = conf_result_animal[2],
                        CI_upper_animal = conf_result_animal[3],
                        P_value_animal = summary_result_animal[4])}
dim(model_1_indexes)
#5 16 (5 QoL measurements; 8 statistics*2 dietary patterns == 16)

model_1_indexes_hpdi <- model_1_indexes %>% select (c(
        QoL_outcomes_hpdi,
        Dietary_features_hpdi,
        adj_OR,
        Standard_error_hpdi,
        Statistic_hpdi,
        low_ci,
        up_ci,
        P_value_hpdi)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_hpdi,
                Dietary_features=Dietary_features_hpdi,
                Odds_ratio=adj_OR,
                Standard_error=Standard_error_hpdi,
                Statistic=Statistic_hpdi,
                CI_lower=low_ci,
                CI_upper=up_ci,
                P_value=P_value_hpdi)

model_1_indexes_animal <- model_1_indexes %>% select (c(
        QoL_outcomes_animal,
        Dietary_features_animal,
        adj_OR.1,
        Standard_error_animal,
        Statistic_animal,
        low_ci.1,
        up_ci.1,
        P_value_animal)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_animal,
                Dietary_features=Dietary_features_animal,
                Odds_ratio=adj_OR.1,
                Standard_error=Standard_error_animal,
                Statistic=Statistic_animal,
                CI_lower=low_ci.1,
                CI_upper=up_ci.1,
                P_value=P_value_animal)

model_1_indexes_final <- bind_rows(model_1_indexes_hpdi,
                                   model_1_indexes_animal)
rm(model_1_indexes_hpdi, model_1_indexes_animal, model_1_indexes)

#dietary patterns
model_1_patterns <- #foreach(i = 1:length(diet_patterns),.combine = rbind) %:%
        foreach(j = 1:ncol(Q15_measurements),.combine = rbind) %do% {
                glm1 = glm(Q15_measurements[,j] ~ 
                                   masterfile$child_ffq_pattern_1_meat_seafood_prepd_meals +
                                   masterfile$child_ffq_pattern_2_cooked_veg_grains_legumes +
                                   masterfile$child_ffq_pattern_3_fruits_raw_veg_cheese +
                                   masterfile$child_ffq_pattern_4_confectioneries_pizza +
                                   masterfile$child_ffq_pattern_5_starchy_foods_sweetened_bev, family='binomial')
                sum1 = summary(glm1)
                conf = or_model_summary(glm1) #package: basecamb
                #child_ffq_pattern_1_meat_seafood_prepd_meals
                summary_result_p1 = sum1$coef[nrow(sum1$coef)-4,]
                conf_result_p1 = conf[nrow(conf)-4,]
                #child_ffq_pattern_2_cooked_veg_grains_legumes
                summary_result_p2 = sum1$coef[nrow(sum1$coef)-3,]
                conf_result_p2 = conf[nrow(conf)-3,]
                #child_ffq_pattern_3_fruits_raw_veg_cheese
                summary_result_p3 = sum1$coef[nrow(sum1$coef)-2,]
                conf_result_p3 = conf[nrow(conf)-2,]
                #child_ffq_pattern_4_confectioneries_pizza
                summary_result_p4 = sum1$coef[nrow(sum1$coef)-1,]
                conf_result_p4 = conf[nrow(conf)-1,]
                #child_ffq_pattern_5_starchy_foods_sweetened_bev
                summary_result_p5 = sum1$coef[nrow(sum1$coef),]
                conf_result_p5 = conf[nrow(conf),]
                data.frame(
                        #child_ffq_pattern_1_meat_seafood_prepd_meals
                        QoL_outcomes_p1 = colnames(Q15_measurements)[j],
                           Dietary_features_p1 = "child_ffq_pattern_1_meat_seafood_prepd_meals",
                           Odds_ratio_p1 = conf_result_p1[1],
                           #estimate_p1 = summary_result_p1[1],
                           Standard_error_p1 = summary_result_p1[2],
                           Statistic_p1 = summary_result_p1[3],
                           CI_lower_p1 = conf_result_p1[2],
                           CI_upper_p1 = conf_result_p1[3],
                           P_value_p1 = summary_result_p1[4],
                        #child_ffq_pattern_2_cooked_veg_grains_legumes
                        QoL_outcomes_p2 = colnames(Q15_measurements)[j],
                        Dietary_features_p2 = "child_ffq_pattern_2_cooked_veg_grains_legumes",
                        Odds_ratio_p2 = conf_result_p2[1],
                        Standard_error_p2 = summary_result_p2[2],
                        Statistic_p2 = summary_result_p2[3],
                        CI_lower_p2 = conf_result_p2[2],
                        CI_upper_p2 = conf_result_p2[3],
                        P_value_p2 = summary_result_p2[4],
                        #child_ffq_pattern_3_fruits_raw_veg_cheese
                        QoL_outcomes_p3 = colnames(Q15_measurements)[j],
                        Dietary_features_p3 = "child_ffq_pattern_3_fruits_raw_veg_cheese",
                        Odds_ratio_p3 = conf_result_p3[1],
                        Standard_error_p3 = summary_result_p3[2],
                        Statistic_p3 = summary_result_p3[3],
                        CI_lower_p3 = conf_result_p3[2],
                        CI_upper_p3 = conf_result_p3[3],
                        P_value_p3 = summary_result_p3[4],
                        #child_ffq_pattern_4_confectioneries_pizza
                        QoL_outcomes_p4 = colnames(Q15_measurements)[j],
                        Dietary_features_p4 = "child_ffq_pattern_4_confectioneries_pizza",
                        Odds_ratio_p4 = conf_result_p4[1],
                        Standard_error_p4 = summary_result_p4[2],
                        Statistic_p4 = summary_result_p4[3],
                        CI_lower_p4 = conf_result_p4[2],
                        CI_upper_p4 = conf_result_p4[3],
                        P_value_p4 = summary_result_p4[4],
                        #child_ffq_pattern_5_starchy_foods_sweetened_bev
                        QoL_outcomes_p5 = colnames(Q15_measurements)[j],
                        Dietary_features_p5 = "child_ffq_pattern_5_starchy_foods_sweetened_bev",
                        Odds_ratio_p5 = conf_result_p5[1],
                        Standard_error_p5 = summary_result_p5[2],
                        Statistic_p5 = summary_result_p5[3],
                        CI_lower_p5 = conf_result_p5[2],
                        CI_upper_p5 = conf_result_p5[3],
                        P_value_p5 = summary_result_p5[4])}
dim(model_1_patterns)
#5 40 (5 QoL measurements,  8 colnames * 5 dietary patterns == 40)

model_1_patterns_1 <- model_1_patterns %>% select (c(
        QoL_outcomes_p1,
        Dietary_features_p1,
        adj_OR,
        Standard_error_p1,
        Statistic_p1,
        low_ci,
        up_ci,
        P_value_p1)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p1,
                Dietary_features=Dietary_features_p1,
                Odds_ratio=adj_OR,
                Standard_error=Standard_error_p1,
                Statistic=Statistic_p1,
                CI_lower=low_ci,
                CI_upper=up_ci,
                P_value=P_value_p1)

model_1_patterns_2 <- model_1_patterns %>% select (c(
        QoL_outcomes_p2,
        Dietary_features_p2,
        adj_OR.1,
        Standard_error_p2,
        Statistic_p2,
        low_ci.1,
        up_ci.1,
        P_value_p2)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p2,
                Dietary_features=Dietary_features_p2,
                Odds_ratio=adj_OR.1,
                Standard_error=Standard_error_p2,
                Statistic=Statistic_p2,
                CI_lower=low_ci.1,
                CI_upper=up_ci.1,
                P_value=P_value_p2)

model_1_patterns_3 <- model_1_patterns %>% select (c(
        QoL_outcomes_p3,
        Dietary_features_p3,
        adj_OR.2,
        Standard_error_p3,
        Statistic_p3,
        low_ci.2,
        up_ci.2,
        P_value_p3)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p3,
                Dietary_features=Dietary_features_p3,
                Odds_ratio=adj_OR.2,
                Standard_error=Standard_error_p3,
                Statistic=Statistic_p3,
                CI_lower=low_ci.2,
                CI_upper=up_ci.2,
                P_value=P_value_p3)

model_1_patterns_4 <- model_1_patterns %>% select (c(
        QoL_outcomes_p4,
        Dietary_features_p4,
        adj_OR.3,
        Standard_error_p4,
        Statistic_p4,
        low_ci.3,
        up_ci.3,
        P_value_p4)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p4,
                Dietary_features=Dietary_features_p4,
                Odds_ratio=adj_OR.3,
                Standard_error=Standard_error_p4,
                Statistic=Statistic_p4,
                CI_lower=low_ci.3,
                CI_upper=up_ci.3,
                P_value=P_value_p4)

model_1_patterns_5 <- model_1_patterns %>% select (c(
        QoL_outcomes_p5,
        Dietary_features_p5,
        adj_OR.4,
        Standard_error_p5,
        Statistic_p5,
        low_ci.4,
        up_ci.4,
        P_value_p5)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p5,
                Dietary_features=Dietary_features_p5,
                Odds_ratio=adj_OR.4,
                Standard_error=Standard_error_p5,
                Statistic=Statistic_p5,
                CI_lower=low_ci.4,
                CI_upper=up_ci.4,
                P_value=P_value_p5)

model_1_patterns_final <- bind_rows(model_1_patterns_1,
                                    model_1_patterns_2,
                                    model_1_patterns_3,
                                    model_1_patterns_4,
                                    model_1_patterns_5)
rm(model_1_patterns_1, model_1_patterns_2, model_1_patterns_3, model_1_patterns_4, model_1_patterns_5, model_1_patterns)

model_1 <- rbind(model_1_indexes_final, model_1_patterns_final)
model_1 <- model_1 %>%
        mutate(Model="1") %>%
        relocate(Model, .before = QoL_outcomes)

                ####1.3 model 2: model 1 + age + sex ####
#dietary indexes
model_2_indexes <- foreach(j = 1:ncol(Q15_measurements),.combine = rbind) %do% {
                glm1 = glm(Q15_measurements[,j] ~ 
                                   masterfile$child_demo_age_years +
                                   masterfile$child_demo_sex +
                                   masterfile$child_ffq_margarine + 
                                   masterfile$child_ffq_hpdi_score +
                                   masterfile$child_ffq_animal_score, family='binomial')
                sum1 = summary(glm1)
                conf = or_model_summary(glm1) #package: basecamb
                #child_ffq_hpdi_score
                summary_result_hpdi = sum1$coef[nrow(sum1$coef)-1,]
                conf_result_hpdi = conf[nrow(conf)-1,]
                #child_ffq_animal_score
                summary_result_animal = sum1$coef[nrow(sum1$coef),]
                conf_result_animal = conf[nrow(conf),]
                data.frame(
                        #child_ffq_hpdi_score
                        QoL_outcomes_hpdi = colnames(Q15_measurements)[j],
                        Dietary_features_hpdi = "child_ffq_hpdi_score",
                        Odds_ratio_hpdi = conf_result_hpdi[1],
                        Standard_error_hpdi = summary_result_hpdi[2],
                        Statistic_hpdi = summary_result_hpdi[3],
                        CI_lower_hpdi = conf_result_hpdi[2],
                        CI_upper_hpdi = conf_result_hpdi[3],
                        P_value_hpdi = summary_result_hpdi[4],
                        #child_ffq_animal_score
                        QoL_outcomes_animal = colnames(Q15_measurements)[j],
                        Dietary_features_animal = "child_ffq_animal_score",
                        Odds_ratio_animal = conf_result_animal[1],
                        Standard_error_animal = summary_result_animal[2],
                        Statistic_animal = summary_result_animal[3],
                        CI_lower_animal = conf_result_animal[2],
                        CI_upper_animal = conf_result_animal[3],
                        P_value_animal = summary_result_animal[4])}
dim(model_2_indexes)
#5 16 (5 QoL measurements; 8 statistics*2 dietary patterns == 16)

model_2_indexes_hpdi <- model_2_indexes %>% select (c(
        QoL_outcomes_hpdi,
        Dietary_features_hpdi,
        adj_OR,
        Standard_error_hpdi,
        Statistic_hpdi,
        low_ci,
        up_ci,
        P_value_hpdi)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_hpdi,
                Dietary_features=Dietary_features_hpdi,
                Odds_ratio=adj_OR,
                Standard_error=Standard_error_hpdi,
                Statistic=Statistic_hpdi,
                CI_lower=low_ci,
                CI_upper=up_ci,
                P_value=P_value_hpdi)

model_2_indexes_animal <- model_2_indexes %>% select (c(
        QoL_outcomes_animal,
        Dietary_features_animal,
        adj_OR.1,
        Standard_error_animal,
        Statistic_animal,
        low_ci.1,
        up_ci.1,
        P_value_animal)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_animal,
                Dietary_features=Dietary_features_animal,
                Odds_ratio=adj_OR.1,
                Standard_error=Standard_error_animal,
                Statistic=Statistic_animal,
                CI_lower=low_ci.1,
                CI_upper=up_ci.1,
                P_value=P_value_animal)

model_2_indexes_final <- bind_rows(model_2_indexes_hpdi,
                                   model_2_indexes_animal)
rm(model_2_indexes_hpdi, model_2_indexes_animal, model_2_indexes)

#dietary patterns
model_2_patterns <- foreach(j = 1:ncol(Q15_measurements),.combine = rbind) %do% {
                glm1 = glm(Q15_measurements[,j] ~ 
                                   masterfile$child_demo_age_years +
                                   masterfile$child_demo_sex +
                                   masterfile$child_ffq_pattern_1_meat_seafood_prepd_meals +
                                   masterfile$child_ffq_pattern_2_cooked_veg_grains_legumes +
                                   masterfile$child_ffq_pattern_3_fruits_raw_veg_cheese +
                                   masterfile$child_ffq_pattern_4_confectioneries_pizza +
                                   masterfile$child_ffq_pattern_5_starchy_foods_sweetened_bev, family='binomial')
                sum1 = summary(glm1)
                conf = or_model_summary(glm1) #package: basecamb
                #child_ffq_pattern_1_meat_seafood_prepd_meals
                summary_result_p1 = sum1$coef[nrow(sum1$coef)-4,]
                conf_result_p1 = conf[nrow(conf)-4,]
                #child_ffq_pattern_2_cooked_veg_grains_legumes
                summary_result_p2 = sum1$coef[nrow(sum1$coef)-3,]
                conf_result_p2 = conf[nrow(conf)-3,]
                #child_ffq_pattern_3_fruits_raw_veg_cheese
                summary_result_p3 = sum1$coef[nrow(sum1$coef)-2,]
                conf_result_p3 = conf[nrow(conf)-2,]
                #child_ffq_pattern_4_confectioneries_pizza
                summary_result_p4 = sum1$coef[nrow(sum1$coef)-1,]
                conf_result_p4 = conf[nrow(conf)-1,]
                #child_ffq_pattern_5_starchy_foods_sweetened_bev
                summary_result_p5 = sum1$coef[nrow(sum1$coef),]
                conf_result_p5 = conf[nrow(conf),]
                data.frame(
                        #child_ffq_pattern_1_meat_seafood_prepd_meals
                        QoL_outcomes_p1 = colnames(Q15_measurements)[j],
                        Dietary_features_p1 = "child_ffq_pattern_1_meat_seafood_prepd_meals",
                        Odds_ratio_p1 = conf_result_p1[1],
                        #estimate_p1 = summary_result_p1[1],
                        Standard_error_p1 = summary_result_p1[2],
                        Statistic_p1 = summary_result_p1[3],
                        CI_lower_p1 = conf_result_p1[2],
                        CI_upper_p1 = conf_result_p1[3],
                        P_value_p1 = summary_result_p1[4],
                        #child_ffq_pattern_2_cooked_veg_grains_legumes
                        QoL_outcomes_p2 = colnames(Q15_measurements)[j],
                        Dietary_features_p2 = "child_ffq_pattern_2_cooked_veg_grains_legumes",
                        Odds_ratio_p2 = conf_result_p2[1],
                        Standard_error_p2 = summary_result_p2[2],
                        Statistic_p2 = summary_result_p2[3],
                        CI_lower_p2 = conf_result_p2[2],
                        CI_upper_p2 = conf_result_p2[3],
                        P_value_p2 = summary_result_p2[4],
                        #child_ffq_pattern_3_fruits_raw_veg_cheese
                        QoL_outcomes_p3 = colnames(Q15_measurements)[j],
                        Dietary_features_p3 = "child_ffq_pattern_3_fruits_raw_veg_cheese",
                        Odds_ratio_p3 = conf_result_p3[1],
                        Standard_error_p3 = summary_result_p3[2],
                        Statistic_p3 = summary_result_p3[3],
                        CI_lower_p3 = conf_result_p3[2],
                        CI_upper_p3 = conf_result_p3[3],
                        P_value_p3 = summary_result_p3[4],
                        #child_ffq_pattern_4_confectioneries_pizza
                        QoL_outcomes_p4 = colnames(Q15_measurements)[j],
                        Dietary_features_p4 = "child_ffq_pattern_4_confectioneries_pizza",
                        Odds_ratio_p4 = conf_result_p4[1],
                        Standard_error_p4 = summary_result_p4[2],
                        Statistic_p4 = summary_result_p4[3],
                        CI_lower_p4 = conf_result_p4[2],
                        CI_upper_p4 = conf_result_p4[3],
                        P_value_p4 = summary_result_p4[4],
                        #child_ffq_pattern_5_starchy_foods_sweetened_bev
                        QoL_outcomes_p5 = colnames(Q15_measurements)[j],
                        Dietary_features_p5 = "child_ffq_pattern_5_starchy_foods_sweetened_bev",
                        Odds_ratio_p5 = conf_result_p5[1],
                        Standard_error_p5 = summary_result_p5[2],
                        Statistic_p5 = summary_result_p5[3],
                        CI_lower_p5 = conf_result_p5[2],
                        CI_upper_p5 = conf_result_p5[3],
                        P_value_p5 = summary_result_p5[4])}
dim(model_2_patterns)
#5 40 

model_2_patterns_1 <- model_2_patterns %>% select (c(
        QoL_outcomes_p1,
        Dietary_features_p1,
        adj_OR,
        Standard_error_p1,
        Statistic_p1,
        low_ci,
        up_ci,
        P_value_p1)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p1,
                Dietary_features=Dietary_features_p1,
                Odds_ratio=adj_OR,
                Standard_error=Standard_error_p1,
                Statistic=Statistic_p1,
                CI_lower=low_ci,
                CI_upper=up_ci,
                P_value=P_value_p1)

model_2_patterns_2 <- model_2_patterns %>% select (c(
        QoL_outcomes_p2,
        Dietary_features_p2,
        adj_OR.1,
        Standard_error_p2,
        Statistic_p2,
        low_ci.1,
        up_ci.1,
        P_value_p2)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p2,
                Dietary_features=Dietary_features_p2,
                Odds_ratio=adj_OR.1,
                Standard_error=Standard_error_p2,
                Statistic=Statistic_p2,
                CI_lower=low_ci.1,
                CI_upper=up_ci.1,
                P_value=P_value_p2)

model_2_patterns_3 <- model_2_patterns %>% select (c(
        QoL_outcomes_p3,
        Dietary_features_p3,
        adj_OR.2,
        Standard_error_p3,
        Statistic_p3,
        low_ci.2,
        up_ci.2,
        P_value_p3)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p3,
                Dietary_features=Dietary_features_p3,
                Odds_ratio=adj_OR.2,
                Standard_error=Standard_error_p3,
                Statistic=Statistic_p3,
                CI_lower=low_ci.2,
                CI_upper=up_ci.2,
                P_value=P_value_p3)

model_2_patterns_4 <- model_2_patterns %>% select (c(
        QoL_outcomes_p4,
        Dietary_features_p4,
        adj_OR.3,
        Standard_error_p4,
        Statistic_p4,
        low_ci.3,
        up_ci.3,
        P_value_p4)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p4,
                Dietary_features=Dietary_features_p4,
                Odds_ratio=adj_OR.3,
                Standard_error=Standard_error_p4,
                Statistic=Statistic_p4,
                CI_lower=low_ci.3,
                CI_upper=up_ci.3,
                P_value=P_value_p4)

model_2_patterns_5 <- model_2_patterns %>% select (c(
        QoL_outcomes_p5,
        Dietary_features_p5,
        adj_OR.4,
        Standard_error_p5,
        Statistic_p5,
        low_ci.4,
        up_ci.4,
        P_value_p5)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p5,
                Dietary_features=Dietary_features_p5,
                Odds_ratio=adj_OR.4,
                Standard_error=Standard_error_p5,
                Statistic=Statistic_p5,
                CI_lower=low_ci.4,
                CI_upper=up_ci.4,
                P_value=P_value_p5)

model_2_patterns_final <- bind_rows(model_2_patterns_1,
                                    model_2_patterns_2,
                                    model_2_patterns_3,
                                    model_2_patterns_4,
                                    model_2_patterns_5)
rm(model_2_patterns_1, model_2_patterns_2, model_2_patterns_3, model_2_patterns_4, model_2_patterns_5, model_2_patterns)

model_2 <- rbind(model_2_indexes_final, model_2_patterns_final)
model_2 <- model_2 %>%
        mutate(Model="2") %>%
        relocate(Model, .before = QoL_outcomes)

                ####1.4 model 3: model 2 + family SES ####
#dietary indexes
model_3_indexes <- foreach(j = 1:ncol(Q15_measurements),.combine = rbind) %do% {
                glm1 = glm(Q15_measurements[,j] ~ 
                                   masterfile$child_demo_age_years +
                                   masterfile$child_demo_sex +
                                   masterfile$family_demo_FAS_group +
                                   masterfile$child_ffq_margarine + 
                                   masterfile$child_ffq_hpdi_score +
                                   masterfile$child_ffq_animal_score, family='binomial')
                sum1 = summary(glm1)
                conf = or_model_summary(glm1) #package: basecamb
                #child_ffq_hpdi_score
                summary_result_hpdi = sum1$coef[nrow(sum1$coef)-1,]
                conf_result_hpdi = conf[nrow(conf)-1,]
                #child_ffq_animal_score
                summary_result_animal = sum1$coef[nrow(sum1$coef),]
                conf_result_animal = conf[nrow(conf),]
                data.frame(
                        #child_ffq_hpdi_score
                        QoL_outcomes_hpdi = colnames(Q15_measurements)[j],
                        Dietary_features_hpdi = "child_ffq_hpdi_score",
                        Odds_ratio_hpdi = conf_result_hpdi[1],
                        Standard_error_hpdi = summary_result_hpdi[2],
                        Statistic_hpdi = summary_result_hpdi[3],
                        CI_lower_hpdi = conf_result_hpdi[2],
                        CI_upper_hpdi = conf_result_hpdi[3],
                        P_value_hpdi = summary_result_hpdi[4],
                        #child_ffq_animal_score
                        QoL_outcomes_animal = colnames(Q15_measurements)[j],
                        Dietary_features_animal = "child_ffq_animal_score",
                        Odds_ratio_animal = conf_result_animal[1],
                        Standard_error_animal = summary_result_animal[2],
                        Statistic_animal = summary_result_animal[3],
                        CI_lower_animal = conf_result_animal[2],
                        CI_upper_animal = conf_result_animal[3],
                        P_value_animal = summary_result_animal[4])}
dim(model_3_indexes)
#5 16 (5 QoL measurements; 8 statistics*2 dietary patterns == 16)

model_3_indexes_hpdi <- model_3_indexes %>% select (c(
        QoL_outcomes_hpdi,
        Dietary_features_hpdi,
        adj_OR,
        Standard_error_hpdi,
        Statistic_hpdi,
        low_ci,
        up_ci,
        P_value_hpdi)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_hpdi,
                Dietary_features=Dietary_features_hpdi,
                Odds_ratio=adj_OR,
                Standard_error=Standard_error_hpdi,
                Statistic=Statistic_hpdi,
                CI_lower=low_ci,
                CI_upper=up_ci,
                P_value=P_value_hpdi)

model_3_indexes_animal <- model_3_indexes %>% select (c(
        QoL_outcomes_animal,
        Dietary_features_animal,
        adj_OR.1,
        Standard_error_animal,
        Statistic_animal,
        low_ci.1,
        up_ci.1,
        P_value_animal)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_animal,
                Dietary_features=Dietary_features_animal,
                Odds_ratio=adj_OR.1,
                Standard_error=Standard_error_animal,
                Statistic=Statistic_animal,
                CI_lower=low_ci.1,
                CI_upper=up_ci.1,
                P_value=P_value_animal)

model_3_indexes_final <- bind_rows(model_3_indexes_hpdi,
                                   model_3_indexes_animal)
rm(model_3_indexes_hpdi, model_3_indexes_animal, model_3_indexes)

#dietary patterns
model_3_patterns <- foreach(j = 1:ncol(Q15_measurements),.combine = rbind) %do% {
                glm1 = glm(Q15_measurements[,j] ~ 
                                   masterfile$child_demo_age_years +
                                   masterfile$child_demo_sex +
                                   masterfile$family_demo_FAS_group +
                                   masterfile$child_ffq_pattern_1_meat_seafood_prepd_meals +
                                   masterfile$child_ffq_pattern_2_cooked_veg_grains_legumes +
                                   masterfile$child_ffq_pattern_3_fruits_raw_veg_cheese +
                                   masterfile$child_ffq_pattern_4_confectioneries_pizza +
                                   masterfile$child_ffq_pattern_5_starchy_foods_sweetened_bev, family='binomial')
                sum1 = summary(glm1)
                conf = or_model_summary(glm1) #package: basecamb
                #child_ffq_pattern_1_meat_seafood_prepd_meals
                summary_result_p1 = sum1$coef[nrow(sum1$coef)-4,]
                conf_result_p1 = conf[nrow(conf)-4,]
                #child_ffq_pattern_2_cooked_veg_grains_legumes
                summary_result_p2 = sum1$coef[nrow(sum1$coef)-3,]
                conf_result_p2 = conf[nrow(conf)-3,]
                #child_ffq_pattern_3_fruits_raw_veg_cheese
                summary_result_p3 = sum1$coef[nrow(sum1$coef)-2,]
                conf_result_p3 = conf[nrow(conf)-2,]
                #child_ffq_pattern_4_confectioneries_pizza
                summary_result_p4 = sum1$coef[nrow(sum1$coef)-1,]
                conf_result_p4 = conf[nrow(conf)-1,]
                #child_ffq_pattern_5_starchy_foods_sweetened_bev
                summary_result_p5 = sum1$coef[nrow(sum1$coef),]
                conf_result_p5 = conf[nrow(conf),]
                data.frame(
                        #child_ffq_pattern_1_meat_seafood_prepd_meals
                        QoL_outcomes_p1 = colnames(Q15_measurements)[j],
                        Dietary_features_p1 = "child_ffq_pattern_1_meat_seafood_prepd_meals",
                        Odds_ratio_p1 = conf_result_p1[1],
                        #estimate_p1 = summary_result_p1[1],
                        Standard_error_p1 = summary_result_p1[2],
                        Statistic_p1 = summary_result_p1[3],
                        CI_lower_p1 = conf_result_p1[2],
                        CI_upper_p1 = conf_result_p1[3],
                        P_value_p1 = summary_result_p1[4],
                        #child_ffq_pattern_2_cooked_veg_grains_legumes
                        QoL_outcomes_p2 = colnames(Q15_measurements)[j],
                        Dietary_features_p2 = "child_ffq_pattern_2_cooked_veg_grains_legumes",
                        Odds_ratio_p2 = conf_result_p2[1],
                        Standard_error_p2 = summary_result_p2[2],
                        Statistic_p2 = summary_result_p2[3],
                        CI_lower_p2 = conf_result_p2[2],
                        CI_upper_p2 = conf_result_p2[3],
                        P_value_p2 = summary_result_p2[4],
                        #child_ffq_pattern_3_fruits_raw_veg_cheese
                        QoL_outcomes_p3 = colnames(Q15_measurements)[j],
                        Dietary_features_p3 = "child_ffq_pattern_3_fruits_raw_veg_cheese",
                        Odds_ratio_p3 = conf_result_p3[1],
                        Standard_error_p3 = summary_result_p3[2],
                        Statistic_p3 = summary_result_p3[3],
                        CI_lower_p3 = conf_result_p3[2],
                        CI_upper_p3 = conf_result_p3[3],
                        P_value_p3 = summary_result_p3[4],
                        #child_ffq_pattern_4_confectioneries_pizza
                        QoL_outcomes_p4 = colnames(Q15_measurements)[j],
                        Dietary_features_p4 = "child_ffq_pattern_4_confectioneries_pizza",
                        Odds_ratio_p4 = conf_result_p4[1],
                        Standard_error_p4 = summary_result_p4[2],
                        Statistic_p4 = summary_result_p4[3],
                        CI_lower_p4 = conf_result_p4[2],
                        CI_upper_p4 = conf_result_p4[3],
                        P_value_p4 = summary_result_p4[4],
                        #child_ffq_pattern_5_starchy_foods_sweetened_bev
                        QoL_outcomes_p5 = colnames(Q15_measurements)[j],
                        Dietary_features_p5 = "child_ffq_pattern_5_starchy_foods_sweetened_bev",
                        Odds_ratio_p5 = conf_result_p5[1],
                        Standard_error_p5 = summary_result_p5[2],
                        Statistic_p5 = summary_result_p5[3],
                        CI_lower_p5 = conf_result_p5[2],
                        CI_upper_p5 = conf_result_p5[3],
                        P_value_p5 = summary_result_p5[4])}
dim(model_3_patterns)
#5 40 

model_3_patterns_1 <- model_3_patterns %>% select (c(
        QoL_outcomes_p1,
        Dietary_features_p1,
        adj_OR,
        Standard_error_p1,
        Statistic_p1,
        low_ci,
        up_ci,
        P_value_p1)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p1,
                Dietary_features=Dietary_features_p1,
                Odds_ratio=adj_OR,
                Standard_error=Standard_error_p1,
                Statistic=Statistic_p1,
                CI_lower=low_ci,
                CI_upper=up_ci,
                P_value=P_value_p1)

model_3_patterns_2 <- model_3_patterns %>% select (c(
        QoL_outcomes_p2,
        Dietary_features_p2,
        adj_OR.1,
        Standard_error_p2,
        Statistic_p2,
        low_ci.1,
        up_ci.1,
        P_value_p2)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p2,
                Dietary_features=Dietary_features_p2,
                Odds_ratio=adj_OR.1,
                Standard_error=Standard_error_p2,
                Statistic=Statistic_p2,
                CI_lower=low_ci.1,
                CI_upper=up_ci.1,
                P_value=P_value_p2)

model_3_patterns_3 <- model_3_patterns %>% select (c(
        QoL_outcomes_p3,
        Dietary_features_p3,
        adj_OR.2,
        Standard_error_p3,
        Statistic_p3,
        low_ci.2,
        up_ci.2,
        P_value_p3)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p3,
                Dietary_features=Dietary_features_p3,
                Odds_ratio=adj_OR.2,
                Standard_error=Standard_error_p3,
                Statistic=Statistic_p3,
                CI_lower=low_ci.2,
                CI_upper=up_ci.2,
                P_value=P_value_p3)

model_3_patterns_4 <- model_3_patterns %>% select (c(
        QoL_outcomes_p4,
        Dietary_features_p4,
        adj_OR.3,
        Standard_error_p4,
        Statistic_p4,
        low_ci.3,
        up_ci.3,
        P_value_p4)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p4,
                Dietary_features=Dietary_features_p4,
                Odds_ratio=adj_OR.3,
                Standard_error=Standard_error_p4,
                Statistic=Statistic_p4,
                CI_lower=low_ci.3,
                CI_upper=up_ci.3,
                P_value=P_value_p4)

model_3_patterns_5 <- model_3_patterns %>% select (c(
        QoL_outcomes_p5,
        Dietary_features_p5,
        adj_OR.4,
        Standard_error_p5,
        Statistic_p5,
        low_ci.4,
        up_ci.4,
        P_value_p5)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p5,
                Dietary_features=Dietary_features_p5,
                Odds_ratio=adj_OR.4,
                Standard_error=Standard_error_p5,
                Statistic=Statistic_p5,
                CI_lower=low_ci.4,
                CI_upper=up_ci.4,
                P_value=P_value_p5)

model_3_patterns_final <- bind_rows(model_3_patterns_1,
                                    model_3_patterns_2,
                                    model_3_patterns_3,
                                    model_3_patterns_4,
                                    model_3_patterns_5)
rm(model_3_patterns_1, model_3_patterns_2, model_3_patterns_3, model_3_patterns_4, model_3_patterns_5, model_3_patterns)

model_3 <- rbind(model_3_indexes_final, model_3_patterns_final)
model_3 <- model_3 %>%
        mutate(Model="3") %>%
        relocate(Model, .before = QoL_outcomes)

                ####1.5 model 4: model 3 + student immigrant status ####
#dietary indexes
model_4_indexes <- foreach(j = 1:ncol(Q15_measurements),.combine = rbind) %do% {
                glm1 = glm(Q15_measurements[,j] ~ 
                                   masterfile$child_demo_age_years +
                                   masterfile$child_demo_sex +
                                   masterfile$family_demo_FAS_group +
                                   masterfile$child_demo_immigrant_status +
                                   masterfile$child_ffq_margarine + 
                                   masterfile$child_ffq_hpdi_score +
                                   masterfile$child_ffq_animal_score, family='binomial')
                sum1 = summary(glm1)
                conf = or_model_summary(glm1) #package: basecamb
                #child_ffq_hpdi_score
                summary_result_hpdi = sum1$coef[nrow(sum1$coef)-1,]
                conf_result_hpdi = conf[nrow(conf)-1,]
                #child_ffq_animal_score
                summary_result_animal = sum1$coef[nrow(sum1$coef),]
                conf_result_animal = conf[nrow(conf),]
                data.frame(
                        #child_ffq_hpdi_score
                        QoL_outcomes_hpdi = colnames(Q15_measurements)[j],
                        Dietary_features_hpdi = "child_ffq_hpdi_score",
                        Odds_ratio_hpdi = conf_result_hpdi[1],
                        Standard_error_hpdi = summary_result_hpdi[2],
                        Statistic_hpdi = summary_result_hpdi[3],
                        CI_lower_hpdi = conf_result_hpdi[2],
                        CI_upper_hpdi = conf_result_hpdi[3],
                        P_value_hpdi = summary_result_hpdi[4],
                        #child_ffq_animal_score
                        QoL_outcomes_animal = colnames(Q15_measurements)[j],
                        Dietary_features_animal = "child_ffq_animal_score",
                        Odds_ratio_animal = conf_result_animal[1],
                        Standard_error_animal = summary_result_animal[2],
                        Statistic_animal = summary_result_animal[3],
                        CI_lower_animal = conf_result_animal[2],
                        CI_upper_animal = conf_result_animal[3],
                        P_value_animal = summary_result_animal[4])}
dim(model_4_indexes)
#5 16 (5 QoL measurements; 8 statistics*2 dietary patterns == 16)

model_4_indexes_hpdi <- model_4_indexes %>% select (c(
        QoL_outcomes_hpdi,
        Dietary_features_hpdi,
        adj_OR,
        Standard_error_hpdi,
        Statistic_hpdi,
        low_ci,
        up_ci,
        P_value_hpdi)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_hpdi,
                Dietary_features=Dietary_features_hpdi,
                Odds_ratio=adj_OR,
                Standard_error=Standard_error_hpdi,
                Statistic=Statistic_hpdi,
                CI_lower=low_ci,
                CI_upper=up_ci,
                P_value=P_value_hpdi)

model_4_indexes_animal <- model_4_indexes %>% select (c(
        QoL_outcomes_animal,
        Dietary_features_animal,
        adj_OR.1,
        Standard_error_animal,
        Statistic_animal,
        low_ci.1,
        up_ci.1,
        P_value_animal)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_animal,
                Dietary_features=Dietary_features_animal,
                Odds_ratio=adj_OR.1,
                Standard_error=Standard_error_animal,
                Statistic=Statistic_animal,
                CI_lower=low_ci.1,
                CI_upper=up_ci.1,
                P_value=P_value_animal)

model_4_indexes_final <- bind_rows(model_4_indexes_hpdi,
                                   model_4_indexes_animal)
rm(model_4_indexes_hpdi, model_4_indexes_animal, model_4_indexes)

#dietary patterns
model_4_patterns <- foreach(j = 1:ncol(Q15_measurements),.combine = rbind) %do% {
                glm1 = glm(Q15_measurements[,j] ~ 
                                   masterfile$child_demo_age_years +
                                   masterfile$child_demo_sex +
                                   masterfile$family_demo_FAS_group +
                                   masterfile$child_demo_immigrant_status +
                                   masterfile$child_ffq_pattern_1_meat_seafood_prepd_meals +
                                   masterfile$child_ffq_pattern_2_cooked_veg_grains_legumes +
                                   masterfile$child_ffq_pattern_3_fruits_raw_veg_cheese +
                                   masterfile$child_ffq_pattern_4_confectioneries_pizza +
                                   masterfile$child_ffq_pattern_5_starchy_foods_sweetened_bev, family='binomial')
                sum1 = summary(glm1)
                conf = or_model_summary(glm1) #package: basecamb
                #child_ffq_pattern_1_meat_seafood_prepd_meals
                summary_result_p1 = sum1$coef[nrow(sum1$coef)-4,]
                conf_result_p1 = conf[nrow(conf)-4,]
                #child_ffq_pattern_2_cooked_veg_grains_legumes
                summary_result_p2 = sum1$coef[nrow(sum1$coef)-3,]
                conf_result_p2 = conf[nrow(conf)-3,]
                #child_ffq_pattern_3_fruits_raw_veg_cheese
                summary_result_p3 = sum1$coef[nrow(sum1$coef)-2,]
                conf_result_p3 = conf[nrow(conf)-2,]
                #child_ffq_pattern_4_confectioneries_pizza
                summary_result_p4 = sum1$coef[nrow(sum1$coef)-1,]
                conf_result_p4 = conf[nrow(conf)-1,]
                #child_ffq_pattern_5_starchy_foods_sweetened_bev
                summary_result_p5 = sum1$coef[nrow(sum1$coef),]
                conf_result_p5 = conf[nrow(conf),]
                data.frame(
                        #child_ffq_pattern_1_meat_seafood_prepd_meals
                        QoL_outcomes_p1 = colnames(Q15_measurements)[j],
                        Dietary_features_p1 = "child_ffq_pattern_1_meat_seafood_prepd_meals",
                        Odds_ratio_p1 = conf_result_p1[1],
                        #estimate_p1 = summary_result_p1[1],
                        Standard_error_p1 = summary_result_p1[2],
                        Statistic_p1 = summary_result_p1[3],
                        CI_lower_p1 = conf_result_p1[2],
                        CI_upper_p1 = conf_result_p1[3],
                        P_value_p1 = summary_result_p1[4],
                        #child_ffq_pattern_2_cooked_veg_grains_legumes
                        QoL_outcomes_p2 = colnames(Q15_measurements)[j],
                        Dietary_features_p2 = "child_ffq_pattern_2_cooked_veg_grains_legumes",
                        Odds_ratio_p2 = conf_result_p2[1],
                        Standard_error_p2 = summary_result_p2[2],
                        Statistic_p2 = summary_result_p2[3],
                        CI_lower_p2 = conf_result_p2[2],
                        CI_upper_p2 = conf_result_p2[3],
                        P_value_p2 = summary_result_p2[4],
                        #child_ffq_pattern_3_fruits_raw_veg_cheese
                        QoL_outcomes_p3 = colnames(Q15_measurements)[j],
                        Dietary_features_p3 = "child_ffq_pattern_3_fruits_raw_veg_cheese",
                        Odds_ratio_p3 = conf_result_p3[1],
                        Standard_error_p3 = summary_result_p3[2],
                        Statistic_p3 = summary_result_p3[3],
                        CI_lower_p3 = conf_result_p3[2],
                        CI_upper_p3 = conf_result_p3[3],
                        P_value_p3 = summary_result_p3[4],
                        #child_ffq_pattern_4_confectioneries_pizza
                        QoL_outcomes_p4 = colnames(Q15_measurements)[j],
                        Dietary_features_p4 = "child_ffq_pattern_4_confectioneries_pizza",
                        Odds_ratio_p4 = conf_result_p4[1],
                        Standard_error_p4 = summary_result_p4[2],
                        Statistic_p4 = summary_result_p4[3],
                        CI_lower_p4 = conf_result_p4[2],
                        CI_upper_p4 = conf_result_p4[3],
                        P_value_p4 = summary_result_p4[4],
                        #child_ffq_pattern_5_starchy_foods_sweetened_bev
                        QoL_outcomes_p5 = colnames(Q15_measurements)[j],
                        Dietary_features_p5 = "child_ffq_pattern_5_starchy_foods_sweetened_bev",
                        Odds_ratio_p5 = conf_result_p5[1],
                        Standard_error_p5 = summary_result_p5[2],
                        Statistic_p5 = summary_result_p5[3],
                        CI_lower_p5 = conf_result_p5[2],
                        CI_upper_p5 = conf_result_p5[3],
                        P_value_p5 = summary_result_p5[4])}
dim(model_4_patterns)
#5 40 

model_4_patterns_1 <- model_4_patterns %>% select (c(
        QoL_outcomes_p1,
        Dietary_features_p1,
        adj_OR,
        Standard_error_p1,
        Statistic_p1,
        low_ci,
        up_ci,
        P_value_p1)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p1,
                Dietary_features=Dietary_features_p1,
                Odds_ratio=adj_OR,
                Standard_error=Standard_error_p1,
                Statistic=Statistic_p1,
                CI_lower=low_ci,
                CI_upper=up_ci,
                P_value=P_value_p1)

model_4_patterns_2 <- model_4_patterns %>% select (c(
        QoL_outcomes_p2,
        Dietary_features_p2,
        adj_OR.1,
        Standard_error_p2,
        Statistic_p2,
        low_ci.1,
        up_ci.1,
        P_value_p2)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p2,
                Dietary_features=Dietary_features_p2,
                Odds_ratio=adj_OR.1,
                Standard_error=Standard_error_p2,
                Statistic=Statistic_p2,
                CI_lower=low_ci.1,
                CI_upper=up_ci.1,
                P_value=P_value_p2)

model_4_patterns_3 <- model_4_patterns %>% select (c(
        QoL_outcomes_p3,
        Dietary_features_p3,
        adj_OR.2,
        Standard_error_p3,
        Statistic_p3,
        low_ci.2,
        up_ci.2,
        P_value_p3)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p3,
                Dietary_features=Dietary_features_p3,
                Odds_ratio=adj_OR.2,
                Standard_error=Standard_error_p3,
                Statistic=Statistic_p3,
                CI_lower=low_ci.2,
                CI_upper=up_ci.2,
                P_value=P_value_p3)

model_4_patterns_4 <- model_4_patterns %>% select (c(
        QoL_outcomes_p4,
        Dietary_features_p4,
        adj_OR.3,
        Standard_error_p4,
        Statistic_p4,
        low_ci.3,
        up_ci.3,
        P_value_p4)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p4,
                Dietary_features=Dietary_features_p4,
                Odds_ratio=adj_OR.3,
                Standard_error=Standard_error_p4,
                Statistic=Statistic_p4,
                CI_lower=low_ci.3,
                CI_upper=up_ci.3,
                P_value=P_value_p4)

model_4_patterns_5 <- model_4_patterns %>% select (c(
        QoL_outcomes_p5,
        Dietary_features_p5,
        adj_OR.4,
        Standard_error_p5,
        Statistic_p5,
        low_ci.4,
        up_ci.4,
        P_value_p5)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p5,
                Dietary_features=Dietary_features_p5,
                Odds_ratio=adj_OR.4,
                Standard_error=Standard_error_p5,
                Statistic=Statistic_p5,
                CI_lower=low_ci.4,
                CI_upper=up_ci.4,
                P_value=P_value_p5)

model_4_patterns_final <- bind_rows(model_4_patterns_1,
                                    model_4_patterns_2,
                                    model_4_patterns_3,
                                    model_4_patterns_4,
                                    model_4_patterns_5)
rm(model_4_patterns_1, model_4_patterns_2, model_4_patterns_3, model_4_patterns_4, model_4_patterns_5, model_4_patterns)

model_4 <- rbind(model_4_indexes_final, model_4_patterns_final)
model_4 <- model_4 %>%
        mutate(Model="4") %>%
        relocate(Model, .before = QoL_outcomes)

                ####1.6 model 5: model 4 + school_region ####
#dietary indexes
model_5_indexes <- foreach(j = 1:ncol(Q15_measurements),.combine = rbind) %do% {
                glm1 = glm(Q15_measurements[,j] ~ 
                                   masterfile$child_demo_age_years +
                                   masterfile$child_demo_sex +
                                   masterfile$family_demo_FAS_group +
                                   masterfile$child_demo_immigrant_status +
                                   masterfile$child_demo_school_region +
                                   masterfile$child_ffq_margarine + 
                                   masterfile$child_ffq_hpdi_score +
                                   masterfile$child_ffq_animal_score, family='binomial')
                sum1 = summary(glm1)
                conf = or_model_summary(glm1) #package: basecamb
                summary_result = sum1$coef[nrow(sum1$coef),]
                conf_result = conf[nrow(conf),]
                #child_ffq_hpdi_score
                summary_result_hpdi = sum1$coef[nrow(sum1$coef)-1,]
                conf_result_hpdi = conf[nrow(conf)-1,]
                #child_ffq_animal_score
                summary_result_animal = sum1$coef[nrow(sum1$coef),]
                conf_result_animal = conf[nrow(conf),]
                data.frame(
                        #child_ffq_hpdi_score
                        QoL_outcomes_hpdi = colnames(Q15_measurements)[j],
                        Dietary_features_hpdi = "child_ffq_hpdi_score",
                        Odds_ratio_hpdi = conf_result_hpdi[1],
                        Standard_error_hpdi = summary_result_hpdi[2],
                        Statistic_hpdi = summary_result_hpdi[3],
                        CI_lower_hpdi = conf_result_hpdi[2],
                        CI_upper_hpdi = conf_result_hpdi[3],
                        P_value_hpdi = summary_result_hpdi[4],
                        #child_ffq_animal_score
                        QoL_outcomes_animal = colnames(Q15_measurements)[j],
                        Dietary_features_animal = "child_ffq_animal_score",
                        Odds_ratio_animal = conf_result_animal[1],
                        Standard_error_animal = summary_result_animal[2],
                        Statistic_animal = summary_result_animal[3],
                        CI_lower_animal = conf_result_animal[2],
                        CI_upper_animal = conf_result_animal[3],
                        P_value_animal = summary_result_animal[4])}
dim(model_5_indexes)
#5 16 (5 QoL measurements; 8 statistics*2 dietary patterns == 16)

model_5_indexes_hpdi <- model_5_indexes %>% select (c(
        QoL_outcomes_hpdi,
        Dietary_features_hpdi,
        adj_OR,
        Standard_error_hpdi,
        Statistic_hpdi,
        low_ci,
        up_ci,
        P_value_hpdi)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_hpdi,
                Dietary_features=Dietary_features_hpdi,
                Odds_ratio=adj_OR,
                Standard_error=Standard_error_hpdi,
                Statistic=Statistic_hpdi,
                CI_lower=low_ci,
                CI_upper=up_ci,
                P_value=P_value_hpdi)

model_5_indexes_animal <- model_5_indexes %>% select (c(
        QoL_outcomes_animal,
        Dietary_features_animal,
        adj_OR.1,
        Standard_error_animal,
        Statistic_animal,
        low_ci.1,
        up_ci.1,
        P_value_animal)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_animal,
                Dietary_features=Dietary_features_animal,
                Odds_ratio=adj_OR.1,
                Standard_error=Standard_error_animal,
                Statistic=Statistic_animal,
                CI_lower=low_ci.1,
                CI_upper=up_ci.1,
                P_value=P_value_animal)

model_5_indexes_final <- bind_rows(model_5_indexes_hpdi,
                                   model_5_indexes_animal)
rm(model_5_indexes_hpdi, model_5_indexes_animal, model_5_indexes)

#dietary patterns
model_5_patterns <- foreach(j = 1:ncol(Q15_measurements),.combine = rbind) %do% {
                glm1 = glm(Q15_measurements[,j] ~ 
                                   masterfile$child_demo_age_years +
                                   masterfile$child_demo_sex +
                                   masterfile$family_demo_FAS_group +
                                   masterfile$child_demo_immigrant_status +
                                   masterfile$child_demo_school_region +
                                   masterfile$child_ffq_pattern_1_meat_seafood_prepd_meals +
                                   masterfile$child_ffq_pattern_2_cooked_veg_grains_legumes +
                                   masterfile$child_ffq_pattern_3_fruits_raw_veg_cheese +
                                   masterfile$child_ffq_pattern_4_confectioneries_pizza +
                                   masterfile$child_ffq_pattern_5_starchy_foods_sweetened_bev, family='binomial')
                sum1 = summary(glm1)
                conf = or_model_summary(glm1) #package: basecamb
                #child_ffq_pattern_1_meat_seafood_prepd_meals
                summary_result_p1 = sum1$coef[nrow(sum1$coef)-4,]
                conf_result_p1 = conf[nrow(conf)-4,]
                #child_ffq_pattern_2_cooked_veg_grains_legumes
                summary_result_p2 = sum1$coef[nrow(sum1$coef)-3,]
                conf_result_p2 = conf[nrow(conf)-3,]
                #child_ffq_pattern_3_fruits_raw_veg_cheese
                summary_result_p3 = sum1$coef[nrow(sum1$coef)-2,]
                conf_result_p3 = conf[nrow(conf)-2,]
                #child_ffq_pattern_4_confectioneries_pizza
                summary_result_p4 = sum1$coef[nrow(sum1$coef)-1,]
                conf_result_p4 = conf[nrow(conf)-1,]
                #child_ffq_pattern_5_starchy_foods_sweetened_bev
                summary_result_p5 = sum1$coef[nrow(sum1$coef),]
                conf_result_p5 = conf[nrow(conf),]
                data.frame(
                        #child_ffq_pattern_1_meat_seafood_prepd_meals
                        QoL_outcomes_p1 = colnames(Q15_measurements)[j],
                        Dietary_features_p1 = "child_ffq_pattern_1_meat_seafood_prepd_meals",
                        Odds_ratio_p1 = conf_result_p1[1],
                        #estimate_p1 = summary_result_p1[1],
                        Standard_error_p1 = summary_result_p1[2],
                        Statistic_p1 = summary_result_p1[3],
                        CI_lower_p1 = conf_result_p1[2],
                        CI_upper_p1 = conf_result_p1[3],
                        P_value_p1 = summary_result_p1[4],
                        #child_ffq_pattern_2_cooked_veg_grains_legumes
                        QoL_outcomes_p2 = colnames(Q15_measurements)[j],
                        Dietary_features_p2 = "child_ffq_pattern_2_cooked_veg_grains_legumes",
                        Odds_ratio_p2 = conf_result_p2[1],
                        Standard_error_p2 = summary_result_p2[2],
                        Statistic_p2 = summary_result_p2[3],
                        CI_lower_p2 = conf_result_p2[2],
                        CI_upper_p2 = conf_result_p2[3],
                        P_value_p2 = summary_result_p2[4],
                        #child_ffq_pattern_3_fruits_raw_veg_cheese
                        QoL_outcomes_p3 = colnames(Q15_measurements)[j],
                        Dietary_features_p3 = "child_ffq_pattern_3_fruits_raw_veg_cheese",
                        Odds_ratio_p3 = conf_result_p3[1],
                        Standard_error_p3 = summary_result_p3[2],
                        Statistic_p3 = summary_result_p3[3],
                        CI_lower_p3 = conf_result_p3[2],
                        CI_upper_p3 = conf_result_p3[3],
                        P_value_p3 = summary_result_p3[4],
                        #child_ffq_pattern_4_confectioneries_pizza
                        QoL_outcomes_p4 = colnames(Q15_measurements)[j],
                        Dietary_features_p4 = "child_ffq_pattern_4_confectioneries_pizza",
                        Odds_ratio_p4 = conf_result_p4[1],
                        Standard_error_p4 = summary_result_p4[2],
                        Statistic_p4 = summary_result_p4[3],
                        CI_lower_p4 = conf_result_p4[2],
                        CI_upper_p4 = conf_result_p4[3],
                        P_value_p4 = summary_result_p4[4],
                        #child_ffq_pattern_5_starchy_foods_sweetened_bev
                        QoL_outcomes_p5 = colnames(Q15_measurements)[j],
                        Dietary_features_p5 = "child_ffq_pattern_5_starchy_foods_sweetened_bev",
                        Odds_ratio_p5 = conf_result_p5[1],
                        Standard_error_p5 = summary_result_p5[2],
                        Statistic_p5 = summary_result_p5[3],
                        CI_lower_p5 = conf_result_p5[2],
                        CI_upper_p5 = conf_result_p5[3],
                        P_value_p5 = summary_result_p5[4])}
dim(model_5_patterns)
#5 40 

model_5_patterns_1 <- model_5_patterns %>% select (c(
        QoL_outcomes_p1,
        Dietary_features_p1,
        adj_OR,
        Standard_error_p1,
        Statistic_p1,
        low_ci,
        up_ci,
        P_value_p1)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p1,
                Dietary_features=Dietary_features_p1,
                Odds_ratio=adj_OR,
                Standard_error=Standard_error_p1,
                Statistic=Statistic_p1,
                CI_lower=low_ci,
                CI_upper=up_ci,
                P_value=P_value_p1)

model_5_patterns_2 <- model_5_patterns %>% select (c(
        QoL_outcomes_p2,
        Dietary_features_p2,
        adj_OR.1,
        Standard_error_p2,
        Statistic_p2,
        low_ci.1,
        up_ci.1,
        P_value_p2)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p2,
                Dietary_features=Dietary_features_p2,
                Odds_ratio=adj_OR.1,
                Standard_error=Standard_error_p2,
                Statistic=Statistic_p2,
                CI_lower=low_ci.1,
                CI_upper=up_ci.1,
                P_value=P_value_p2)

model_5_patterns_3 <- model_5_patterns %>% select (c(
        QoL_outcomes_p3,
        Dietary_features_p3,
        adj_OR.2,
        Standard_error_p3,
        Statistic_p3,
        low_ci.2,
        up_ci.2,
        P_value_p3)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p3,
                Dietary_features=Dietary_features_p3,
                Odds_ratio=adj_OR.2,
                Standard_error=Standard_error_p3,
                Statistic=Statistic_p3,
                CI_lower=low_ci.2,
                CI_upper=up_ci.2,
                P_value=P_value_p3)

model_5_patterns_4 <- model_5_patterns %>% select (c(
        QoL_outcomes_p4,
        Dietary_features_p4,
        adj_OR.3,
        Standard_error_p4,
        Statistic_p4,
        low_ci.3,
        up_ci.3,
        P_value_p4)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p4,
                Dietary_features=Dietary_features_p4,
                Odds_ratio=adj_OR.3,
                Standard_error=Standard_error_p4,
                Statistic=Statistic_p4,
                CI_lower=low_ci.3,
                CI_upper=up_ci.3,
                P_value=P_value_p4)

model_5_patterns_5 <- model_5_patterns %>% select (c(
        QoL_outcomes_p5,
        Dietary_features_p5,
        adj_OR.4,
        Standard_error_p5,
        Statistic_p5,
        low_ci.4,
        up_ci.4,
        P_value_p5)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p5,
                Dietary_features=Dietary_features_p5,
                Odds_ratio=adj_OR.4,
                Standard_error=Standard_error_p5,
                Statistic=Statistic_p5,
                CI_lower=low_ci.4,
                CI_upper=up_ci.4,
                P_value=P_value_p5)

model_5_patterns_final <- bind_rows(model_5_patterns_1,
                                    model_5_patterns_2,
                                    model_5_patterns_3,
                                    model_5_patterns_4,
                                    model_5_patterns_5)
rm(model_5_patterns_1, model_5_patterns_2, model_5_patterns_3, model_5_patterns_4, model_5_patterns_5, model_5_patterns)

model_5 <- rbind(model_5_indexes_final, model_5_patterns_final)
model_5 <- model_5 %>%
        mutate(Model="5") %>%
        relocate(Model, .before = QoL_outcomes)

                ####1.7 model 6: model 5 + BMI_SDS ####
#dietary indexes
model_6_indexes <- foreach(j = 1:ncol(Q15_measurements),.combine = rbind) %do% {
                glm1 = glm(Q15_measurements[,j] ~ 
                                   masterfile$child_demo_age_years +
                                   masterfile$child_demo_sex +
                                   masterfile$family_demo_FAS_group +
                                   masterfile$child_demo_immigrant_status +
                                   masterfile$child_demo_school_region +
                                   masterfile$child_anthro_bmi_WHO_sds +
                                   masterfile$child_ffq_margarine + 
                                   masterfile$child_ffq_hpdi_score +
                                   masterfile$child_ffq_animal_score, family='binomial')
                sum1 = summary(glm1)
                conf = or_model_summary(glm1) #package: basecamb
                summary_result = sum1$coef[nrow(sum1$coef),]
                #child_ffq_hpdi_score
                summary_result_hpdi = sum1$coef[nrow(sum1$coef)-1,]
                conf_result_hpdi = conf[nrow(conf)-1,]
                #child_ffq_animal_score
                summary_result_animal = sum1$coef[nrow(sum1$coef),]
                conf_result_animal = conf[nrow(conf),]
                data.frame(
                        #child_ffq_hpdi_score
                        QoL_outcomes_hpdi = colnames(Q15_measurements)[j],
                        Dietary_features_hpdi = "child_ffq_hpdi_score",
                        Odds_ratio_hpdi = conf_result_hpdi[1],
                        Standard_error_hpdi = summary_result_hpdi[2],
                        Statistic_hpdi = summary_result_hpdi[3],
                        CI_lower_hpdi = conf_result_hpdi[2],
                        CI_upper_hpdi = conf_result_hpdi[3],
                        P_value_hpdi = summary_result_hpdi[4],
                        #child_ffq_animal_score
                        QoL_outcomes_animal = colnames(Q15_measurements)[j],
                        Dietary_features_animal = "child_ffq_animal_score",
                        Odds_ratio_animal = conf_result_animal[1],
                        Standard_error_animal = summary_result_animal[2],
                        Statistic_animal = summary_result_animal[3],
                        CI_lower_animal = conf_result_animal[2],
                        CI_upper_animal = conf_result_animal[3],
                        P_value_animal = summary_result_animal[4])}
dim(model_6_indexes)
#5 16 (5 QoL measurements; 8 statistics*2 dietary patterns == 16)

model_6_indexes_hpdi <- model_6_indexes %>% select (c(
        QoL_outcomes_hpdi,
        Dietary_features_hpdi,
        adj_OR,
        Standard_error_hpdi,
        Statistic_hpdi,
        low_ci,
        up_ci,
        P_value_hpdi)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_hpdi,
                Dietary_features=Dietary_features_hpdi,
                Odds_ratio=adj_OR,
                Standard_error=Standard_error_hpdi,
                Statistic=Statistic_hpdi,
                CI_lower=low_ci,
                CI_upper=up_ci,
                P_value=P_value_hpdi)

model_6_indexes_animal <- model_6_indexes %>% select (c(
        QoL_outcomes_animal,
        Dietary_features_animal,
        adj_OR.1,
        Standard_error_animal,
        Statistic_animal,
        low_ci.1,
        up_ci.1,
        P_value_animal)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_animal,
                Dietary_features=Dietary_features_animal,
                Odds_ratio=adj_OR.1,
                Standard_error=Standard_error_animal,
                Statistic=Statistic_animal,
                CI_lower=low_ci.1,
                CI_upper=up_ci.1,
                P_value=P_value_animal)

model_6_indexes_final <- bind_rows(model_6_indexes_hpdi,
                                   model_6_indexes_animal)
rm(model_6_indexes_hpdi, model_6_indexes_animal, model_6_indexes)

#dietary patterns
model_6_patterns <- foreach(j = 1:ncol(Q15_measurements),.combine = rbind) %do% {
                glm1 = glm(Q15_measurements[,j] ~ 
                                   masterfile$child_demo_age_years +
                                   masterfile$child_demo_sex +
                                   masterfile$family_demo_FAS_group +
                                   masterfile$child_demo_immigrant_status +
                                   masterfile$child_demo_school_region +
                                   masterfile$child_anthro_bmi_WHO_sds +
                                   masterfile$child_ffq_pattern_1_meat_seafood_prepd_meals +
                                   masterfile$child_ffq_pattern_2_cooked_veg_grains_legumes +
                                   masterfile$child_ffq_pattern_3_fruits_raw_veg_cheese +
                                   masterfile$child_ffq_pattern_4_confectioneries_pizza +
                                   masterfile$child_ffq_pattern_5_starchy_foods_sweetened_bev, family='binomial')
                sum1 = summary(glm1)
                conf = or_model_summary(glm1) #package: basecamb
                #child_ffq_pattern_1_meat_seafood_prepd_meals
                summary_result_p1 = sum1$coef[nrow(sum1$coef)-4,]
                conf_result_p1 = conf[nrow(conf)-4,]
                #child_ffq_pattern_2_cooked_veg_grains_legumes
                summary_result_p2 = sum1$coef[nrow(sum1$coef)-3,]
                conf_result_p2 = conf[nrow(conf)-3,]
                #child_ffq_pattern_3_fruits_raw_veg_cheese
                summary_result_p3 = sum1$coef[nrow(sum1$coef)-2,]
                conf_result_p3 = conf[nrow(conf)-2,]
                #child_ffq_pattern_4_confectioneries_pizza
                summary_result_p4 = sum1$coef[nrow(sum1$coef)-1,]
                conf_result_p4 = conf[nrow(conf)-1,]
                #child_ffq_pattern_5_starchy_foods_sweetened_bev
                summary_result_p5 = sum1$coef[nrow(sum1$coef),]
                conf_result_p5 = conf[nrow(conf),]
                data.frame(
                        #child_ffq_pattern_1_meat_seafood_prepd_meals
                        QoL_outcomes_p1 = colnames(Q15_measurements)[j],
                        Dietary_features_p1 = "child_ffq_pattern_1_meat_seafood_prepd_meals",
                        Odds_ratio_p1 = conf_result_p1[1],
                        #estimate_p1 = summary_result_p1[1],
                        Standard_error_p1 = summary_result_p1[2],
                        Statistic_p1 = summary_result_p1[3],
                        CI_lower_p1 = conf_result_p1[2],
                        CI_upper_p1 = conf_result_p1[3],
                        P_value_p1 = summary_result_p1[4],
                        #child_ffq_pattern_2_cooked_veg_grains_legumes
                        QoL_outcomes_p2 = colnames(Q15_measurements)[j],
                        Dietary_features_p2 = "child_ffq_pattern_2_cooked_veg_grains_legumes",
                        Odds_ratio_p2 = conf_result_p2[1],
                        Standard_error_p2 = summary_result_p2[2],
                        Statistic_p2 = summary_result_p2[3],
                        CI_lower_p2 = conf_result_p2[2],
                        CI_upper_p2 = conf_result_p2[3],
                        P_value_p2 = summary_result_p2[4],
                        #child_ffq_pattern_3_fruits_raw_veg_cheese
                        QoL_outcomes_p3 = colnames(Q15_measurements)[j],
                        Dietary_features_p3 = "child_ffq_pattern_3_fruits_raw_veg_cheese",
                        Odds_ratio_p3 = conf_result_p3[1],
                        Standard_error_p3 = summary_result_p3[2],
                        Statistic_p3 = summary_result_p3[3],
                        CI_lower_p3 = conf_result_p3[2],
                        CI_upper_p3 = conf_result_p3[3],
                        P_value_p3 = summary_result_p3[4],
                        #child_ffq_pattern_4_confectioneries_pizza
                        QoL_outcomes_p4 = colnames(Q15_measurements)[j],
                        Dietary_features_p4 = "child_ffq_pattern_4_confectioneries_pizza",
                        Odds_ratio_p4 = conf_result_p4[1],
                        Standard_error_p4 = summary_result_p4[2],
                        Statistic_p4 = summary_result_p4[3],
                        CI_lower_p4 = conf_result_p4[2],
                        CI_upper_p4 = conf_result_p4[3],
                        P_value_p4 = summary_result_p4[4],
                        #child_ffq_pattern_5_starchy_foods_sweetened_bev
                        QoL_outcomes_p5 = colnames(Q15_measurements)[j],
                        Dietary_features_p5 = "child_ffq_pattern_5_starchy_foods_sweetened_bev",
                        Odds_ratio_p5 = conf_result_p5[1],
                        Standard_error_p5 = summary_result_p5[2],
                        Statistic_p5 = summary_result_p5[3],
                        CI_lower_p5 = conf_result_p5[2],
                        CI_upper_p5 = conf_result_p5[3],
                        P_value_p5 = summary_result_p5[4])}
dim(model_6_patterns)
#5 40 

model_6_patterns_1 <- model_6_patterns %>% select (c(
        QoL_outcomes_p1,
        Dietary_features_p1,
        adj_OR,
        Standard_error_p1,
        Statistic_p1,
        low_ci,
        up_ci,
        P_value_p1)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p1,
                Dietary_features=Dietary_features_p1,
                Odds_ratio=adj_OR,
                Standard_error=Standard_error_p1,
                Statistic=Statistic_p1,
                CI_lower=low_ci,
                CI_upper=up_ci,
                P_value=P_value_p1)

model_6_patterns_2 <- model_6_patterns %>% select (c(
        QoL_outcomes_p2,
        Dietary_features_p2,
        adj_OR.1,
        Standard_error_p2,
        Statistic_p2,
        low_ci.1,
        up_ci.1,
        P_value_p2)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p2,
                Dietary_features=Dietary_features_p2,
                Odds_ratio=adj_OR.1,
                Standard_error=Standard_error_p2,
                Statistic=Statistic_p2,
                CI_lower=low_ci.1,
                CI_upper=up_ci.1,
                P_value=P_value_p2)

model_6_patterns_3 <- model_6_patterns %>% select (c(
        QoL_outcomes_p3,
        Dietary_features_p3,
        adj_OR.2,
        Standard_error_p3,
        Statistic_p3,
        low_ci.2,
        up_ci.2,
        P_value_p3)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p3,
                Dietary_features=Dietary_features_p3,
                Odds_ratio=adj_OR.2,
                Standard_error=Standard_error_p3,
                Statistic=Statistic_p3,
                CI_lower=low_ci.2,
                CI_upper=up_ci.2,
                P_value=P_value_p3)

model_6_patterns_4 <- model_6_patterns %>% select (c(
        QoL_outcomes_p4,
        Dietary_features_p4,
        adj_OR.3,
        Standard_error_p4,
        Statistic_p4,
        low_ci.3,
        up_ci.3,
        P_value_p4)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p4,
                Dietary_features=Dietary_features_p4,
                Odds_ratio=adj_OR.3,
                Standard_error=Standard_error_p4,
                Statistic=Statistic_p4,
                CI_lower=low_ci.3,
                CI_upper=up_ci.3,
                P_value=P_value_p4)

model_6_patterns_5 <- model_6_patterns %>% select (c(
        QoL_outcomes_p5,
        Dietary_features_p5,
        adj_OR.4,
        Standard_error_p5,
        Statistic_p5,
        low_ci.4,
        up_ci.4,
        P_value_p5)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p5,
                Dietary_features=Dietary_features_p5,
                Odds_ratio=adj_OR.4,
                Standard_error=Standard_error_p5,
                Statistic=Statistic_p5,
                CI_lower=low_ci.4,
                CI_upper=up_ci.4,
                P_value=P_value_p5)

model_6_patterns_final <- bind_rows(model_6_patterns_1,
                                    model_6_patterns_2,
                                    model_6_patterns_3,
                                    model_6_patterns_4,
                                    model_6_patterns_5)
rm(model_6_patterns_1, model_6_patterns_2, model_6_patterns_3, model_6_patterns_4, model_6_patterns_5, model_6_patterns)

model_6 <- rbind(model_6_indexes_final, model_6_patterns_final)
model_6 <- model_6 %>%
        mutate(Model="6") %>%
        relocate(Model, .before = QoL_outcomes)

                ####1.8 model 7: model 5 + predictors with greater importance than dietary features ####

#dietary indexes
model_7_indexes <- foreach(j = 1:ncol(Q15_measurements),.combine = rbind) %do% {
                glm1 = glm(Q15_measurements[,j] ~ 
                                   masterfile$child_demo_age_years +
                                   masterfile$child_demo_sex +
                                   masterfile$family_demo_FAS_group +
                                   masterfile$child_demo_immigrant_status +
                                   masterfile$child_demo_school_region +
                                   masterfile$child_phyact_screen_time_hrs_per_wk +
                                   masterfile$child_demo_school_year + 
                                   masterfile$child_anthro_bmi_WHO_sds +
                                   masterfile$child_ffq_margarine + 
                                   masterfile$child_ffq_hpdi_score +
                                   masterfile$child_ffq_animal_score, family='binomial')
                sum1 = summary(glm1)
                conf = or_model_summary(glm1) #package: basecamb
                #child_ffq_hpdi_score
                summary_result_hpdi = sum1$coef[nrow(sum1$coef)-1,]
                conf_result_hpdi = conf[nrow(conf)-1,]
                #child_ffq_animal_score
                summary_result_animal = sum1$coef[nrow(sum1$coef),]
                conf_result_animal = conf[nrow(conf),]
                data.frame(
                        #child_ffq_hpdi_score
                        QoL_outcomes_hpdi = colnames(Q15_measurements)[j],
                        Dietary_features_hpdi = "child_ffq_hpdi_score",
                        Odds_ratio_hpdi = conf_result_hpdi[1],
                        Standard_error_hpdi = summary_result_hpdi[2],
                        Statistic_hpdi = summary_result_hpdi[3],
                        CI_lower_hpdi = conf_result_hpdi[2],
                        CI_upper_hpdi = conf_result_hpdi[3],
                        P_value_hpdi = summary_result_hpdi[4],
                        #child_ffq_animal_score
                        QoL_outcomes_animal = colnames(Q15_measurements)[j],
                        Dietary_features_animal = "child_ffq_animal_score",
                        Odds_ratio_animal = conf_result_animal[1],
                        Standard_error_animal = summary_result_animal[2],
                        Statistic_animal = summary_result_animal[3],
                        CI_lower_animal = conf_result_animal[2],
                        CI_upper_animal = conf_result_animal[3],
                        P_value_animal = summary_result_animal[4])}
dim(model_7_indexes)
#5 16 (5 QoL measurements; 8 statistics*2 dietary patterns == 16)

model_7_indexes_hpdi <- model_7_indexes %>% select (c(
        QoL_outcomes_hpdi,
        Dietary_features_hpdi,
        adj_OR,
        Standard_error_hpdi,
        Statistic_hpdi,
        low_ci,
        up_ci,
        P_value_hpdi)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_hpdi,
                Dietary_features=Dietary_features_hpdi,
                Odds_ratio=adj_OR,
                Standard_error=Standard_error_hpdi,
                Statistic=Statistic_hpdi,
                CI_lower=low_ci,
                CI_upper=up_ci,
                P_value=P_value_hpdi)

model_7_indexes_animal <- model_7_indexes %>% select (c(
        QoL_outcomes_animal,
        Dietary_features_animal,
        adj_OR.1,
        Standard_error_animal,
        Statistic_animal,
        low_ci.1,
        up_ci.1,
        P_value_animal)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_animal,
                Dietary_features=Dietary_features_animal,
                Odds_ratio=adj_OR.1,
                Standard_error=Standard_error_animal,
                Statistic=Statistic_animal,
                CI_lower=low_ci.1,
                CI_upper=up_ci.1,
                P_value=P_value_animal)

model_7_indexes_final <- bind_rows(model_7_indexes_hpdi,
                                   model_7_indexes_animal)
rm(model_7_indexes_hpdi, model_7_indexes_animal, model_7_indexes)

#dietary patterns
model_7_patterns <- foreach(j = 1:ncol(Q15_measurements),.combine = rbind) %do% {
                glm1 = glm(Q15_measurements[,j] ~ 
                                   masterfile$child_demo_age_years +
                                   masterfile$child_demo_sex +
                                   masterfile$family_demo_FAS_group +
                                   masterfile$child_demo_immigrant_status +
                                   masterfile$child_demo_school_region +
                                   masterfile$child_phyact_screen_time_hrs_per_wk +
                                   masterfile$child_demo_school_year + 
                                   masterfile$child_anthro_bmi_WHO_sds +
                                   masterfile$child_ffq_pattern_1_meat_seafood_prepd_meals +
                                   masterfile$child_ffq_pattern_2_cooked_veg_grains_legumes +
                                   masterfile$child_ffq_pattern_3_fruits_raw_veg_cheese +
                                   masterfile$child_ffq_pattern_4_confectioneries_pizza +
                                   masterfile$child_ffq_pattern_5_starchy_foods_sweetened_bev, family='binomial')
                sum1 = summary(glm1)
                conf = or_model_summary(glm1) #package: basecamb
                #child_ffq_pattern_1_meat_seafood_prepd_meals
                summary_result_p1 = sum1$coef[nrow(sum1$coef)-4,]
                conf_result_p1 = conf[nrow(conf)-4,]
                #child_ffq_pattern_2_cooked_veg_grains_legumes
                summary_result_p2 = sum1$coef[nrow(sum1$coef)-3,]
                conf_result_p2 = conf[nrow(conf)-3,]
                #child_ffq_pattern_3_fruits_raw_veg_cheese
                summary_result_p3 = sum1$coef[nrow(sum1$coef)-2,]
                conf_result_p3 = conf[nrow(conf)-2,]
                #child_ffq_pattern_4_confectioneries_pizza
                summary_result_p4 = sum1$coef[nrow(sum1$coef)-1,]
                conf_result_p4 = conf[nrow(conf)-1,]
                #child_ffq_pattern_5_starchy_foods_sweetened_bev
                summary_result_p5 = sum1$coef[nrow(sum1$coef),]
                conf_result_p5 = conf[nrow(conf),]
                data.frame(
                        #child_ffq_pattern_1_meat_seafood_prepd_meals
                        QoL_outcomes_p1 = colnames(Q15_measurements)[j],
                        Dietary_features_p1 = "child_ffq_pattern_1_meat_seafood_prepd_meals",
                        Odds_ratio_p1 = conf_result_p1[1],
                        #estimate_p1 = summary_result_p1[1],
                        Standard_error_p1 = summary_result_p1[2],
                        Statistic_p1 = summary_result_p1[3],
                        CI_lower_p1 = conf_result_p1[2],
                        CI_upper_p1 = conf_result_p1[3],
                        P_value_p1 = summary_result_p1[4],
                        #child_ffq_pattern_2_cooked_veg_grains_legumes
                        QoL_outcomes_p2 = colnames(Q15_measurements)[j],
                        Dietary_features_p2 = "child_ffq_pattern_2_cooked_veg_grains_legumes",
                        Odds_ratio_p2 = conf_result_p2[1],
                        Standard_error_p2 = summary_result_p2[2],
                        Statistic_p2 = summary_result_p2[3],
                        CI_lower_p2 = conf_result_p2[2],
                        CI_upper_p2 = conf_result_p2[3],
                        P_value_p2 = summary_result_p2[4],
                        #child_ffq_pattern_3_fruits_raw_veg_cheese
                        QoL_outcomes_p3 = colnames(Q15_measurements)[j],
                        Dietary_features_p3 = "child_ffq_pattern_3_fruits_raw_veg_cheese",
                        Odds_ratio_p3 = conf_result_p3[1],
                        Standard_error_p3 = summary_result_p3[2],
                        Statistic_p3 = summary_result_p3[3],
                        CI_lower_p3 = conf_result_p3[2],
                        CI_upper_p3 = conf_result_p3[3],
                        P_value_p3 = summary_result_p3[4],
                        #child_ffq_pattern_4_confectioneries_pizza
                        QoL_outcomes_p4 = colnames(Q15_measurements)[j],
                        Dietary_features_p4 = "child_ffq_pattern_4_confectioneries_pizza",
                        Odds_ratio_p4 = conf_result_p4[1],
                        Standard_error_p4 = summary_result_p4[2],
                        Statistic_p4 = summary_result_p4[3],
                        CI_lower_p4 = conf_result_p4[2],
                        CI_upper_p4 = conf_result_p4[3],
                        P_value_p4 = summary_result_p4[4],
                        #child_ffq_pattern_5_starchy_foods_sweetened_bev
                        QoL_outcomes_p5 = colnames(Q15_measurements)[j],
                        Dietary_features_p5 = "child_ffq_pattern_5_starchy_foods_sweetened_bev",
                        Odds_ratio_p5 = conf_result_p5[1],
                        Standard_error_p5 = summary_result_p5[2],
                        Statistic_p5 = summary_result_p5[3],
                        CI_lower_p5 = conf_result_p5[2],
                        CI_upper_p5 = conf_result_p5[3],
                        P_value_p5 = summary_result_p5[4])}
dim(model_7_patterns)
#5 40 

model_7_patterns_1 <- model_7_patterns %>% select (c(
        QoL_outcomes_p1,
        Dietary_features_p1,
        adj_OR,
        Standard_error_p1,
        Statistic_p1,
        low_ci,
        up_ci,
        P_value_p1)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p1,
                Dietary_features=Dietary_features_p1,
                Odds_ratio=adj_OR,
                Standard_error=Standard_error_p1,
                Statistic=Statistic_p1,
                CI_lower=low_ci,
                CI_upper=up_ci,
                P_value=P_value_p1)

model_7_patterns_2 <- model_7_patterns %>% select (c(
        QoL_outcomes_p2,
        Dietary_features_p2,
        adj_OR.1,
        Standard_error_p2,
        Statistic_p2,
        low_ci.1,
        up_ci.1,
        P_value_p2)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p2,
                Dietary_features=Dietary_features_p2,
                Odds_ratio=adj_OR.1,
                Standard_error=Standard_error_p2,
                Statistic=Statistic_p2,
                CI_lower=low_ci.1,
                CI_upper=up_ci.1,
                P_value=P_value_p2)

model_7_patterns_3 <- model_7_patterns %>% select (c(
        QoL_outcomes_p3,
        Dietary_features_p3,
        adj_OR.2,
        Standard_error_p3,
        Statistic_p3,
        low_ci.2,
        up_ci.2,
        P_value_p3)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p3,
                Dietary_features=Dietary_features_p3,
                Odds_ratio=adj_OR.2,
                Standard_error=Standard_error_p3,
                Statistic=Statistic_p3,
                CI_lower=low_ci.2,
                CI_upper=up_ci.2,
                P_value=P_value_p3)

model_7_patterns_4 <- model_7_patterns %>% select (c(
        QoL_outcomes_p4,
        Dietary_features_p4,
        adj_OR.3,
        Standard_error_p4,
        Statistic_p4,
        low_ci.3,
        up_ci.3,
        P_value_p4)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p4,
                Dietary_features=Dietary_features_p4,
                Odds_ratio=adj_OR.3,
                Standard_error=Standard_error_p4,
                Statistic=Statistic_p4,
                CI_lower=low_ci.3,
                CI_upper=up_ci.3,
                P_value=P_value_p4)

model_7_patterns_5 <- model_7_patterns %>% select (c(
        QoL_outcomes_p5,
        Dietary_features_p5,
        adj_OR.4,
        Standard_error_p5,
        Statistic_p5,
        low_ci.4,
        up_ci.4,
        P_value_p5)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p5,
                Dietary_features=Dietary_features_p5,
                Odds_ratio=adj_OR.4,
                Standard_error=Standard_error_p5,
                Statistic=Statistic_p5,
                CI_lower=low_ci.4,
                CI_upper=up_ci.4,
                P_value=P_value_p5)

model_7_patterns_final <- bind_rows(model_7_patterns_1,
                                    model_7_patterns_2,
                                    model_7_patterns_3,
                                    model_7_patterns_4,
                                    model_7_patterns_5)
rm(model_7_patterns_1, model_7_patterns_2, model_7_patterns_3, model_7_patterns_4, model_7_patterns_5, model_7_patterns)

model_7 <- rbind(model_7_indexes_final, model_7_patterns_final)
model_7 <- model_7 %>%
        mutate(Model="7") %>%
        relocate(Model, .before = QoL_outcomes)

                ####1.9 bind all models to ensure model did not collapse due to over-correction ####
qol_diet_log <- bind_rows(model_1, model_2, model_3, model_4, model_5, model_6, model_7)
str(qol_diet_log)
#245  9
qol_diet_log <- qol_diet_log %>% mutate_if(is.character, as.factor)
qol_diet_log <- qol_diet_log %>%
        mutate(QoL_outcomes = fct_relevel(QoL_outcomes, "child_pedsql_health_related_qol_score_total_catg_Q15", 
                                          "child_pedsql_phy_score_catg_Q15", 
                                          "child_pedsql_emo_score_catg_Q15",
                                          "child_pedsql_soc_score_catg_Q15",
                                          "child_pedsql_sch_score_catg_Q15")) %>%
        mutate(Dietary_features = fct_relevel(Dietary_features, "child_ffq_hpdi_score", 
                                              "child_ffq_animal_score", 
                                              "child_ffq_pattern_1_meat_seafood_prepd_meals",
                                              "child_ffq_pattern_2_cooked_veg_grains_legumes",
                                              "child_ffq_pattern_3_fruits_raw_veg_cheese",
                                              "child_ffq_pattern_4_confectioneries_pizza",
                                              "child_ffq_pattern_5_starchy_foods_sweetened_bev"))

#order by QoL_outcomes and 'Dietary_features'
qol_diet_log <- qol_diet_log[order(qol_diet_log$QoL_outcomes, qol_diet_log$Dietary_features), ]

setwd("~/OneDrive - UMCG/Prolepsis_internship/merged_2015_2018/masterfile/4_logistic_regression/")
setwd("C:/Users/Siobhan Brushett/OneDrive - UMCG/Prolepsis_internship/merged_2015_2018/masterfile/4_logistic_regression/")

write.table(qol_diet_log, "2023_07_05_qol_diet_log_all_models_updated.txt", sep="\t", row.names=F, quote = F) 

#models did not collapse from over-correction or collinearity
#continue with models: 1, 2, 4 and 7


                #### =============== 2. MODELS FOR MANUSCRIPT - FDR CORRECTION   ================ ####
final_log_models <- bind_rows(model_1, #crude model
                              model_2, #standard correction for age and sex
                              model_4, #anthropometric model
                              model_7) #anthropometric and prediction model
str(final_log_models)
final_log_models <- final_log_models %>% mutate_if(is.character, as.factor)

final_log_models <- final_log_models %>%
        mutate(QoL_outcomes = fct_relevel(QoL_outcomes, "child_pedsql_health_related_qol_score_total_catg_Q15", 
                                          "child_pedsql_phy_score_catg_Q15", 
                                          "child_pedsql_emo_score_catg_Q15",
                                          "child_pedsql_soc_score_catg_Q15",
                                          "child_pedsql_sch_score_catg_Q15")) %>%
        mutate(Dietary_features = fct_relevel(Dietary_features, "child_ffq_hpdi_score", 
                                              "child_ffq_animal_score", 
                                              "child_ffq_pattern_1_meat_seafood_prepd_meals",
                                              "child_ffq_pattern_2_cooked_veg_grains_legumes",
                                              "child_ffq_pattern_3_fruits_raw_veg_cheese",
                                              "child_ffq_pattern_4_confectioneries_pizza",
                                              "child_ffq_pattern_5_starchy_foods_sweetened_bev"))

final_log_models$FDR = p.adjust(final_log_models$P_value, method = "BH")
sum(final_log_models$FDR<0.05)
#65
sum(final_log_models$P_value<0.05)
#70

final_log_models <- final_log_models[order(final_log_models$Model, final_log_models$Dietary_features), ]

setwd("~/OneDrive - UMCG/Prolepsis_internship/merged_2015_2018/masterfile/4_logistic_regression/")
setwd("C:/Users/Siobhan Brushett//OneDrive - UMCG/Prolepsis_internship/merged_2015_2018/masterfile/4_logistic_regression/")
str(final_log_models)
write.table(final_log_models, "2023_07_05_qol_diet_log_manuscript_models_FDR_updated.txt", sep="\t", row.names=F, quote = F) 
#140 10

