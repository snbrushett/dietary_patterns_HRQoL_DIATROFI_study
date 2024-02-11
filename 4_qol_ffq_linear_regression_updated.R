###     SCRIPT: LINEAR REGRESSION OF QOL (CONTINOUS PERCENTILE) AND DIETARY FEATURES 
###     AUTHOR(S): SIOBHAN BRUSHETT 
###     DESCRIPTION: FITTING LINEAR REGRESSION OF CONTINOUS PERCENTILE HRQOL OUTCOMES AND TOTAL DIETARY INDEX - SENSITIVITY ANALYSIS
###     PROJECT: DIATROFI_BASELINE_2015_2018
###     NOTE(S): ONLY MODELS SELECTED FOR MANUSCRIPT (AS DEPICTED IN LOGISTIC REGRESSION WERE SELECTED FOR LINEAR REGRESSION)

#libraries
library(tidyverse)
library(foreach)

## 0. LOAD DATA
## 1. FIT LINEAR MODELS
## 2. MODELS FOR MANUSCRIPT - FDR CORRECTION

                #### =============== 0. LOAD DATA ================ ####
masterfile <- read.delim("2023_05_09_diatrofi_2015_2018_masterfile_updated.txt", sep = "\t", header = T, stringsAsFactors = T)
#6583 66
str(masterfile)
masterfile$child_demo_school_code <- as.factor(masterfile$child_demo_school_code)

                #### =============== 1. FIT LINEAR MODELS   ================ ####
#NOTE: for dietary indexes AND NOT DIETARY PATTERNS, models need to be corrected for margarine intake as it was NOT included in the indexes (as in previous studies)

                ####1.1 prepare vectors for linear models ####
#prepare different HRQoL measurements
QoL_perc_measurements <- masterfile [, c("pseudo_ids", 
                                         "child_pedsql_health_related_qol_score_total_perc",
                                         #"child_pedsql_psysoc_score_perc",
                                         "child_pedsql_phy_score_perc",
                                         "child_pedsql_emo_score_perc",
                                         "child_pedsql_soc_score_perc",
                                         "child_pedsql_sch_score_perc")]
rownames(QoL_perc_measurements) <- masterfile$pseudo_ids
head(QoL_perc_measurements)
QoL_perc_measurements$pseudo_ids <- NULL

                ####1.2 model 1: crude model of QoL and diet ####
#dietary indexes
model_1_indexes <- foreach(j = 1:ncol(QoL_perc_measurements),.combine = rbind) %do% {
                lm1 = lm(QoL_perc_measurements[,j] ~ 
                                 masterfile$child_ffq_margarine + 
                                 masterfile$child_ffq_hpdi_score +
                                 masterfile$child_ffq_animal_score)
                sum1 = summary(lm1)
                #child_ffq_hpdi_score
                summary_result_hpdi = sum1$coef[nrow(sum1$coef)-1,]
                #child_ffq_animal_score
                summary_result_animal = sum1$coef[nrow(sum1$coef),]
                data.frame(#child_ffq_hpdi_score
                           QoL_outcomes_hpdi = colnames(QoL_perc_measurements)[j],
                           Dietary_features_hpdi = "child_ffq_hpdi_score",
                           Beta_coefficient_hpdi = summary_result_hpdi[1],
                           Standard_error_hpdi = summary_result_hpdi[2],
                           Statistic_hpdi = summary_result_hpdi[3],
                           P_value_hpdi = summary_result_hpdi[4],
                           #child_ffq_animal_score
                           QoL_outcomes_animal = colnames(QoL_perc_measurements)[j],
                           Dietary_features_animal = "child_ffq_animal_score",
                           Beta_coefficient_animal = summary_result_animal[1],
                           Standard_error_animal = summary_result_animal[2],
                           Statistic_animal = summary_result_animal[3],
                           P_value_animal = summary_result_animal[4])
        }
dim(model_1_indexes)
#5 12 

model_1_indexes_hpdi <- model_1_indexes %>% select (c(
        QoL_outcomes_hpdi,
        Dietary_features_hpdi,
        Beta_coefficient_hpdi,
        Standard_error_hpdi,
        Statistic_hpdi,
        P_value_hpdi)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_hpdi,
                Dietary_features=Dietary_features_hpdi,
                Beta_coefficient=Beta_coefficient_hpdi,
                Standard_error=Standard_error_hpdi,
                Statistic=Statistic_hpdi,
                P_value=P_value_hpdi)

model_1_indexes_animal <- model_1_indexes %>% select (c(
        QoL_outcomes_animal,
        Dietary_features_animal,
        Beta_coefficient_animal,
        Standard_error_animal,
        Statistic_animal,
        P_value_animal)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_animal,
                Dietary_features=Dietary_features_animal,
                Beta_coefficient=Beta_coefficient_animal,
                Standard_error=Standard_error_animal,
                Statistic=Statistic_animal,
                P_value=P_value_animal)

model_1_indexes_final <- bind_rows(model_1_indexes_hpdi,
                                   model_1_indexes_animal)
rm(model_1_indexes_hpdi, model_1_indexes_animal, model_1_indexes)

#dietary patterns
model_1_patterns <- foreach(j = 1:ncol(QoL_perc_measurements),.combine = rbind) %do% {
                lm1 = lm(QoL_perc_measurements[,j] ~ 
                                 masterfile$child_ffq_pattern_1_meat_seafood_prepd_meals +
                                 masterfile$child_ffq_pattern_2_cooked_veg_grains_legumes +
                                 masterfile$child_ffq_pattern_3_fruits_raw_veg_cheese +
                                 masterfile$child_ffq_pattern_4_confectioneries_pizza +
                                 masterfile$child_ffq_pattern_5_starchy_foods_sweetened_bev)
                sum1 = summary(lm1)
                #child_ffq_pattern_1_meat_seafood_prepd_meals
                summary_result_p1 = sum1$coef[nrow(sum1$coef)-4,]
                #child_ffq_pattern_2_cooked_veg_grains_legumes
                summary_result_p2 = sum1$coef[nrow(sum1$coef)-3,]
                #child_ffq_pattern_3_fruits_raw_veg_cheese
                summary_result_p3 = sum1$coef[nrow(sum1$coef)-2,]
                #child_ffq_pattern_4_confectioneries_pizza
                summary_result_p4 = sum1$coef[nrow(sum1$coef)-1,]
                #child_ffq_pattern_5_starchy_foods_sweetened_bev
                summary_result_p5 = sum1$coef[nrow(sum1$coef),]
                data.frame(
                        #child_ffq_pattern_1_meat_seafood_prepd_meals
                        QoL_outcomes_p1 = colnames(QoL_perc_measurements)[j],
                           Dietary_features_p1 = "child_ffq_pattern_1_meat_seafood_prepd_meals",
                           Beta_coefficient_p1 = summary_result_p1[1],
                           Standard_error_p1 = summary_result_p1[2],
                           Statistic_p1 = summary_result_p1[3],
                           P_value_p1 = summary_result_p1[4],
                        #child_ffq_pattern_2_cooked_veg_grains_legumes
                        QoL_outcomes_p2 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p2 = "child_ffq_pattern_2_cooked_veg_grains_legumes",
                        Beta_coefficient_p2 = summary_result_p2[1],
                        Standard_error_p2 = summary_result_p2[2],
                        Statistic_p2 = summary_result_p2[3],
                        P_value_p2 = summary_result_p2[4],
                        #child_ffq_pattern_3_fruits_raw_veg_cheese
                        QoL_outcomes_p3 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p3 = "child_ffq_pattern_3_fruits_raw_veg_cheese",
                        Beta_coefficient_p3 = summary_result_p3[1],
                        Standard_error_p3 = summary_result_p3[2],
                        Statistic_p3 = summary_result_p3[3],
                        P_value_p3 = summary_result_p3[4],
                        #child_ffq_pattern_4_confectioneries_pizza
                        QoL_outcomes_p4 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p4 = "child_ffq_pattern_4_confectioneries_pizza",
                        Beta_coefficient_p4 = summary_result_p4[1],
                        Standard_error_p4 = summary_result_p4[2],
                        Statistic_p4 = summary_result_p4[3],
                        P_value_p4 = summary_result_p4[4],
                        #child_ffq_pattern_5_starchy_foods_sweetened_bev
                        QoL_outcomes_p5 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p5 = "child_ffq_pattern_5_starchy_foods_sweetened_bev",
                        Beta_coefficient_p5 = summary_result_p5[1],
                        Standard_error_p5 = summary_result_p5[2],
                        Statistic_p5 = summary_result_p5[3],
                        P_value_p5 = summary_result_p5[4])
        }
dim(model_1_patterns)
#5 30

model_1_patterns_1 <- model_1_patterns %>% select (c(
        QoL_outcomes_p1,
        Dietary_features_p1,
        Beta_coefficient_p1,
        Standard_error_p1,
        Statistic_p1,
        P_value_p1)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p1,
                Dietary_features=Dietary_features_p1,
                Beta_coefficient=Beta_coefficient_p1,
                Standard_error=Standard_error_p1,
                Statistic=Statistic_p1,
                P_value=P_value_p1)

model_1_patterns_2 <- model_1_patterns %>% select (c(
        QoL_outcomes_p2,
        Dietary_features_p2,
        Beta_coefficient_p2,
        Standard_error_p2,
        Statistic_p2,
        P_value_p2)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p2,
                Dietary_features=Dietary_features_p2,
                Beta_coefficient=Beta_coefficient_p2,
                Standard_error=Standard_error_p2,
                Statistic=Statistic_p2,
                P_value=P_value_p2)

model_1_patterns_3 <- model_1_patterns %>% select (c(
        QoL_outcomes_p3,
        Dietary_features_p3,
        Beta_coefficient_p3,
        Standard_error_p3,
        Statistic_p3,
        P_value_p3)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p3,
                Dietary_features=Dietary_features_p3,
                Beta_coefficient=Beta_coefficient_p3,
                Standard_error=Standard_error_p3,
                Statistic=Statistic_p3,
                P_value=P_value_p3)

model_1_patterns_4 <- model_1_patterns %>% select (c(
        QoL_outcomes_p4,
        Dietary_features_p4,
        Beta_coefficient_p4,
        Standard_error_p4,
        Statistic_p4,
        P_value_p4)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p4,
                Dietary_features=Dietary_features_p4,
                Beta_coefficient=Beta_coefficient_p4,
                Standard_error=Standard_error_p4,
                Statistic=Statistic_p4,
                P_value=P_value_p4)

model_1_patterns_5 <- model_1_patterns %>% select (c(
        QoL_outcomes_p5,
        Dietary_features_p5,
        Beta_coefficient_p5,
        Standard_error_p5,
        Statistic_p5,
        P_value_p5)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p5,
                Dietary_features=Dietary_features_p5,
                Beta_coefficient=Beta_coefficient_p5,
                Standard_error=Standard_error_p5,
                Statistic=Statistic_p5,
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
model_2_indexes <- foreach(j = 1:ncol(QoL_perc_measurements),.combine = rbind) %do% {
                lm1 = lm(QoL_perc_measurements[,j] ~ 
                                 masterfile$child_demo_age_years +
                                 masterfile$child_demo_sex +
                                 masterfile$child_ffq_margarine + 
                                 masterfile$child_ffq_hpdi_score +
                                 masterfile$child_ffq_animal_score)
                sum1 = summary(lm1)
                #child_ffq_hpdi_score
                summary_result_hpdi = sum1$coef[nrow(sum1$coef)-1,]
                #child_ffq_animal_score
                summary_result_animal = sum1$coef[nrow(sum1$coef),]
                data.frame(#child_ffq_hpdi_score
                        QoL_outcomes_hpdi = colnames(QoL_perc_measurements)[j],
                        Dietary_features_hpdi = "child_ffq_hpdi_score",
                        Beta_coefficient_hpdi = summary_result_hpdi[1],
                        Standard_error_hpdi = summary_result_hpdi[2],
                        Statistic_hpdi = summary_result_hpdi[3],
                        P_value_hpdi = summary_result_hpdi[4],
                        #child_ffq_animal_score
                        QoL_outcomes_animal = colnames(QoL_perc_measurements)[j],
                        Dietary_features_animal = "child_ffq_animal_score",
                        Beta_coefficient_animal = summary_result_animal[1],
                        Standard_error_animal = summary_result_animal[2],
                        Statistic_animal = summary_result_animal[3],
                        P_value_animal = summary_result_animal[4])
}
dim(model_2_indexes)
#5 12 

model_2_indexes_hpdi <- model_2_indexes %>% select (c(
        QoL_outcomes_hpdi,
        Dietary_features_hpdi,
        Beta_coefficient_hpdi,
        Standard_error_hpdi,
        Statistic_hpdi,
        P_value_hpdi)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_hpdi,
                Dietary_features=Dietary_features_hpdi,
                Beta_coefficient=Beta_coefficient_hpdi,
                Standard_error=Standard_error_hpdi,
                Statistic=Statistic_hpdi,
                P_value=P_value_hpdi)

model_2_indexes_animal <- model_2_indexes %>% select (c(
        QoL_outcomes_animal,
        Dietary_features_animal,
        Beta_coefficient_animal,
        Standard_error_animal,
        Statistic_animal,
        P_value_animal)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_animal,
                Dietary_features=Dietary_features_animal,
                Beta_coefficient=Beta_coefficient_animal,
                Standard_error=Standard_error_animal,
                Statistic=Statistic_animal,
                P_value=P_value_animal)

model_2_indexes_final <- bind_rows(model_2_indexes_hpdi,
                                   model_2_indexes_animal)
rm(model_2_indexes_hpdi, model_2_indexes_animal, model_2_indexes)

#dietary patterns
model_2_patterns <- foreach(j = 1:ncol(QoL_perc_measurements),.combine = rbind) %do% {
                lm1 = lm(QoL_perc_measurements[,j] ~ 
                                 masterfile$child_demo_age_years +
                                 masterfile$child_demo_sex +
                                 masterfile$child_ffq_pattern_1_meat_seafood_prepd_meals +
                                 masterfile$child_ffq_pattern_2_cooked_veg_grains_legumes +
                                 masterfile$child_ffq_pattern_3_fruits_raw_veg_cheese +
                                 masterfile$child_ffq_pattern_4_confectioneries_pizza +
                                 masterfile$child_ffq_pattern_5_starchy_foods_sweetened_bev)
                sum1 = summary(lm1)
                #child_ffq_pattern_1_meat_seafood_prepd_meals
                summary_result_p1 = sum1$coef[nrow(sum1$coef)-4,]
                #child_ffq_pattern_2_cooked_veg_grains_legumes
                summary_result_p2 = sum1$coef[nrow(sum1$coef)-3,]
                #child_ffq_pattern_3_fruits_raw_veg_cheese
                summary_result_p3 = sum1$coef[nrow(sum1$coef)-2,]
                #child_ffq_pattern_4_confectioneries_pizza
                summary_result_p4 = sum1$coef[nrow(sum1$coef)-1,]
                #child_ffq_pattern_5_starchy_foods_sweetened_bev
                summary_result_p5 = sum1$coef[nrow(sum1$coef),]
                data.frame(
                        #child_ffq_pattern_1_meat_seafood_prepd_meals
                        QoL_outcomes_p1 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p1 = "child_ffq_pattern_1_meat_seafood_prepd_meals",
                        Beta_coefficient_p1 = summary_result_p1[1],
                        Standard_error_p1 = summary_result_p1[2],
                        Statistic_p1 = summary_result_p1[3],
                        P_value_p1 = summary_result_p1[4],
                        #child_ffq_pattern_2_cooked_veg_grains_legumes
                        QoL_outcomes_p2 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p2 = "child_ffq_pattern_2_cooked_veg_grains_legumes",
                        Beta_coefficient_p2 = summary_result_p2[1],
                        Standard_error_p2 = summary_result_p2[2],
                        Statistic_p2 = summary_result_p2[3],
                        P_value_p2 = summary_result_p2[4],
                        #child_ffq_pattern_3_fruits_raw_veg_cheese
                        QoL_outcomes_p3 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p3 = "child_ffq_pattern_3_fruits_raw_veg_cheese",
                        Beta_coefficient_p3 = summary_result_p3[1],
                        Standard_error_p3 = summary_result_p3[2],
                        Statistic_p3 = summary_result_p3[3],
                        P_value_p3 = summary_result_p3[4],
                        #child_ffq_pattern_4_confectioneries_pizza
                        QoL_outcomes_p4 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p4 = "child_ffq_pattern_4_confectioneries_pizza",
                        Beta_coefficient_p4 = summary_result_p4[1],
                        Standard_error_p4 = summary_result_p4[2],
                        Statistic_p4 = summary_result_p4[3],
                        P_value_p4 = summary_result_p4[4],
                        #child_ffq_pattern_5_starchy_foods_sweetened_bev
                        QoL_outcomes_p5 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p5 = "child_ffq_pattern_5_starchy_foods_sweetened_bev",
                        Beta_coefficient_p5 = summary_result_p5[1],
                        Standard_error_p5 = summary_result_p5[2],
                        Statistic_p5 = summary_result_p5[3],
                        P_value_p5 = summary_result_p5[4])
}
dim(model_2_patterns)
#5 30

model_2_patterns_1 <- model_2_patterns %>% select (c(
        QoL_outcomes_p1,
        Dietary_features_p1,
        Beta_coefficient_p1,
        Standard_error_p1,
        Statistic_p1,
        P_value_p1)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p1,
                Dietary_features=Dietary_features_p1,
                Beta_coefficient=Beta_coefficient_p1,
                Standard_error=Standard_error_p1,
                Statistic=Statistic_p1,
                P_value=P_value_p1)

model_2_patterns_2 <- model_2_patterns %>% select (c(
        QoL_outcomes_p2,
        Dietary_features_p2,
        Beta_coefficient_p2,
        Standard_error_p2,
        Statistic_p2,
        P_value_p2)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p2,
                Dietary_features=Dietary_features_p2,
                Beta_coefficient=Beta_coefficient_p2,
                Standard_error=Standard_error_p2,
                Statistic=Statistic_p2,
                P_value=P_value_p2)

model_2_patterns_3 <- model_2_patterns %>% select (c(
        QoL_outcomes_p3,
        Dietary_features_p3,
        Beta_coefficient_p3,
        Standard_error_p3,
        Statistic_p3,
        P_value_p3)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p3,
                Dietary_features=Dietary_features_p3,
                Beta_coefficient=Beta_coefficient_p3,
                Standard_error=Standard_error_p3,
                Statistic=Statistic_p3,
                P_value=P_value_p3)

model_2_patterns_4 <- model_2_patterns %>% select (c(
        QoL_outcomes_p4,
        Dietary_features_p4,
        Beta_coefficient_p4,
        Standard_error_p4,
        Statistic_p4,
        P_value_p4)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p4,
                Dietary_features=Dietary_features_p4,
                Beta_coefficient=Beta_coefficient_p4,
                Standard_error=Standard_error_p4,
                Statistic=Statistic_p4,
                P_value=P_value_p4)

model_2_patterns_5 <- model_2_patterns %>% select (c(
        QoL_outcomes_p5,
        Dietary_features_p5,
        Beta_coefficient_p5,
        Standard_error_p5,
        Statistic_p5,
        P_value_p5)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p5,
                Dietary_features=Dietary_features_p5,
                Beta_coefficient=Beta_coefficient_p5,
                Standard_error=Standard_error_p5,
                Statistic=Statistic_p5,
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

                ####1.4 model 4: model 1 + SES + student_immigrant_status ####
#dietary indexes
model_4_indexes <- foreach(j = 1:ncol(QoL_perc_measurements),.combine = rbind) %do% {
                lm1 = lm(QoL_perc_measurements[,j] ~ 
                                 masterfile$child_demo_age_years +
                                 masterfile$child_demo_sex +
                                 masterfile$family_demo_FAS_group +
                                 masterfile$child_demo_immigrant_status +
                                 masterfile$child_ffq_margarine + 
                                 masterfile$child_ffq_hpdi_score +
                                 masterfile$child_ffq_animal_score)
                sum1 = summary(lm1)
                #child_ffq_hpdi_score
                summary_result_hpdi = sum1$coef[nrow(sum1$coef)-1,]
                #child_ffq_animal_score
                summary_result_animal = sum1$coef[nrow(sum1$coef),]
                data.frame(#child_ffq_hpdi_score
                        QoL_outcomes_hpdi = colnames(QoL_perc_measurements)[j],
                        Dietary_features_hpdi = "child_ffq_hpdi_score",
                        Beta_coefficient_hpdi = summary_result_hpdi[1],
                        Standard_error_hpdi = summary_result_hpdi[2],
                        Statistic_hpdi = summary_result_hpdi[3],
                        P_value_hpdi = summary_result_hpdi[4],
                        #child_ffq_animal_score
                        QoL_outcomes_animal = colnames(QoL_perc_measurements)[j],
                        Dietary_features_animal = "child_ffq_animal_score",
                        Beta_coefficient_animal = summary_result_animal[1],
                        Standard_error_animal = summary_result_animal[2],
                        Statistic_animal = summary_result_animal[3],
                        P_value_animal = summary_result_animal[4])
}
dim(model_4_indexes)
#5 12 

model_4_indexes_hpdi <- model_4_indexes %>% select (c(
        QoL_outcomes_hpdi,
        Dietary_features_hpdi,
        Beta_coefficient_hpdi,
        Standard_error_hpdi,
        Statistic_hpdi,
        P_value_hpdi)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_hpdi,
                Dietary_features=Dietary_features_hpdi,
                Beta_coefficient=Beta_coefficient_hpdi,
                Standard_error=Standard_error_hpdi,
                Statistic=Statistic_hpdi,
                P_value=P_value_hpdi)

model_4_indexes_animal <- model_4_indexes %>% select (c(
        QoL_outcomes_animal,
        Dietary_features_animal,
        Beta_coefficient_animal,
        Standard_error_animal,
        Statistic_animal,
        P_value_animal)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_animal,
                Dietary_features=Dietary_features_animal,
                Beta_coefficient=Beta_coefficient_animal,
                Standard_error=Standard_error_animal,
                Statistic=Statistic_animal,
                P_value=P_value_animal)

model_4_indexes_final <- bind_rows(model_4_indexes_hpdi,
                                   model_4_indexes_animal)
rm(model_4_indexes_hpdi, model_4_indexes_animal, model_4_indexes)

#dietary patterns
model_4_patterns <- foreach(j = 1:ncol(QoL_perc_measurements),.combine = rbind) %do% {
                lm1 = lm(QoL_perc_measurements[,j] ~ 
                                 masterfile$child_demo_age_years +
                                 masterfile$child_demo_sex +
                                 masterfile$family_demo_FAS_group +
                                 masterfile$child_demo_immigrant_status +
                                 masterfile$child_demo_school_region +
                                 masterfile$child_ffq_pattern_1_meat_seafood_prepd_meals +
                                 masterfile$child_ffq_pattern_2_cooked_veg_grains_legumes +
                                 masterfile$child_ffq_pattern_3_fruits_raw_veg_cheese +
                                 masterfile$child_ffq_pattern_4_confectioneries_pizza +
                                 masterfile$child_ffq_pattern_5_starchy_foods_sweetened_bev)
                sum1 = summary(lm1)
                #child_ffq_pattern_1_meat_seafood_prepd_meals
                summary_result_p1 = sum1$coef[nrow(sum1$coef)-4,]
                #child_ffq_pattern_2_cooked_veg_grains_legumes
                summary_result_p2 = sum1$coef[nrow(sum1$coef)-3,]
                #child_ffq_pattern_3_fruits_raw_veg_cheese
                summary_result_p3 = sum1$coef[nrow(sum1$coef)-2,]
                #child_ffq_pattern_4_confectioneries_pizza
                summary_result_p4 = sum1$coef[nrow(sum1$coef)-1,]
                #child_ffq_pattern_5_starchy_foods_sweetened_bev
                summary_result_p5 = sum1$coef[nrow(sum1$coef),]
                data.frame(
                        #child_ffq_pattern_1_meat_seafood_prepd_meals
                        QoL_outcomes_p1 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p1 = "child_ffq_pattern_1_meat_seafood_prepd_meals",
                        Beta_coefficient_p1 = summary_result_p1[1],
                        Standard_error_p1 = summary_result_p1[2],
                        Statistic_p1 = summary_result_p1[3],
                        P_value_p1 = summary_result_p1[4],
                        #child_ffq_pattern_2_cooked_veg_grains_legumes
                        QoL_outcomes_p2 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p2 = "child_ffq_pattern_2_cooked_veg_grains_legumes",
                        Beta_coefficient_p2 = summary_result_p2[1],
                        Standard_error_p2 = summary_result_p2[2],
                        Statistic_p2 = summary_result_p2[3],
                        P_value_p2 = summary_result_p2[4],
                        #child_ffq_pattern_3_fruits_raw_veg_cheese
                        QoL_outcomes_p3 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p3 = "child_ffq_pattern_3_fruits_raw_veg_cheese",
                        Beta_coefficient_p3 = summary_result_p3[1],
                        Standard_error_p3 = summary_result_p3[2],
                        Statistic_p3 = summary_result_p3[3],
                        P_value_p3 = summary_result_p3[4],
                        #child_ffq_pattern_4_confectioneries_pizza
                        QoL_outcomes_p4 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p4 = "child_ffq_pattern_4_confectioneries_pizza",
                        Beta_coefficient_p4 = summary_result_p4[1],
                        Standard_error_p4 = summary_result_p4[2],
                        Statistic_p4 = summary_result_p4[3],
                        P_value_p4 = summary_result_p4[4],
                        #child_ffq_pattern_5_starchy_foods_sweetened_bev
                        QoL_outcomes_p5 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p5 = "child_ffq_pattern_5_starchy_foods_sweetened_bev",
                        Beta_coefficient_p5 = summary_result_p5[1],
                        Standard_error_p5 = summary_result_p5[2],
                        Statistic_p5 = summary_result_p5[3],
                        P_value_p5 = summary_result_p5[4])
}
dim(model_4_patterns)
#5 30

model_4_patterns_1 <- model_4_patterns %>% select (c(
        QoL_outcomes_p1,
        Dietary_features_p1,
        Beta_coefficient_p1,
        Standard_error_p1,
        Statistic_p1,
        P_value_p1)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p1,
                Dietary_features=Dietary_features_p1,
                Beta_coefficient=Beta_coefficient_p1,
                Standard_error=Standard_error_p1,
                Statistic=Statistic_p1,
                P_value=P_value_p1)

model_4_patterns_2 <- model_4_patterns %>% select (c(
        QoL_outcomes_p2,
        Dietary_features_p2,
        Beta_coefficient_p2,
        Standard_error_p2,
        Statistic_p2,
        P_value_p2)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p2,
                Dietary_features=Dietary_features_p2,
                Beta_coefficient=Beta_coefficient_p2,
                Standard_error=Standard_error_p2,
                Statistic=Statistic_p2,
                P_value=P_value_p2)

model_4_patterns_3 <- model_4_patterns %>% select (c(
        QoL_outcomes_p3,
        Dietary_features_p3,
        Beta_coefficient_p3,
        Standard_error_p3,
        Statistic_p3,
        P_value_p3)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p3,
                Dietary_features=Dietary_features_p3,
                Beta_coefficient=Beta_coefficient_p3,
                Standard_error=Standard_error_p3,
                Statistic=Statistic_p3,
                P_value=P_value_p3)

model_4_patterns_4 <- model_4_patterns %>% select (c(
        QoL_outcomes_p4,
        Dietary_features_p4,
        Beta_coefficient_p4,
        Standard_error_p4,
        Statistic_p4,
        P_value_p4)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p4,
                Dietary_features=Dietary_features_p4,
                Beta_coefficient=Beta_coefficient_p4,
                Standard_error=Standard_error_p4,
                Statistic=Statistic_p4,
                P_value=P_value_p4)

model_4_patterns_5 <- model_4_patterns %>% select (c(
        QoL_outcomes_p5,
        Dietary_features_p5,
        Beta_coefficient_p5,
        Standard_error_p5,
        Statistic_p5,
        P_value_p5)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p5,
                Dietary_features=Dietary_features_p5,
                Beta_coefficient=Beta_coefficient_p5,
                Standard_error=Standard_error_p5,
                Statistic=Statistic_p5,
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

                ####1.5 model 7: model 5 + predictors with greater importance than dietary features ####

#dietary indexes
model_7_indexes <- foreach(j = 1:ncol(QoL_perc_measurements),.combine = rbind) %do% {
                lm1 = lm(QoL_perc_measurements[,j] ~ 
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
                                 masterfile$child_ffq_animal_score)
                sum1 = summary(lm1)
                #child_ffq_hpdi_score
                summary_result_hpdi = sum1$coef[nrow(sum1$coef)-1,]
                #child_ffq_animal_score
                summary_result_animal = sum1$coef[nrow(sum1$coef),]
                data.frame(#child_ffq_hpdi_score
                        QoL_outcomes_hpdi = colnames(QoL_perc_measurements)[j],
                        Dietary_features_hpdi = "child_ffq_hpdi_score",
                        Beta_coefficient_hpdi = summary_result_hpdi[1],
                        Standard_error_hpdi = summary_result_hpdi[2],
                        Statistic_hpdi = summary_result_hpdi[3],
                        P_value_hpdi = summary_result_hpdi[4],
                        #child_ffq_animal_score
                        QoL_outcomes_animal = colnames(QoL_perc_measurements)[j],
                        Dietary_features_animal = "child_ffq_animal_score",
                        Beta_coefficient_animal = summary_result_animal[1],
                        Standard_error_animal = summary_result_animal[2],
                        Statistic_animal = summary_result_animal[3],
                        P_value_animal = summary_result_animal[4])
}
dim(model_7_indexes)
#5 12 

model_7_indexes_hpdi <- model_7_indexes %>% select (c(
        QoL_outcomes_hpdi,
        Dietary_features_hpdi,
        Beta_coefficient_hpdi,
        Standard_error_hpdi,
        Statistic_hpdi,
        P_value_hpdi)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_hpdi,
                Dietary_features=Dietary_features_hpdi,
                Beta_coefficient=Beta_coefficient_hpdi,
                Standard_error=Standard_error_hpdi,
                Statistic=Statistic_hpdi,
                P_value=P_value_hpdi)

model_7_indexes_animal <- model_7_indexes %>% select (c(
        QoL_outcomes_animal,
        Dietary_features_animal,
        Beta_coefficient_animal,
        Standard_error_animal,
        Statistic_animal,
        P_value_animal)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_animal,
                Dietary_features=Dietary_features_animal,
                Beta_coefficient=Beta_coefficient_animal,
                Standard_error=Standard_error_animal,
                Statistic=Statistic_animal,
                P_value=P_value_animal)

model_7_indexes_final <- bind_rows(model_7_indexes_hpdi,
                                   model_7_indexes_animal)
rm(model_7_indexes_hpdi, model_7_indexes_animal, model_7_indexes)

#dietary patterns
model_7_patterns <- foreach(j = 1:ncol(QoL_perc_measurements),.combine = rbind) %do% {
                lm1 = lm(QoL_perc_measurements[,j] ~ 
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
                                 masterfile$child_ffq_pattern_5_starchy_foods_sweetened_bev)
                sum1 = summary(lm1)
                #child_ffq_pattern_1_meat_seafood_prepd_meals
                summary_result_p1 = sum1$coef[nrow(sum1$coef)-4,]
                #child_ffq_pattern_2_cooked_veg_grains_legumes
                summary_result_p2 = sum1$coef[nrow(sum1$coef)-3,]
                #child_ffq_pattern_3_fruits_raw_veg_cheese
                summary_result_p3 = sum1$coef[nrow(sum1$coef)-2,]
                #child_ffq_pattern_4_confectioneries_pizza
                summary_result_p4 = sum1$coef[nrow(sum1$coef)-1,]
                #child_ffq_pattern_5_starchy_foods_sweetened_bev
                summary_result_p5 = sum1$coef[nrow(sum1$coef),]
                data.frame(
                        #child_ffq_pattern_1_meat_seafood_prepd_meals
                        QoL_outcomes_p1 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p1 = "child_ffq_pattern_1_meat_seafood_prepd_meals",
                        Beta_coefficient_p1 = summary_result_p1[1],
                        Standard_error_p1 = summary_result_p1[2],
                        Statistic_p1 = summary_result_p1[3],
                        P_value_p1 = summary_result_p1[4],
                        #child_ffq_pattern_2_cooked_veg_grains_legumes
                        QoL_outcomes_p2 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p2 = "child_ffq_pattern_2_cooked_veg_grains_legumes",
                        Beta_coefficient_p2 = summary_result_p2[1],
                        Standard_error_p2 = summary_result_p2[2],
                        Statistic_p2 = summary_result_p2[3],
                        P_value_p2 = summary_result_p2[4],
                        #child_ffq_pattern_3_fruits_raw_veg_cheese
                        QoL_outcomes_p3 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p3 = "child_ffq_pattern_3_fruits_raw_veg_cheese",
                        Beta_coefficient_p3 = summary_result_p3[1],
                        Standard_error_p3 = summary_result_p3[2],
                        Statistic_p3 = summary_result_p3[3],
                        P_value_p3 = summary_result_p3[4],
                        #child_ffq_pattern_4_confectioneries_pizza
                        QoL_outcomes_p4 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p4 = "child_ffq_pattern_4_confectioneries_pizza",
                        Beta_coefficient_p4 = summary_result_p4[1],
                        Standard_error_p4 = summary_result_p4[2],
                        Statistic_p4 = summary_result_p4[3],
                        P_value_p4 = summary_result_p4[4],
                        #child_ffq_pattern_5_starchy_foods_sweetened_bev
                        QoL_outcomes_p5 = colnames(QoL_perc_measurements)[j],
                        Dietary_features_p5 = "child_ffq_pattern_5_starchy_foods_sweetened_bev",
                        Beta_coefficient_p5 = summary_result_p5[1],
                        Standard_error_p5 = summary_result_p5[2],
                        Statistic_p5 = summary_result_p5[3],
                        P_value_p5 = summary_result_p5[4])
}
dim(model_7_patterns)
#5 30

model_7_patterns_1 <- model_7_patterns %>% select (c(
        QoL_outcomes_p1,
        Dietary_features_p1,
        Beta_coefficient_p1,
        Standard_error_p1,
        Statistic_p1,
        P_value_p1)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p1,
                Dietary_features=Dietary_features_p1,
                Beta_coefficient=Beta_coefficient_p1,
                Standard_error=Standard_error_p1,
                Statistic=Statistic_p1,
                P_value=P_value_p1)

model_7_patterns_2 <- model_7_patterns %>% select (c(
        QoL_outcomes_p2,
        Dietary_features_p2,
        Beta_coefficient_p2,
        Standard_error_p2,
        Statistic_p2,
        P_value_p2)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p2,
                Dietary_features=Dietary_features_p2,
                Beta_coefficient=Beta_coefficient_p2,
                Standard_error=Standard_error_p2,
                Statistic=Statistic_p2,
                P_value=P_value_p2)

model_7_patterns_3 <- model_7_patterns %>% select (c(
        QoL_outcomes_p3,
        Dietary_features_p3,
        Beta_coefficient_p3,
        Standard_error_p3,
        Statistic_p3,
        P_value_p3)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p3,
                Dietary_features=Dietary_features_p3,
                Beta_coefficient=Beta_coefficient_p3,
                Standard_error=Standard_error_p3,
                Statistic=Statistic_p3,
                P_value=P_value_p3)

model_7_patterns_4 <- model_7_patterns %>% select (c(
        QoL_outcomes_p4,
        Dietary_features_p4,
        Beta_coefficient_p4,
        Standard_error_p4,
        Statistic_p4,
        P_value_p4)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p4,
                Dietary_features=Dietary_features_p4,
                Beta_coefficient=Beta_coefficient_p4,
                Standard_error=Standard_error_p4,
                Statistic=Statistic_p4,
                P_value=P_value_p4)

model_7_patterns_5 <- model_7_patterns %>% select (c(
        QoL_outcomes_p5,
        Dietary_features_p5,
        Beta_coefficient_p5,
        Standard_error_p5,
        Statistic_p5,
        P_value_p5)) %>% 
        rename(
                QoL_outcomes=QoL_outcomes_p5,
                Dietary_features=Dietary_features_p5,
                Beta_coefficient=Beta_coefficient_p5,
                Standard_error=Standard_error_p5,
                Statistic=Statistic_p5,
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

                #### =============== 1. MODELS FOR MANUSCRIPT - FDR CORRECTION   ================ ####
final_linear_models <- bind_rows(model_1, #crude model
                                 model_2, #standard correction for age and sex
                                 model_4, #anthropometric model
                                 model_7) #anthropometric and prediction model
str(final_linear_models)
final_linear_models <- final_linear_models %>% mutate_if(is.character, as.factor)

final_linear_models <- final_linear_models %>%
        mutate(QoL_outcomes = fct_relevel(QoL_outcomes, "child_pedsql_health_related_qol_score_total_perc", 
                                          "child_pedsql_phy_score_perc", 
                                          "child_pedsql_emo_score_perc",
                                          "child_pedsql_soc_score_perc",
                                          "child_pedsql_sch_score_perc")) %>%
        mutate(Dietary_features = fct_relevel(Dietary_features, "child_ffq_hpdi_score", 
                                              "child_ffq_animal_score", 
                                              "child_ffq_pattern_1_meat_seafood_prepd_meals",
                                              "child_ffq_pattern_2_cooked_veg_grains_legumes",
                                              "child_ffq_pattern_3_fruits_raw_veg_cheese",
                                              "child_ffq_pattern_4_confectioneries_pizza",
                                              "child_ffq_pattern_5_starchy_foods_sweetened_bev"))

final_linear_models$FDR = p.adjust(final_linear_models$P_value, method = "BH")
sum(final_linear_models$FDR<0.05)
#77
sum(final_linear_models$P_value<0.05)
#82

#final_linear_models <- final_linear_models[order(final_linear_models$Model, final_linear_models$Dietary_features), ]
final_linear_models <- final_linear_models[order(final_linear_models$QoL_outcomes, final_linear_models$Dietary_features), ]
#140  8

write.table(final_linear_models, "2023_06_21_qol_diet_linear_manuscript_models_FDR_updated.txt", sep="\t", row.names=F, quote = F) 
