###     SCRIPT: PREDICTION MODELLING USING REGULARIZED LOGISTIC REGRESSION
###     AUTHOR(S): SIOBHAN BRUSHETT 
###     DESCRIPTION: PREDICT WHICH FEATURES ARE ASSOCIATED WITH QOL USING REGULARIZED LOGISTIC REGRESSION OF TIDY MODELS
###     PROJECT: DIATROFI_BASELINE_2015_2018
###     NOTE(S): https://juliasilge.com/blog/project-feederwatch/ #downsampling and regularized logistic regression
###              https://juliasilge.com/blog/ikea-prices/ #feature engineering and recipes
###              https://recipes.tidymodels.org/reference/step_unknown.html #new levels when using step_dummy with missings in categorical data
###              https://www.kaggle.com/code/dansbecker/what-is-log-loss 
###              https://rpubs.com/eR_ic/regression #linear regression prediction modelling

###     LAST UPDATED: 15 MAY 2023

#libraries
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
library(skimr)
library(reshape2)
library(tidymodels)
library(themis) #test down-sampling approach
library(vip) #which predictors are driving the classification model

#functions
#spearman correlation (Alex and Trishla)
spearman <- function(x,y) {
  matchID <- intersect(rownames(x), rownames(y))
  x1 <- x[matchID,]
  y1 <- y[matchID,]
  result_cor <- matrix(nrow=ncol(x1), ncol=ncol(y1))
  rownames(result_cor) <- colnames(x1)
  colnames(result_cor) <- colnames(y1)
  result_pvalue <- matrix(nrow=ncol(x1), ncol=ncol(y1))
  for (i in 1:ncol(y1)) {
    for (j in 1:ncol(x1)) {
      cor1<-try(cor.test(x1[,j], y1[,i], method = "spearman"))
      if(class(cor1)[1] != "try-error") {
        result_cor[j,i]= cor1$estimate
        result_pvalue[j,i]= cor1$p.value
      } else {
        result_cor[j,i]= NA
        result_pvalue[j,i]= NA
      }
    }}
  result = list()
  result$p.val= result_pvalue
  result$cor= result_cor
  return(result)
}
result <- function(x,y){
  correlation <- spearman(x,y)
  a<- melt(correlation$cor)
  a<- cbind(a, melt(correlation$p.val)[,"value"])
  result= a[order(a[,4]),]
  colnames(result)=c("factor1", "factor2", "CorCoefficient","pvalue")
  return(result)
}

## 0. LOAD DATA
## 1. EXPLORATORY PLOTS
## 2. VARIABLE SELECTION AND SPEARMAN CORRELATION PLOTS
## 3. FEATURE ENGINEERING: HRQoL Q15 PREDICTION MODELS
## 4. FEATURE ENGINEERING: HRQoL CONTINOUS PREDICTION MODELS

        #### =============== 0. LOAD DATA   ================ ####
masterfile <- read.delim("2023_05_09_diatrofi_2015_2018_masterfile_updated.txt", sep = "\t", header = T, stringsAsFactors = T)
#6583 66
str(masterfile)
masterfile$child_demo_school_code <- as.factor(masterfile$child_demo_school_code)

        #### =============== 1. EXPLORATORY PLOTS   ================ ####
#healthy plant dietary index (hPDI) vs total HRQoL by several other possible predictors
pdf(file = "2023_05_09_diet_by_hrqol_q15_covariates.pdf", useDingbats = F, onefile = T)
masterfile %>% select(child_pedsql_health_related_qol_score_total_catg_Q15,
                      child_ffq_hpdi_score) %>%
        na.omit() %>%
        group_by(child_pedsql_health_related_qol_score_total_catg_Q15) %>%
        summarise(child_ffq_hpdi_score=mean(child_ffq_hpdi_score)) %>%
        ggplot(aes(child_ffq_hpdi_score, child_pedsql_health_related_qol_score_total_catg_Q15, fill=child_pedsql_health_related_qol_score_total_catg_Q15))+
        geom_col(position="dodge") + 
        labs(x="hPDI score", y="Total HRQoL score", fill=NULL) +
        theme_bw() +
        scale_colour_brewer(palette = "Set2") + 
        scale_fill_brewer(palette = "Set2") +
        theme(legend.position="none")

masterfile %>% select(child_pedsql_health_related_qol_score_total_catg_Q15,
                      child_ffq_hpdi_score, child_demo_school_region) %>%
        na.omit() %>%
        group_by(child_pedsql_health_related_qol_score_total_catg_Q15, child_demo_school_region) %>%
        summarise(child_ffq_hpdi_score=mean(child_ffq_hpdi_score)) %>%
        ggplot(aes(child_ffq_hpdi_score, child_demo_school_region, fill=child_pedsql_health_related_qol_score_total_catg_Q15)) +
        geom_col(position="dodge") + 
        labs(x="hPDI score", y="School region", fill="Total HRQoL score") +
        theme_bw() +
        scale_colour_brewer(palette = "Set2") + 
        scale_fill_brewer(palette = "Set2")

masterfile %>% select(child_pedsql_health_related_qol_score_total_catg_Q15,
                      child_ffq_hpdi_score, family_demo_FAS_group_dic) %>%
        na.omit() %>%
        group_by(child_pedsql_health_related_qol_score_total_catg_Q15, family_demo_FAS_group_dic) %>%
        summarise(child_ffq_hpdi_score=mean(child_ffq_hpdi_score)) %>%
        ggplot(aes(child_ffq_hpdi_score, family_demo_FAS_group_dic, fill=child_pedsql_health_related_qol_score_total_catg_Q15)) +
        geom_col(position="dodge") + 
        labs(x="hPDI score", y="Family Affluent Scale (FAS)", fill="Total HRQoL score") +
        theme_bw() +
        scale_colour_brewer(palette = "Set2") + 
        scale_fill_brewer(palette = "Set2")

masterfile %>% select(child_pedsql_health_related_qol_score_total_catg_Q15,
                      child_ffq_hpdi_score, child_demo_immigrant_status) %>%
        na.omit() %>%
        group_by(child_pedsql_health_related_qol_score_total_catg_Q15, child_demo_immigrant_status) %>%
        summarise(child_ffq_hpdi_score=mean(child_ffq_hpdi_score)) %>%
        ggplot(aes(child_ffq_hpdi_score, child_demo_immigrant_status, fill=child_pedsql_health_related_qol_score_total_catg_Q15)) +
        geom_col(position="dodge") + 
        labs(x="hPDI score", y="Child immigrant status", fill="Total HRQoL score") +
        theme_bw() +
        scale_colour_brewer(palette = "Set2") + 
        scale_fill_brewer(palette = "Set2")

masterfile %>% select(child_pedsql_health_related_qol_score_total_catg_Q15,
                      child_ffq_hpdi_score, child_anthro_bmi_WHO_sds_grouped) %>%
        na.omit() %>%
        group_by(child_pedsql_health_related_qol_score_total_catg_Q15, child_anthro_bmi_WHO_sds_grouped) %>%
        summarise(child_ffq_hpdi_score=mean(child_ffq_hpdi_score)) %>%
        ggplot(aes(child_ffq_hpdi_score, child_anthro_bmi_WHO_sds_grouped, fill=child_pedsql_health_related_qol_score_total_catg_Q15)) +
        geom_col(position="dodge") + 
        labs(x="hPDI score", y="Child BMI WHO thresholds", fill="Total HRQoL score") +
        theme_bw() +
        scale_colour_brewer(palette = "Set2") + 
        scale_fill_brewer(palette = "Set2")

#Dietary patterns 3 (fruits, raw veg, cheese) vs total HRQoL by several other possible predictors
masterfile %>% select(child_pedsql_health_related_qol_score_total_catg_Q15,
                      child_ffq_pattern_3_fruits_raw_veg_cheese) %>%
        na.omit() %>%
        group_by(child_pedsql_health_related_qol_score_total_catg_Q15) %>%
        summarise(child_ffq_pattern_3_fruits_raw_veg_cheese=median(child_ffq_pattern_3_fruits_raw_veg_cheese)) %>%
        ggplot(aes(child_ffq_pattern_3_fruits_raw_veg_cheese, child_pedsql_health_related_qol_score_total_catg_Q15, fill=child_pedsql_health_related_qol_score_total_catg_Q15))+
        geom_col(position="dodge") + 
        labs(x="Pattern 3: fruit, raw veg and cheese", y="Total HRQoL score", fill=NULL) +
        theme_bw() +
        scale_colour_brewer(palette = "Set2") + 
        scale_fill_brewer(palette = "Set2") +
        theme(legend.position="none")

masterfile %>% select(child_pedsql_health_related_qol_score_total_catg_Q15,
                      child_ffq_pattern_3_fruits_raw_veg_cheese, child_demo_school_region) %>%
        na.omit() %>%
        group_by(child_pedsql_health_related_qol_score_total_catg_Q15, child_demo_school_region) %>%
        summarise(child_ffq_pattern_3_fruits_raw_veg_cheese=median(child_ffq_pattern_3_fruits_raw_veg_cheese)) %>%
        ggplot(aes(child_ffq_pattern_3_fruits_raw_veg_cheese, child_demo_school_region, fill=child_pedsql_health_related_qol_score_total_catg_Q15)) +
        geom_col(position="dodge") + 
        labs(x="Pattern 3: fruit, raw veg and cheese", y="School region", fill="Total HRQoL score") +
        theme_bw() +
        scale_colour_brewer(palette = "Set2") + 
        scale_fill_brewer(palette = "Set2")

masterfile %>% select(child_pedsql_health_related_qol_score_total_catg_Q15,
                      child_ffq_pattern_3_fruits_raw_veg_cheese, family_demo_FAS_group_dic) %>%
        na.omit() %>%
        group_by(child_pedsql_health_related_qol_score_total_catg_Q15, family_demo_FAS_group_dic) %>%
        summarise(child_ffq_pattern_3_fruits_raw_veg_cheese=median(child_ffq_pattern_3_fruits_raw_veg_cheese)) %>%
        ggplot(aes(child_ffq_pattern_3_fruits_raw_veg_cheese, family_demo_FAS_group_dic, fill=child_pedsql_health_related_qol_score_total_catg_Q15)) +
        geom_col(position="dodge") + 
        labs(x="Pattern 3: fruit, raw veg and cheese", y="Family Affluent Scale (FAS)", fill="Total HRQoL score") +
        theme_bw() +
        scale_colour_brewer(palette = "Set2") + 
        scale_fill_brewer(palette = "Set2")

masterfile %>% select(child_pedsql_health_related_qol_score_total_catg_Q15,
                      child_ffq_pattern_3_fruits_raw_veg_cheese, child_demo_immigrant_status) %>%
        na.omit() %>%
        group_by(child_pedsql_health_related_qol_score_total_catg_Q15, child_demo_immigrant_status) %>%
        summarise(child_ffq_pattern_3_fruits_raw_veg_cheese=median(child_ffq_pattern_3_fruits_raw_veg_cheese)) %>%
        ggplot(aes(child_ffq_pattern_3_fruits_raw_veg_cheese, child_demo_immigrant_status, fill=child_pedsql_health_related_qol_score_total_catg_Q15)) +
        geom_col(position="dodge") + 
        labs(x="Pattern 3: fruit, raw veg and cheese", y="Child immigrant status", fill="Total HRQoL score") +
        theme_bw() +
        scale_colour_brewer(palette = "Set2") + 
        scale_fill_brewer(palette = "Set2")

masterfile %>% select(child_pedsql_health_related_qol_score_total_catg_Q15,
                      child_ffq_pattern_3_fruits_raw_veg_cheese, child_anthro_bmi_WHO_sds_grouped) %>%
        na.omit() %>%
        group_by(child_pedsql_health_related_qol_score_total_catg_Q15, child_anthro_bmi_WHO_sds_grouped) %>%
        summarise(child_ffq_pattern_3_fruits_raw_veg_cheese=median(child_ffq_pattern_3_fruits_raw_veg_cheese)) %>%
        ggplot(aes(child_ffq_pattern_3_fruits_raw_veg_cheese, child_anthro_bmi_WHO_sds_grouped, fill=child_pedsql_health_related_qol_score_total_catg_Q15)) +
        geom_col(position="dodge") + 
        labs(x="Pattern 3: fruit, raw veg and cheese", y="Child BMI WHO thresholds", fill="Total HRQoL score") +
        theme_bw() +
        scale_colour_brewer(palette = "Set2") + 
        scale_fill_brewer(palette = "Set2")               

#Dietary patterns 5 (starchy foods and sweetend beverages) vs total HRQoL by several other possible predictors
masterfile %>% select(child_pedsql_health_related_qol_score_total_catg_Q15,
                      child_ffq_pattern_5_starchy_foods_sweetened_bev) %>%
        na.omit() %>%
        group_by(child_pedsql_health_related_qol_score_total_catg_Q15) %>%
        summarise(child_ffq_pattern_5_starchy_foods_sweetened_bev=median(child_ffq_pattern_5_starchy_foods_sweetened_bev)) %>%
        ggplot(aes(child_ffq_pattern_5_starchy_foods_sweetened_bev, child_pedsql_health_related_qol_score_total_catg_Q15, fill=child_pedsql_health_related_qol_score_total_catg_Q15))+
        geom_col(position="dodge") + 
        labs(x="Pattern 5: starchy foods and sweetend beverages", y="Total HRQoL score", fill=NULL) +
        theme_bw() +
        scale_colour_brewer(palette = "Set2") + 
        scale_fill_brewer(palette = "Set2") +
        theme(legend.position="none")

masterfile %>% select(child_pedsql_health_related_qol_score_total_catg_Q15,
                      child_ffq_pattern_5_starchy_foods_sweetened_bev, child_demo_school_region) %>%
        na.omit() %>%
        group_by(child_pedsql_health_related_qol_score_total_catg_Q15, child_demo_school_region) %>%
        summarise(child_ffq_pattern_5_starchy_foods_sweetened_bev=median(child_ffq_pattern_5_starchy_foods_sweetened_bev)) %>%
        ggplot(aes(child_ffq_pattern_5_starchy_foods_sweetened_bev, child_demo_school_region, fill=child_pedsql_health_related_qol_score_total_catg_Q15)) +
        geom_col(position="dodge") + 
        labs(x="Pattern 5: starchy foods and sweetend beverages", y="School region", fill="Total HRQoL score") +
        theme_bw() +
        scale_colour_brewer(palette = "Set2") + 
        scale_fill_brewer(palette = "Set2")

masterfile %>% select(child_pedsql_health_related_qol_score_total_catg_Q15,
                      child_ffq_pattern_5_starchy_foods_sweetened_bev, family_demo_FAS_group_dic) %>%
        na.omit() %>%
        group_by(child_pedsql_health_related_qol_score_total_catg_Q15, family_demo_FAS_group_dic) %>%
        summarise(child_ffq_pattern_5_starchy_foods_sweetened_bev=median(child_ffq_pattern_5_starchy_foods_sweetened_bev)) %>%
        ggplot(aes(child_ffq_pattern_5_starchy_foods_sweetened_bev, family_demo_FAS_group_dic, fill=child_pedsql_health_related_qol_score_total_catg_Q15)) +
        geom_col(position="dodge") + 
        labs(x="Pattern 5: starchy foods and sweetend beverages", y="Family Affluent Scale (FAS)", fill="Total HRQoL score") +
        theme_bw() +
        scale_colour_brewer(palette = "Set2") + 
        scale_fill_brewer(palette = "Set2")

masterfile %>% select(child_pedsql_health_related_qol_score_total_catg_Q15,
                      child_ffq_pattern_5_starchy_foods_sweetened_bev, child_demo_immigrant_status) %>%
        na.omit() %>%
        group_by(child_pedsql_health_related_qol_score_total_catg_Q15, child_demo_immigrant_status) %>%
        summarise(child_ffq_pattern_5_starchy_foods_sweetened_bev=median(child_ffq_pattern_5_starchy_foods_sweetened_bev)) %>%
        ggplot(aes(child_ffq_pattern_5_starchy_foods_sweetened_bev, child_demo_immigrant_status, fill=child_pedsql_health_related_qol_score_total_catg_Q15)) +
        geom_col(position="dodge") + 
        labs(x="Pattern 5: starchy foods and sweetend beverages", y="Child immigrant status", fill="Total HRQoL score") +
        theme_bw() +
        scale_colour_brewer(palette = "Set2") + 
        scale_fill_brewer(palette = "Set2")

masterfile %>% select(child_pedsql_health_related_qol_score_total_catg_Q15,
                      child_ffq_pattern_5_starchy_foods_sweetened_bev, child_anthro_bmi_WHO_sds_grouped) %>%
        na.omit() %>%
        group_by(child_pedsql_health_related_qol_score_total_catg_Q15, child_anthro_bmi_WHO_sds_grouped) %>%
        summarise(child_ffq_pattern_5_starchy_foods_sweetened_bev=median(child_ffq_pattern_5_starchy_foods_sweetened_bev)) %>%
        ggplot(aes(child_ffq_pattern_5_starchy_foods_sweetened_bev, child_anthro_bmi_WHO_sds_grouped, fill=child_pedsql_health_related_qol_score_total_catg_Q15)) +
        geom_col(position="dodge") + 
        labs(x="Pattern 5: starchy foods and sweetend beverages", y="Child BMI WHO thresholds", fill="Total HRQoL score") +
        theme_bw() +
        scale_colour_brewer(palette = "Set2") + 
        scale_fill_brewer(palette = "Set2")               

dev.off()               

        #### =============== 2. VARIABLE SELECTION AND SPEARMAN CORRELATION PLOTS ================ ####
str(masterfile)

#filter masterfile for variables of interest
colnames(masterfile)
#child_demo_school_year
#child_demo_school_region
#child_demo_age_years
#child_demo_sex
#child_anthro_bmi_WHO_sds
#child_demo_immigrant_status (this encompasses child_demo_country_birth and parent_demo_country_birth)
#family_demo_FAS_group (I expect this to be associated with child_demo_lives_with, family_demo_children_in_family_n, 
#                       family_demo_household_members_n and child_demo_insurance, parent_demo_employment,
#                       mother_demo_age_years, mother_demo_edu, father_demo_age_years, father_demo_edu)
correlations <- masterfile
table(correlations$family_demo_FAS_group, correlations$child_demo_lives_with) #NOTE: (84% of families live with married parents/guardians) 
table(correlations$family_demo_FAS_group, correlations$family_demo_children_in_family_n) #low-middle FAS have bigger families (as expected)
table(correlations$family_demo_FAS_group, correlations$family_demo_household_members_n) #low-middle FAS have larger households (as expected)
table(correlations$family_demo_FAS_group, correlations$child_demo_insurance) #low FAS are less covered (18%) in comparison to middle-high FAS (~3-6%)
table(correlations$family_demo_FAS_group, correlations$parent_demo_employment) #unemployment is higher in low FAS (16%) in comparison to middle (4.6%) and high (1.2%) FAS
table(correlations$family_demo_FAS_group, correlations$mother_demo_age_years)
gghistogram(data=subset(correlations, !is.na(family_demo_FAS_group) & !is.na(mother_demo_age_years)),
            x = "mother_demo_age_years", bins = 15, 
            fill = "family_demo_FAS_group", color = "family_demo_FAS_group",
            palette = "Set1", 
            add = "mean", mean.color = "black",
            mean.size = 1, mean.linetype = "dashed",
            alpha = 0.5) + 
  labs(title = "Histogram of mother age by FAS group", 
       x = "mother_demo_age_years", y = "Frequency", 
       fill = "family_demo_FAS_group") + 
  theme_pubclean()
#average age of the mother is lower in low FAS, then middle FAS and higher age for high FAS

table(correlations$family_demo_FAS_group, correlations$mother_demo_edu)
#lower FAS, lower education; higher FAS, middle-high education

table(correlations$family_demo_FAS_group, correlations$father_demo_age_years)
gghistogram(data=subset(correlations, !is.na(family_demo_FAS_group) & !is.na(father_demo_age_years)),
            x = "father_demo_age_years", bins = 15, 
            fill = "family_demo_FAS_group", color = "family_demo_FAS_group",
            palette = "Set1", 
            add = "mean", mean.color = "black",
            mean.size = 1, mean.linetype = "dashed",
            alpha = 0.5) + 
  labs(title = "Histogram of mother age by FAS group", 
       x = "father_demo_age_years", y = "Frequency", 
       fill = "family_demo_FAS_group") + 
  theme_pubclean()
#average age of the father is lower in low FAS, then middle FAS and higher age for high FAS

table(correlations$family_demo_FAS_group, correlations$father_demo_edu)
#lower FAS, lower education; higher FAS, middle-high education

#child_phyact_exercise_hrs_per_wk (correlation with child_phyact_sport_outside_school and child_phyact_play_outside_days_per_wk
correlations$child_phyact_sport_outside_school <- as.numeric(correlations$child_phyact_sport_outside_school)
cor(correlations$child_phyact_exercise_hrs_per_wk, correlations$child_phyact_sport_outside_school, use="complete.obs")
#0.4053047
table(correlations$child_phyact_exercise_hrs_per_wk, correlations$child_phyact_sport_outside_school)

cor(correlations$child_phyact_exercise_hrs_per_wk, correlations$child_phyact_play_outside_days_per_wk, use="complete.obs")
#0.1147949
table(correlations$child_phyact_exercise_hrs_per_wk, correlations$child_phyact_play_outside_days_per_wk)
#weak correlation, but expected because the weight differs, i.e. hrs per week vs days per week

#use child_phyact_exercise_hrs_per_wk to represent physical activity
#child_phyact_screen_time_hrs_per_wk

#include all dietary features

rm(correlations)

#Thus, the list of predictors include (n=17):
spearman_var <- masterfile %>% 
  select(c(#predictor variables
           child_demo_school_year,
           child_demo_school_region,
           child_demo_age_years,
           child_demo_sex, 
           child_anthro_bmi_WHO_sds,
           child_demo_immigrant_status,
           family_demo_FAS_group,  
           #child_phyact_exercise_hrs_per_wk, #one of the components of QoL is Physical Functioning
           child_phyact_screen_time_hrs_per_wk, 
           child_ffq_hpdi_score, 
           child_ffq_animal_score,
           child_ffq_diet_quality_score,
           child_ffq_pattern_1_meat_seafood_prepd_meals,
           child_ffq_pattern_2_cooked_veg_grains_legumes,
           child_ffq_pattern_3_fruits_raw_veg_cheese,
           child_ffq_pattern_4_confectioneries_pizza,
           child_ffq_pattern_5_starchy_foods_sweetened_bev,
           #outcome variables
           child_pedsql_phy_score_catg_Q15,
           child_pedsql_emo_score_catg_Q15,
           child_pedsql_soc_score_catg_Q15,
           child_pedsql_psysoc_score_catg_Q15, 
           child_pedsql_health_related_qol_score_total_catg_Q15,
           child_pedsql_phy_score_perc,
           child_pedsql_emo_score_perc,
           child_pedsql_soc_score_perc,
           child_pedsql_sch_score_perc,
           child_pedsql_psysoc_score_perc,
           child_pedsql_health_related_qol_score_total_perc)) 

str(spearman_var)
#NOTE: retrospectively removed child_phyact_exercise_hrs_per_wk as it landed up having an importance of 22.3%
#this makes sense given that one of the components of QoL is Physical Functioning

for (i in 1:ncol(spearman_var)){spearman_var[,i] <- as.numeric(spearman_var[,i])}
str(spearman_var)
#6583 27

spearman_corr <- result(spearman_var, spearman_var)
#729  4

#remove 1:1 correlations
spearman_corr_unique <- spearman_corr %>% 
  as_tibble() %>% 
  mutate(duplicates = if_else(factor1 == factor2,
                              TRUE,
                              FALSE)) %>% 
  filter(duplicates == FALSE)

spearman_corr_unique$duplicates <- NULL 
#702  4

spearman_FDR <- spearman_corr_unique
spearman_FDR$FDR<-p.adjust(spearman_FDR$pvalue, method = "BH")
#702  5
spearman_FDR <- spearman_FDR %>% filter(
  FDR<0.05
)
#558  5

pdf(file = "2023_05_15_hrqol_predictors_outcomes_spearman.pdf", useDingbats = F, onefile = T, width = 15, height=18)
ggplot(data = spearman_corr, aes(factor1, factor2, fill = CorCoefficient))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ 
  #theme(axis.text.x=element_blank(), #remove x axis labels
  #      axis.ticks.x=element_blank(), #remove x axis ticks
  #      axis.text.y=element_blank(),  #remove y axis labels
  #      axis.ticks.y=element_blank()  #remove y axis ticks
  #)+
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text=element_text(size=8))

ggplot(data = spearman_FDR, aes(factor1, factor2, fill = CorCoefficient))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ 
  #theme(axis.text.x=element_blank(), #remove x axis labels
  #      axis.ticks.x=element_blank(), #remove x axis ticks
  #      axis.text.y=element_blank(),  #remove y axis labels
  #      axis.ticks.y=element_blank()  #remove y axis ticks
  #)+
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text=element_text(size=8))

dev.off()

#NOTE: The Psycho-Social Functioning Composite Score will not be included as an outcome for linear and logistic regression as it is
#      highly correlated with Emotional, Social and School Functioning (as these scales make up this composite score)

        #### =============== 3. FEATURE ENGINEERING: HRQoL Q15 PREDICTION MODELS   ================ ####
hrqol_Q15_split <- masterfile %>% 
        filter(!is.na(child_pedsql_health_related_qol_score_total_catg_Q15)) %>%
        select(c(child_demo_school_year,
                 child_demo_school_region,
                 child_demo_age_years,
                 child_demo_sex, 
                 child_anthro_bmi_WHO_sds,
                 child_demo_immigrant_status,
                 family_demo_FAS_group,  
                 #child_phyact_exercise_hrs_per_wk, #one of the components of QoL is Physical Functioning
                 child_phyact_screen_time_hrs_per_wk, 
                 child_ffq_hpdi_score, 
                 child_ffq_animal_score,
                 child_ffq_diet_quality_score,
                 child_ffq_pattern_1_meat_seafood_prepd_meals,
                 child_ffq_pattern_2_cooked_veg_grains_legumes,
                 child_ffq_pattern_3_fruits_raw_veg_cheese,
                 child_ffq_pattern_4_confectioneries_pizza,
                 child_ffq_pattern_5_starchy_foods_sweetened_bev,
                 child_pedsql_health_related_qol_score_total_catg_Q15)) 
#6481 17
skim(hrqol_Q15_split)

set.seed(123)
#split the data into training and test sets, and make sure that the proportion of poor_avg and good hrqol are equally distributed between these sets
hrqol_Q15_split <- hrqol_Q15_split %>%
        initial_split(strata = child_pedsql_health_related_qol_score_total_catg_Q15)

hrqol_Q15_split
#<Training/Testing/Total>
#<4860/1621/6481>

hrqol_Q15_train <- training (hrqol_Q15_split)
hrqol_Q15_test <- testing (hrqol_Q15_split)

#confirm that the proportion of poor_avg and good hrqol are equally distributed between these sets
hrqol_Q15_train %>%
        count(child_pedsql_health_related_qol_score_total_catg_Q15) %>%
        mutate(prop = n/sum(n)*100)
#child_pedsql_health_related_qol_score_total_catg_Q15    n     prop
#1                                           0_poor_avg  740 15.22634
#2                                               1_good 4120 84.77366

hrqol_Q15_test %>%
        count(child_pedsql_health_related_qol_score_total_catg_Q15) %>%
        mutate(prop = n/sum(n)*100)
#child_pedsql_health_related_qol_score_total_catg_Q15    n     prop
#1                                           0_poor_avg  247 15.23751
#2                                               1_good 1374 84.76249

#hrqol is equally represented in both the training and test set

#for later resampling
set.seed(123)
hrqol_folds <- vfold_cv(hrqol_Q15_train, strata = child_pedsql_health_related_qol_score_total_catg_Q15)
hrqol_folds
#<split [4374/486]> Fold01
#etc..

summary(hrqol_Q15_train)

#make recipe for tidy models
hrqol_recipe <- 
        recipe(child_pedsql_health_related_qol_score_total_catg_Q15 ~ ., data = hrqol_Q15_train) %>% #relate hrqol to all other variables in the dataframe
        step_impute_mean(all_numeric_predictors()) %>% #impute numeric variables (NOTE: imputation is needed to run the downstream logistic regression)
        step_impute_bag(all_nominal(), -all_outcomes()) %>% #impute categorical variables
        step_dummy(all_nominal(), -all_outcomes()) %>% #generate dummy variables for all categorical variables, excluding the outcome, hrqol
        step_nzv(all_predictors()) #remove near zero variance for all predictor (numerical and nominal) variables

#confirm that the recipe works
prep(hrqol_recipe)
#NOTE: Sparse, unbalanced variable filter removed: child_demo_immigrant_status_X1_1st_gen_immigrant

#specify the regularized logistic regression model
glmnet_spec <- 
        logistic_reg(penalty = tune(), mixture=1) %>% #lasso
        set_engine("glmnet")

#workflow: to compare the model with and without down-sampling (for the unbalanced outcome)
wf_set <- 
        workflow_set(
                list(basic = hrqol_recipe,
                     downsampling = hrqol_recipe %>% step_downsample(child_pedsql_health_related_qol_score_total_catg_Q15)),
                list(glmnet = glmnet_spec)
        )
wf_set

#use tuning to evaluate different possible penalty values for each option, include several metrics to understand model performance
narrower_penality <- penalty(range=c(-3,0))

doParallel::registerDoParallel()
set.seed(123)
tune_rs <- 
        workflow_map(
                wf_set,
                "tune_grid",
                resamples = hrqol_folds,
                grid = 15,
                metrics = metric_set(accuracy, mn_log_loss, sensitivity, specificity),
                param_info = parameters(narrower_penality)
        )

#evaluate models
autoplot(tune_rs) + theme(legend.position = "none")
dev.off()

rank_results(tune_rs, rank_metric = "sensitivity")
rank_results(tune_rs, rank_metric = "mn_log_loss")

##basic glment (no down-sampling)
basic_rs <-
        tune_rs %>%
        extract_workflow_set_result("basic_glmnet")

autoplot(basic_rs)

best_penalty <- 
        basic_rs %>%
        select_by_one_std_err(-penalty, metric = "mn_log_loss")

best_penalty

final_fit <-  
        wf_set %>% 
        extract_workflow("basic_glmnet") %>%
        finalize_workflow(best_penalty) %>% #penalty is based on mn_log_loss
        last_fit(hrqol_Q15_split) #fit the final best model to the training set and evaluate the test set

final_fit

#evaluate model
collect_metrics(final_fit)
#metric  .estimator .estimate .config             
#<chr>    <chr>          <dbl> <chr>               
#1 accuracy binary         0.848 Preprocessor1_Model1
#2 roc_auc  binary         0.615 Preprocessor1_Model1

#model performance using a confusion matrix
collect_predictions(final_fit) %>%
        conf_mat(child_pedsql_health_related_qol_score_total_catg_Q15, .pred_class)
#                Truth
#Prediction   0_poor_avg 1_good
# 0_poor_avg          0      0
# 1_good            247   1374
#0/(0+246)*100: predicts the poor to average group with 0% accuracy (sensitivity)
#1374/(1374+0)*100: predicts the good group with 100% accuracy (specificity)
#CONCLUSION: very poor performance on predicting the poor_avg group; thus low sensitivity 

#down-sampling glmnet
downsample_rs <-
        tune_rs %>%
        extract_workflow_set_result("downsampling_glmnet")

autoplot(downsample_rs)
dev.off()

best_penalty <- 
        downsample_rs %>%
        select_by_one_std_err(-penalty, metric = "mn_log_loss")

best_penalty

final_fit <-  
        wf_set %>% 
        extract_workflow("downsampling_glmnet") %>%
        finalize_workflow(best_penalty) %>%
        last_fit(hrqol_Q15_split)

final_fit

#evaluate model
collect_metrics(final_fit)
#1 accuracy binary         0.611 Preprocessor1_Model1
#2 roc_auc  binary         0.612 Preprocessor1_Model1

#model performance
collect_predictions(final_fit) %>%
        conf_mat(child_pedsql_health_related_qol_score_total_catg_Q15, .pred_class)
`               #Truth
#Prediction   0_poor_avg 1_good
#0_poor_avg        123    506
#1_good            124    868

#123/(123+124)*100: predicts the poor to average group with 49.8% accuracy (sensitivity)
#868/(868+506)*100: predicts the good group with 63.2% accuracy (specificity)
#CONCLUSION: predictions are better with down-sampling

#which predictors are driving the classification model?
hrqol_vip <-
        extract_fit_engine(final_fit) %>%
        vi()

hrqol_vip
#Variable                                              Importance Sign 
#<chr>                                                      <dbl> <chr>
#1 child_demo_school_region_eastern_macedonia_and_thrace    0.498   POS  
#2 child_phyact_screen_time_hrs_per_wk                      0.446   NEG  
#3 child_demo_sex_X1_female                                 0.364   POS  
#4 family_demo_FAS_group_X2_high                            0.250   POS  
#5 child_ffq_pattern_5_starchy_foods_sweetened_bev          0.240   NEG  
#6 child_demo_school_year_X2017_2018                        0.225   NEG  
#7 child_ffq_pattern_3_fruits_raw_veg_cheese                0.179   POS  
#8 child_demo_immigrant_status_X2_2nd_gen_immigrant         0.128   NEG  
#9 child_demo_school_region_central_macedonia               0.119   NEG  
#10 family_demo_FAS_group_X1_middle                          0.0953  POS  
#11 child_anthro_bmi_WHO_sds                                 0.0935  NEG  
#12 child_ffq_pattern_4_confectioneries_pizza                0.0778  NEG  
#13 child_demo_age_years                                     0.0748  NEG  
#14 child_ffq_pattern_1_meat_seafood_prepd_meals             0.0736  POS  
#15 child_ffq_pattern_2_cooked_veg_grains_legumes            0.0707  NEG  
#16 child_ffq_animal_score                                   0.0619  POS  
#17 child_ffq_hpdi_score                                     0.00268 NEG  
#18 child_ffq_diet_quality_score                             0       NEG

#variables do not have a very high importance (%) in relation to QoL - we are already investigating child_demo_school_region

pdf(file = "2023_05_09_predictors_of_hrqol.pdf", useDingbats = F, onefile = T, width = 20)
hrqol_vip %>%
        group_by(Sign) %>%
        slice_max(Importance, n = 15) %>%
        ungroup() %>%
        ggplot(aes(Importance, fct_reorder(Variable, Importance), fill = Sign)) + 
        geom_col() +
        facet_wrap(vars(Sign), scales = "free_y") +
        labs(y = NULL) +
        theme(legend.position = "none") +
        scale_colour_brewer(palette = "Set2") + 
        scale_fill_brewer(palette = "Set2")

        #child_demo_school_region
masterfile %>%
        filter(!is.na(child_pedsql_health_related_qol_score_total_catg_Q15) &
                       !is.na(child_demo_school_region)) %>%
        count(child_demo_school_region)

test <- masterfile %>%
        filter(!is.na(child_pedsql_health_related_qol_score_total_catg_Q15) &
                       !is.na(child_demo_school_region)) %>%
        select(c(child_pedsql_health_related_qol_score_total_catg_Q15, child_demo_school_region)) %>%
        group_by(child_pedsql_health_related_qol_score_total_catg_Q15,child_demo_school_region) %>%
        summarize(n=n()) %>%
        group_by(child_demo_school_region) %>% mutate(percent = prop.table(n)*100)

test %>% ggplot(aes(fill=child_pedsql_health_related_qol_score_total_catg_Q15, y=percent, x=child_demo_school_region)) + 
        geom_bar(position="stack", stat="identity") +
        geom_text(aes(label = round(percent, digits = 1 )), vjust=-0.5) +
        theme_bw() +
        labs(x ="child_demo_school_region", y = "percentage (%)", fill = "Total HRQoL") +
        scale_colour_brewer(palette = "Set2") + 
        scale_fill_brewer(palette = "Set2") +
        scale_x_discrete(labels=c('attica_and_central_greece' = 'attica_and_central_greece\n N = 4089',
                                  'central_macedonia' = 'central_macedonia \n N = 1857',
                                  'eastern_macedonia_and_thrace' = 'eastern_macedonia_and_thrace\n N = 382'))

#use to investigate numbers only (do not save)
#test %>%
#        ggplot(aes(x=child_demo_school_region, y=n, fill = child_pedsql_health_related_qol_score_total_catg_Q15, group = child_pedsql_health_related_qol_score_total_catg_Q15)) + 
#        geom_bar(stat="identity", position = "dodge") +
#        geom_text(aes(label = n), vjust=-0.5) +
#        theme_bw() +
#        labs(x ="child_demo_school_region", y = "count (n)", fill = "Total HRQoL") +
#        scale_colour_brewer(palette = "Set2") + 
#        scale_fill_brewer(palette = "Set2")

#child_phyact_screen_time_hrs_per_wk
masterfile$child_phyact_screen_time_hrs_per_wk <- as.factor(masterfile$child_phyact_screen_time_hrs_per_wk)
masterfile %>%
        filter(!is.na(child_pedsql_health_related_qol_score_total_catg_Q15) &
                       !is.na(child_phyact_screen_time_hrs_per_wk)) %>%
        count(child_phyact_screen_time_hrs_per_wk)

test <- masterfile %>%
        filter(!is.na(child_pedsql_health_related_qol_score_total_catg_Q15) &
                       !is.na(child_phyact_screen_time_hrs_per_wk)) %>%
        select(c(child_pedsql_health_related_qol_score_total_catg_Q15, child_phyact_screen_time_hrs_per_wk)) %>%
        group_by(child_pedsql_health_related_qol_score_total_catg_Q15,child_phyact_screen_time_hrs_per_wk) %>%
        summarize(n=n()) %>%
        group_by(child_phyact_screen_time_hrs_per_wk) %>% mutate(percent = prop.table(n)*100)

test %>% ggplot(aes(fill=child_pedsql_health_related_qol_score_total_catg_Q15, y=percent, x=child_phyact_screen_time_hrs_per_wk)) + 
        geom_bar(position="stack", stat="identity") +
        geom_text(aes(label = round(percent, digits = 1 )), vjust=-0.5) +
        theme_bw() +
        labs(x ="child_phyact_screen_time_hrs_per_wk", y = "percentage (%)", fill = "Total HRQoL") +
        scale_colour_brewer(palette = "Set2") + 
        scale_fill_brewer(palette = "Set2") +
        scale_x_discrete(labels=c('0.003' = '0-1 hrs week\n N = 357',
                                  '0.009' = '1-2 hrs week \n N = 690',
                                  '0.02' = '2-4 hrs week\n N = 1034',
                                  '0.03' = '4-6 hrs week\n N = 1168',
                                  '0.04' = '6-8 hrs week \n N = 981',
                                  '0.05' = '8-10 hrs week\n N = 829',
                                  '0.06' = '>10 hrs week\n N = 1136'))

#use to investigate numbers only (do not save)
#test %>%
#        ggplot(aes(x=child_phyact_screen_time_hrs_per_wk, y=n, fill = child_pedsql_health_related_qol_score_total_catg_Q15, group = child_pedsql_health_related_qol_score_total_catg_Q15)) + 
#        geom_bar(stat="identity", position = "dodge") +
#        geom_text(aes(label = n), vjust=-0.5) +
#        theme_bw() +
#        labs(x ="child_phyact_screen_time_hrs_per_wk", y = "count (n)", fill = "Total HRQoL") +
#        scale_colour_brewer(palette = "Set2") + 
#        scale_fill_brewer(palette = "Set2")

dev.off()

rm(basic_rs, best_penalty, downsample_rs, final_fit, glmnet_spec, 
   hrqol_Q15_test, hrqol_folds, hrqol_Q15_split, hrqol_Q15_train, hrqol_recipe, hrqol_vip, 
   narrower_penality, test, tune_rs, wf_set)

        #### =============== 4. FEATURE ENGINEERING: HRQoL CONTINOUS PREDICTION MODELS   ================ ####
hrqol_perc_split <- masterfile %>% 
  filter(!is.na(child_pedsql_health_related_qol_score_total_perc)) %>%
  select(c(child_demo_school_year,
           child_demo_school_region,
           child_demo_age_years,
           child_demo_sex, 
           child_anthro_bmi_WHO_sds,
           child_demo_immigrant_status,
           family_demo_FAS_group,  
           #child_phyact_exercise_hrs_per_wk, #one of the components of QoL is Physical Functioning
           child_phyact_screen_time_hrs_per_wk, 
           child_ffq_hpdi_score, 
           child_ffq_animal_score,
           child_ffq_diet_quality_score,
           child_ffq_pattern_1_meat_seafood_prepd_meals,
           child_ffq_pattern_2_cooked_veg_grains_legumes,
           child_ffq_pattern_3_fruits_raw_veg_cheese,
           child_ffq_pattern_4_confectioneries_pizza,
           child_ffq_pattern_5_starchy_foods_sweetened_bev,
           child_pedsql_health_related_qol_score_total_perc)) 
#6481 17
str(hrqol_perc_split)

set.seed(123)
#split the data into training and test sets, and make sure that the proportion of poor_avg and good hrqol are equally distributed between these sets
hrqol_perc_split <- hrqol_perc_split %>%
  initial_split(strata = child_pedsql_health_related_qol_score_total_perc)

hrqol_perc_split
#<Training/Testing/Total>
#<4859/1622/6481>

hrqol_perc_train <- training (hrqol_perc_split)
hrqol_perc_test <- testing (hrqol_perc_split)

#make recipe for tidy models
hrqol_recipe <- 
  recipe(child_pedsql_health_related_qol_score_total_perc ~ ., data = hrqol_perc_train) %>% #relate hrqol to all other variables in the dataframe
  step_impute_mean(all_numeric_predictors()) %>% #impute numeric variables (NOTE: imputation is needed to run the downstream logistic regression)
  step_impute_bag(all_nominal(), -all_outcomes()) %>% #impute categorical variables
  step_dummy(all_nominal(), -all_outcomes()) %>% #generate dummy variables for all categorical variables, excluding the outcome, hrqol
  step_nzv(all_predictors()) #remove near zero variance for all predictor (numerical and nominal) variables

#confirm that the recipe works
prep(hrqol_recipe)
#NOTE: Sparse, unbalanced variable filter removed: child_demo_immigrant_status_X1_1st_gen_immigrant

#specify linear regression model
lm_model <- linear_reg() %>% set_engine("lm")

#workflow:
lm_workflow <- 
  workflow() %>%
  add_model(lm_model) %>%
  add_recipe(hrqol_recipe)

#train linear regression model
lm_fit <- fit(lm_workflow, hrqol_perc_train)

#predict QoL for the test set and bind it to lm_fit
results <- hrqol_perc_test %>% 
  bind_cols(lm_fit %>% 
              #predict QoL
              predict(new_data = hrqol_perc_test) %>% 
              rename(predictions = .pred))

#compare predictions
results %>% 
  select(c(child_pedsql_health_related_qol_score_total_perc, predictions)) %>% 
  slice_head(n = 10)

#visualise the results of the prediction
results %>% 
  ggplot(mapping = aes(x = child_pedsql_health_related_qol_score_total_perc, y = predictions)) +
  geom_point(size = 1.6) +
  #overlay a regression line
  geom_smooth(method = "lm", se = F) +
  ggtitle("Total HRQoL (percentiles) Predictions") +
  xlab("Observed") +
  ylab("Predicted") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#evaluate linear regression model
#multiple regression metrics
eval_metrics <- metric_set(rmse, rsq)

#evaluate RMSE and R2 based on the results
eval_metrics(data = results,
             truth = child_pedsql_health_related_qol_score_total_perc,
             estimate = predictions) %>% 
  select(-2)

#A tibble: 2 × 2
#.metric .estimate
#<chr>       <dbl>
#1 rmse      29.2  #on average, predictions are wrong by about 30 percentiles 
#2 rsq        0.0473 #only 4% of the variation in HRQoL percentile is explained by the predictors
