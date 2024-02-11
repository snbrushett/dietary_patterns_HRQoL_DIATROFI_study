###     SCRIPT: TABLE 1: DESCRIPTIVE STATISTICS
###     AUTHOR(S): SIOBHAN BRUSHETT 
###     DESCRIPTION: DESCRIPTIVE STATISTICS OF COVARIATES AND OUTCOMES
###     PROJECT: DIATROFI_BASELINE_2015_2018

###     LAST UPDATED: 16 MAY 2023

#libraries
library(skimr)
library(tidyverse)
library(DataExplorer)

## 0. LOAD DATA
## 1. TABLE 1: DESCRIPTIVE STATISTICS

        #### =============== 0. LOAD DATA ================ ####
#MAC-book location
setwd("~/OneDrive - UMCG/Prolepsis_internship/merged_2015_2018/masterfile/2_tidy_masterfile_incl_marg_add_hrqol_percentiles/")
masterfile <- read.delim("2023_05_09_diatrofi_2015_2018_masterfile_updated.txt", sep = "\t", header = T, stringsAsFactors = T)
#6583 66
str(masterfile)
masterfile$child_demo_school_code <- as.factor(masterfile$child_demo_school_code)

masterfile_covariates <- masterfile %>% select(
  c(child_demo_age_years, #(mostly) normally distributed
    child_demo_sex,
    family_demo_FAS_group,
    child_demo_immigrant_status,
    child_demo_school_region,
    child_phyact_screen_time_hrs_per_wk, #not normally distributed
    child_demo_school_year,
    child_anthro_bmi_WHO_sds)) #normally distributed

#overall summary stats (using skimr package)
skim(masterfile_covariates)

#investigate distributions (using DataExlporer package)
masterfile_covariates %>% plot_histogram(nrow = 4L, ncol=3L)
masterfile_covariates %>% plot_bar(nrow = 4L, ncol=2L)
dev.off()

        #### =============== 1. TABLE 1: DESCRIPTIVE STATISTICS   ================ ####
#NOTE: screen time was recoded as follows:
#there are 168 hrs in a week, thus:
#"0_hrs_week_0_1"="0.003", #0.5/168
#"1_hrs_week_1_2"="0.009", #1.5/168
#"2_hrs_week_2_4"="0.02", #3/168
#"3_hrs_week_4_6"="0.03", #5/168
#"4_hrs_week_6_8"="0.04", #7/168
#"5_hrs_week_8_10"="0.05", #9/168
#"6_hrs_week_more_than_10"="0.06") #10/168

        ####1.1 Total summary statistics####
IQR(masterfile$child_pedsql_health_related_qol_score_total, na.rm = T)

summary(masterfile)
prop.table(table(masterfile$child_demo_school_year))*100
prop.table(table(masterfile$child_demo_school_region))*100
sd(masterfile$child_demo_age_years)
prop.table(table(masterfile$child_demo_sex))*100
prop.table(table(masterfile$child_anthro_bmi_WHO_sds_grouped))*100
prop.table(table(masterfile$child_demo_immigrant_status))*100

#median screen time: 0.03 --> 5/168; thus 5 hrs a week
IQR(masterfile$child_phyact_screen_time_hrs_per_wk, na.rm = T)
#0.03 --> 5/168; thus around +/- 5 hrs

prop.table(table(masterfile$family_demo_FAS_group))*100

        ####1.2 HRQoL ####
#given that the outcome is not normally distributed, non-parametric tests will be used.
#For independent variables with 2 categorical levels, the Wilcoxon rank-sum test was applied (the non-parametric alternative of t-test)
#For independent variables with 3 or more categorical levels, the Kruskal-Wallis test was applied (the non-parametric alternative of one-way anova)

#school year
masterfile %>%
  #filter for complete cases
  filter(is.na(child_pedsql_health_related_qol_score_total) == F) %>%
  #group by child_pedsql_health_related_qol_score_total_catg_Q15
  group_by(child_demo_school_year) %>%
  #select numeric variables
  select (child_pedsql_health_related_qol_score_total) %>%
  #obtain summary statistis
  summarise_all(list(median = ~ median (., na.rm=T), 
                     IQR = ~ IQR (., na.rm = T)))

wilcox.test(masterfile$child_pedsql_health_related_qol_score_total ~ masterfile$child_demo_school_year)
#3.42e-12

#school region
masterfile %>%
  #filter for complete cases
  filter(is.na(child_pedsql_health_related_qol_score_total) == F) %>%
  #group by child_pedsql_health_related_qol_score_total_catg_Q15
  group_by(child_demo_school_region) %>%
  #select numeric variables
  select (child_pedsql_health_related_qol_score_total) %>%
  #obtain summary statistis
  summarise_all(list(median = ~ median (., na.rm=T), 
                     IQR = ~ IQR (., na.rm = T)))

kruskal.test(masterfile$child_pedsql_health_related_qol_score_total ~ masterfile$child_demo_school_region)
#0.04076

#sex
masterfile %>%
  #filter for complete cases
  filter(is.na(child_pedsql_health_related_qol_score_total) == F) %>%
  #group by child_pedsql_health_related_qol_score_total_catg_Q15
  group_by(child_demo_sex) %>%
  #select numeric variables
  select (child_pedsql_health_related_qol_score_total) %>%
  #obtain summary statistis
  summarise_all(list(median = ~ median (., na.rm=T), 
                     IQR = ~ IQR (., na.rm = T)))

wilcox.test(masterfile$child_pedsql_health_related_qol_score_total ~ masterfile$child_demo_sex)
#6.939e-08

#BMI grouped
masterfile %>%
  #filter for complete cases
  filter(is.na(child_pedsql_health_related_qol_score_total) == F) %>%
  #group by child_pedsql_health_related_qol_score_total_catg_Q15
  group_by(child_anthro_bmi_WHO_sds_grouped) %>%
  #select numeric variables
  select (child_pedsql_health_related_qol_score_total) %>%
  #obtain summary statistis
  summarise_all(list(median = ~ median (., na.rm=T), 
                     IQR = ~ IQR (., na.rm = T)))

kruskal.test(masterfile$child_pedsql_health_related_qol_score_total ~ masterfile$child_anthro_bmi_WHO_sds_grouped)
#1.768e-05

#immigrant status
masterfile %>%
  #filter for complete cases
  filter(is.na(child_pedsql_health_related_qol_score_total) == F) %>%
  #group by child_pedsql_health_related_qol_score_total_catg_Q15
  group_by(child_demo_immigrant_status) %>%
  #select numeric variables
  select (child_pedsql_health_related_qol_score_total) %>%
  #obtain summary statistis
  summarise_all(list(median = ~ median (., na.rm=T), 
                     IQR = ~ IQR (., na.rm = T)))

kruskal.test(masterfile$child_pedsql_health_related_qol_score_total ~ masterfile$child_demo_immigrant_status)
#6.182e-09

#FAS
masterfile %>%
  #filter for complete cases
  filter(is.na(child_pedsql_health_related_qol_score_total) == F) %>%
  #group by child_pedsql_health_related_qol_score_total_catg_Q15
  group_by(family_demo_FAS_group) %>%
  #select numeric variables
  select (child_pedsql_health_related_qol_score_total) %>%
  #obtain summary statistis
  summarise_all(list(median = ~ median (., na.rm=T), 
                     IQR = ~ IQR (., na.rm = T)))

kruskal.test(masterfile$child_pedsql_health_related_qol_score_total ~ masterfile$family_demo_FAS_group)
#5.946e-09