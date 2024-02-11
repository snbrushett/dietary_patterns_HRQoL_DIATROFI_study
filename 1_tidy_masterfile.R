###     AUTHOR(S): SIOBHAN BRUSHETT 
###     DESCRIPTION: FURTHER TIDY MASTERFILE AND ADD HRQOL PERCENTILES
###     PROJECT: DIATROFI_BASELINE_2015_2018

###     LAST UPDATED: 9 MAY 2023

#libraries
library(tidyverse)
library(skimr)

## 0. LOAD DATA
## 1. ADD MARGARINE INTAKE TO MASTERFILE (FOR ADJUSTMENT IN REGRESSION ANALYSES)
## 2. TIDY MASTERFILE
## 3. GENERATE CONTINUOUS QOL USING PERCENTILES
## 4. WRITE FILE

        #### =============== 0. LOAD AND UPDATE DATA FRAME   ================ ####
masterfile <- read.delim("2023_03_20_diatrofi_2015_2018_masterfile.txt", sep = "\t", header = T, stringsAsFactors = T)
str(masterfile)
#6583 59

skim(masterfile)
summary(masterfile)

        #### =============== 1. ADD MARGARINE INTAKE TO MASTERFILE (FOR ADJUSTMENT IN REGRESSION ANALYSES)   ================ ####
#NOTE: because margarine was not included in the indexes, margarine needs to be adjusted for in regression analyses
dia_merged <- read.delim("diatrofi_merged_2015_2018_baseline.txt", sep = "\t", header = T, stringsAsFactors = T)
#6583  117
row.names(dia_merged) <- dia_merged$pseudo_ids

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

#clean ffq by exclusion criteria based on the following study: https://doi.org/10.1186/s12966-016-0353-2 
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
sum(dia_merged_ffq$count_na>67.5) #participants with more than 90% of the items indicated recorded as '7_times_per_day_more_than_2'
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

#derive the weighted frequencies of the ffq items
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

summary(dia_merged_ffq)

#subset only for margarine
dia_merged_ffq$pseudo_ids <- rownames(dia_merged_ffq)
dia_merged_ffq <- dia_merged_ffq[,c("pseudo_ids",
                                    "child_ffq_margarine")]

masterfile_updated <- left_join(masterfile, dia_merged_ffq) #joined by pseudo_ids
masterfile_updated <- masterfile_updated %>%
  relocate(child_ffq_margarine, .before = child_pedsql_phy_score) 
str(masterfile_updated)

rm(dia_merged, dia_merged_ffq, masterfile)

        #### =============== 2. TIDY MASTERFILE   ================ ####
##child_demo_school_year##
table(masterfile_updated$child_demo_school_year, useNA="ifany")

#combine central_greece and attica into one category 
masterfile_updated$child_demo_school_year <- as.character(masterfile_updated$child_demo_school_year)
masterfile_updated$child_demo_school_year[masterfile_updated$child_demo_school_year=="2015-2016"] <- "2015_2016"
masterfile_updated$child_demo_school_year[masterfile_updated$child_demo_school_year=="2017-2018"] <- "2017_2018" 
masterfile_updated$child_demo_school_year <- as.factor(masterfile_updated$child_demo_school_year)

##child_demo_school_region##
table(masterfile_updated$child_demo_school_region, useNA="ifany")

#combine central_greece and attica into one category 
masterfile_updated$child_demo_school_region <- as.character(masterfile_updated$child_demo_school_region)
masterfile_updated$child_demo_school_region[masterfile_updated$child_demo_school_region=="attica" |
                                              masterfile_updated$child_demo_school_region=="central_greece"] <- "attica_and_central_greece"
#to avoid skewing data
masterfile_updated$child_demo_school_region[masterfile_updated$child_demo_school_region=="peloponnese" | #only 27 participants
                                              masterfile_updated$child_demo_school_region=="thessaly" | #only 67 participants
                                              masterfile_updated$child_demo_school_region=="western_greece"] <- NA #only 60 participants
masterfile_updated$child_demo_school_region <- as.factor(masterfile_updated$child_demo_school_region)

##child_demo_insurance##
table(masterfile_updated$child_demo_insurance, useNA="ifany")
masterfile_updated$child_demo_insurance <- as.character(masterfile_updated$child_demo_insurance)
masterfile_updated$child_demo_insurance[masterfile_updated$child_demo_insurance=="unknown"] <- NA
masterfile_updated$child_demo_insurance[masterfile_updated$child_demo_insurance!="not_covered"] <- "1_covered"
masterfile_updated$child_demo_insurance[masterfile_updated$child_demo_insurance=="not_covered"] <- "0_not_covered"
masterfile_updated$child_demo_insurance <- as.factor(masterfile_updated$child_demo_insurance)

##child_phyact_exercise_hrs_per_wk##
table(masterfile_updated$child_phyact_exercise_hrs_per_wk, useNA="ifany")
#there are 168 hrs in a week

masterfile_updated$child_phyact_exercise_hrs_per_wk <- recode_factor(masterfile_updated$child_phyact_exercise_hrs_per_wk, 
                                                                     "0_hrs_week_0_1"="0.003", #0.5/168
                                                                     "1_hrs_week_1_2"="0.009", #1.5/168
                                                                     "2_hrs_week_2_4"="0.02", #3/168
                                                                     "3_hrs_week_4_6"="0.03", #5/168
                                                                     "4_hrs_week_6_8"="0.04", #7/168
                                                                     "5_hrs_week_8_10"="0.05", #9/168
                                                                     "6_hrs_week_more_than_10"="0.06") #10/168
masterfile_updated$child_phyact_exercise_hrs_per_wk <- as.numeric(as.character(masterfile_updated$child_phyact_exercise_hrs_per_wk))

##child_phyact_screen_time_hrs_per_wk##
table(masterfile_updated$child_phyact_screen_time_hrs_per_wk, useNA="ifany")
#there are 168 hrs in a week

masterfile_updated$child_phyact_screen_time_hrs_per_wk <- recode_factor(masterfile_updated$child_phyact_screen_time_hrs_per_wk, 
                                                                        "0_hrs_week_0_1"="0.003", #0.5/168
                                                                        "1_hrs_week_1_2"="0.009", #1.5/168
                                                                        "2_hrs_week_2_4"="0.02", #3/168
                                                                        "3_hrs_week_4_6"="0.03", #5/168
                                                                        "4_hrs_week_6_8"="0.04", #7/168
                                                                        "5_hrs_week_8_10"="0.05", #9/168
                                                                        "6_hrs_week_more_than_10"="0.06") #10/168
masterfile_updated$child_phyact_screen_time_hrs_per_wk <- as.numeric(as.character(masterfile_updated$child_phyact_screen_time_hrs_per_wk))

##child_phyact_screen_time_hrs_per_wk##
table(masterfile_updated$child_phyact_play_outside_days_per_wk, useNA="ifany")

masterfile_updated$child_phyact_play_outside_days_per_wk <- recode_factor(masterfile_updated$child_phyact_play_outside_days_per_wk, 
                                                                          "0_never_rarely"="0.5", 
                                                                          "1_times_week_1_2"="1.5", 
                                                                          "2_times_week_3_4"="3.5", 
                                                                          "3_times_week_5-6"="5.5", 
                                                                          "4_daily"="7") 
masterfile_updated$child_phyact_play_outside_days_per_wk <- as.numeric(as.character(masterfile_updated$child_phyact_play_outside_days_per_wk))

str(masterfile_updated)

        #### =============== 3. GENERATE CONTINUOUS QOL USING PERCENTILES   ================ ####
#calculate the percentile for each value of QoL (using the cumulative distribution function)
masterfile_updated$child_pedsql_phy_score_perc <- ecdf(masterfile_updated$child_pedsql_phy_score)(masterfile_updated$child_pedsql_phy_score)*100
masterfile_updated$child_pedsql_emo_score_perc <- ecdf(masterfile_updated$child_pedsql_emo_score)(masterfile_updated$child_pedsql_emo_score)*100
masterfile_updated$child_pedsql_soc_score_perc <- ecdf(masterfile_updated$child_pedsql_soc_score)(masterfile_updated$child_pedsql_soc_score)*100
masterfile_updated$child_pedsql_sch_score_perc <- ecdf(masterfile_updated$child_pedsql_sch_score)(masterfile_updated$child_pedsql_sch_score)*100
masterfile_updated$child_pedsql_psysoc_score_perc <- ecdf(masterfile_updated$child_pedsql_psysoc_score)(masterfile_updated$child_pedsql_psysoc_score)*100
masterfile_updated$child_pedsql_health_related_qol_score_total_perc <- ecdf(masterfile_updated$child_pedsql_health_related_qol_score_total)(masterfile_updated$child_pedsql_health_related_qol_score_total)*100


pdf(file = "2023_05_09_hrqol_as_percentile.pdf", useDingbats = F, onefile = T, width = 20, height=12)
ggplot(masterfile_updated, aes(child_pedsql_phy_score)) + stat_ecdf(geom = "step")
hist(masterfile_updated$child_pedsql_phy_score_perc)

ggplot(masterfile_updated, aes(child_pedsql_emo_score)) + stat_ecdf(geom = "step")
hist(masterfile_updated$child_pedsql_emo_score_perc)

ggplot(masterfile_updated, aes(child_pedsql_soc_score)) + stat_ecdf(geom = "step")
hist(masterfile_updated$child_pedsql_soc_score_perc)

ggplot(masterfile_updated, aes(child_pedsql_sch_score)) + stat_ecdf(geom = "step")
hist(masterfile_updated$child_pedsql_sch_score_perc)

ggplot(masterfile_updated, aes(child_pedsql_psysoc_score)) + stat_ecdf(geom = "step")
hist(masterfile_updated$child_pedsql_psysoc_score_perc)

ggplot(masterfile_updated, aes(child_pedsql_health_related_qol_score_total)) + stat_ecdf(geom = "step")
hist(masterfile_updated$child_pedsql_health_related_qol_score_total_perc)
dev.off()

percentile_qc <- masterfile_updated %>% select (c(
  child_pedsql_phy_score, #over 2000 participants with a score of 100, thus, the percentile conversion depicts this
  child_pedsql_phy_score_perc,
  child_pedsql_emo_score,
  child_pedsql_emo_score_perc,
  child_pedsql_soc_score,
  child_pedsql_soc_score_perc,
  child_pedsql_sch_score,
  child_pedsql_sch_score_perc,
  child_pedsql_psysoc_score,
  child_pedsql_psysoc_score_perc,
  child_pedsql_health_related_qol_score_total,
  child_pedsql_health_related_qol_score_total_perc))
#coincides with the ggplot above, however, total scores (psysoc and hrqol) fair better as percentiles than individual domains as there is quite a distinction
#in the scores of the individual domains - dichotomization is probably better for these domains
rm(percentile_qc)

        #### =============== 4. WRITE FILE   ================ ####
str(masterfile_updated)
#6583  66
masterfile_updated$child_demo_school_code <- as.factor(masterfile_updated$child_demo_school_code)
summary(masterfile_updated)

write.table(masterfile_updated, "2023_05_09_diatrofi_2015_2018_masterfile_updated.txt", sep="\t", row.names=F, quote = F)
