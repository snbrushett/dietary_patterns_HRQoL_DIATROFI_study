###     SCRIPT: DIATROFI DERIVE MICROBIOME RELATED FOOD SCORE
###     AUTHOR(S): SIOBHAN BRUSHETT 
###     DESCRIPTION: DERIVE DIETARY MICROBIOME RELATED FOOD SCORE FROM hPDI AND THE ANIMAL SCORE INDECES USED IN LITERATURE ON ASSOCIATIONS OF FOOD ITEMS WITH THE GUT MICROBIOME
###     PROJECT: MERGED DIATROFI 2015-2018
###     NOTE(S): NB!! Use in combination with '2023_02_22_dietary_diversity_microbiome_related.xlsx' spreadsheet
###              doi:10.1016/j.jacc.2017.05.047. (Plant Diversity Index paper)
###              https://www.nature.com/articles/s41591-020-01183-8 (microbiome related paper)

#libraries
library(tidyverse)
library(DataExplorer)
library(ggpubr)
library(lattice)
library(reshape2)

##Contents of this file

## 0. LOAD DATA
## 1. CLEAN FFQ ITEMS BY EXCLUSION CRITERIA
## 2. DERIVE FFQ WEIGHTED FREQUENCIES
## 3. DERIVE DIETARY DIVERSITY INDECES

#functions
#pearson correlation (Alex and Trishla)
pearson <- function(x,y) {
        matchID <- intersect(rownames(x), rownames(y))
        x1 <- x[matchID,]
        y1 <- y[matchID,]
        result_cor <- matrix(nrow=ncol(x1), ncol=ncol(y1))
        rownames(result_cor) <- colnames(x1)
        colnames(result_cor) <- colnames(y1)
        result_pvalue <- matrix(nrow=ncol(x1), ncol=ncol(y1))
        for (i in 1:ncol(y1)) {
                for (j in 1:ncol(x1)) {
                        cor1<-try(cor.test(x1[,j], y1[,i], method = "pearson"))
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
        correlation <- pearson(x,y)
        a<- melt(correlation$cor)
        a<- cbind(a, melt(correlation$p.val)[,"value"])
        result= a[order(a[,4]),]
        colnames(result)=c("factor1", "factor2", "CorCoefficient","pvalue")
        return(result)
}

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

summary(dia_merged_ffq)

        #### =============== 3. DERIVE FOOD GROUPS AND MICROBIOME RELATED FOOD SCORE   ================ ####

#NOTE: The only item that is not included is margarine.
#"margarine’s fatty acid composition has changed over time from high trans to high unsaturated fats, 
#we did not include these foods in the indices, but adjusted for them in the analysis" - article: doi:10.1016/j.jacc.2017.05.047.

food_groups <- dia_merged_ffq

        ####3.1. hPDI: more healthful ####
#whole-grains
food_groups$whole_grains <- (food_groups$child_ffq_bread_wholegrain + 
                                    food_groups$child_ffq_rice_wholegrain +
                                    food_groups$child_ffq_pasta_wholegrain)/3
#fruits
food_groups$fruits <- (food_groups$child_ffq_fruit_orange +
                               food_groups$child_ffq_fruit_apple_pear +
                               food_groups$child_ffq_fruit_winter_other +
                               food_groups$child_ffq_fruit_banana +
                               food_groups$child_ffq_fruit_summer_other +
                               food_groups$child_ffq_fruit_dried +
                               food_groups$child_ffq_olives)/7
#vegetables
food_groups$vegetables <- (food_groups$child_ffq_tomato_cucumber_carrot_pepper +
                                   food_groups$child_ffq_lettuce_cabbage_spinach_rocket +
                                   food_groups$child_ffq_broccoli_pumpkin_cauliflower +
                                   food_groups$child_ffq_leeks_spinach_celery +
                                   food_groups$child_ffq_grean_beans_peas_okra)/5
#nuts
food_groups$nuts <- food_groups$child_ffq_nuts

#legumes
food_groups$legumes <- food_groups$child_ffq_lentils_beans_chickpeas

#vegetable oils
food_groups$plant_oils <- (food_groups$child_ffq_olive_oil +
                                   food_groups$child_ffq_seed_oil)/2

#tea and coffee
food_groups$tea_coffee <- (food_groups$child_ffq_herbal_teas +
                                   food_groups$child_ffq_coffee)/2

#healthful plant-based (for individual hPDI score)
food_groups$misc_plant_healthful <- food_groups$child_ffq_veg_risotto_stuff_veg

        ####3.2. hPDI: less healthful ####
#fruit juices
food_groups$fruit_juices <- (food_groups$child_ffq_fruit_juice_fresh +
                                     food_groups$child_ffq_fruit_juice_standard)/2

#refined grains
food_groups$refined_grains <- (food_groups$child_ffq_bread_white +
                                       food_groups$child_ffq_bread_bagel_burger +
                                       food_groups$child_ffq_bread_crumbs +
                                       food_groups$child_ffq_toast_sandwich +
                                       food_groups$child_ffq_cereal_products +
                                       food_groups$child_ffq_rice_white +
                                       food_groups$child_ffq_pasta)/7

#potatoes
food_groups$potatoes <- (food_groups$child_ffq_potatoes +
                                 food_groups$child_ffq_french_fries)/2

#sugar sweetened beverages
food_groups$sugary_beverages <- (food_groups$child_ffq_soft_drinks +
                                         food_groups$child_ffq_soft_drinks_light+
                                         food_groups$child_ffq_energy_drinks)/3

#sweets and desserts
food_groups$sweets_desserts <- (food_groups$child_ffq_pastries +
                                        food_groups$child_ffq_dessert +
                                        food_groups$child_ffq_jelly_composte +
                                        food_groups$child_ffq_cakes +
                                        food_groups$child_ffq_chocolate_milk +
                                        food_groups$child_ffq_chocolate_dark +
                                        food_groups$child_ffq_honey_marmalade_sugar)/7

#less healthful plant-based (for individual hPDI score)
food_groups$misc_plant_less_healthful <- (food_groups$child_ffq_homemade_pies +
                                                  food_groups$child_ffq_standard_pies +
                                                  food_groups$child_ffq_crisps_popcorn)/3

        ####3.3. Animal Score: more healthful ####

#dairy
food_groups$dairy_healthful <- (food_groups$child_ffq_milk_yoghurt_full_fat +
                                        food_groups$child_ffq_milk_yoghurt_low_fat +
                                        food_groups$child_ffq_cottage_light_cheese_low_fat +
                                        food_groups$child_ffq_feta_anthotyros_cheese)/4
#meat
food_groups$meat_healthful <- food_groups$child_ffq_chicken_turkey

#eggs
food_groups$eggs <- food_groups$child_ffq_egg

#seafood
food_groups$seafood <- (food_groups$child_ffq_fish_small +
                                food_groups$child_ffq_fish_big +
                                food_groups$child_ffq_seafood)/3

#healthful animal-based (for individual animal score)
food_groups$misc_animal_healthful <- food_groups$child_ffq_mousaka

        ####3.4. Animal Score: less healthful ####

#animal fats
food_groups$animal_fats <- food_groups$child_ffq_butter

#meat
food_groups$meat_less_healthful <- (food_groups$child_ffq_beef +
                                            food_groups$child_ffq_burger_meatballs_minced_meat +
                                            food_groups$child_ffq_pork +
                                            food_groups$child_ffq_hamburger +
                                            food_groups$child_ffq_hotdog +
                                            food_groups$child_ffq_lamb_goat_ribs +
                                            food_groups$child_ffq_meat_cold_cuts +
                                            food_groups$child_ffq_sausage_bacon +
                                            food_groups$child_ffq_meat_cold_cuts_low_no_fat)/9
#dairy
food_groups$dairy_less_healthful <- (food_groups$child_ffq_yellow_or_cream_cheese +
                                             food_groups$child_ffq_ice_cream_milkshake)/2

#less healthful animal-based (for individual animal score)
food_groups$misc_animal_less_healthful <- (food_groups$child_ffq_pork_souvlaki_pita +
                                              food_groups$child_ffq_chicken_souvlaki_pita +
                                              food_groups$child_ffq_pizza +
                                              food_groups$child_ffq_chocolate_milk_drink +
                                              food_groups$child_ffq_mayo + 
                                              food_groups$child_ffq_mayo_light)/6

        ####3.5. Miscellaneous plant and animal based foods (for total score) ####

#more healthy
food_groups$misc_healthful <- (food_groups$child_ffq_veg_risotto_stuff_veg +
                                       food_groups$child_ffq_mousaka)/2

#less healthy
food_groups$misc_less_healthful <- (food_groups$child_ffq_pork_souvlaki_pita +
                                            food_groups$child_ffq_chicken_souvlaki_pita +
                                            food_groups$child_ffq_homemade_pies +
                                            food_groups$child_ffq_standard_pies +
                                            food_groups$child_ffq_pizza +
                                            food_groups$child_ffq_crisps_popcorn +
                                            food_groups$child_ffq_chocolate_milk_drink + 
                                            food_groups$child_ffq_mayo +
                                            food_groups$child_ffq_mayo_light)/9

        ####3.6. Dietary quality scores (hPDI, Animal score and Total score) ####

        #3.6.1 food group distributions
summary(food_groups)
food_groups <- as.data.frame(food_groups[,c("whole_grains",
                                            "fruits",
                                            "vegetables",
                                            "nuts",
                                            "legumes",
                                            "plant_oils",
                                            "tea_coffee",
                                            "fruit_juices",
                                            "misc_plant_healthful",
                                            "misc_animal_healthful",
                                            "misc_healthful",
                                            "refined_grains",
                                            "potatoes",
                                            "sugary_beverages",
                                            "sweets_desserts",
                                            "dairy_healthful",
                                            "meat_healthful",
                                            "eggs",
                                            "seafood",
                                            "animal_fats",
                                            "meat_less_healthful",
                                            "dairy_less_healthful",
                                            "misc_plant_less_healthful",
                                            "misc_animal_less_healthful",
                                            "misc_less_healthful")])
str(food_groups)
#6074  25

pdf(file = "2023_02_23_food_groups_hist_plots.pdf", useDingbats = F, onefile = T, width = 20, height=12)
food_groups %>% plot_histogram(nrow = 3L, ncol=2L)
dev.off()

food_groups <- as.data.frame(food_groups)
food_groups$pseudo_ids <- rownames(food_groups)
food_groups <- food_groups %>%                  
        mutate_if(is.numeric,
                  round,
                  digits = 2)
summary(food_groups)
head(food_groups)

food_groups_hpdi <- as.data.frame(food_groups[,c("pseudo_ids",
                                   "whole_grains",
                                   "fruits",
                                   "vegetables",
                                   "nuts",
                                   "legumes",
                                   "plant_oils",
                                   "tea_coffee",
                                   "misc_plant_healthful",
                                   "fruit_juices",
                                   "refined_grains", 
                                   "potatoes", 
                                   "sugary_beverages",
                                   "sweets_desserts",
                                   "misc_plant_less_healthful")])
#6074  15
str(food_groups_hpdi)
rownames(food_groups_hpdi) <- food_groups_hpdi$pseudo_ids
food_groups_hpdi$pseudo_ids <- NULL
head(food_groups_hpdi)

food_groups_animal <- as.data.frame(food_groups[,c("pseudo_ids",
                                                   "dairy_healthful",
                                                   "meat_healthful",
                                                   "eggs", 
                                                   "seafood",
                                                   "misc_animal_healthful",
                                                   "animal_fats",
                                                   "meat_less_healthful",
                                                   "dairy_less_healthful",
                                                   "misc_animal_less_healthful")])
str(food_groups_animal)
rownames(food_groups_animal) <- food_groups_animal$pseudo_ids
food_groups_animal$pseudo_ids <- NULL
head(food_groups_animal)

food_groups_total <- as.data.frame(food_groups[,c("pseudo_ids",
                                                  "whole_grains",
                                                  "fruits",
                                                  "vegetables",
                                                  "nuts",
                                                  "legumes",
                                                  "plant_oils",
                                                  "tea_coffee",
                                                  "dairy_healthful",
                                                  "meat_healthful",
                                                  "eggs",
                                                  "seafood",
                                                  "misc_healthful",
                                                  "fruit_juices",
                                                  "refined_grains",
                                                  "potatoes",
                                                  "sugary_beverages",
                                                  "sweets_desserts",
                                                  "animal_fats",
                                                  "meat_less_healthful",
                                                  "dairy_less_healthful",
                                                  "misc_less_healthful")])
rownames(food_groups_total) <- food_groups_total$pseudo_ids
food_groups_total$pseudo_ids <- NULL
head(food_groups_total)

        #3.6.2 hPDI SCORE: 

        #3.6.2.1 healthful products (n=8)
food_groups_hpdi_quant <- food_groups_hpdi 
healthful <- c("whole_grains",
               "fruits",
               "vegetables",
               "nuts",
               "legumes",
               "plant_oils",
               "tea_coffee",
               "misc_plant_healthful")

healthful <- which(colnames(food_groups_hpdi_quant) %in% healthful)

#generate quantiles
for (i in healthful){
        food_groups_hpdi_quant[,i] [food_groups_hpdi[,i] < quantile(food_groups_hpdi[,i], prob=0.20, type=1, na.rm=T)] <- 1 #Q1]
        food_groups_hpdi_quant[,i] [food_groups_hpdi[,i] >= quantile(food_groups_hpdi[,i], prob=0.20, type=1, na.rm=T) &
                                            food_groups_hpdi[,i] < quantile(food_groups_hpdi[,i], prob=0.40, type=1, na.rm=T)] <- 2 #Q1-Q2
        food_groups_hpdi_quant[,i] [food_groups_hpdi[,i] >= quantile(food_groups_hpdi[,i], prob=0.40, type=1, na.rm=T) &
                                            food_groups_hpdi[,i] <= quantile(food_groups_hpdi[,i], prob=0.60, type=1, na.rm=T)] <- 3 #Q2-Q3
        food_groups_hpdi_quant[,i] [food_groups_hpdi[,i] > quantile(food_groups_hpdi[,i], prob=0.60, type=1, na.rm=T) &
                                            food_groups_hpdi[,i]<=quantile(food_groups_hpdi[,i], prob=0.80, type=1, na.rm=T)] <- 4 #Q3-Q4
        food_groups_hpdi_quant[,i] [food_groups_hpdi[,i] > quantile(food_groups_hpdi[,i], prob=0.80, type=1, na.rm=T)] <- 5 #Q5
        
        print(table(food_groups_hpdi_quant[,i], useNA="ifany"))}

        #3.6.2.2 less healthful products (n=6)
less_healthful <- c("fruit_juices",
                    "refined_grains", 
                    "potatoes", 
                    "sugary_beverages",
                    "sweets_desserts",
                    "misc_plant_less_healthful")

less_healthful <- which(colnames(food_groups_hpdi_quant) %in% less_healthful)

#generate quantiles
for (i in less_healthful){
        food_groups_hpdi_quant[,i] [food_groups_hpdi[,i] < quantile(food_groups_hpdi[,i], prob=0.20, type=1, na.rm=T)] <- 5 #Q1]
        food_groups_hpdi_quant[,i] [food_groups_hpdi[,i] >= quantile(food_groups_hpdi[,i], prob=0.20, type=1, na.rm=T) &
                                            food_groups_hpdi[,i] < quantile(food_groups_hpdi[,i], prob=0.40, type=1, na.rm=T)] <- 4 #Q1-Q2
        food_groups_hpdi_quant[,i] [food_groups_hpdi[,i] >= quantile(food_groups_hpdi[,i], prob=0.40, type=1, na.rm=T) &
                                            food_groups_hpdi[,i] <= quantile(food_groups_hpdi[,i], prob=0.60, type=1, na.rm=T)] <- 3 #Q2-Q3
        food_groups_hpdi_quant[,i] [food_groups_hpdi[,i] > quantile(food_groups_hpdi[,i], prob=0.60, type=1, na.rm=T) &
                                            food_groups_hpdi[,i]<=quantile(food_groups_hpdi[,i], prob=0.80, type=1, na.rm=T)] <- 2 #Q3-Q4
        food_groups_hpdi_quant[,i] [food_groups_hpdi[,i] > quantile(food_groups_hpdi[,i], prob=0.80, type=1, na.rm=T)] <- 1 #Q5
        
        print(table(food_groups_hpdi_quant[,i], useNA="ifany"))}

head(food_groups_hpdi_quant)
food_groups_hpdi_quant$child_ffq_hpdi_score <- rowSums(food_groups_hpdi_quant) #score should range between 14-70
summary(food_groups_hpdi_quant$child_ffq_hpdi_score)

        #3.6.3 ANIMAL SCORE : 
        #3.6.3.1 healthful products (n=5)
food_groups_animal_quant <- food_groups_animal 
healthful <- c("dairy_healthful",
               "meat_healthful",
               "eggs", 
               "seafood",
               "misc_animal_healthful")

healthful <- which(colnames(food_groups_animal_quant) %in% healthful)

#generate quantiles
for (i in healthful){
        food_groups_animal_quant[,i] [food_groups_animal[,i] < quantile(food_groups_animal[,i], prob=0.20, type=1, na.rm=T)] <- 1 #Q1]
        food_groups_animal_quant[,i] [food_groups_animal[,i] >= quantile(food_groups_animal[,i], prob=0.20, type=1, na.rm=T) &
                                              food_groups_animal[,i] < quantile(food_groups_animal[,i], prob=0.40, type=1, na.rm=T)] <- 2 #Q1-Q2
        food_groups_animal_quant[,i] [food_groups_animal[,i] >= quantile(food_groups_animal[,i], prob=0.40, type=1, na.rm=T) &
                                              food_groups_animal[,i] <= quantile(food_groups_animal[,i], prob=0.60, type=1, na.rm=T)] <- 3 #Q2-Q3
        food_groups_animal_quant[,i] [food_groups_animal[,i] > quantile(food_groups_animal[,i], prob=0.60, type=1, na.rm=T) &
                                              food_groups_animal[,i]<=quantile(food_groups_animal[,i], prob=0.80, type=1, na.rm=T)] <- 4 #Q3-Q4
        food_groups_animal_quant[,i] [food_groups_animal[,i] > quantile(food_groups_animal[,i], prob=0.80, type=1, na.rm=T)] <- 5 #Q5
        
        print(table(food_groups_animal_quant[,i], useNA="ifany"))}

#3.6.3.2 less healthful products (n=4)
less_healthful <- c("animal_fats",
                    "meat_less_healthful",
                    "dairy_less_healthful",
                    "misc_animal_less_healthful")

less_healthful <- which(colnames(food_groups_animal_quant) %in% less_healthful)

#generate quantiles
for (i in less_healthful){
        food_groups_animal_quant[,i] [food_groups_animal[,i] < quantile(food_groups_animal[,i], prob=0.20, type=1, na.rm=T)] <- 5 #Q1]
        food_groups_animal_quant[,i] [food_groups_animal[,i] >= quantile(food_groups_animal[,i], prob=0.20, type=1, na.rm=T) &
                                              food_groups_animal[,i] < quantile(food_groups_animal[,i], prob=0.40, type=1, na.rm=T)] <- 4 #Q1-Q2
        food_groups_animal_quant[,i] [food_groups_animal[,i] >= quantile(food_groups_animal[,i], prob=0.40, type=1, na.rm=T) &
                                              food_groups_animal[,i] <= quantile(food_groups_animal[,i], prob=0.60, type=1, na.rm=T)] <- 3 #Q2-Q3
        food_groups_animal_quant[,i] [food_groups_animal[,i] > quantile(food_groups_animal[,i], prob=0.60, type=1, na.rm=T) &
                                              food_groups_animal[,i]<=quantile(food_groups_animal[,i], prob=0.80, type=1, na.rm=T)] <- 2 #Q3-Q4
        food_groups_animal_quant[,i] [food_groups_animal[,i] > quantile(food_groups_animal[,i], prob=0.80, type=1, na.rm=T)] <- 1 #Q5
        
        print(table(food_groups_animal_quant[,i], useNA="ifany"))}

head(food_groups_animal_quant)
food_groups_animal_quant$child_ffq_animal_score <- rowSums(food_groups_animal_quant) #score should range between 9-45
summary(food_groups_animal_quant$child_ffq_animal_score)

        #3.6.4 TOTAL hPDI AND ANIMAL SCORE DIET QUALITY INDEX: 
        #3.6.4.1 healthful products (n=12)
food_groups_total_quant <- food_groups_total 
healthful <- c("whole_grains",
               "fruits",
               "vegetables",
               "nuts",
               "legumes",
               "plant_oils",
               "tea_coffee",
               "dairy_healthful",
               "meat_healthful",
               "eggs",
               "seafood",
               "misc_healthful")

healthful <- which(colnames(food_groups_total_quant) %in% healthful)

#generate quantiles
for (i in healthful){
        food_groups_total_quant[,i] [food_groups_total[,i] < quantile(food_groups_total[,i], prob=0.20, type=1, na.rm=T)] <- 1 #Q1]
        food_groups_total_quant[,i] [food_groups_total[,i] >= quantile(food_groups_total[,i], prob=0.20, type=1, na.rm=T) &
                                             food_groups_total[,i] < quantile(food_groups_total[,i], prob=0.40, type=1, na.rm=T)] <- 2 #Q1-Q2
        food_groups_total_quant[,i] [food_groups_total[,i] >= quantile(food_groups_total[,i], prob=0.40, type=1, na.rm=T) &
                                             food_groups_total[,i] <= quantile(food_groups_total[,i], prob=0.60, type=1, na.rm=T)] <- 3 #Q2-Q3
        food_groups_total_quant[,i] [food_groups_total[,i] > quantile(food_groups_total[,i], prob=0.60, type=1, na.rm=T) &
                                             food_groups_total[,i]<=quantile(food_groups_total[,i], prob=0.80, type=1, na.rm=T)] <- 4 #Q3-Q4
        food_groups_total_quant[,i] [food_groups_total[,i] > quantile(food_groups_total[,i], prob=0.80, type=1, na.rm=T)] <- 5 #Q5
        
        print(table(food_groups_total_quant[,i], useNA="ifany"))}

#QC using fruits as an example
quantile(food_groups_total$fruits, prob=0.20, type=1, na.rm=T) #change prob= argument for different intervals
test <- food_groups_total
test$id <- rownames(food_groups_total)
test_2 <- food_groups_total_quant
test_2$id <- rownames(food_groups_total_quant)
test_2 <- test_2[, c("id", "fruits")]
test_2 <- as.data.frame(test_2 %>% rename(fruit_raw = fruits))
test <- as.data.frame(test[, c("id", "fruits")])
test <- left_join (test, test_2, by="id")
rm(test, test_2)

        #3.6.4.2 less healthful products (n=9)
less_healthful <- c("fruit_juices",
                    "refined_grains",
                    "potatoes",
                    "sugary_beverages",
                    "sweets_desserts",
                    "animal_fats",
                    "meat_less_healthful",
                    "dairy_less_healthful",
                    "misc_less_healthful")

less_healthful <- which(colnames(food_groups_total_quant) %in% less_healthful)

#generate quantiles
for (i in less_healthful){
        food_groups_total_quant[,i] [food_groups_total[,i] < quantile(food_groups_total[,i], prob=0.20, type=1, na.rm=T)] <- 5 #Q1]
        food_groups_total_quant[,i] [food_groups_total[,i] >= quantile(food_groups_total[,i], prob=0.20, type=1, na.rm=T) &
                                             food_groups_total[,i] < quantile(food_groups_total[,i], prob=0.40, type=1, na.rm=T)] <- 4 #Q1-Q2
        food_groups_total_quant[,i] [food_groups_total[,i] >= quantile(food_groups_total[,i], prob=0.40, type=1, na.rm=T) &
                                             food_groups_total[,i] <= quantile(food_groups_total[,i], prob=0.60, type=1, na.rm=T)] <- 3 #Q2-Q3
        food_groups_total_quant[,i] [food_groups_total[,i] > quantile(food_groups_total[,i], prob=0.60, type=1, na.rm=T) &
                                             food_groups_total[,i]<=quantile(food_groups_total[,i], prob=0.80, type=1, na.rm=T)] <- 2 #Q3-Q4
        food_groups_total_quant[,i] [food_groups_total[,i] > quantile(food_groups_total[,i], prob=0.80, type=1, na.rm=T)] <- 1 #Q5
        
        print(table(food_groups_total_quant[,i], useNA="ifany"))}

#QC using refined_grains as an example
quantile(food_groups_total$refined_grains, prob=0.20, type=1, na.rm=T) #change prob= argument for different intervals
test <- food_groups_total_quant
test$id <- rownames(food_groups_total_quant)
test_2 <- food_groups_total
test_2$id <- rownames(food_groups_total)
test_2 <- test_2[, c("id", "refined_grains")]
test_2 <- as.data.frame(test_2 %>% rename(refined_grains_raw = refined_grains))
test <- as.data.frame(test[, c("id", "refined_grains")])
test <- left_join (test, test_2, by="id")
rm(test, test_2)

head(food_groups_total_quant)
food_groups_total_quant$diet_quality_score <- rowSums(food_groups_total_quant) #score should range between 21-105
summary(food_groups_total_quant$diet_quality_score)

food_groups_hpdi_quant$pseudo_ids <- rownames(food_groups_hpdi_quant)
food_groups_hpdi_quant <- food_groups_hpdi_quant[,c("pseudo_ids",
                                                    "child_ffq_hpdi_score")]

food_groups_animal_quant$pseudo_ids <- rownames(food_groups_animal_quant)
food_groups_animal_quant <- food_groups_animal_quant[, c("pseudo_ids",
                                                         "child_ffq_animal_score")]

food_groups_total_quant$pseudo_ids <- rownames(food_groups_total_quant)
food_groups_total_quant <- food_groups_total_quant[,c("pseudo_ids",
                                                      "diet_quality_score")]
food_groups_total_quant <- food_groups_total_quant %>% rename(child_ffq_quality_score=diet_quality_score)

food_scores_all <- list(food_groups_hpdi_quant, food_groups_animal_quant,food_groups_total_quant)
food_scores_all <- food_scores_all %>% reduce(full_join)

child_age <- dia_merged[,c("pseudo_ids", "child_demo_age_years")]
child_age$child_demo_age_years_cat[child_age$child_demo_age_years<13] <- "0_children"
child_age$child_demo_age_years_cat[child_age$child_demo_age_years>12] <- "1_adolescents"
child_age$child_demo_age_years_cat <- as.factor(child_age$child_demo_age_years_cat)
table(child_age$child_demo_age_years_cat, useNA="ifany")


dia_merged_ffq$pseudo_ids <- row.names(dia_merged_ffq)
ffq_scores_age <- list(food_scores_all, dia_merged_ffq,child_age)
ffq_scores_age <- ffq_scores_age %>% reduce(full_join)
head(ffq_scores_age)

pdf(file = "2023_03_20_diet_quality_scores_hist_plots.pdf", useDingbats = F, onefile = T, width = 20, height=12)
hist(ffq_scores_age$child_ffq_hpdi_score)
gghistogram(ffq_scores_age, x = "child_ffq_hpdi_score",
            add = "mean", rug = TRUE,
            fill = "child_demo_age_years_cat", palette = c("#00AFBB", "#E7B800"),
            add_density = TRUE)
hist(ffq_scores_age$child_ffq_animal_score)
gghistogram(ffq_scores_age, x = "child_ffq_animal_score",
            add = "mean", rug = TRUE,
            fill = "child_demo_age_years_cat", palette = c("#00AFBB", "#E7B800"),
            add_density = TRUE)
hist(ffq_scores_age$child_ffq_quality_score)
gghistogram(ffq_scores_age, x = "child_ffq_quality_score",
            add = "mean", rug = TRUE,
            fill = "child_demo_age_years_cat", palette = c("#00AFBB", "#E7B800"),
            add_density = TRUE)
dev.off()

#this score should be positively associated/correlated with healthful food items and negatively correlated with less healthful food items
#thus, check pearson correlation

str(food_scores_all)
food_scores_all <- food_scores_all %>% rename(child_ffq_diet_quality_score=child_ffq_quality_score)
food_scores_all$pseudo_ids <- as.factor(food_scores_all$pseudo_ids)
#6074  4

write.table(food_scores_all, "2023_03_20_ffq_diet_quality_score.txt", sep="\t", row.names=F, quote = F)

        #### =============== 4. PEARSON CORRELATION OF SCORE WITH FFQ GROUPS AND ITEMS   ================ ####

        #FFQ GROUPS
#food_groups_quant$pseudo_ids <- NULL
#pearson_corr <- result(food_groups_quant, food_groups_quant)
#484 4

#pearson_FDR <- pearson_corr
#pearson_FDR$FDR<-p.adjust(pearson_FDR$pvalue, method = "BH")
#pearson_FDR <- pearson_FDR %>% filter(
#        FDR<0.05
#)
#474 5

#pdf(file = "2023_02_23_ffq_groups_micro_score_pearson_corr.pdf", useDingbats = F, onefile = T, width = 15, height=18)
#ggplot(data = pearson_corr, aes(factor1, factor2, fill = CorCoefficient))+
#        geom_tile(color = "white")+
#        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                             midpoint = 0, limit = c(-1,1), space = "Lab", 
#                             name="Pearson\nCorrelation") +
#        theme_minimal()+ 
#        #theme(axis.text.x=element_blank(), #remove x axis labels
#        #      axis.ticks.x=element_blank(), #remove x axis ticks
#        #      axis.text.y=element_blank(),  #remove y axis labels
#        #      axis.ticks.y=element_blank()  #remove y axis ticks
#        #)+
#        coord_fixed() +
#        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#       theme(axis.text=element_text(size=8))

#dev.off()

##remove 1:1 correlations
#pearson_corr_FDR <- pearson_corr
#pearson_corr_FDR$FDR<-p.adjust(pearson_corr_FDR$pvalue, method = "BH")
#pearson_corr_FDR <- pearson_corr_FDR %>% 
#        as_tibble() %>% 
#        mutate(duplicates = if_else(factor1 == factor2,
#                                    TRUE,
#                                    FALSE)) %>% 
#        filter(duplicates == FALSE)

#pearson_corr_FDR$duplicates <- NULL 

#write.table(pearson_corr_FDR, "2023_02_23_ffq_item_micro_score_pearson_corr.txt", sep="\t", row.names=F, quote = F)
#5700  4

        #FFQ ITEMS
str(ffq_scores_age)
ffq_scores_corr <- ffq_scores_age %>% select (-c(pseudo_ids,
                                                 child_demo_age_years,
                                                 child_demo_age_years_cat))
str(ffq_scores_corr)
pearson_corr <- result(ffq_scores_corr, ffq_scores_corr)
#6084 4

pearson_FDR <- pearson_corr
pearson_FDR$FDR<-p.adjust(pearson_FDR$pvalue, method = "BH")
pearson_FDR <- pearson_FDR %>% filter(
        FDR<0.05
)
#5808 5

pearson_corr_plot <- pearson_corr
pearson_corr_plot$factor1 <- sub("child_ffq_", "", pearson_corr_plot$factor1)
pearson_corr_plot$factor2 <- sub("child_ffq_", "", pearson_corr_plot$factor2)

table(pearson_corr_plot$factor1, useNA="ifany")
pearson_corr_plot$factor1 <- as.character(pearson_corr_plot$factor1)
pearson_corr_plot$factor1 <- recode_factor(pearson_corr_plot$factor1,
                                           "hpdi_score"="0_hpdi_score",
                                           "animal_score"="0_animal_score",
                                           "quality_score"="0_total_dietary_quality_score",
                                           "bread_white"="1_bread_white",
                                           "bread_bagel_burger"="1_bread_bagel_burger",
                                           "bread_crumbs"="1_bread_crumbs",
                                           "toast_sandwich"="1_toast_sandwich",
                                           "cereal_products"="1_cereal_products",
                                           "rice_white"="1_rice_white",
                                           "pasta"="1_pasta",
                                           "potatoes"="1_potatoes",
                                           "french_fries"="1_french_fries",
                                           "yellow_or_cream_cheese"="1_yellow_or_cream_cheese",
                                           "ice_cream_milkshake"="1_ice_cream_milkshake",
                                           "beef"="1_beef",
                                           "burger_meatballs_minced_meat"="1_burger_meatballs_minced_meat",
                                           "pork"="1_pork",
                                           "pork_souvlaki_pita"="1_pork_souvlaki_pita",
                                           "chicken_souvlaki_pita"="1_chicken_souvlaki_pita",
                                           "hamburger"="1_hamburger",
                                           "hotdog"="1_hotdog",
                                           "lamb_goat_ribs"="1_lamb_goat_ribs",
                                           "meat_cold_cuts"="1_meat_cold_cuts",
                                           "sausage_bacon"="1_sausage_bacon",
                                           "meat_cold_cuts_low_no_fat"="1_meat_cold_cuts_low_no_fat",
                                           "homemade_pies"="1_homemade_pies",
                                           "standard_pies"="1_standard_pies",
                                           "pizza"="1_pizza",
                                           "pastries"="1_pastries",
                                           "dessert"="1_dessert",
                                           "jelly_composte"="1_jelly_composte",
                                           "cakes"="1_cakes",
                                           "chocolate_milk"="1_chocolate_milk",
                                           "chocolate_dark"="1_chocolate_dark",
                                           "honey_marmalade_sugar"="1_honey_marmalade_sugar",
                                           "crisps_popcorn"="1_crisps_popcorn",
                                           "fruit_juice_fresh"="1_fruit_juice_fresh",
                                           "fruit_juice_standard"="1_fruit_juice_standard",
                                           "soft_drinks"="1_soft_drinks",
                                           "soft_drinks_light"="1_soft_drinks_light",
                                           "energy_drinks"="1_energy_drinks",
                                           "chocolate_milk_drink"="1_chocolate_milk_drink",
                                           "mayo"="1_mayo",
                                           "mayo_light"="1_mayo_light",
                                           "butter"="1_butter",
                                           "bread_wholegrain"="2_bread_wholegrain",
                                           "rice_wholegrain"="2_rice_wholegrain",
                                           "pasta_wholegrain"="2_pasta_wholegrain",
                                           "tomato_cucumber_carrot_pepper"="2_tomato_cucumber_carrot_pepper",
                                           "lettuce_cabbage_spinach_rocket"="2_lettuce_cabbage_spinach_rocket",
                                           "broccoli_pumpkin_cauliflower"="2_broccoli_pumpkin_cauliflower",
                                           "leeks_spinach_celery"="2_leeks_spinach_celery",
                                           "grean_beans_peas_okra"="2_grean_beans_peas_okra",
                                           "fruit_orange"="2_fruit_orange",
                                           "fruit_apple_pear"="2_fruit_apple_pear",
                                           "fruit_winter_other"="2_fruit_winter_other",
                                           "fruit_banana"="2_fruit_banana",
                                           "fruit_summer_other"="2_fruit_summer_other",
                                           "fruit_dried"="2_fruit_dried",
                                           "olives"="2_olives",
                                           "lentils_beans_chickpeas"="2_lentils_beans_chickpeas",
                                           "nuts"="2_nuts",
                                           "milk_yoghurt_full_fat"="2_milk_yoghurt_full_fat",
                                           "milk_yoghurt_low_fat"="2_milk_yoghurt_low_fat",
                                           "cottage_light_cheese_low_fat"="2_cottage_light_cheese_low_fat",
                                           "feta_anthotyros_cheese"="2_feta_anthotyros_cheese",
                                           "egg"="2_egg",
                                           "chicken_turkey"="2_chicken_turkey",
                                           "fish_small"="2_fish_small",
                                           "fish_big"="2_fish_big",
                                           "seafood"="2_seafood",
                                           "veg_risotto_stuff_veg"="2_veg_risotto_stuff_veg",
                                           "mousaka"="2_mousaka",
                                           "coffee"="2_coffee",
                                           "herbal_teas"="2_herbal_teas",
                                           "olive_oil"="2_olive_oil",
                                           "seed_oil"="2_seed_oil",
                                           "margarine"="3_NA_margarine",
                                                     .ordered=T)
table(pearson_corr_plot$factor1, useNA="ifany")

table(pearson_corr_plot$factor2, useNA="ifany")
pearson_corr_plot$factor2 <- as.character(pearson_corr_plot$factor2)
pearson_corr_plot$factor2 <- recode_factor(pearson_corr_plot$factor2,
                                           "hpdi_score"="0_hpdi_score",
                                           "animal_score"="0_animal_score",
                                           "quality_score"="0_total_dietary_quality_score",
                                           "bread_white"="1_bread_white",
                                           "bread_bagel_burger"="1_bread_bagel_burger",
                                           "bread_crumbs"="1_bread_crumbs",
                                           "toast_sandwich"="1_toast_sandwich",
                                           "cereal_products"="1_cereal_products",
                                           "rice_white"="1_rice_white",
                                           "pasta"="1_pasta",
                                           "potatoes"="1_potatoes",
                                           "french_fries"="1_french_fries",
                                           "yellow_or_cream_cheese"="1_yellow_or_cream_cheese",
                                           "ice_cream_milkshake"="1_ice_cream_milkshake",
                                           "beef"="1_beef",
                                           "burger_meatballs_minced_meat"="1_burger_meatballs_minced_meat",
                                           "pork"="1_pork",
                                           "pork_souvlaki_pita"="1_pork_souvlaki_pita",
                                           "chicken_souvlaki_pita"="1_chicken_souvlaki_pita",
                                           "hamburger"="1_hamburger",
                                           "hotdog"="1_hotdog",
                                           "lamb_goat_ribs"="1_lamb_goat_ribs",
                                           "meat_cold_cuts"="1_meat_cold_cuts",
                                           "sausage_bacon"="1_sausage_bacon",
                                           "meat_cold_cuts_low_no_fat"="1_meat_cold_cuts_low_no_fat",
                                           "homemade_pies"="1_homemade_pies",
                                           "standard_pies"="1_standard_pies",
                                           "pizza"="1_pizza",
                                           "pastries"="1_pastries",
                                           "dessert"="1_dessert",
                                           "jelly_composte"="1_jelly_composte",
                                           "cakes"="1_cakes",
                                           "chocolate_milk"="1_chocolate_milk",
                                           "chocolate_dark"="1_chocolate_dark",
                                           "honey_marmalade_sugar"="1_honey_marmalade_sugar",
                                           "crisps_popcorn"="1_crisps_popcorn",
                                           "fruit_juice_fresh"="1_fruit_juice_fresh",
                                           "fruit_juice_standard"="1_fruit_juice_standard",
                                           "soft_drinks"="1_soft_drinks",
                                           "soft_drinks_light"="1_soft_drinks_light",
                                           "energy_drinks"="1_energy_drinks",
                                           "chocolate_milk_drink"="1_chocolate_milk_drink",
                                           "mayo"="1_mayo",
                                           "mayo_light"="1_mayo_light",
                                           "butter"="1_butter",
                                           "bread_wholegrain"="2_bread_wholegrain",
                                           "rice_wholegrain"="2_rice_wholegrain",
                                           "pasta_wholegrain"="2_pasta_wholegrain",
                                           "tomato_cucumber_carrot_pepper"="2_tomato_cucumber_carrot_pepper",
                                           "lettuce_cabbage_spinach_rocket"="2_lettuce_cabbage_spinach_rocket",
                                           "broccoli_pumpkin_cauliflower"="2_broccoli_pumpkin_cauliflower",
                                           "leeks_spinach_celery"="2_leeks_spinach_celery",
                                           "grean_beans_peas_okra"="2_grean_beans_peas_okra",
                                           "fruit_orange"="2_fruit_orange",
                                           "fruit_apple_pear"="2_fruit_apple_pear",
                                           "fruit_winter_other"="2_fruit_winter_other",
                                           "fruit_banana"="2_fruit_banana",
                                           "fruit_summer_other"="2_fruit_summer_other",
                                           "fruit_dried"="2_fruit_dried",
                                           "olives"="2_olives",
                                           "lentils_beans_chickpeas"="2_lentils_beans_chickpeas",
                                           "nuts"="2_nuts",
                                           "milk_yoghurt_full_fat"="2_milk_yoghurt_full_fat",
                                           "milk_yoghurt_low_fat"="2_milk_yoghurt_low_fat",
                                           "cottage_light_cheese_low_fat"="2_cottage_light_cheese_low_fat",
                                           "feta_anthotyros_cheese"="2_feta_anthotyros_cheese",
                                           "egg"="2_egg",
                                           "chicken_turkey"="2_chicken_turkey",
                                           "fish_small"="2_fish_small",
                                           "fish_big"="2_fish_big",
                                           "seafood"="2_seafood",
                                           "veg_risotto_stuff_veg"="2_veg_risotto_stuff_veg",
                                           "mousaka"="2_mousaka",
                                           "coffee"="2_coffee",
                                           "herbal_teas"="2_herbal_teas",
                                           "olive_oil"="2_olive_oil",
                                           "seed_oil"="2_seed_oil",
                                           "margarine"="3_NA_margarine",
                                           .ordered=T)
table(pearson_corr_plot$factor2, useNA="ifany")

pdf(file = "2023_03_20_ffq_items_dietary_scores_pearson_corr.pdf", useDingbats = F, onefile = T, width = 15, height=18)
ggplot(data = pearson_corr_plot, aes(factor1, factor2, fill = CorCoefficient))+
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                             midpoint = 0, limit = c(-1,1), space = "Lab", 
                             name="Pearson\nCorrelation") +
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

#remove 1:1 correlations
pearson_corr_FDR <- pearson_corr
pearson_corr_FDR$FDR<-p.adjust(pearson_corr_FDR$pvalue, method = "BH")
pearson_corr_FDR <- pearson_corr_FDR %>% 
        as_tibble() %>% 
        mutate(duplicates = if_else(factor1 == factor2,
                                    TRUE,
                                    FALSE)) %>% 
        filter(duplicates == FALSE)

pearson_corr_FDR$duplicates <- NULL 

write.table(pearson_corr_FDR, "2023_03_20_ffq_item_micro_score_pearson_corr.txt", sep="\t", row.names=F, quote = F)
#5700  4
