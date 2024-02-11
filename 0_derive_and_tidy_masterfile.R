###     SCRIPT: MASTERFILE 2015-2018 BASELINE DATA AND SUMMARY STATISTICS
###     AUTHOR(S): SIOBHAN BRUSHETT 
###     DESCRIPTION: TIDY AND COMBINE CLEANED DATAFRAMES TO GENERATE FINAL MASTERFILE OF DIATROFI 2015-2018 BASELINE DATA AND RUN SUMMARY STATISTICS
###     FOR FURTHER ANALYSIS
###     PROJECT: DIATROFI_BASELINE_2015_2018
###     NOTE(S): https://rstudio-pubs-static.s3.amazonaws.com/308486_16ebecee006a414a9de941fd104ac292.html
###              http://rpkgs.datanovia.com/ggpubr/reference/stat_compare_means.html

#libraries
library(tidyverse)
library(DataExplorer)
library(reshape2)
library(ggpubr)
library(patchwork)

##Contents of this file

## 0. LOAD DATA
## 1. DERIVE VARIABLES OF INTEREST
## 2. TIDY MASTERFILE 
## 3. SUMMARY STATISTICS AND CORRELATIONS
## 4. TABLE 2: QOL BY DIETARY INDEX AND PATTERNS WITH P-VALUES
## 5. PLOTS: QOL BY DIETARY INDEX AND PATTERNS WITH P-VALUES 


#functions
#summary stats function on GitHub: https://github.com/GRONINGEN-MICROBIOME-CENTRE/Lifelines_NEXT/blob/main/Create_summary_statistics_metadata.R
#authors: Paula Sureda and Arnau Vich
summary_statistics_metadata <- function (metadata_input, category_table) {
        
        # Packages needed       
        library (psych)  #describe r function
        
        # Create other functions to calculate the different parameters
        
        ## Categorical values - create function to calculate the counts and the percentage for categorical variables
        tblFun <- function(x) {
                # Create a table
                tbl <- table(x)
                # Combine columnes/rows to get the counts and percentage (creates new table -> res)
                res <- cbind(tbl,round(prop.table(tbl)*100,2))
                # Give names to the columns
                colnames(res) <- c('Count','Percentage')
                res
        }
        
        ## NA sum function - counts the number of NA
        nzsum <- function(x) {
                sum (is.na(x))
        }
        
        if (missing(category_table)) {
                
                ## Calculate table1 with the whole data:
                
                my_results = matrix(ncol = 9, nrow = ncol(metadata_input))
                
                for (k in 1:ncol(metadata_input)){
                        
                        if (is.numeric(metadata_input[,k])) {
                                # Keep in "x" the result from describe function (done in the columns) - for each factor  
                                x = describe(metadata_input[,k])
                                z = nzsum(metadata_input[,k])
                                # In the new table ("x"): keep different values in the different columns
                                my_results[k,1] = "numerical"
                                my_results[k,2] = x$median
                                my_results[k,3] = x$mean
                                my_results[k,4] = x$sd
                                my_results[k,5] = x$n
                                my_results[k,6] = z
                                my_results[k,7] = x$min
                                my_results[k,8] = x$max
                                my_results[k,9] = x$range
                        }
                        
                        # Condition: if the column values are categorical  
                        else {
                                # Keep in "x" the result from tblFun function (done in the columns) - for each factor
                                x = tblFun(metadata_input[,k])
                                z = nzsum(metadata_input[,k])
                                # In the new table ("x"): keep different values in the different columns 
                                my_results[k,1]="categorical"
                                # toString to keep the possible different values/categories in the same vector/column
                                my_results[k,2]=toString(rownames(x))
                                # First column table x = 'Count'
                                my_results[k,3]=toString(x[,1]) 
                                # Second column table x = 'Percentage'
                                my_results[k,4]=toString(x[,2])
                                # Sum of the values on column1 ("x")
                                my_results[k,5]=sum(x[,1])
                                my_results[k,6]= z
                                my_results[k,7] = NA
                                my_results[k,8] = NA
                                my_results[k,9] = NA
                        }
                }
                
                
                # The column names from the original table = row names from the new table 
                rownames(my_results) = colnames(metadata_input)
                # Give names to the columns of the new table 
                colnames(my_results) = c("Type", "Categories/Median", "Counts/Mean", "%/SD", "Number_non_zeros", "Number_NA",
                                         "Min", "Max", "Range") 
                
                # Export the new table
                write.table (my_results, file = "./meta_data_summary_stats.txt" , quote = F, sep = "\t")  
        }
        
}

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

        #### =============== 0. LOAD AND MERGE DATA   ================ ####

        #general data
dia_merged <- read.delim("diatrofi_merged_2015_2018_baseline.txt", sep = "\t", header = T, stringsAsFactors = T)
str(dia_merged)
#6583  117

        #general data: BMI SDS 
bmi_sds <- read.delim("2023_03_07_bmi_WHO_sds_with_cutoffs.txt", sep = "\t", header = T, stringsAsFactors = T)
str(bmi_sds)
#5056  3

        #predictors: dietary diversity 
#ffq_diversity <- read.delim("2023_03_02_ffq_shannon_diversity_entire_RA_grouped.txt", sep = "\t", header = T, stringsAsFactors = T)
#6074  4
#str(ffq_diversity)
#it was decided to not use dietary diversity because 1) the index does not work well with frequencies, for example, a
#frequency of 0 and a frequency of 1 for all items would result in the same diversity index, and 2) this concept is reflected
#better in the (validated) dietary index

        #predictors: dietary patterns 
ffq_patterns <- read.delim("2023_02_20_ffq_pca_varimax_regression_loadings.txt", sep = "\t", header = T, stringsAsFactors = T)
#6047  5
ffq_patterns$pseudo_ids <- rownames(ffq_patterns)
str(ffq_patterns)

        #predictors: dietary index/score
ffq_diet_score <- read.delim("2023_03_20_ffq_diet_quality_score.txt", sep = "\t", header = T, stringsAsFactors = T)
#6074 4
str(ffq_diet_score)

        #outcome: Q15 and Q25
qol_dich <- read.delim("2023_03_07_pedsql_dich_Q15_Q25.txt", sep = "\t", header = T, stringsAsFactors = T)
#6583  13

masterfile <- list(dia_merged, bmi_sds, ffq_patterns, ffq_diet_score, qol_dich) #ffq_diversity retrospectively removed
masterfile <- masterfile %>% reduce(full_join) #joined by pseudo_ids
#6583  139

        #### =============== 1. DERIVE VARIABLES OF INTEREST   ================ ####

        #parent employment
table(masterfile$mother_demo_occupation, useNA="ifany")
table(masterfile$father_demo_occupation, useNA="ifany")
masterfile$parent_demo_employment[masterfile$mother_demo_occupation=="employed" & masterfile$father_demo_occupation=="employed"] <- "0_both_parents_employed"
masterfile$parent_demo_employment[is.na(masterfile$parent_demo_employment) & 
                                      masterfile$mother_demo_occupation=="employed" | 
                                     is.na(masterfile$parent_demo_employment) & masterfile$father_demo_occupation=="employed"] <- "1_one_parent_employed"
masterfile$parent_demo_employment[is.na(masterfile$parent_demo_employment) & 
                                     masterfile$mother_demo_occupation=="unemployed_or_retired" &
                                                         masterfile$father_demo_occupation=="unemployed_or_retired"] <- "2_both_parents_unemployed"

masterfile$parent_demo_employment[is.na(masterfile$mother_demo_occupation) | 
                                     is.na(masterfile$father_demo_occupation)] <- NA

masterfile$parent_demo_employment <- as.factor(masterfile$parent_demo_employment)
table(masterfile$parent_demo_employment, useNA="ifany")
#0_both_parents_employed     1_one_parent_employed 2_both_parents_unemployed    <NA> 
#4006                      1521                       460                       596


#QC:
#parent_emp <- masterfile[,c("mother_demo_occupation", "father_demo_occupation", "parent_demo_employment")]
##variable as expected
#rm(parent_emp)

        #parent country of birth
table(masterfile$mother_demo_country_birth, useNA = "ifany")
table(masterfile$father_demo_country_birth, useNA = "ifany")

masterfile$parent_demo_country_birth[masterfile$mother_demo_country_birth=="Greece" & masterfile$father_demo_country_birth=="Greece"] <- "0_both_parents_Greek"
masterfile$parent_demo_country_birth[is.na(masterfile$parent_demo_country_birth) & 
                                     masterfile$mother_demo_country_birth=="Greece" | 
                                     is.na(masterfile$parent_demo_country_birth) & masterfile$father_demo_country_birth=="Greece"] <- "1_one_parent_Greek"
masterfile$parent_demo_country_birth[is.na(masterfile$parent_demo_country_birth) & 
                                     masterfile$mother_demo_country_birth!="Greece" &
                                     masterfile$father_demo_country_birth!="Greece"] <- "2_both_parents_not_Greek"

masterfile$parent_demo_country_birth[is.na(masterfile$mother_demo_country_birth) | 
                                     is.na(masterfile$father_demo_country_birth)] <- NA

masterfile$parent_demo_country_birth <- as.factor(masterfile$parent_demo_country_birth)
table(masterfile$parent_demo_country_birth, useNA="ifany")
#0_both_parents_Greek       1_one_parent_Greek 2_both_parents_not_Greek      <NA> 
#4159                      690                     1438                      296 

#QC
#parent_birth_country <- masterfile[,c("mother_demo_country_birth", "father_demo_country_birth", "parent_demo_country_birth")]
##variable as expected
#rm(parent_birth_country)

        #child sex
table(masterfile$child_demo_sex, useNA="ifany")
masterfile$child_demo_sex <- as.character(masterfile$child_demo_sex)
masterfile$child_demo_sex[masterfile$child_demo_sex==0] <- "0_male"
masterfile$child_demo_sex[masterfile$child_demo_sex==1] <- "1_female"
masterfile$child_demo_sex <- as.factor(masterfile$child_demo_sex)
table(masterfile$child_demo_sex)

        #child grade
masterfile$child_demo_grade <- recode_factor(masterfile$child_demo_grade, 
                                             "1_pre_kindergarten"="1_pre_primary_school",
                                             "2_kindergarten"="1_pre_primary_school",
                                             "3_1st_cl_primary"="2_primary_school",
                                             "4_2nd_cl_primary"="2_primary_school",
                                             "5_3rd_cl_primary"="2_primary_school",
                                             "6_4th_cl_primary"="2_primary_school",
                                             "7_5th_cl_primary"="2_primary_school",
                                             "8_6th_cl_primary"="2_primary_school",
                                             "9_1st_cl_lower_secondary"="3_secondary_school",
                                             "10_2nd_cl_lower_secondary"="3_secondary_school",
                                             "11_3rd_cl_lower_secondary"="3_secondary_school",
                                             "12_1st_cl_upper_secondary"="3_secondary_school",
                                             "13_2nd_cl_upper_secondary"="3_secondary_school",
                                             "14_3rd_cl_upper_secondary"="3_secondary_school",
                                             "15_4th_cl_upper_secondary_evening"="3_secondary_school",
                                             .ordered=T)
table(masterfile$child_demo_grade, useNA="ifany")

        #child immigration status
table(masterfile$child_demo_country_birth, useNA="ifany")
table(masterfile$parent_demo_country_birth, useNA="ifany")

masterfile$child_demo_immigrant_status[masterfile$child_demo_country_birth=="Greece" & masterfile$parent_demo_country_birth=="0_both_parents_Greek"] <- "0_non_immigrants"
masterfile$child_demo_immigrant_status[is.na(masterfile$child_demo_immigrant_status) & 
                                               masterfile$child_demo_country_birth!="Greece" &
                                               masterfile$parent_demo_country_birth!="0_both_parents_Greek"] <- "1_1st_gen_immigrant"
masterfile$child_demo_immigrant_status[is.na(masterfile$child_demo_immigrant_status) & 
                                               masterfile$child_demo_country_birth=="Greece" &
                                               masterfile$parent_demo_country_birth!="0_both_parents_Greek"] <- "2_2nd_gen_immigrant"
masterfile$child_demo_immigrant_status[is.na(masterfile$child_demo_country_birth) | 
                                               is.na(masterfile$parent_demo_country_birth)] <- NA

masterfile$child_demo_immigrant_status <- as.factor(masterfile$child_demo_immigrant_status)
table(masterfile$child_demo_immigrant_status, useNA="ifany")
#0_non_immigrants 1_1st_gen_immigrant 2_2nd_gen_immigrant      <NA> 
#4137                 312                1805                 329 

#QC
#child_immg <- masterfile[,c("child_demo_country_birth", "parent_demo_country_birth", "child_demo_immigrant_status")]
##variable as expected
#rm(child_immg)

        #child country of birth
table(masterfile$child_demo_country_birth, useNA = "ifany")
masterfile$child_demo_country_birth <- as.character(masterfile$child_demo_country_birth)

masterfile$child_demo_country_birth[masterfile$child_demo_country_birth!="Greece"] <- "1_not_Greece"
masterfile$child_demo_country_birth[masterfile$child_demo_country_birth=="Greece"] <- "0_Greece"

masterfile$child_demo_country_birth <- as.factor(masterfile$child_demo_country_birth)
table(masterfile$child_demo_country_birth, useNA="ifany")
#0_Greece 1_not_Greece        <NA> 
#6196         360              27

prop.table(table(masterfile$child_demo_country_birth))*100
#0_Greek 1_not_Greek 
#94.508847    5.491153

        #dichotomize FAS
table(masterfile$family_FAS_group, useNA="ifany")
#0_low 1_middle   2_high     <NA> 
#2005     3326      989      263
prop.table(table(masterfile$family_FAS_group, useNA = "ifany"))*100
#0_low          1_middle        2_high     <NA>
#30.457238      50.524077       15.023545  3.995139 

masterfile$family_FAS_group_dic[masterfile$family_FAS_group=="0_low"] <- "0_low"
masterfile$family_FAS_group_dic[masterfile$family_FAS_group=="1_middle" |
                                        masterfile$family_FAS_group=="2_high"] <- "1_mod_high"
masterfile$family_FAS_group_dic <- as.factor(masterfile$family_FAS_group_dic)
table(masterfile$family_FAS_group_dic, useNA="ifany")
#0_low 1_mod_high       <NA> 
#2005       4315        263

prop.table(table(masterfile$family_FAS_group_dic, useNA="ifany"))*100

        #### =============== 2. TIDY MASTERFILE   ================ ####
#col_tidy <-  colnames(masterfile)
#write.table(col_tidy, "colnames_masterfile_to_be_cleaned.txt", sep="\t", row.names=F, quote = F)
#rm(col_tidy)

#rename variables
head(masterfile)

masterfile <- masterfile %>% rename (
        family_demo_FAS_score = family_FAS_score,
        family_demo_FAS_group = family_FAS_group,
        family_demo_FAS_group_dic = family_FAS_group_dic,
        child_ffq_pattern_1_meat_seafood_prepd_meals = Factor1,
        child_ffq_pattern_2_cooked_veg_grains_legumes = Factor2,
        child_ffq_pattern_3_fruits_raw_veg_cheese = Factor3,
        child_ffq_pattern_4_confectioneries_pizza = Factor4,
        child_ffq_pattern_5_starchy_foods_sweetened_bev = Factor5
)
#6583  143

#reorder masterfile
masterfile <- masterfile[,c("pseudo_ids",
                            "child_demo_school_year",
                            "child_demo_school_code",
                            "child_demo_school_region",
                            "child_demo_grade",
                            "child_demo_age_months",
                            "child_demo_age_years",
                            "child_demo_sex",
                            "child_anthro_bmi_WHO_sds",
                            "child_anthro_bmi_WHO_sds_grouped",
                            "child_anthro_bmi_grouped_7",
                            "child_anthro_bmi_grouped_4",
                            "child_demo_country_birth",
                            "child_demo_immigrant_status",
                            "child_demo_lives_with",
                            "child_demo_insurance",
                            "child_phyact_sport_outside_school",
                            "child_phyact_exercise_hrs_per_wk",
                            "child_phyact_screen_time_hrs_per_wk",
                            "child_phyact_play_outside_days_per_wk",
                            "parent_demo_employment",
                            "parent_demo_country_birth",
                            "mother_demo_age_years",
                            "mother_demo_edu",
                            "father_demo_age_years",
                            "father_demo_edu",
                            "family_demo_FAS_score",
                            "family_demo_FAS_group",
                            "family_demo_FAS_group_dic",
                            "family_demo_children_in_family_n",
                            "family_demo_household_members_n",
                            "child_kidmed_score",
                            "child_kidmed_adherence_level",
                            "child_ffq_pattern_1_meat_seafood_prepd_meals",
                            "child_ffq_pattern_2_cooked_veg_grains_legumes",
                            "child_ffq_pattern_3_fruits_raw_veg_cheese",
                            "child_ffq_pattern_4_confectioneries_pizza",
                            "child_ffq_pattern_5_starchy_foods_sweetened_bev",
                            "child_ffq_hpdi_score",
                            "child_ffq_animal_score",
                            "child_ffq_diet_quality_score",
                            "child_pedsql_phy_score",
                            "child_pedsql_emo_score",
                            "child_pedsql_soc_score",
                            "child_pedsql_sch_score",
                            "child_pedsql_psysoc_score",
                            "child_pedsql_health_related_qol_score_total",
                            "child_pedsql_phy_score_catg_Q15",
                            "child_pedsql_emo_score_catg_Q15",
                            "child_pedsql_soc_score_catg_Q15",
                            "child_pedsql_sch_score_catg_Q15",
                            "child_pedsql_psysoc_score_catg_Q15",
                            "child_pedsql_health_related_qol_score_total_catg_Q15",
                            "child_pedsql_phy_score_catg_Q25",
                            "child_pedsql_emo_score_catg_Q25",
                            "child_pedsql_soc_score_catg_Q25",
                            "child_pedsql_sch_score_catg_Q25",
                            "child_pedsql_psysoc_score_catg_Q25",
                            "child_pedsql_health_related_qol_score_total_catg_Q25")]
#6583  59
masterfile$pseudo_ids <- as.factor(masterfile$pseudo_ids)
str(masterfile)

write.table(masterfile, "2023_03_20_diatrofi_2015_2018_masterfile.txt", sep="\t", row.names=T, quote = F)

        #### =============== 3. SUMMARY STATISTICS AND CORRELATIONS   ================ ####
#summary statistics
metadata_input <- masterfile
metadata_input$pseudo_ids <- NULL
dim(metadata_input)
#6583   58

summary_statistics_metadata(metadata_input)

#summary plots
pdf(file = "2023_03_20_diatrofi_2015_2018_masterfile_plots.pdf", useDingbats = F, onefile = T, width = 20, height=12)
metadata_input %>% plot_histogram(nrow = 4L, ncol=3L)
metadata_input %>% plot_bar(nrow = 4L, ncol=2L)
dev.off()

#spearman correlations
masterfile_corr <- masterfile[,c("child_pedsql_phy_score",
                                 "child_pedsql_emo_score",
                                 "child_pedsql_soc_score",
                                 "child_pedsql_sch_score",
                                 "child_pedsql_psysoc_score",
                                 "child_pedsql_health_related_qol_score_total",
                                 "child_pedsql_phy_score_catg_Q15",
                                 "child_pedsql_emo_score_catg_Q15",
                                 "child_pedsql_soc_score_catg_Q15",
                                 "child_pedsql_sch_score_catg_Q15",
                                 "child_pedsql_psysoc_score_catg_Q15",
                                 "child_pedsql_health_related_qol_score_total_catg_Q15",
                                 "child_pedsql_phy_score_catg_Q25",
                                 "child_pedsql_emo_score_catg_Q25",
                                 "child_pedsql_soc_score_catg_Q25",
                                 "child_pedsql_sch_score_catg_Q25",
                                 "child_pedsql_psysoc_score_catg_Q25",
                                 "child_pedsql_health_related_qol_score_total_catg_Q25",
                                 "child_ffq_pattern_1_meat_seafood_prepd_meals",
                                 "child_ffq_pattern_2_cooked_veg_grains_legumes",
                                 "child_ffq_pattern_3_fruits_raw_veg_cheese",
                                 "child_ffq_pattern_4_confectioneries_pizza",
                                 "child_ffq_pattern_5_starchy_foods_sweetened_bev",
                                 "child_ffq_hpdi_score",
                                 "child_ffq_animal_score",
                                 "child_ffq_diet_quality_score",
                                 "child_kidmed_score",
                                 "child_kidmed_adherence_level",
                                 "child_demo_school_year",
                                 "child_demo_school_code",
                                 "child_demo_school_region",
                                 "child_demo_grade",
                                 "child_demo_age_months",
                                 "child_demo_age_years",
                                 "child_demo_sex",
                                 "child_anthro_bmi_WHO_sds",
                                 "child_anthro_bmi_WHO_sds_grouped",
                                 "child_anthro_bmi_grouped_7",
                                 "child_anthro_bmi_grouped_4",
                                 "child_demo_country_birth",
                                 "child_demo_immigrant_status",
                                 "child_demo_lives_with",
                                 "child_demo_insurance",
                                 "child_phyact_sport_outside_school",
                                 "child_phyact_exercise_hrs_per_wk",
                                 "child_phyact_screen_time_hrs_per_wk",
                                 "child_phyact_play_outside_days_per_wk",
                                 "parent_demo_employment",
                                 "parent_demo_country_birth",
                                 "mother_demo_age_years",
                                 "mother_demo_edu",
                                 "father_demo_age_years",
                                 "father_demo_edu",
                                 "family_demo_FAS_score",
                                 "family_demo_FAS_group",
                                 "family_demo_FAS_group_dic",
                                 "family_demo_children_in_family_n",
                                 "family_demo_household_members_n")]

for (i in 1:ncol(masterfile_corr)){masterfile_corr[,i] <- as.numeric(masterfile_corr[,i])}
str(masterfile_corr)
#6583 58

spearman_corr <- result(masterfile_corr, masterfile_corr)
#3364  4

#remove 1:1 correlations
spearman_corr_unique <- spearman_corr %>% 
        as_tibble() %>% 
        mutate(duplicates = if_else(factor1 == factor2,
                                    TRUE,
                                    FALSE)) %>% 
        filter(duplicates == FALSE)

spearman_corr_unique$duplicates <- NULL 

#remove NA's (for incomplete cases for correlations)
spearman_corr_unique <- spearman_corr_unique[complete.cases(spearman_corr_unique), ]

spearman_FDR <- spearman_corr_unique
spearman_FDR$FDR<-p.adjust(spearman_FDR$pvalue, method = "BH")
#3306  5
spearman_FDR_corr <- spearman_FDR %>% filter(
        FDR<0.05
)
#2254  5

write.table(spearman_FDR, "2023_03_20_diatrofi_2015_2018_masterfile_corrlns.txt", sep="\t", row.names=F, quote = F)

pdf(file = "2023_03_20_diatrofi_2015_2018_masterfile_corrlns.pdf", useDingbats = F, onefile = T, width = 15, height=18)
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

ggplot(data = spearman_FDR_corr, aes(factor1, factor2, fill = CorCoefficient))+
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

        #### =============== 4. TABLE 2: QOL BY DIETARY INDEX AND PATTERNS WITH P-VALUES   ================ ####

masterfile <- read.delim("2023_03_20_diatrofi_2015_2018_masterfile.txt", sep = "\t", header = T, stringsAsFactors = T)
str(masterfile)

masterfile_sumstats <- masterfile
str(masterfile_sumstats)

#general summary stats for whole cohort/study group
mean(masterfile_sumstats$child_ffq_hpdi_score, na.rm=T)
sd(masterfile_sumstats$child_ffq_hpdi_score, na.rm=T)

mean(masterfile_sumstats$child_ffq_animal_score, na.rm=T)
sd(masterfile_sumstats$child_ffq_animal_score, na.rm=T)

mean(masterfile_sumstats$child_ffq_diet_quality_score, na.rm=T)
sd(masterfile_sumstats$child_ffq_diet_quality_score, na.rm=T)

median(masterfile_sumstats$child_ffq_pattern_1_meat_seafood_prepd_meals, na.rm=T)
IQR(masterfile_sumstats$child_ffq_pattern_1_meat_seafood_prepd_meals, na.rm=T)

median(masterfile_sumstats$child_ffq_pattern_2_cooked_veg_grains_legumes, na.rm=T)
IQR(masterfile_sumstats$child_ffq_pattern_2_cooked_veg_grains_legumes, na.rm=T)

median(masterfile_sumstats$child_ffq_pattern_3_fruits_raw_veg_cheese, na.rm=T)
IQR(masterfile_sumstats$child_ffq_pattern_3_fruits_raw_veg_cheese, na.rm=T)

median(masterfile_sumstats$child_ffq_pattern_4_confectioneries_pizza, na.rm=T)
IQR(masterfile_sumstats$child_ffq_pattern_4_confectioneries_pizza, na.rm=T)

median(masterfile_sumstats$child_ffq_pattern_5_starchy_foods_sweetened_bev, na.rm=T)
IQR(masterfile_sumstats$child_ffq_pattern_5_starchy_foods_sweetened_bev, na.rm=T)

        ####4.1. QoL (Q15) by dietary quality index ####

#HRQoL: Q15 summary stats
hrqol_numeric <- masterfile_sumstats %>%
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
x <- t.test(child_ffq_hpdi_score ~ child_pedsql_health_related_qol_score_total_catg_Q15, data=masterfile_sumstats)
x$p.value

x <- t.test(child_ffq_animal_score ~ child_pedsql_health_related_qol_score_total_catg_Q15, data=masterfile_sumstats)
x$p.value

x <- t.test(child_ffq_diet_quality_score ~ child_pedsql_health_related_qol_score_total_catg_Q15, data=masterfile_sumstats)
x$p.value

wilcox.test(child_ffq_pattern_1_meat_seafood_prepd_meals ~ child_pedsql_health_related_qol_score_total_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_2_cooked_veg_grains_legumes ~ child_pedsql_health_related_qol_score_total_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_3_fruits_raw_veg_cheese ~ child_pedsql_health_related_qol_score_total_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_4_confectioneries_pizza ~ child_pedsql_health_related_qol_score_total_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_5_starchy_foods_sweetened_bev ~ child_pedsql_health_related_qol_score_total_catg_Q15, data=masterfile_sumstats)


#psysoc: Q15 summary stats
psysoc_numeric <- masterfile_sumstats %>%
        filter(is.na(child_pedsql_psysoc_score_catg_Q15) == F) %>%
        group_by(child_pedsql_psysoc_score_catg_Q15) %>%
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

#psysoc: Q15 p-values
t.test(child_ffq_hpdi_score ~ child_pedsql_psysoc_score_catg_Q15, data=masterfile_sumstats)
t.test(child_ffq_animal_score ~ child_pedsql_psysoc_score_catg_Q15, data=masterfile_sumstats)
t.test(child_ffq_diet_quality_score ~ child_pedsql_psysoc_score_catg_Q15, data=masterfile_sumstats)

wilcox.test(child_ffq_pattern_1_meat_seafood_prepd_meals ~ child_pedsql_psysoc_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_2_cooked_veg_grains_legumes ~ child_pedsql_psysoc_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_3_fruits_raw_veg_cheese ~ child_pedsql_psysoc_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_4_confectioneries_pizza ~ child_pedsql_psysoc_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_5_starchy_foods_sweetened_bev ~ child_pedsql_psysoc_score_catg_Q15, data=masterfile_sumstats)

#phy: Q15 summary stats
phy_numeric <- masterfile_sumstats %>%
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
t.test(child_ffq_hpdi_score ~ child_pedsql_phy_score_catg_Q15, data=masterfile_sumstats)
t.test(child_ffq_animal_score ~ child_pedsql_phy_score_catg_Q15, data=masterfile_sumstats)
t.test(child_ffq_diet_quality_score ~ child_pedsql_phy_score_catg_Q15, data=masterfile_sumstats)

wilcox.test(child_ffq_pattern_1_meat_seafood_prepd_meals ~ child_pedsql_phy_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_2_cooked_veg_grains_legumes ~ child_pedsql_phy_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_3_fruits_raw_veg_cheese ~ child_pedsql_phy_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_4_confectioneries_pizza ~ child_pedsql_phy_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_5_starchy_foods_sweetened_bev ~ child_pedsql_phy_score_catg_Q15, data=masterfile_sumstats)

#emo: Q15 summary stats
emo_numeric <- masterfile_sumstats %>%
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
t.test(child_ffq_hpdi_score ~ child_pedsql_emo_score_catg_Q15, data=masterfile_sumstats)
t.test(child_ffq_animal_score ~ child_pedsql_emo_score_catg_Q15, data=masterfile_sumstats)
t.test(child_ffq_diet_quality_score ~ child_pedsql_emo_score_catg_Q15, data=masterfile_sumstats)

wilcox.test(child_ffq_pattern_1_meat_seafood_prepd_meals ~ child_pedsql_emo_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_2_cooked_veg_grains_legumes ~ child_pedsql_emo_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_3_fruits_raw_veg_cheese ~ child_pedsql_emo_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_4_confectioneries_pizza ~ child_pedsql_emo_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_5_starchy_foods_sweetened_bev ~ child_pedsql_emo_score_catg_Q15, data=masterfile_sumstats)

#soc: Q15 summary stats
soc_numeric <- masterfile_sumstats %>%
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
t.test(child_ffq_hpdi_score ~ child_pedsql_soc_score_catg_Q15, data=masterfile_sumstats)
t.test(child_ffq_animal_score ~ child_pedsql_soc_score_catg_Q15, data=masterfile_sumstats)
t.test(child_ffq_diet_quality_score ~ child_pedsql_soc_score_catg_Q15, data=masterfile_sumstats)

wilcox.test(child_ffq_pattern_1_meat_seafood_prepd_meals ~ child_pedsql_soc_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_2_cooked_veg_grains_legumes ~ child_pedsql_soc_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_3_fruits_raw_veg_cheese ~ child_pedsql_soc_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_4_confectioneries_pizza ~ child_pedsql_soc_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_5_starchy_foods_sweetened_bev ~ child_pedsql_soc_score_catg_Q15, data=masterfile_sumstats)

#sch: Q15 summary stats
sch_numeric <- masterfile_sumstats %>%
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
t.test(child_ffq_hpdi_score ~ child_pedsql_sch_score_catg_Q15, data=masterfile_sumstats)
t.test(child_ffq_animal_score ~ child_pedsql_sch_score_catg_Q15, data=masterfile_sumstats)
t.test(child_ffq_diet_quality_score ~ child_pedsql_sch_score_catg_Q15, data=masterfile_sumstats)

wilcox.test(child_ffq_pattern_1_meat_seafood_prepd_meals ~ child_pedsql_sch_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_2_cooked_veg_grains_legumes ~ child_pedsql_sch_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_3_fruits_raw_veg_cheese ~ child_pedsql_sch_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_4_confectioneries_pizza ~ child_pedsql_sch_score_catg_Q15, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_5_starchy_foods_sweetened_bev ~ child_pedsql_sch_score_catg_Q15, data=masterfile_sumstats)

        ####4.2. QoL (Q25) by dietary quality index ####
#HRQoL: Q25 summary stats
hrqol_numeric <- masterfile_sumstats %>%
        filter(is.na(child_pedsql_health_related_qol_score_total_catg_Q25) == F) %>%
        group_by(child_pedsql_health_related_qol_score_total_catg_Q25) %>%
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

#HRQoL: Q25 p-values
x <- t.test(child_ffq_hpdi_score ~ child_pedsql_health_related_qol_score_total_catg_Q25, data=masterfile_sumstats)
x$p.value

x <- t.test(child_ffq_animal_score ~ child_pedsql_health_related_qol_score_total_catg_Q25, data=masterfile_sumstats)
x$p.value

x <- t.test(child_ffq_diet_quality_score ~ child_pedsql_health_related_qol_score_total_catg_Q25, data=masterfile_sumstats)
x$p.value

wilcox.test(child_ffq_pattern_1_meat_seafood_prepd_meals ~ child_pedsql_health_related_qol_score_total_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_2_cooked_veg_grains_legumes ~ child_pedsql_health_related_qol_score_total_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_3_fruits_raw_veg_cheese ~ child_pedsql_health_related_qol_score_total_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_4_confectioneries_pizza ~ child_pedsql_health_related_qol_score_total_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_5_starchy_foods_sweetened_bev ~ child_pedsql_health_related_qol_score_total_catg_Q25, data=masterfile_sumstats)


#psysoc: Q25 summary stats
psysoc_numeric <- masterfile_sumstats %>%
        filter(is.na(child_pedsql_psysoc_score_catg_Q25) == F) %>%
        group_by(child_pedsql_psysoc_score_catg_Q25) %>%
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

#psysoc: Q25 p-values
t.test(child_ffq_hpdi_score ~ child_pedsql_psysoc_score_catg_Q25, data=masterfile_sumstats)
t.test(child_ffq_animal_score ~ child_pedsql_psysoc_score_catg_Q25, data=masterfile_sumstats)
t.test(child_ffq_diet_quality_score ~ child_pedsql_psysoc_score_catg_Q25, data=masterfile_sumstats)

wilcox.test(child_ffq_pattern_1_meat_seafood_prepd_meals ~ child_pedsql_psysoc_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_2_cooked_veg_grains_legumes ~ child_pedsql_psysoc_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_3_fruits_raw_veg_cheese ~ child_pedsql_psysoc_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_4_confectioneries_pizza ~ child_pedsql_psysoc_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_5_starchy_foods_sweetened_bev ~ child_pedsql_psysoc_score_catg_Q25, data=masterfile_sumstats)

#phy: Q25 summary stats
phy_numeric <- masterfile_sumstats %>%
        filter(is.na(child_pedsql_phy_score_catg_Q25) == F) %>%
        group_by(child_pedsql_phy_score_catg_Q25) %>%
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

#phy: Q25 p-values
t.test(child_ffq_hpdi_score ~ child_pedsql_phy_score_catg_Q25, data=masterfile_sumstats)
t.test(child_ffq_animal_score ~ child_pedsql_phy_score_catg_Q25, data=masterfile_sumstats)
t.test(child_ffq_diet_quality_score ~ child_pedsql_phy_score_catg_Q25, data=masterfile_sumstats)

wilcox.test(child_ffq_pattern_1_meat_seafood_prepd_meals ~ child_pedsql_phy_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_2_cooked_veg_grains_legumes ~ child_pedsql_phy_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_3_fruits_raw_veg_cheese ~ child_pedsql_phy_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_4_confectioneries_pizza ~ child_pedsql_phy_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_5_starchy_foods_sweetened_bev ~ child_pedsql_phy_score_catg_Q25, data=masterfile_sumstats)

#emo: Q25 summary stats
emo_numeric <- masterfile_sumstats %>%
        filter(is.na(child_pedsql_emo_score_catg_Q25) == F) %>%
        group_by(child_pedsql_emo_score_catg_Q25) %>%
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

#emo: Q25 p-values
t.test(child_ffq_hpdi_score ~ child_pedsql_emo_score_catg_Q25, data=masterfile_sumstats)
t.test(child_ffq_animal_score ~ child_pedsql_emo_score_catg_Q25, data=masterfile_sumstats)
t.test(child_ffq_diet_quality_score ~ child_pedsql_emo_score_catg_Q25, data=masterfile_sumstats)

wilcox.test(child_ffq_pattern_1_meat_seafood_prepd_meals ~ child_pedsql_emo_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_2_cooked_veg_grains_legumes ~ child_pedsql_emo_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_3_fruits_raw_veg_cheese ~ child_pedsql_emo_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_4_confectioneries_pizza ~ child_pedsql_emo_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_5_starchy_foods_sweetened_bev ~ child_pedsql_emo_score_catg_Q25, data=masterfile_sumstats)

#soc: Q25 summary stats
soc_numeric <- masterfile_sumstats %>%
        filter(is.na(child_pedsql_soc_score_catg_Q25) == F) %>%
        group_by(child_pedsql_soc_score_catg_Q25) %>%
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

#soc: Q25 p-values
t.test(child_ffq_hpdi_score ~ child_pedsql_soc_score_catg_Q25, data=masterfile_sumstats)
t.test(child_ffq_animal_score ~ child_pedsql_soc_score_catg_Q25, data=masterfile_sumstats)
t.test(child_ffq_diet_quality_score ~ child_pedsql_soc_score_catg_Q25, data=masterfile_sumstats)

wilcox.test(child_ffq_pattern_1_meat_seafood_prepd_meals ~ child_pedsql_soc_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_2_cooked_veg_grains_legumes ~ child_pedsql_soc_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_3_fruits_raw_veg_cheese ~ child_pedsql_soc_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_4_confectioneries_pizza ~ child_pedsql_soc_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_5_starchy_foods_sweetened_bev ~ child_pedsql_soc_score_catg_Q25, data=masterfile_sumstats)

#sch: Q25 summary stats
sch_numeric <- masterfile_sumstats %>%
        filter(is.na(child_pedsql_sch_score_catg_Q25) == F) %>%
        group_by(child_pedsql_sch_score_catg_Q25) %>%
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

#sch: Q25 p-values
t.test(child_ffq_hpdi_score ~ child_pedsql_sch_score_catg_Q25, data=masterfile_sumstats)
t.test(child_ffq_animal_score ~ child_pedsql_sch_score_catg_Q25, data=masterfile_sumstats)
t.test(child_ffq_diet_quality_score ~ child_pedsql_sch_score_catg_Q25, data=masterfile_sumstats)

wilcox.test(child_ffq_pattern_1_meat_seafood_prepd_meals ~ child_pedsql_sch_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_2_cooked_veg_grains_legumes ~ child_pedsql_sch_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_3_fruits_raw_veg_cheese ~ child_pedsql_sch_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_4_confectioneries_pizza ~ child_pedsql_sch_score_catg_Q25, data=masterfile_sumstats)
wilcox.test(child_ffq_pattern_5_starchy_foods_sweetened_bev ~ child_pedsql_sch_score_catg_Q25, data=masterfile_sumstats)

        #### =============== 5. PLOTS: QOL BY DIETARY INDEX AND PATTERNS WITH P-VALUES   ================ ####

###NOTE: plots dietary quality index and Q25, as well as all other QoL outcomes (for both Q15 and Q25) vs dietary patterns still need to be coded
###      use current basis below to adjust and add code accordingly (when needed)

        ####5.1. plots: QoL (Q15) by dietary quality index ####

#total score                
diet_score_plot <- masterfile_sumstats[, c("child_ffq_diet_quality_score",
                                       "child_pedsql_health_related_qol_score_total_catg_Q15")]
diet_score_plot <- diet_score_plot %>% rename(qol_total_Q15=child_pedsql_health_related_qol_score_total_catg_Q15,
                                              diet_score=child_ffq_diet_quality_score)
diet_score_plot <- diet_score_plot[complete.cases(diet_score_plot), ]
summary(diet_score_plot)
#6019  2

p1 <- ggline(diet_score_plot,
               x = "qol_total_Q15", y = "diet_score", 
               color = "qol_total_Q15", 
               order = c("0_poor_avg", "1_good"),
               add = c("mean_sd", "violin"),
               main = "Diet quality score by Total QoL (Q15)",
               ylab = "mean diet quality score (SE)", xlab = "",
               legend = "none")
p1 <- ggpar(p1, 
            palette = "Dark2",
            legend = "none")

p1 <- p1 + scale_x_discrete(labels=c('0_poor_avg' = 'poor/avg \n N=914', 
                                              '1_good' = 'good\n N=5105'))

p1 <- p1 + stat_compare_means(method = "t.test", label.y=95) 

#psysoc
diet_score_plot <- masterfile_sumstats[, c("child_ffq_diet_quality_score",
                                        "child_pedsql_psysoc_score_catg_Q15")]
diet_score_plot <- diet_score_plot %>% rename(qol_psysoc_Q15=child_pedsql_psysoc_score_catg_Q15,
                                              diet_score=child_ffq_diet_quality_score)
diet_score_plot <- diet_score_plot[complete.cases(diet_score_plot), ]
summary(diet_score_plot)
#5999  2

p2 <- ggline(diet_score_plot,
               x = "qol_psysoc_Q15", y = "diet_score", 
               color = "qol_psysoc_Q15", 
               order = c("0_poor_avg", "1_good"),
             add = c("mean_sd", "violin"),
               main = "Diet quality score by Psychosocial QoL (Q15)",
               ylab = "mean diet quality score (SE)", xlab = "",
               legend = "none")
p2 <- ggpar(p2, 
            palette = "Dark2",
            legend = "none")

p2 <- p2 + scale_x_discrete(labels=c('0_poor_avg' = 'poor/avg \n N=978', 
                                     '1_good' = 'good\n N=5021'))

p2 <- p2 + stat_compare_means(method = "t.test", label.y=95) 

#phys
diet_score_plot <- masterfile_sumstats[, c("child_ffq_diet_quality_score",
                                           "child_pedsql_phy_score_catg_Q15")]
diet_score_plot <- diet_score_plot %>% rename(qol_phy_Q15=child_pedsql_phy_score_catg_Q15,
                                              diet_score=child_ffq_diet_quality_score)

diet_score_plot <- diet_score_plot[complete.cases(diet_score_plot), ]
summary(diet_score_plot)
#5996  2

p3 <- ggline(diet_score_plot,
               x = "qol_phy_Q15", y = "diet_score", 
               color = "qol_phy_Q15", 
               order = c("0_poor_avg", "1_good"),
             add = c("mean_sd", "violin"),
               main = "Diet quality score by Physical QoL (Q15)",
               ylab = "mean diet quality score (SE)", xlab = "",
               legend = "none")
p3 <- ggpar(p3, 
            palette = "Dark2",
            legend = "none")

p3 <- p3 + scale_x_discrete(labels=c('0_poor_avg' = 'poor/avg \n N=906', 
                                     '1_good' = 'good\n N=5090'))

p3 <- p3 + stat_compare_means(method = "t.test", label.y=95) 

#emo
diet_score_plot <- masterfile_sumstats[, c("child_ffq_diet_quality_score",
                                           "child_pedsql_emo_score_catg_Q15")]
diet_score_plot <- diet_score_plot %>% rename(qol_emo_Q15=child_pedsql_emo_score_catg_Q15,
                                              diet_score=child_ffq_diet_quality_score)

diet_score_plot <- diet_score_plot[complete.cases(diet_score_plot), ]
summary(diet_score_plot)
#6008  2

p4 <- ggline(diet_score_plot,
               x = "qol_emo_Q15", y = "diet_score", 
               color = "qol_emo_Q15", 
               order = c("0_poor_avg", "1_good"),
             add = c("mean_sd", "violin"),
               main = "Diet quality score by Emotional QoL (Q15)",
               ylab = "mean diet quality score (SE)", xlab = "",
               legend = "none")
p4 <- ggpar(p4, 
            palette = "Dark2",
            legend = "none")

p4 <- p4 + scale_x_discrete(labels=c('0_poor_avg' = 'poor/avg \n N=906', 
                                     '1_good' = 'good\n N=5102'))

p4 <- p4 + stat_compare_means(method = "t.test", label.y=95)

#soc
diet_score_plot <- masterfile_sumstats[, c("child_ffq_diet_quality_score",
                                           "child_pedsql_soc_score_catg_Q15")]
diet_score_plot <- diet_score_plot %>% rename(qol_soc_Q15=child_pedsql_soc_score_catg_Q15,
                                              diet_score=child_ffq_diet_quality_score)

diet_score_plot <- diet_score_plot[complete.cases(diet_score_plot), ]
summary(diet_score_plot)
#5996  2

p5 <- ggline(diet_score_plot,
               x = "qol_soc_Q15", y = "diet_score", 
               color = "qol_soc_Q15", 
               order = c("0_poor_avg", "1_good"),
             add = c("mean_sd", "violin"),
               main = "Diet quality score by Social QoL (Q15)",
               ylab = "mean diet quality score (SE)", xlab = "",
               legend = "none")
p5 <- ggpar(p5, 
            palette = "Dark2",
            legend = "none")

p5 <- p5 + scale_x_discrete(labels=c('0_poor_avg' = 'poor/avg \n N=902', 
                                     '1_good' = 'good\n N=5094'))

p5 <- p5 + stat_compare_means(method = "t.test", label.y=95)

#sch
diet_score_plot <- masterfile_sumstats[, c("child_ffq_diet_quality_score",
                                           "child_pedsql_sch_score_catg_Q15")]
diet_score_plot <- diet_score_plot %>% rename(qol_sch_Q15=child_pedsql_sch_score_catg_Q15,
                                              diet_score=child_ffq_diet_quality_score)

diet_score_plot <- diet_score_plot[complete.cases(diet_score_plot), ]
summary(diet_score_plot)
#5938  2

p6 <- ggline(diet_score_plot,
               x = "qol_sch_Q15", y = "diet_score", 
               color = "qol_sch_Q15", 
               order = c("0_poor_avg", "1_good"),
             add = c("mean_sd", "violin"),
               main = "Diet quality score by School Functioning QoL (Q15)",
               ylab = "mean diet quality score (SE)", xlab = "",
               legend = "none")
p6 <- ggpar(p6, 
            palette = "Dark2",
            legend = "none")

p6 <- p6 + scale_x_discrete(labels=c('0_poor_avg' = 'poor/avg \n N=906', 
                                     '1_good' = 'good\n N=5102'))

p6 <- p6 + stat_compare_means(method = "t.test", label.y=95)


pdf(file = "2023_03_16_diet_score_qol_Q15_all.pdf", useDingbats = F, onefile = T, width = 20, height=12)
p <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3) + 
        plot_annotation(title = "Diet Quality Index Score by QoL outcomes (Q15)", tag_levels = "a") &
        theme(plot.title = element_text(hjust = 0.5, size = 12))

print(p)
dev.off()

        ####5.2. plots: QoL (Q15) by dietary patterns ####
                ####6.2.1 total QoL (Q15) by dietary patterns####
hrqol_patterns <- masterfile_sumstats[,c("child_pedsql_health_related_qol_score_total_catg_Q15",
                                         "child_ffq_pattern_1_meat_seafood_prepd_meals",
                                         "child_ffq_pattern_2_cooked_veg_grains_legumes",
                                         "child_ffq_pattern_3_fruits_raw_veg_cheese",
                                         "child_ffq_pattern_4_confectioneries_pizza",
                                         "child_ffq_pattern_5_starchy_foods_sweetened_bev")]
str(hrqol_patterns)
hrqol_patterns <- hrqol_patterns[complete.cases(hrqol_patterns), ]
#6019 6
summary(hrqol_patterns)

hrqol_patterns <- hrqol_patterns %>% rename(qol_total_Q15=child_pedsql_health_related_qol_score_total_catg_Q15,
                                            pattern_1=child_ffq_pattern_1_meat_seafood_prepd_meals,
                                            pattern_2=child_ffq_pattern_2_cooked_veg_grains_legumes,
                                            pattern_3=child_ffq_pattern_3_fruits_raw_veg_cheese,
                                            pattern_4=child_ffq_pattern_4_confectioneries_pizza,
                                            pattern_5=child_ffq_pattern_5_starchy_foods_sweetened_bev)

#pattern 1
pat1 <- ggline(hrqol_patterns,
             x = "qol_total_Q15", y = "pattern_1", 
             color = "qol_total_Q15", 
             order = c("0_poor_avg", "1_good"),
             add = c("median_iqr", "violin"),
             main = "Dietary Pattern 1 by Total QoL (Q15)",
             ylab = "median dietay pattern 1 (IQR)", xlab = "",
             legend = "none")
pat1 <- ggpar(pat1, 
            palette = "Dark2",
            legend = "none")

pat1 <- pat1 + scale_x_discrete(labels=c('0_poor_avg' = 'poor/avg \n N=914', 
                                     '1_good' = 'good\n N=5105'))

pat1 <- pat1 + stat_compare_means(label.y=28) 

#pattern 2
pat2 <- ggline(hrqol_patterns,
               x = "qol_total_Q15", y = "pattern_2", 
               color = "qol_total_Q15", 
               order = c("0_poor_avg", "1_good"),
               add = c("median_iqr", "violin"),
               main = "Dietary Pattern 2 by Total QoL (Q15)",
               ylab = "median dietay pattern 2 (IQR)", xlab = "",
               legend = "none")
pat2 <- ggpar(pat2, 
              palette = "Dark2",
              legend = "none")

pat2 <- pat2 + scale_x_discrete(labels=c('0_poor_avg' = 'poor/avg \n N=914', 
                                         '1_good' = 'good\n N=5105'))

pat2 <- pat2 + stat_compare_means(label.y=28) 

#pattern 3
pat3 <- ggline(hrqol_patterns,
               x = "qol_total_Q15", y = "pattern_3", 
               color = "qol_total_Q15", 
               order = c("0_poor_avg", "1_good"),
               add = c("median_iqr", "violin"),
               main = "Dietary Pattern 3 by Total QoL (Q15)",
               ylab = "median dietay pattern 3 (IQR)", xlab = "",
               legend = "none")
pat3 <- ggpar(pat3, 
              palette = "Dark2",
              legend = "none")

pat3 <- pat3 + scale_x_discrete(labels=c('0_poor_avg' = 'poor/avg \n N=914', 
                                         '1_good' = 'good\n N=5105'))

pat3 <- pat3 + stat_compare_means(label.y=28) 

#pattern 4
pat4 <- ggline(hrqol_patterns,
               x = "qol_total_Q15", y = "pattern_4", 
               color = "qol_total_Q15", 
               order = c("0_poor_avg", "1_good"),
               add = c("median_iqr", "violin"),
               main = "Dietary Pattern 4 by Total QoL (Q15)",
               ylab = "median dietay pattern 4 (IQR)", xlab = "",
               legend = "none")
pat4 <- ggpar(pat4, 
              palette = "Dark2",
              legend = "none")

pat4 <- pat4 + scale_x_discrete(labels=c('0_poor_avg' = 'poor/avg \n N=914', 
                                         '1_good' = 'good\n N=5105'))

pat4 <- pat4 + stat_compare_means(label.y=28) 

#pattern 5
pat5 <- ggline(hrqol_patterns,
               x = "qol_total_Q15", y = "pattern_5", 
               color = "qol_total_Q15", 
               order = c("0_poor_avg", "1_good"),
               add = c("median_iqr", "violin"),
               main = "Dietary Pattern 5 by Total QoL (Q15)",
               ylab = "median dietay pattern 5 (IQR)", xlab = "",
               legend = "none")
pat5 <- ggpar(pat5, 
              palette = "Dark2",
              legend = "none")

pat5 <- pat5 + scale_x_discrete(labels=c('0_poor_avg' = 'poor/avg \n N=914', 
                                         '1_good' = 'good\n N=5105'))

pat5 <- pat5 + stat_compare_means(label.y=28) 

pdf(file = "2023_03_16_diet_patterns_qol_total_Q15_all.pdf", useDingbats = F, onefile = T, width = 20, height=12)
p <- pat1 + pat2 + pat3 + pat4 + pat5 + plot_layout(ncol = 3) + 
        plot_annotation(title = "Dietary Patterns HRQoL (Q15)", tag_levels = "a") &
        theme(plot.title = element_text(hjust = 0.5, size = 12))

print(p)
dev.off()
