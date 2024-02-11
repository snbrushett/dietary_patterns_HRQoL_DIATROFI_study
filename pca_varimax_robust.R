###     SCRIPT: DIATROFI FFQ DIETARY PATTERNS VIA (ROBUST) PCA VARIMAX ROTATION
###     AUTHOR(S): SIOBHAN BRUSHETT 
###     DESCRIPTION: VARIMAC PCA AND ROBUST PCA ON FFQ ITEMS FROM THE BASELINE 2015-2016 AND 2017-2018 DIATROFI DATA
###     PROJECT: MERGED DIATROFI 2015-2018
###     NOTE(S): PCA: https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome/blob/master/Projects/IBD_onset_Diet/05_PCA.R
###              PAPERS: https://doi.org/10.1093/ecco-jcc/jjab219 (Genetics Department, Groningen); https://doi.org/10.1093/eurpub/ckaa178 (Prolepsis)
###              TESTS: https://www.statology.org/bartletts-test-of-sphericity/; https://www.statisticshowto.com/kaiser-meyer-olkin/

setwd("C:/Users/Siobhan Brushett/OneDrive - UMCG/Prolepsis_internship/merged_2015_2018/")

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
 
#functions
#cluster analysis (Laura and Arnau) #updated_SB: removed scaling as all variables have the same value range
do.clustering <- function(diet_data,
                          dist.method="euclidean",        #ordinary straight line distance between 2 points in the euclidean space. With this distance, euclidean space becomes a metric space. root of (x1-x2)^2 + (y1-y2)^2
                          cluster.method = "complete",    
                          h= 0.5){                        #Clustering height at which the tree is cut into groups. The higher this number, the more stringent and the less groups.
        if (dist.method == "euclidean"){         
                data.dist = dist(t(diet_data))                #t: dist() calculating the distance matrix expects vectors to be horizontal. Transform data so that cols=id, rows=foods 
        } else {
                data.dist = as.dist(1-cor(diet_data)^2)               #else: correlation based distance, considers 2 objects similar if their features are highly correlated, even though the observed values mau be far apart in euclidean distance. The lower the distance e.g. 0 the more correlated. 
        }
        hclust.object=hclust(data.dist, method = cluster.method) #hierarchical clustering of similar foods into groups
        plot(hang.dendrogram(as.dendrogram(hclust.object), cex=0.5, hang=-1)) #horiz=T, type="triangle", hang=-1, error: "hang" is not a graphical parameter
        
        cutree_returned = cutree(hclust.object,h=h)                      #Cut the tree resulting from hclust into several group by specifying h 
        cutree_returned
}

##Contents of this file

## 0. LOAD DATA
## 1. CLEAN FFQ ITEMS BY EXCLUSION CRITERIA
## 2. DERIVE FFQ WEIGHTED FREQUENCIES
## 3. TESTS TO VALIDATE PC FACTOR ANALYSIS
## 4. PCA USING VARIMAX ROTATION
## 5. PCA USING OBLIQUE TRANSFORMATION
## 6. ROBUST PCA AND COMPARISON WITH STANDARD PCA
## 7. PCA USING HIERARCHICAL CLUSTERING 
## 8. CLUSTER ANALYSIS (EUCLIDEAN AND BRAY CURTIS DISTANCES)


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

        #### =============== 5. PCA USING OBLIQUE TRANSFORMATION   ================ ####
#pca_oblique <- psych::principal(dia_merged_ffq, rotate="oblimin", nfactors=5, scores=TRUE) #first 5 components
#summary(pca_oblique)

#pdf(file = "2023_02_16_diatrofi_ffq_pca_oblique_correlation.pdf", useDingbats = F, onefile = T, width = 20, height=12)
#corrplot(pca_oblique$loadings,is.corr=F)
##fa.diagram(pca_oblique)
#dev.off()

#pc_obl_loadings <- pca_oblique$loadings

#ffq_pc_obl <- as.matrix(dia_merged_ffq)
#ffq_pc_obl <- ffq_pc_obl %*% pc_obl_loadings

#write.table(ffq_pc_obl, "2023_02_20_ffq_pca_oblique_loadings.txt", sep="\t", row.names=T, quote = F)

#defining key items of loadings
#oblique_item_loadings <- as.data.frame(loadings(pca_oblique)[])
#summary(oblique_item_loadings)

#oblique_item_loadings$items <- rownames(oblique_item_loadings)
#TC1_loadings <- oblique_item_loadings %>%
#        filter(TC1 > 0.45) %>%
#        select (c(items, TC1)) %>%
#        na.omit(TC1_loadings) %>%
#        arrange(desc(TC1))
#19 items

#TC2_loadings <- oblique_item_loadings %>% 
#        filter(TC2 > 0.45) %>%
#        select (c(items, TC2)) %>%
#        na.omit(TC2_loadings) %>%
#        arrange(desc(TC2))
#10 items

#TC3_loadings <- oblique_item_loadings %>% 
#        filter(TC3 > 0.45) %>%
#        select (c(items, TC3)) %>%
#        na.omit(TC3_loadings) %>%
#        arrange(desc(TC3))
#2 items

#TC4_loadings <- oblique_item_loadings %>% 
#        filter(TC4 > 0.45) %>%
#        select (c(items, TC4)) %>%
#        na.omit(TC4_loadings) %>%
#        arrange(desc(TC4))
#5 items

#TC5_loadings <- oblique_item_loadings %>% 
#        filter(TC5 > 0.5) %>%
#        select (c(items, TC5)) %>%
#        na.omit(TC5_loadings) %>%
#        arrange(desc(TC5))
#3 items

#clean RStudio
#rm(pca_oblique, ffq_pc_obl, oblique_item_loadings, TC1_loadings, TC2_loadings, TC3_loadings, TC4_loadings, TC5_loadings)

        #### =============== 6. ROBUST PCA AND COMPARISON WITH STANDARD PCA   ================ ####

#ffq_robust_pca <- PcaHubert(dia_merged_ffq, k=5, alpha =0.75, scale=FALSE, center = FALSE)
#sdev <- sqrt(ffq_robust_pca$eigenvalues) #compute SD of PCs from eigenvalues
#raw_loadings_robust_pca <- ffq_robust_pca$loadings[,1:5] %*% diag(sdev, 5, 5)
#robust_pca_loadings_varimax <- varimax(raw_loadings_robust_pca, normalize = FALSE)$loadings
#inv_loadings <- t(pracma::pinv(robust_pca_loadings_varimax))
#res_pca_rot_scores <- scale(dia_merged_ffq) %*% inv_loadings
#res_pca_rot_scores <- as_tibble(res_pca_rot_scores) %>% mutate(pseudo_ids = rownames(dia_merged_ffq)) 

#write.table(res_pca_rot_scores,file='2023_02_15_ffq_robust_pca_varimax_loadings.txt', quote=FALSE, sep='\t',row.names = F)

#print(inv_loadings)

        #### =============== 7. PCA USING HIERARCHICAL CLUSTERING   ================ ####

#pdf(file = "2023_02_16_diatrofi_ffq_pca_for_hcpc.pdf", useDingbats = F, onefile = T, width = 20, height=12)
#visualizing pcas
#res.pca <- PCA(dia_merged_ffq)
#fviz_pca_biplot(res.pca)
#dev.off()

#hcpc
#res.hcpc <- HCPC(res.pca, graph=F)

#pdf(file = "2023_02_16_diatrofi_ffq_hcpc_dendogram_participants.pdf", useDingbats = F, onefile = T, width = 20, height=12)
#fviz_dend(res.hcpc)
#dev.off()
#data is too large for clustering participants by denogram; instead, continue with follow-on approach 

#pdf(file = "2023_02_16_diatrofi_ffq_hcpc_participant_clusters.pdf", useDingbats = F, onefile = T, width = 20, height=12)
##options(ggrepel.max.overlaps = Inf)
#fviz_cluster(res.hcpc,
#             repel=T, #avoid label overlapping
#             show.clust.cent=T, #show cluster centers
#             elipse.type="norm",
#             palette="jco", #color palette ?ggpubr::ggpar
#             ggtheme=theme_minimal(),
#             main="Factor map")
#dev.off()

#compare with previous cluster analysis 
#clean Rstudio
#rm(res.hcpc, res.pca)

        #### =============== 8. CLUSTER ANALYSIS (EUCLIDEAN AND BRAY CURTIS DISTANCES)   ================ ####
        
        ####8.1: EUCLIDEAN DISTANCE####
#NOTE: values were not scaled given that the values are already standardized 
clus_analysis <- dia_merged_ffq
head(rownames(dia_merged_ffq))
str(dia_merged_ffq)
#to avoid missings in rowVars function below, na.rm=T argument was added, however, missings have already been accounted for

clus_analysis <- as.matrix(clus_analysis)

#generate an index based on sample variance, as hierarchial clustering is sensitive to the order of cases
rv <- rowVars(clus_analysis) #calculating row variance for each participant: sample variance of food portions by participant
idx <- order(-rv) #index, which is based on the row variance (sorted from smallest to largest variance; high variance indicates differences in food portions by participant,
#i.e. 1 serving of fruits vs 5 servings of vegetables per day, vs low variance where serving sizes are more or less the same)

do.clustering(clus_analysis[idx,])->s1       #investigate possible h

#investigating possible clusters and cuts
hc <- hclust(dist(t(clus_analysis[idx,])))
#cutree(hc, k=1:4) #cut tree based on clustering
#hc_8_12 <- cutree(hc, k = c(8,12))
#table(grp8 = hc_8_12[,"8"], grp25 = hc_8_12[,"12"])

cutree(hc, h=48) #cut tree based on height

#generate circularized plot
my_colour <- colorRampPalette(brewer.pal(8, "Dark2"))
pdf(file = "2023_02_16_diatrofi_ffq_items_CA_eucl.pdf", useDingbats = F, onefile = T, width = 25, height=25)
hc <- as.dendrogram(hclust(dist(t(clus_analysis[idx,]))))
hc <- hc %>%
        color_branches(h = 48, col=my_colour(8)) %>%
        color_labels(h = 48, col=my_colour(8)) %>%
        set("branches_lwd", 3)
circlize_dendrogram(hc,
                    dend_track_height = 0.8)
dev.off()
#cut the tree at height 48

##NOTE: Generate Centroid ONLY if this is the agreed method for dietary pattern analysis

        ####8.2: BRAY CURTIS DISTANCE####
#Bray Curtis is best for categorical or binary data (Euclidean is best for continous data)
#Thus use the unweighted frequencies as input

#str(dia_merged_ffq_bray_curtis)
#head(rownames(dia_merged_ffq_bray_curtis))
#clus_analysis <- dia_merged_ffq_bray_curtis
##to avoid missings in rowVars function below, na.rm=T argument was added, however, rows with all NA's were filtered out first

#for (i in 1:ncol(clus_analysis)){
#        print(table(clus_analysis[,i], useNA = "ifany"))
#        clus_analysis[,i] <- gsub("_.*", "", clus_analysis[,i])
#        clus_analysis[,i] <- as.numeric(as.character(clus_analysis[,i]))
#        print(table(clus_analysis[,i], useNA = "ifany"))
#}

#clus_analysis <- as.matrix(clus_analysis)

#rv <- rowVars(clus_analysis) 
#idx <- order(-rv) 

#manually apply clustering to use Bray Curtis distance
#generate bray curtis matrix
#data.bray=vegdist(t(clus_analysis[idx,]))
#data.bray=as.matrix(data.bray)
#data.dist <- as.dist(data.bray)

#hclust.object=hclust(data.dist, method="complete") 
#plot(hang.dendrogram(as.dendrogram(hclust.object), cex=0.5, hang=-1)) #horiz=T, type="triangle", hang=-1, error: "hang" is not a graphical parameter

#hc <- hclust.object

#cutree(hc, h=0.65) #cut tree based on height
#clustering_overview <- as.data.frame(cutree(hc, h=0.65))

#generate circularized plot
#my_colour <- colorRampPalette(brewer.pal(8, "Dark2"))
#pdf(file = "2023_02_16_diatrofi_ffq_items_CA_BC.pdf", useDingbats = F, onefile = T, width = 25, height=25)
#hc <- as.dendrogram(hc)
#hc <- hc %>%
#        color_branches(h = 0.65, col=my_colour(10)) %>%
#        color_labels(h = 0.65, col=my_colour(10)) %>%
#        set("branches_lwd", 3)
#circlize_dendrogram(hc,
#                    dend_track_height = 0.8)
#dev.off()
##cut the tree at height 400