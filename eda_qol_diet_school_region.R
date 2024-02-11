library(ggpubr)
library(rstatix) #tukey test
library(ggmisc)

#MAC-book location
setwd("~/OneDrive - UMCG/Prolepsis_internship/merged_2015_2018/masterfile/")
masterfile <- read.delim("2023_04_17_diatrofi_2015_2018_masterfile_incl_marg.txt", sep = "\t", header = T, stringsAsFactors = T)
#6583 60
str(masterfile)
masterfile$child_demo_school_code <- as.factor(masterfile$child_demo_school_code)

test <- masterfile
#investigate normalacy of the hPDI score in relation to school_region

  ####hPDI and QoL####
setwd("~/OneDrive - UMCG/Prolepsis_internship/merged_2015_2018/masterfile/eda_diet_school_region/")
pdf(file = "2023_05_02_diet_school_region_hist_tukey.pdf", useDingbats = F, onefile = T, width = 20, height=12)
gghistogram(data=subset(test, !is.na(child_demo_school_region) & !is.na(child_ffq_hpdi_score)),
            x = "child_ffq_hpdi_score", bins = 15, 
            fill = "child_demo_school_region", color = "child_demo_school_region",
            palette = "Set1", 
            add = "mean", mean.color = "black",
            mean.size = 1, mean.linetype = "dashed",
            alpha = 0.5) + 
  #scale_fill_manual(values = c("attica" = "#66C2A5", 
  #                             "central_greece" = "#FC8D62", 
  #                             "central_macedonia" = "#8DA0CB", 
  #                             "eastern_macedonia_and_thrace" = "#E78AC3",
  #                             "peloponnese" = "#A6D854",
  #                             "thessaly"= "#FFD92F",
  #                             "western_greece"= "#E5C494")) + 
  labs(title = "Histogram of hPDI by School Region", 
       x = "hPDI", y = "Frequency", 
       fill = "School Region") + 
  theme_pubclean()

#distribution of hPDI across school region looks to be normally distributed, thus we can perform a one-way ANOVA test

# Perform one-way ANOVA test
model <- aov(child_ffq_hpdi_score ~ child_demo_school_region, data = test)
summary(model)
p_value <- summary(model)[[1]][["Pr(>F)"]][[1]]
p_value#1.154168e-17

#calculate Tukey's HSD test
tukey <- aov(child_ffq_hpdi_score ~ child_demo_school_region, data = test) %>%
  tukey_hsd()
tukey_sig <- tukey %>% filter(
  p.adj.signif != "ns"
)

# Plot using ggline and add violin plot
#significant according to Tukey
my_comparisions=list(c("attica", "central_macedonia"),
                     c("attica", "eastern_macedonia_and_thrace"),
                     c("central_greece", "central_macedonia"),
                     c("central_macedonia", "eastern_macedonia_and_thrace"),
                     c("central_macedonia", "thessaly"),
                     c("central_macedonia", "western_greece"))

p1 <- ggline(data=subset(test, !is.na(child_demo_school_region) & !is.na(child_ffq_hpdi_score)),
       x = "child_demo_school_region", y = "child_ffq_hpdi_score", 
       color = "child_demo_school_region", 
       add = c("mean_se", "violin"),
       main = "Mean Diet quality score (SE) by School Region",
       ylab = "Mean diet quality score (SE)", xlab = "",
       legend = "none") +
  stat_compare_means(comparisons = my_comparisions) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1 <- ggpar(p1, 
            palette = "Dark2",
            legend = "none")
print(p1)
dev.off()

    ####HRQoL and school region####
####generate the continuous percentile of HRQoL ####
pdf(file = "2023_05_02_hrqol_as_percentile.pdf", useDingbats = F, onefile = T, width = 20, height=12)
ggplot(test, aes(child_pedsql_health_related_qol_score_total)) + stat_ecdf(geom = "step")

#calculate the cumulative distribution function
hqrol_ecdf <- ecdf(test$child_pedsql_health_related_qol_score_total)

# Calculate the percentile for each value of QoL
test$hrqol_percentile <- ecdf(test$child_pedsql_health_related_qol_score_total)(test$child_pedsql_health_related_qol_score_total)*100
hist(test$hrqol_percentile)

hqrol_test <- test %>% select (c(
  child_pedsql_health_related_qol_score_total,
  hrqol_percentile))
#coincides with the ggplot above

pdf(file = "2023_05_02_QoL_school_region_hist_tukey.pdf", useDingbats = F, onefile = T, width = 20, height=12)
gghistogram(data=subset(test, !is.na(child_demo_school_region) & !is.na(hrqol_percentile)),
            x = "hrqol_percentile", bins = 15, 
            fill = "child_demo_school_region", color = "child_demo_school_region",
            palette = "Set1", 
            add = "mean", mean.color = "black",
            mean.size = 1, mean.linetype = "dashed",
            alpha = 0.5) + 
  #scale_fill_manual(values = c("attica" = "#66C2A5", 
  #                             "central_greece" = "#FC8D62", 
  #                             "central_macedonia" = "#8DA0CB", 
  #                             "eastern_macedonia_and_thrace" = "#E78AC3",
  #                             "peloponnese" = "#A6D854",
  #                             "thessaly"= "#FFD92F",
  #                             "western_greece"= "#E5C494")) + 
  labs(title = "Histogram of HRQoL by School Region", 
       x = "HRQoL", y = "Frequency", 
       fill = "School Region") + 
  theme_pubclean()

#distribution of hrqol (in percentiles) across school region looks to be evenly distributed

# Perform one-way ANOVA test
model <- aov(hrqol_percentile ~ child_demo_school_region, data = test)
summary(model)
#1.77e-11 ***

#calculate Tukey's HSD test
tukey <- aov(hrqol_percentile ~ child_demo_school_region, data = test) %>%
  tukey_hsd()

tukey_sig <- tukey %>% filter(
  p.adj.signif != "ns")

# Plot using ggline and add violin plot
#significant according to Tukey
my_comparisions=list(c("attica", "central_greece"),
                     c("attica", "central_macedonia"),
                     c("attica", "eastern_macedonia_and_thrace"),
                     c("central_greece", "thessaly"),
                     c("central_greece", "western_greece"),
                     c("central_macedonia", "thessaly"),
                     c("central_macedonia", "western_greece"),
                     c("eastern_macedonia_and_thrace", "thessaly"),
                     c("eastern_macedonia_and_thrace", "western_greece"),
                     c("peloponnese", "thessaly"),
                     c("peloponnese", "western_greece"))

p1 <- ggline(data=subset(test, !is.na(child_demo_school_region) & !is.na(hrqol_percentile)),
             x = "child_demo_school_region", y = "hrqol_percentile", 
             color = "child_demo_school_region", 
             add = c("mean_se", "violin"),
             main = "Mean HRQoL (in percentiles; SE) by School Region",
             ylab = "Mean HRQoL (in percentiles; SE)", xlab = "",
             legend = "none") +
  stat_compare_means(comparisons = my_comparisions) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1 <- ggpar(p1, 
            palette = "Dark2",
            legend = "none")
print(p1)
dev.off()


####investigate the interaction of diet and school region on quality of life####
pdf(file = "2023_05_02_interactions_school_region_diet_qol_incl_SES.pdf", useDingbats = F, onefile = T, width = 20, height=12)
ggplot(data=subset(test, !is.na(child_demo_school_region) & !is.na(hrqol_percentile) & !is.na(child_ffq_hpdi_score)),
       aes(x = child_ffq_hpdi_score, y = hrqol_percentile, colour = child_demo_school_region)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw()

#I need to remove the areas with little variance (i.e. peloponnese, thessaly and western greece)

test$school_region_zv <- test$child_demo_school_region
test$school_region_zv <- as.character(test$school_region_zv)
test$school_region_zv[test$school_region_zv=="peloponnese" |
                        test$school_region_zv=="thessaly" |
                        test$school_region_zv=="western_greece"] <- NA
table(test$school_region_zv, useNA = "ifany")

ggplot(data=subset(test, !is.na(school_region_zv) & !is.na(hrqol_percentile) & !is.na(child_ffq_hpdi_score)),
       aes(x = child_ffq_hpdi_score, y = hrqol_percentile, colour = school_region_zv)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw()

#what is the correlation between school region and SES
x <- chisq.test(test$school_region_zv, test$family_demo_FAS_group, correct=FALSE)
x$p.value
#school regions and FAS groups are significantly different

school_FAS <- test %>%
  with(gmodels::CrossTable(school_region_zv,
                           family_demo_FAS_group,
                           chisq=T,
                           prop.r=T,
                           prop.t=T,
                           prop.chisq=T,
                           expected=T))

#there is more low SES in eastern macedonia and thrace

ggplot(data=subset(test, !is.na(school_region_zv) & !is.na(hrqol_percentile) 
                   & !is.na(child_ffq_hpdi_score) & !is.na(family_demo_FAS_group)),
       aes(x = child_ffq_hpdi_score, y = hrqol_percentile, colour = school_region_zv)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~family_demo_FAS_group) +
  theme_bw()

ggplot(data=subset(test, !is.na(school_region_zv) & !is.na(hrqol_percentile) 
                   & !is.na(child_ffq_hpdi_score) & !is.na(child_demo_immigrant_status)),
       aes(x = child_ffq_hpdi_score, y = hrqol_percentile, colour = school_region_zv)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~child_demo_immigrant_status) +
  theme_bw()
dev.off()