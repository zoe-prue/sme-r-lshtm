# SME Final Assessment
# 9/2/2026
# ‘Does being employed reduce your risk of starting antiretroviral therapy 
# one month after testing positive for HIV 
# and does this vary by level of education?’

###################################
# library calls

library(haven)
# install.packages("gtsummary")
library(gtsummary)    # For regression tables
library(marginaleffects)  # For predictions and contrasts
library(lmtest)       # For likelihood ratio tests
library(tidyverse)    # For data manipulation
library(here)
library(dplyr)
library(tidyr)
library(DescTools)
library(gmodels)
# install.packages("table1")
library(table1)

###################################
# global variables

# global variables
setwd("~/Desktop/sme-r-lshtm/")

# Read the hib dataset
hiv <- read_dta("raw-data/SMEassessment2026.dta") |>
  mutate(across(where(is.labelled), as.factor))

# sapply(mortality, class)

source(here("scripts", "or_function.R"))
source(here("scripts", "mh_functions_updated.R"))

options(digits = 3, scipen = 999)

###################################
# script

head(hiv)

# 3.define analysis plan
  # a.minimal adjustment set
  # ANSWER: age, sex, education, rural

# 4.data reduction (recoding of variables)
  # a.ordered categorical? save as new column. justify why in methods
    # i.rationale
    # ii.number of categories
    # iii.was this before data analysis or guided by data analysis?

# checking class
class(hiv$empstatus)
class(hiv$ARTstart)
class(hiv$agegp6)
class(hiv$sex)
class(hiv$educ)
class(hiv$ruralsite)

# only 3 observations in 55+ age group -> collapse into other age group
hiv <- hiv %>%
  mutate(
    agegp6 = fct_collapse(
      agegp6,
      "4" = c("4", "5")
    )
  )

# 5.cross tabulations

CrossTable(hiv$empstatus, hiv$ARTstart, prop.r = TRUE, chisq = TRUE)

  # a.data missingness check

colSums(is.na(hiv))

# 46 missing values in cd4grp6, none in the rest
# but remember, 12% of participants unreachable for follow-up

  # b.data missing will be removed from specific models, not the entire dataset!

# 6.descriptive analyses
  # will write up

# 7.summary tables

# Cross-tab with row percentages

# Cross-tabulation with row percentages

hiv |>
  count(ARTstart)

hiv |>
   count(empstatus, ARTstart) |>
   group_by(empstatus) |>
   mutate(percentage = n / sum(n) * 100,
          n_pct = sprintf("%d (%.1f%%)", n, percentage)) |>
   select(empstatus, ARTstart, n_pct) |>
   pivot_wider(names_from = ARTstart, values_from = n_pct)

tab1 <-
  hiv |>
  mutate(empstatus = recode(empstatus,
                      '0' = "Unemployed",
                      '1' = "Employed")) |>
  mutate(agegp6 = recode(agegp6,
                          '0' = '<25',
                          '1' = '25-29',
                          '2' = '30-34',
                          '3' = '35-44',
                          '4' = '≥45')) |>
  mutate(sex = recode(sex,
                      '0' = "Female",
                      '1' = "Male")) |>
  mutate(ruralsite = recode(ruralsite,
                      '0' = "Non-Rural",
                      '1' = "Rural")) |>
  mutate(educ = recode(educ,
                       '0' = "None/Primary",
                       '1' = "Secondary/Tertiary")) |>
  rename("Employment" = empstatus,
         "Age Group" = agegp6,
         "Sex" = sex,
         "Education" = educ,
         "Lives in Rural Area" = ruralsite) |>
  tbl_summary(
    by = ARTstart,
    statistic = all_categorical() ~ "{n} ({p}%)",
    percent = "row",
    missing = "no"
  ) |>
  add_overall(last = FALSE) |>
  bold_labels() |>
  modify_header(
    stat_1 ~ "**Did Not ART Within 1 Month**",
    stat_2 ~ "**Started ART Within 1 Month**"
  )

tab1

# 8.power/sample size calculation?

# ?

# 9.univariable analysis: crude measure of effect (calculate_or function)

# tab_crude_or <- table(hiv$empstatus, hiv$ARTstart)
# tab_crude_or
# calculate_or(tab_crude_or)

# PRIMARY logistic regression model 
mod_null <- glm(ARTstart ~ 1, family = binomial, data = hiv)

# crude logistic regression models
mod_empstatus_crude <- glm(ARTstart ~ empstatus, family = binomial, data = hiv)
mod_empstatus_crude |> 
  tbl_regression(
    exponentiate = TRUE, 
    label = list(empstatus = "Employment Status")) |>
  add_global_p(method="lrt")

mod_educ_crude <- glm(ARTstart ~ educ, family = binomial, data = hiv)
mod_educ_crude |> 
  tbl_regression(
    exponentiate = TRUE, 
    label = list(empstatus = "Education")) |>
  add_global_p(method="lrt")

mod_agegp6_crude <- glm(ARTstart ~ agegp6, family = binomial, data = hiv)
mod_agegp6_crude |> 
  tbl_regression(
    exponentiate = TRUE, 
    label = list(empstatus = "Age Group")) |>
  add_global_p(method="lrt")

mod_sex_crude <- glm(ARTstart ~ sex, family = binomial, data = hiv)
mod_sex_crude |> 
  tbl_regression(
    exponentiate = TRUE, 
    label = list(empstatus = "Sex")) |>
  add_global_p(method="lrt")

mod_rural_crude <- glm(ARTstart ~ ruralsite, family = binomial, data = hiv)
mod_rural_crude |> 
  tbl_regression(
    exponentiate = TRUE, 
    label = list(empstatus = "Ruralness of Home")) |>
  add_global_p(method="lrt")

# LRT tests
lrtest(mod_null, mod_empstatus_crude)
lrtest(mod_null, mod_educ_crude)
lrtest(mod_null, mod_agegp6_crude)
lrtest(mod_null, mod_sex_crude)
lrtest(mod_null, mod_rural_crude)

# reminder: 0 is unemployed
# Odds Ratio: 0.77 
# (95% CI: 0.585 - 1.014)
# P-value: 0.0622
# there was a 33% decrease in odds of starting ART within 1 month of HIV diagnoses 
# in people unemployed compared to those employed

# 10.stratified analysis (M-H on a couple of potential confounders? mhor() function)
  # a.compare crude vs adjusted estimates

mhor(hiv, outcome = "ARTstart", exposure = "empstatus", strata = "agegp6") # maybe confounder? poor evidence
mhor(hiv, outcome = "ARTstart", exposure = "empstatus", strata = "sex") # maybe not confounder
mhor(hiv, outcome = "ARTstart", exposure = "empstatus", strata = "educ") # maybe not confounder
mhor(hiv, outcome = "ARTstart", exposure = "empstatus", strata = "ruralsite") # maybe not confounder

# 11.multivariable analysis using regression

# PRIMARY logistic regression model 
mod_hiv <- glm(ARTstart ~ empstatus, family = binomial, data = hiv)
# Display with odds ratios and 95% CIs
mod_hiv |> 
  tbl_regression(
    exponentiate = TRUE, 
    label = list(empstatus = "Employment Status")) |>
  add_global_p(method="lrt")

# logistic regression with EACH CONFOUNDER individually

# primary logistic regression model + agegp6
mod_adj_agegp6 <- glm(ARTstart ~ empstatus + agegp6, family = binomial, data = hiv)
# Display with odds ratios and 95% CIs
mod_adj_agegp6 |> 
  tbl_regression(
    exponentiate = TRUE, 
    label = list(empstatus = "Employment Status", agegp6 = "Age Group")) |>
  add_global_p(method="lrt")

# Compare unadjusted and adjusted models side by side -> just a tbale to let me see
# tbl_merge(
#   list(
#     tbl_regression(mod_hiv, 
#                    exponentiate = TRUE, 
#                    include = empstatus,
#                    label = list(empstatus = "Employment Status")),
#     tbl_regression(mod_adj_agegp6, 
#                    exponentiate = TRUE, 
#                    include = c(empstatus, agegp6),
#                    label = list(empstatus = "Employment Status", agegp6 = "Age Group"))
#   ),
#   tab_spanner = c("**Unadjusted**", "**Adjusted for age**"))

# primary logistic regression model + sex
mod_adj_sex <- glm(ARTstart ~ empstatus + sex, family = binomial, data = hiv)
# Display with odds ratios and 95% CIs
mod_adj_sex |> 
  tbl_regression(
    exponentiate = TRUE, 
    label = list(empstatus = "Employment Status", sex = "Sex")) |>
  add_global_p(method="lrt")

# primary logistic regression model + educ
mod_adj_educ <- glm(ARTstart ~ empstatus + educ, family = binomial, data = hiv)
# Display with odds ratios and 95% CIs
mod_adj_educ |> 
  tbl_regression(
    exponentiate = TRUE, 
    label = list(empstatus = "Employment Status", educ = "Educational Status")) |>
  add_global_p(method="lrt")

# primary logistic regression model + ruralsite
mod_adj_rural <- glm(ARTstart ~ empstatus + ruralsite, family = binomial, data = hiv)
# Display with odds ratios and 95% CIs
mod_adj_rural |> 
  tbl_regression(
    exponentiate = TRUE, 
    label = list(empstatus = "Employment Status", ruralsite = "Rural Area")) |>
  add_global_p(method="lrt")

# a.fully adjusted measure of effect (3 categorical confounders) (LRT)

# adjusting for confounding for all

model_adjust_all <- glm(ARTstart ~ empstatus + agegp6 + sex + educ + ruralsite, 
                    data = hiv,
                    family = binomial()) 

tbl_regression(model_adjust_all, exponentiate = TRUE)

# Logistic regression with all confounders AND WITHOUT THE EXPOSURE

model_no_emp <- glm(ARTstart ~  agegp6 + sex + educ + ruralsite, # two potential confounders/effect modifiers are mfpos and agegrp
                     data = hiv,
                     family = binomial())

tbl_regression(model_no_emp, exponentiate = TRUE)

# LRT between model with 3 categorical confounders and the exposure 
# vs model with confounders and w/out exposure

lrtest(model_adjust_all, model_no_emp)

# p-value = 0.016
# this provides evidence that the model with the exposure 
# is a better fit than the model excluding the exposure

# Compare unadjusted and adjusted models side by side

tab2 <- tbl_merge(
  list(
    tbl_regression(mod_hiv, exponentiate = TRUE, include = "empstatus", 
                   label = list(empstatus ~ "Employment Status")),
    tbl_regression(model_adjust_all, exponentiate = TRUE, include = "empstatus",
                   label = list(empstatus ~ "Employment Status"))),
  tab_spanner = c("**Unadjusted**", "**Adjusted For Covariates**")
) 

tab2

as_flex_table(tab2)

  # d.modifying effect of a covariate with 'fully adjusted' (if necessary)
    # i.testing for effect modification? (LRT? Wald?)

# interaction logistic regression model + agegp6
mod_effect_agegp6 <- glm(ARTstart ~ empstatus * agegp6, family = binomial, data = hiv)
# Display with odds ratios and 95% CIs
mod_effect_agegp6 |> 
  tbl_regression(
    exponentiate = TRUE, 
    label = list(empstatus = "Employment Status", agegp6 = "Age Group")) |>
  add_global_p(method="lrt")
# lrtest(mod_adj_agegp6, mod_effect_agegp6) # just checking that theyre the same p-values
# p=0.022

# interaction logistic regression model + sex
mod_effect_sex <- glm(ARTstart ~ empstatus * sex, family = binomial, data = hiv)
# Display with odds ratios and 95% CIs
mod_effect_sex |> 
  tbl_regression(
    exponentiate = TRUE, 
    label = list(empstatus = "Employment Status", sex = "Sex")) |>
  add_global_p(method="lrt")
# p=0.4

# interaction logistic regression model + educ
mod_effect_educ <- glm(ARTstart ~ empstatus * educ, family = binomial, data = hiv)
# Display with odds ratios and 95% CIs
mod_effect_educ |> 
  tbl_regression(
    exponentiate = TRUE, 
    label = list(empstatus = "Employment Status", educ = "Educational Status")) |>
  add_global_p(method="lrt")
# 0.016

# effect logistic regression model + ruralsite
mod_effect_rural <- glm(ARTstart ~ empstatus * ruralsite, family = binomial, data = hiv)
# Display with odds ratios and 95% CIs
mod_effect_rural |> 
  tbl_regression(
    exponentiate = TRUE, 
    label = list(empstatus = "Employment Status", ruralsite = "Rural Area")) |>
  add_global_p(method="lrt")
# p=0.4

# model in totality?

model_adjust_all <- glm(ARTstart ~ empstatus + agegp6 + sex + educ + ruralsite, 
                        data = hiv,
                        family = binomial()) 

tbl_regression(model_adjust_all, exponentiate = TRUE)

# model with one effect modifier (educ)
model_chosen_effect <- glm(ARTstart ~ empstatus * educ + agegp6 + sex + ruralsite, 
                        data = hiv,
                        family = binomial())

# ADD CONF.INT = TRUE???

tbl_regression(model_chosen_effect, exponentiate = TRUE) 

# TEST HETEROGENEITY OF THESE ODDS RATIOS?????? REPORT STRATUM SPECIFIC

# LRT

lrtest(model_adjust_all, model_chosen_effect)

# stat specific ORs effect mod

library(emmeans)

strat_spec_effect <- emmeans(model_chosen_effect, revpairwise ~ empstatus | educ, type = 'response') |> 
  pluck('contrasts') |>
  confint()

strat_spec_effect

# DO NOT REPORT THE POOLED OR FOR EFFECT MODIFICATION
# tab3 <- tbl_merge(
#   list(
#     tbl_regression(model_adjust_all, exponentiate = TRUE, include = "empstatus", 
#                    label = list(empstatus ~ "Employment Status")),
#     tbl_regression(model_chosen_effect, exponentiate = TRUE, include = "empstatus",
#                    label = list(empstatus ~ "Employment Status"))),
#   tab_spanner = c("**Adjusted for Covariates**", "**Adjusted for Covariates and Educational Interaction**")
# ) 
# 
# tab3
# 
# as_flex_table(tab3)

unique(hiv$reastest2)
levels(hiv$reastest2)
summary(hiv)



