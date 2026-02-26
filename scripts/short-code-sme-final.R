# SME Final Assessment
# 9/2/2026
# ‘Does being employed reduce your risk of starting antiretroviral therapy 
# one month after testing positive for HIV 
# and does this vary by level of education?’

###################################
# library calls

library(haven)
# install.packages("gtsummary")
# install.packages("gtregression")
library(gtregression)
library(gtsummary)    # For regression tables
library(marginaleffects)  # For predictions and contrasts
library(lmtest)       # For likelihood ratio tests
library(tidyverse) 
library(emmeans)# For data manipulation
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

# Read the Whitehall dataset
hiv <- read_dta("raw-data/SMEassessment2026.dta") |>
  mutate(across(where(is.labelled), as.factor))

# sapply(mortality, class)

source(here("scripts", "or_function.R"))
source(here("scripts", "mh_functions_updated.R"))

options(digits = 3, scipen = 999)

###################################
# script

# only 3 observations in 55+ age group -> collapse into other age group
hiv <- hiv %>%
  mutate(
    agegp6 = fct_collapse(
      agegp6,
      "4" = c("4", "5")))

#data missingness check

colSums(is.na(hiv))

# tab1 <-
#   hiv |>
#   select(ARTstart, empstatus, agegp6, sex, educ, ruralsite) |>
#   mutate(empstatus = recode(empstatus,
#                             '0' = "None/Primary Education",
#                             '1' = "Secondary/Tertiary Education")) |>
#   mutate(agegp6 = recode(agegp6,
#                          '0' = '<25',
#                          '1' = '25-29',
#                          '2' = '30-34',
#                          '3' = '35-44',
#                          '4' = '≥45')) |>
#   mutate(sex = recode(sex,
#                       '0' = "Female",
#                       '1' = "Male")) |>
#   mutate(ruralsite = recode(ruralsite,
#                             '0' = "Non-Rural",
#                             '1' = "Rural")) |>
#   mutate(educ = recode(educ,
#                        '0' = "None/Primary",
#                        '1' = "Secondary/Tertiary")) |>
#   rename("Employment" = empstatus,
#          "Age Group" = agegp6,
#          "Sex" = sex,
#          "Education" = educ,
#          "Lives in Rural Area" = ruralsite) |>
#   tbl_summary(
#     by = ARTstart,
#     statistic = all_categorical() ~ "{n} ({p}%)",
#     percent = "row",
#     missing = "no"
#   ) |>
#   add_overall(last = FALSE) |>
#   bold_labels() |>
#   modify_header(
#     stat_1 ~ "**Did Not ART Within 1 Month**",
#     stat_2 ~ "**Started ART Within 1 Month**"
#   )
# 
# tab1
# 
# as_flex_table((tab1))

# # global variables
# setwd("~/Desktop/sme-r-lshtm/outputs")
# 
# tab1 |>
#   as_gt() |>
#   gt::gtsave(filename = "tab1.png") # use extensions .png, .html, .docx, .rtf, .tex, .ltx

# options(gtsummary.print_engine = "tab1") 

# crude measure of effect (calculate_or function)

tab_crude_or <- table(hiv$empstatus, hiv$ARTstart)
tab_crude_or
calculate_or(tab_crude_or)

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
# mod_hiv |> 
#   tbl_regression(
#     exponentiate = TRUE, 
#     label = list(empstatus = "Employment Status")) |>
#   add_global_p(method="lrt")

# mod_hiv |>
#   tbl_regression(
#     exponentiate = TRUE,
#     label = list(
#       empstatus = "Employment Status"
#     )
#   ) |>
#   modify_caption("**Unadjusted Odds Ratios of Beginning ART Within 1 Month**") |>
#   add_global_p(method = "lrt")

# test <- mod_hiv |>
#   tbl_regression(
#     exponentiate = TRUE,
#     label = list(empstatus = "Employment Status")) |>
#   modify_caption("**Unadjusted Odds Ratios of Beginning ART Within 1 Month**") |>
#   add_global_p(method = "lrt") 
# 
# test
# 
# # change what prints in the rows
# test |>
#   modify_table_body(
#     ~ .x |>
#       dplyr::mutate(
#         label = dplyr::case_when(
#           variable == "empstatus" & label == "0" ~ "Unemployed",
#           variable == "empstatus" & label == "1" ~ "Employed",
#           TRUE ~ label
#         )
#       )
#   )

# adjusting for confounding for all

model_adjust_all <- glm(ARTstart ~ empstatus + agegp6 + sex + educ + ruralsite, 
                        data = hiv,
                        family = binomial()) 

tbl_regression(model_adjust_all, exponentiate = TRUE)

# # Logistic regression with all confounders AND WITHOUT THE EXPOSURE
# 
# model_no_emp <- glm(ARTstart ~  agegp6 + sex + educ + ruralsite, # two potential confounders/effect modifiers are mfpos and agegrp
#                     data = hiv,
#                     family = binomial())
# 
# tbl_regression(model_no_emp, exponentiate = TRUE)
# 
# # LRT between model with 3 categorical confounders and the exposure 
# # vs model with confounders and w/out exposure
# 
# lrtest(model_adjust_all, model_no_emp)

# merge the adjusted and unadjusted models into a table
# 
# tab2 <- tbl_merge(
#   list(
#     tbl_regression(mod_hiv,
#                    exponentiate = TRUE,
#                    include = "empstatus",
#                    label = list(empstatus ~ "Employment Status")),
#     tbl_regression(model_adjust_all, exponentiate = TRUE, include = "empstatus",
#                    label = list(empstatus ~ "Employment Status"))),
#   tab_spanner = c("**Unadjusted**", "**Adjusted For Covariates**")
# )
# 
# tab2

# # change what prints in the rows
# tab2 |>
#   modify_table_body(
#     ~ .x |>
#       dplyr::mutate(
#         label = dplyr::case_when(
#           variable == "empstatus" & label == "0" ~ "Unemployed",
#           variable == "empstatus" & label == "1" ~ "Employed",
#           TRUE ~ label
#         )
#       )
#   )
# 
# tab2 |>
#   as_gt() |>
#   gt::gtsave(filename = "tab2.png") # use extensions .png, .html, .docx, .rtf, .tex, .ltx
# 
# options(gtsummary.print_engine = "tab1")

# asterick for what covariates -> footnote

model_adjust_all <- glm(ARTstart ~ empstatus + agegp6 + sex + educ + ruralsite, 
                        data = hiv,
                        family = binomial()) 

tbl_regression(model_adjust_all, exponentiate = TRUE)

# model with one effect modifier (educ)
model_chosen_effect <- glm(ARTstart ~ empstatus * educ + agegp6 + sex + ruralsite, 
                           data = hiv,
                           family = binomial()) 

tbl_regression(model_chosen_effect, exponentiate = TRUE) 

lrtest(model_adjust_all, model_chosen_effect)

# stat specific ORs effect mod

strat_spec_effect <- emmeans(model_chosen_effect, ~ empstatus | educ) |>
  confint()

strat_spec_effect


# # change what prints in the rows
# tab3 |>
#   modify_table_body(
#     ~ .x |>
#       dplyr::mutate(
#         label = dplyr::case_when(
#           variable == "empstatus" & label == "0" ~ "Unemployed",
#           variable == "empstatus" & label == "1" ~ "Employed",
#           TRUE ~ label
#         )
#       )
#   )

# lrtest(model_adjust_all, model_chosen_effect)

# see HETEROGENEITY OF THESE ODDS RATIOS

# strat_spec_effect <- emmeans(model_chosen_effect, ~ empstatus | educ)

# emmeans(model_chosen_effect,
#         revpairwise ~ empstatus | educ,
#         at = list(educ = unique(hiv$educ)),
#         type = "response") |>
#   pluck("contrasts")

# Extract contrasts and convert to data frame
# ADD CONFIDENCE INTERVAL TO THIS FUNCTION??
# contrasts_tbl <- emmeans(model_chosen_effect,
#                          revpairwise ~ empstatus | educ,
#                          at = list(educ = unique(hiv$educ)),
#                          type = "response") |>
#   pluck("contrasts") |>
#   as.data.frame()
# 
# library(dplyr)
# library(tidyr)
# library(stringr)
# 
# # Create baseline rows first
# baseline <- tibble(
#   Education = c("None/Primary Education", "Secondary/Tertiary Education"),
#   Emp_status = "Unemployed",
#   OR = 1.0,
#   `95% CI` = "-",
#   P_value = NA
# )
# 
# table_effmod <- contrasts_tbl |>
#   separate(contrast, into = c("empstatus_ref", "empstatus_cmp"), 
#            sep = " / ", remove = FALSE, extra = "merge") |>
#   mutate(
#     Education = ifelse(educ == 0, "None/Primary Education", "Secondary/Tertiary Education"),
#     Emp_status = "Employed",
#     OR = round(odds.ratio, 2),
#     logOR_se = SE,
#     lower_CI = round(exp(log(odds.ratio) - 1.96 * SE), 2),
#     upper_CI = round(exp(log(odds.ratio) + 1.96 * SE), 2),
#     `95% CI` = sprintf("(%.2f, %.2f)", lower_CI, upper_CI),
#     P_value = p.value
#   ) |>
#   select(Education, Emp_status, OR, `95% CI`, P_value) |>
#   # Bind with baseline rows and arrange
#   bind_rows(baseline) |>
#   arrange(Education, desc(Emp_status == "Unemployed"))
# 
# table_effmod
# 
# # pkgs <- c("flextable", "broom", "report", "effectsize", "rempsyc")
# # install.packages(pkgs)
# library(rempsyc)
# 
# table_effmod
# 
# table_effmod_nice <- nice_table(table_effmod,
#            note = c(
#   "LRT p-value",
#   "* p < .05, ** p < .01, *** p < .001"))
# 
# table_effmod_nice
# 
# 
# # REPORT STRATUM SPECIFIC
# # refer to example
# 
# # LRT
# 
# lrtest(model_adjust_all, model_chosen_effect)
# 
# # collinearity matrix
# 
# install.packages("olsrr")
# library(olsrr)
# 
# collinear <- lm(ARTstart ~ empstatus + agegp6 + sex + educ + ruralsite, data = hiv)
# ols_vif_tol(collinear)
# 
# 
# # ## traing emmeans table another way
# # 
# # library(emmeans)
# # 
# # # Get EMMs
# # emm <- emmeans(model_chosen_effect, ~ empstatus | educ)
# # 
# # # Convert to data frame
# # emm_df <- as.data.frame(emm)
# # print(emm_df)
# # 
# # # For contrasts (pairwise comparisons)
# # pairs_df <- emmeans(model_chosen_effect, pairwise ~ empstatus | educ)$contrasts %>%
# #   as.data.frame() %>%
# #   mutate(
# #     estimate = round(estimate, 2),
# #     SE = round(SE, 2),
# #     p.value = round(p.value, 3)
# #   )
# # 
# # print(pairs_df[, c("contrast", "estimate", "SE", "p.value")])
# # 
# # library(gt)
# # library(dplyr)
# # 
# # contrasts_gt <- pairs_df %>%
# #   gt() %>%
# #   cols_label(contrast = "Contrast") %>%
# #   fmt_number(columns = c(estimate, SE), decimals = 2) %>%
# #   fmt_number(columns = p.value, decimals = 3) %>%
# #   tab_style(
# #     style = cell_text(weight = "bold"),
# #     locations = cells_body(columns = p.value, rows = p.value < 0.05)
# #   )
# # 
# # contrasts_gt
# # 
# # conf_df <- confint(emm) %>% as.data.frame()
# # print(conf_df)
# # 
# # # Assuming 'emm_pairs' is your emmeans contrasts object
# # # emm_pairs <- pairs(emm, type = "response", adjust = "none")  # or your existing object
# # 
# # # Convert and format table
# # table_df <- as.data.frame(pairs_df) %>%
# #   select(contrast, estimate, lower.CL, upper.CL) %>%
# #   mutate(
# #     Education = case_when(
# #       grepl("No education", contrast) ~ "No education",
# #       grepl("Some education", contrast) ~ "Some education",
# #       TRUE ~ NA_character_
# #     ),
# #     `Visual acuity` = case_when(
# #       grepl("Unimpaired", contrast) ~ "Unimpaired",
# #       grepl("Impaired", contrast) ~ "Impaired",
# #       TRUE ~ NA_character_
# #     ),
# #     OR = sprintf("%.1f", odds.ratio),
# #     `95% CI` = sprintf("(%s,%s)", 
# #                        sprintf("%.1f", asymp.L), 
# #                        sprintf("%.1f", asymp.U))
# #   ) %>%
# #   select(Model = Education, `Employment`, OR, `95% CI`) %>%
# #   filter(!is.na(Education))
# # 
# # print(table_df)
# # 
# # 
# # 
# # 
# # 
