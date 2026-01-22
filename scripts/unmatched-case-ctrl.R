# Unmatched Case-Control Analysis Practical 6
# 22/1/2026
# # In this practical, you will learn how to analyze unmatched 
# case-control studies using logistic regression in R. 
# You'll work with data from the Mwanza, 
# Tanzania study on HIV infection conducted in the 1990s. 
# By the end of this session, you should be able to:
# 
# - Calculate crude odds ratios and confidence intervals for case-control studies
# - Assess confounding using stratified analysis
# - Test for effect modification (using statistical interaction term)
# - Handle missing data appropriately in epidemiological analyses
# - Explore dose-response relationships


###################################
# library calls

# Load required packages
library(haven)            
library(here)
library(tidyverse)
library(gtsummary)
install.packages("data.table")
library(data.table)

# Source the Mantel-Haenszel functions
# Loads the function so that they are ready to use
setwd("~/Desktop/sme-r-lshtm/scripts")
source(here("scripts", "mh_functions.R")) 
source(here("scripts", "or_function.R"))

# Set options for cleaner output
options(digits = 3, scipen = 999)

###################################
# global variables

# Import the dataset
mwanza <- read_stata(here("raw-data", "mwanza.dta")) |> 
  mutate(across(where(is.labelled), as_factor))

###################################
# script

# Data Import and Preparation

## Reading the Mwanza dataset

# Examine the structure
glimpse(mwanza)
summary(mwanza)

# Question 1: Generate Appropriate Variables #

# Recategorize the variables "ed" and "age1" into two new variables:
#   
# - `ed2`: 1 = none/adult only, 2 = 1 or more years of formal education
# - `age2`: 1 = 15-24, 2 = 25-34, 3 = 35+ years
# 
# Turn these variables into factors and check they were created correctly.

# First examine the original variables
mwanza |> count(ed)
mwanza |> count(age1)

# Recategorize and label education level
mwanza <- mwanza |> 
  mutate(ed2 = factor(
    case_when(
      ed == 1 ~ "None/adult only",
      ed %in% c(2, 3, 4) ~ "≥1 years"),
    levels = c("None/adult only", "≥1 years")
  ))

# Recategorize and label age
mwanza <- mwanza |> 
  mutate(age2 = factor(
    case_when(
      age1 %in% c(1, 2) ~ "15-24",
      age1 %in% c(3, 4) ~ "25-34",
      age1 %in% c(5, 6) ~ "35+"),
    levels = c("15-24", "25-34", "35+")
  ))

# Check the recoding worked correctly
mwanza |> count(ed, ed2)
mwanza |> count(age1, age2)

# Make "case" a factor
mwanza <- mwanza |> 
  mutate(case = factor(case,
                       levels = c(0, 1),
                       labels = c("Control", "HIV case")
  ))

# Also make the case variable into a factor:

# Check
mwanza |> count(case)

# Question 2: Obtain the Unadjusted Odds Ratio for Education #

# Reproduce the 2×2 table of education (ed2) by HIV status showing counts and percentages.
# 
# Calculate the unadjusted odds ratio and 95% confidence interval. 

# Create 2x2 table with counts
tab_case <- table(mwanza$case, mwanza$ed2)
tab_case

# With proportions (by row - showing % within cases and controls)

mwanza |> 
  
  # counts cases by education level
  count(case, ed2) |>  
  
  # group calculation going forward by case status
  group_by(case) |>
  
  # calculate % within each case group and create variable which 
  # combines n and % and prints in this particular format
  mutate(percentage = n / sum(n) * 100, 
         n_pct = sprintf("%d (%.1f%%)", n, percentage)) |> 
  
  # keep only these selected variables 
  select(case, ed2, n_pct) |> 
  
  # transform dataframe from long format to wide format
  pivot_wider(names_from = ed2, values_from = n_pct) 

# An important exercise is to run the code line by line so that you understand 
# what each bit of code is doing 
# - you will have to make an edit to the code above to run it line by.


# Why do you think we use row percentages here?
# ANSWER: we can look at the odds ratio as odds of exposure within the outcome group
# odds of HIV in educated vs non-educated group
# only view this way bc it is a case-control study?

# Now calculate with the custom calculate_or function 
# (created in the or_function.R script)

# Get OR with custom function

calculate_or(tab_case)

# Calculate unadjusted OR with M-H

mhor(mwanza, "case", "ed2")

# How does the output compare with the lecture notes?
# ANSWER: 

# Understanding the output:
#   
# - The estimate is the odds ratio comparing the odds of being a case 
# for those with "≥1 years" education to "None/adult only"
# - An OR > 1 means higher odds of HIV with more education; OR < 1 means lower odds
# 
# Test for association using chi-squared test:

# Chi-squared test
chisq.test(table(mwanza$case, mwanza$ed2))

# Fisher's exact test (for small cell counts)
# for when less than 5 people in the 2x2 table
fisher.test(table(mwanza$case, mwanza$ed2))

# The null hypothesis of the test is that there is no association 
# between education and case status. 
# The p-value indicates the strength of evidence against this null hypothesis. 
# Fisher's exact test is particularly useful for small cell counts, 
# though with this sample size the chi-squared approximation is adequate. 
# In case-control studies, the p-value from the logistic regression 
# (shown in tbl_regression) tests the same hypothesis: is the OR = 1?
# ANSWER: there is strong evidence to rejec tthe null hypothesis that 
# the OR is equal to 1 between the two categorical exposure groups

# Optional - practical 9 will go in more detail about logistic regression
# note the similar code to Poisson regression from practical 2

m1 <- glm(case ~ ed2, family = binomial, data = mwanza)
summary(m1)

tbl_regression(m1,
               exponentiate = TRUE,
               label = list(ed2 = "Education"))

tbl_regression

# Question 3: Assessing Whether Age is a Confounder or Effect Modifier #

# A variable may act as a confounder 
# (a common cause of exposure and outcome) 
# or an effect modifier 
#(where the effect of the exposure on the outcome differs across strata). 
#To assess effect modification, we can examine whether the education-HIV association 
# differs across age groups using stratified analysis and a hypothesis test 
# for effect modification.

## Exploratory stratified analysis

# First, let's look at descriptive 2×2 tables within each age stratum:

# Cross-tabulations within each age group

mwanza |> # stratified analysis between each age group
  filter(age2 == "15-24") |> 
  count(case, ed2) |> 
  group_by(case) |> 
  mutate(percentage = n / sum(n) * 100,
         n_pct = sprintf("%d (%.1f%%)", n, percentage)) |> 
  select(case, ed2, n_pct) |> 
  pivot_wider(names_from = ed2, values_from = n_pct) 


mwanza |> # stratified analysis between each age group
  filter(age2 == "25-34") |> 
  count(case, ed2) |> 
  group_by(case) |> 
  mutate(percentage = n / sum(n) * 100) |> 
  mutate(n_pct = sprintf("%d (%.1f%%)", n, percentage)) |> 
  select(case, ed2, n_pct) |> 
  pivot_wider(names_from = ed2, values_from = n_pct)


mwanza |> # stratified analysis between each age group
  filter(age2 == "35+") |> 
  count(case, ed2) |> 
  group_by(case) |> 
  mutate(percentage = n / sum(n) * 100) |> 
  mutate(n_pct = sprintf("%d (%.1f%%)", n, percentage)) |> 
  select(case, ed2, n_pct) |> 
  pivot_wider(names_from = ed2, values_from = n_pct)

## Adjusted analysis

# Now let's use the Mantel-Haenszel method to get odds ratios, both unadjusted and adjusted for age:
# 
# Mantel-Haenszel analysis stratified by age
mhor(mwanza, "case", "ed2")
mhor(mwanza, "case", "ed2", strata = "age2")

# The adjusted OR is like a weighted average of the age-specific ORs 
# (this is the Mantel-Haenszel estimator). Compare the crude and adjusted estimates:
#   
# - If the adjusted OR moves substantially toward 1.0, age *explains* part of the association
# - If the adjusted OR moves away from 1.0, age was *masking* part of the association
# - If they're similar, age is *not distorting* the observed relationship
# 
# Note that this does not necessarily mean that age is not a confounder 
# **in the population** even if it does not distort any observed associations 
# in the sample: confounding is not a property of the sample or a property of the data.
# 
# Questions:
# 
# - What is the M-H adjusted odds ratio for education?
# ANSWER: 2.331

# - Does it provide evidence against the null hypothesis 
# of no association after adjusting for age?
# ANSWER: Yes, because the 95% confidence interval 
# (95% CI: 1.549 - 3.508) doesn't contain the null value 


# - Is there evidence of effect modification (do the stratum-specific ORs differ)?
# ANSWER: The stratum specific OR's do vary, 
# especially without their 95% confidence intervals all overlapping
# showing evidence for differing effect modification between groups
# and the p-value for homogeneity is 0.0051, 
# providing evidence to reject the null hypothesis

# Question 4: Does Religion Confound the Association? #

# Investigate whether religion (rel) confounds the association 
# between education (ed2) and HIV infection (case):
#   
# First, handle missing values in the religion variable (coded as 9):

# Recode religion and handle missing values
mwanza <- mwanza |> 
  mutate(rel = factor(
    case_when(
      rel == 1 ~ "Muslim",
      rel == 2 ~ "Catholic",
      rel == 3 ~ "Protestant",
      rel == 4 ~ "Other",
      rel == 9 ~ NA_character_
    )
  ))

# Check the distribution
mwanza |> count(rel)

## check associations

# Religion and HIV status
mwanza |> 
  filter(!is.na(rel)) |> 
  count(case, rel) |> 
  group_by(case) |> 
  mutate(percentage = n / sum(n) * 100) |> 
  mutate(n_pct = sprintf("%d (%.1f%%)", n, percentage)) |> 
  select(case, rel, n_pct) |> 
  pivot_wider(names_from = rel, values_from = n_pct)


# # Now examine association between religion and education 
mwanza |> 
  filter(!is.na(rel) & case == "Control") |> 
  count(ed2, rel) |> 
  group_by(ed2) |> 
  mutate(percentage = n / sum(n) * 100) |> 
  mutate(n_pct = sprintf("%d (%.1f%%)", n, percentage)) |> 
  select(ed2, rel, n_pct) |> 
  pivot_wider(names_from = rel, values_from = n_pct)

# Why did we restrict this analysis of the association 
# between religion and education to controls only?
# ANSWER: 

## Unadjusted and adjusted analysis

# When assessing confounding, the unadjusted and adjusted ORs must be based on 
# the same individuals. We use `filter(!is.na(rel))` for the unadjusted model 
# to match the adjusted model (which automatically drops missing values). 
# The _cc suffix stands for 'complete cases'

# Create complete cases dataset
mwanza_rel_cc <- mwanza |> filter(!is.na(rel))

# Unadjusted OR on complete cases
mhor(mwanza_rel_cc, "case", "ed2")

# Mantel-Haenszel analysis stratified by religion
mhor_religion <- mhor(mwanza_rel_cc, "case", "ed2", strata = "rel")

## Test for effect modification with religion

# The test for homogeneity from the mhor() output above tests whether the stratum-specific ORs differ across religion categories.
# 
# Questions:
#   
# - What are the stratum-specific odds ratios?
# ANSWER:

# - What is the religion-adjusted odds ratio (from the model without interaction)?
# ANSWER:
# Stratum-specific odds ratios:
# rel                 OR       95% CI
# -------------------------------------------------- 
# Catholic           2.252  (1.259 - 4.030)
# Muslim             2.022  (0.602 - 6.788)
# Other              2.020  (0.780 - 5.227)
# Protestant         1.394  (0.670 - 2.899)

# - Does the adjusted estimate provide evidence against the null hypothesis of no association?
# ANSWER: Yes, strong evidence, because the OR is 1.914 (95% CI: 1.294 - 2.833)

# - Is there evidence of effect modification between religion and education?
# ANSWER: No, even though the OR are difference between strata
# the test of homogeneity yields P-value = 0.7901
# which is no evidnece against the null hypothesis of homoegeneity among ORs

# - Draft a table summarizing your findings 
# about whether religion confounds the education-HIV relationship
# ANSWER:

## MAKE A TABLE HERE PLEASE ZOE #

# Display as table
mhor_table_religion <- mhor_religion |> 
  as_tibble() |> 
  select(mhor_or)

# Question 5: Dealing with Missing Values for a Potential Confounder #

# Replace "9" with NA using na_if()
mwanza <- mwanza |> 
  mutate(npa = na_if(npa, 9))

# Check the recoding
mwanza |> count(npa)

## Comparing unadjusted ORs with and without restricting to complete cases

# Unadjusted OR including all observations (with missing npa)
mhor(mwanza, "case", "ed2")
# using all of the cases, even people without values for npa ###!!

# Unadjusted OR restricted to complete cases
mwanza_npa_cc <- mwanza |> filter(!is.na(npa))
mhor(mwanza_npa_cc, "case", "ed2")
# not using all of the cases, removing people without values for npa ###!!

## Adjusted for number of partners

# First convert npa to factor with labels
mwanza <- mwanza |> 
  mutate(npa_f = factor(npa,
                        levels = c(1, 2, 3, 4),
                        labels = c("0-1", "2-4", "5-9", "10-19")))

# Remake complete cases with factor
mwanza_npa_cc <- mwanza |> filter(!is.na(npa))

# Mantel-Haenszel analysis stratified by npa
mhor(mwanza_npa_cc, "case", "ed2", strata = "npa_f") #this removes people missing npa values

# Questions to consider:
#   
# - Compare the unadjusted ORs with and without missing npa values.
# Are they different?
#ANSWER: unadjusted with na values: Odds Ratio: 2.416, 95% CI: 1.678 - 3.478
# unadjusted with without na values: Odds Ratio: 2.311, 95% CI: 1.597 - 3.345

# - Compare the unadjusted and adjusted ORs (both on complete cases). 
# ANSWER: unadjusted: OR: 2.311 (95% CI: 1.597 - 3.345)
# adjusted: M-H OR: 2.417 (95% CI: 1.640 - 3.561)

# Does npa confound the relationship?
#ANSEWR: there is a small amount of confounding by number of sexual partners (npa)

#  - If we had compared the unadjusted OR (including missing) to the 
#adjusted OR (excluding missing), what might we have incorrectly concluded?
#ANSWER: the unadjusted OR with the missing npa values could have overestimated (?)
# the relationship between the exposure and the outcome
# which the adjusted with the excluded values would more accurately represent the association

# Question 6 (Optional): Exploring a Dose-Response Relationship# 

# Test for a dose-response effect of number of sexual partners on HIV infection:
# 
# Create a new variable `npa2` with midpoint scores for each category:

# Create npa2 with midpoint scores
mwanza <- mwanza |> 
  mutate(npa2 = case_when(
    npa == 1 ~ 0,   # Midpoint of 0-1
    npa == 2 ~ 3,   # Midpoint of 2-4
    npa == 3 ~ 7,   # Midpoint of 5-9
    npa == 4 ~ 15   # Midpoint of 10-19
  ))

# Check the recoding
mwanza |> count(npa, npa2)

## Odds ratios treating partners as categorical

# Calculate ORs for each partner group vs reference using tabodds
tabodds(mwanza, "case", "npa2")

## Test for departure from linear trend

# We test whether the categorical approach (3 df) 
# fits the data better than the linear trend approach (1 df). 

# To test for departure from linear trend, we need to compare:
# Chi-squared for categorical association (3 df)
# vs Chi-squared for linear trend (1 df, from tabodds output above)

# Chi-squared test for categorical association
npa_table <- table(mwanza$case, mwanza$npa_f)
cat_chi2 <- chisq.test(npa_table)$statistic
cat_chi2

# Get trend chi-squared from tabodds
trend_result <- tabodds(mwanza, "case", "npa2")
trend_chi2 <- trend_result$trend_chi2

# Departure test
departure_chi2 <- cat_chi2 - trend_chi2
departure_df <- 2  # 3 df categorical - 1 df trend
departure_p <- 1 - pchisq(departure_chi2, departure_df)

cat("Test for departure from linear trend:\n")
cat(sprintf("Chi-squared(%d) = %.3f\n", departure_df, departure_chi2))
cat(sprintf("P-value = %.4f\n", departure_p))

# The p-value 0.59, so there is insufficient evidence to conclude the relationship 
# between the npa2 and log-odds of HIV departs from linearity.
# meaning that the linear association is best representative of the dose-dependent relationship

# # Key R Functions Used
# 
# - `mhor()`: Mantel-Haenszel odds ratios and test for homogeneity (custom function)
# - `tabodds()`: Odds ratios for ordered exposures with score test for trend (custom function) 
# - `na_if()`: Convert specific values to missing
# - `chisq.test()`: Chi-squared test for association
# - `fisher.test()`: Fisher's exact test for small samples




