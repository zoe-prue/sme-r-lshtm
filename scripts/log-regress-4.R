# Logistic regression 4 practical 12
# 30/1/2025
# PURPOSE

###################################
# library calls

library(haven)
library(tidyverse)
library(gtsummary)
library(here)
library(lmtest)
library(emmeans)

###################################
# global variables

# Source custom functions
source(here("scripts/mh_functions_updated.R"))

## Data import, exploration, management
mortality <- read_dta(here("raw-data/mortality.dta"))
mortality <- mortality |> mutate(across(where(is.labelled), as_factor))

###################################
# script

glimpse(mortality)
summary(mortality)
#View(mortality)

# 1a. Explore the data. 
# What type of variable (in R) are "systolic" and "died"? 
# To which type of variable (in the data sense) 
# does this type of variable (in R) correspond?

# Check variable type
class(mortality$systolic) # numeric
class(mortality$died) # numeric

mortality <- mortality |>
  mutate(died = factor(
    died,
    levels = c(0, 1),
    labels = c("Alive", "Dead")
  )) # need to turn them into factors to be like categorical ordered variables

# Plot systolic BP against death. What does the plot show?

# Scatterplot
mortality |> 
  ggplot(aes(x = systolic, y = died)) +
  geom_point() +
  theme_bw() +
  labs(title = "Systolic blood pressure and death",
       x = "Systolic blood pressure",
       y = "Death")

# Box-and-whiskers plot
mortality |> 
  ggplot(aes(x = died, y = systolic)) +
  geom_boxplot() +
  labs(title = "Systolic blood pressure and death",
       y = "Systolic blood pressure",
       x = "Death") +
  theme_bw() +
  coord_flip()

# ANSWER: the boxplot shows that the median systolic bp
# in the group that survived there is a lowr systolic bp
# compared to the group that died
# the range of data is similar, but with fewer outliers in the died group
# but also a wider interquartile range in the died group

# 1b. Group systolic into three levels and add labels.
# - <120: normal
# - 120-139: pre-hypertension
# - 140+: hypertension

# Could also do this with case_when!
mortality <- mortality |>
  mutate(systolic_grp = cut(
    systolic,
    breaks = c(0, 119, 139, +Inf),
    labels = c("Normal", "Pre-hypertension", "Hypertension")
  ))

# Check it worked
table(mortality$systolic, mortality$systolic_grp)

# Explore the relation between systolic_grp and died.

mortality |> 
  tbl_cross(row = systolic_grp, col = died, percent = "row")

# simple table with percentages

# higher perfcentage of normal in the alive group compared to in the dead group
# in the normal group
# but this pattern shifts in the hypertensive groups and pre-hypertensive groups
# shifting towards more people being in the died group vs the alive

# 1c + 1d. Calculate the OR of death in each sBP group compared to baseline, 
# and test for linear trend. The tabodds() function calculates ORs for each level 
# versus the reference group and provides a score test for trend.

# OR and test for trend
tabodds(mortality, outcome = "died", exposure = "systolic_grp")

# this is the crude OR? for each group?
# increased odds of death outcome in pre-hypertensive (OR = 1.323) 
# and hypertensive group (OR = 3.633)
# with a very low p-value for a chi-squared test for trend in odds
# this is a part of the m-h function

# 1e. Fit a logistic regression model for death, 
# with exposure being grouped sBP treated as a continuous variable.

# Recode systolic_grp as numeric
mortality <- mortality |>
  mutate(systolic_quant = as.numeric(systolic_grp))

# Check
table(mortality$systolic_grp, mortality$systolic_quant)

# Logistic regression with linear relation - log odds scale
model_cont <- glm(died ~ systolic_quant,
                  data = mortality,
                  family = binomial()) 
# OR scale
tbl_regression(model_cont, exponentiate = TRUE)

# What is the value of the slope of this regression line? 
# What is the meaning of this slope? 
# How does its value compare to your answer above? 
# What is the interpretation of the estimate 
# of log(odds) for systolic blood pressure-grouped?

# ANSWER: the slope of the regression line is 1.86
# each increase in "unit" or the category of the systolic bp
# normal to pre-hypertensive to hypertensive groupings
# results in an 86% increase in odds of death
# compared to the test for linear trend above using the mh_function
# the difference is assumed to be the same increase between groups

# 1g. Fit a model with grouped sBP assuming factors 
# and perform a LRT on the following hypotheses:
# - H0: the effect of grouped sBP on log(odds) of death is linear (continuous regression)
# - H1: the relationship is not linear (more complicated) (categorical regression)

model_cat <- glm(died ~ systolic_grp, # this is the categorical model, comparing to the other continuous model from before
                 data = mortality,
                 family = binomial())

tbl_regression(model_cat, exponentiate = TRUE)

# LRT
lrtest(model_cont, model_cat)

# 1h. Summarise your findings
# p = 0.05506 
# meaning there is some evidence to suggest 
# that the non-linear model is a better fit
# than the linear model, suggesting that a more complex model is better
# worded differently: 
# the p-value gives some evidence to reject the null hypothesis 
# that the linear regression model is a better fit
# than the more complex non-linear model

# 1i. Recode the sBP groups so that it's a continuous variable, 
# with each group containing a mid-point value in mmHg. 
# Notice that we need to specify the data type to tell R 
# that it's a continuous variable using as.numeric()

# Recode the variable with (arbitrary) midpoints
mortality <- mortality |>
  mutate(systolic_grp_mid = as.numeric(
    case_when(
      systolic_grp == "Normal" ~ 110,
      systolic_grp == "Pre-hypertension" ~ 130,
      systolic_grp == "Hypertension" ~ 150)
  ))

# Check
table(mortality$systolic_grp, mortality$systolic_grp_mid)

# Now fit the equivalent model to 1f. What changes?

glm(died ~ systolic_grp_mid,
    data = mortality,
    family = binomial()) |>
  tbl_regression(exponentiate = TRUE)

# the answer is showing that for every singluar sBP increase
# the odds of death increase by 3%
# the odds increase by a much smaller amount
# because the change in the x-var is much smaller

# 2a. Consider a model with:
# - outcome: death
# - exposure: visual impairment
# - confounder: age (agegrp)
# 
# Does our estimate of effect of visual impairment differ 
# if we treat age as continuous rather than categorical? 
# Fit a model using each specification of the variable.

# Transform agegrp into a continuous variable
mortality <- mortality |>
  mutate(agegrp_cont = as.numeric(agegrp))

table(mortality$agegrp, mortality$agegrp_cont)

# Continuous model
model_cont <- glm(died ~ vimp + agegrp_cont,
                  data = mortality,
                  family = binomial())
tbl_regression(model_cont, exponentiate = TRUE)

# Categorical model
model_cat <- glm(died ~ vimp + agegrp,
                 data = mortality,
                 family = binomial())
tbl_regression(model_cat, exponentiate = TRUE)

# LRT
lrtest(model_cont, model_cat)

# ANSWER: The adjusted odds ratio for vimp is 2.197. 
# When agegrp was specified as a categorical variable, 
# the adjusted odds ratio for vimp is 2.202. 
# The crude odds ratio is 5.57. 
# In this example there is a similar effect of vimp on 
# log odds of death when treating age group as a 
# categorical variable or assuming a linear trend for age group.

# 2b. Perform an LRT for the interaction 
# between vimp and agegrp treated as continuous. 
# Is there evidence that the association of visual impairment 
# on death differs by age group, where age group is treated as continuous?

# Logistic regression with interaction
model_int <- glm(died ~ vimp * agegrp_cont,
                 data = mortality,
                 family = binomial())
tbl_regression(model_int, exponentiate = TRUE)

# Model without interaction for LRT
model_noint <- glm(died ~ vimp + agegrp_cont,
                   data = mortality,
                   family = binomial())

# LRT
lrtest(model_noint, model_int)

# ANSWER: the LRT is showing that there is not evidence
# the model with interaction is better than the one without

# Calculate stratum-specific ORs for visual impairment 
# for the four age group strata using the continuous model
# The at argument specifies which values of agegrp_cont to evaluate 
# (so here, all unique values) and it must be passed as a list.

emmeans(model_int,
        revpairwise ~ vimp | agegrp_cont,
        at = list(agegrp_cont = unique(mortality$agegrp_cont)),
        type = "response") |>
  pluck("contrasts")

# 3. Using the Mwanza dataset, 
# investigate the association between HIV infection 
# and number of injections in the past year.

# Import the dataset
mwanza <- read_dta(here("raw-data/mwanza.dta"))

View(mwanza)
glimpse(mwanza)
summary(mwanza)

# The variable "inj" is coded as follows:
# 1: does not inject
# 2-4: increasing number of injections per year
# 9: missing value

mwanza <- mwanza |>
  mutate(hiv = as.factor(case_when(case == 1 ~ "HIV+",
                                   case == 0 ~ "HIV-")))

# Replace inj "9" with "NA"
mwanza$inj <- na_if(mwanza$inj, 9)
table(mwanza$inj)
mwanza |> 
  tbl_cross(row = hiv, col = inj)

# Logistic regression including the zero group (coded as 1...)
glm(hiv ~ inj,
    data = mwanza,
    family = binomial()) |>
  tbl_regression(exponentiate = TRUE)

# Again, the odds shown should not be interpreted literally 
# as this is a case control study, but the hypothesis test results are valid. 
# There is strong evidence that the odds of HIV infection 
# increase with increasing numbers of injections. (OR = 1.26)

# Remember that the zero group should be excluded 
# to confirm that a trend with increasing numbers 
# is not simply induced by a difference between 
# the zero category and the rest.

# Logistic regression excluding the zero group
mwanza |>
  filter(inj != 1) |>
  glm(hiv ~ inj,
      data = _,
      family = binomial()) |>
  tbl_regression(exponentiate = TRUE)

# exclude the zero group to mroe accurately 
# represent the dose-dependent relationship

# Although the z statistic (Wald test) is somewhat reduced 
# after excluding the zero group, there is  still strong evidence 
# of a trend for increasing odds of HIV infection with 
# increasing number of injections in the past year





