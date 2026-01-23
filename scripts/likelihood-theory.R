# Lieklooh Theory 1 and 2
# 23/1/2026
# # Point of this practical is not for use analyzing data, 
# but rather to understand these concepts and visualize them
# Below is an R Markdown setup chunk to set default parameters for display 
# - for more info on setup chunks, see intro_06

# READ INTO WHY POISSON MATCHES RATE AND BINOMIAL MATCHES RISK

# use 4 functions to explore aexact and approximate likelihoods
# -   `blik()` - binomial likelihood
# -   `plik()` - Poisson likelihood
# -   `bloglik()` - binomial log-likelihood
# -   `ploglik()` - Poisson log-likelihood

# don't use these functions to analyze data
# we are sourcing these functions from other files
# if you are rewriting functions >3 times, you should be sourcing them


knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.width = 8,
  fig.height = 6
)

###################################
# library calls

# Load required packages
library(tidyverse)  
install.packages("patchwork")
library(patchwork)  # For combining multiple plots - might need to install!
library(here)

###################################
# global variables

# Source the likelihood functions
# Adjust these paths to where you saved the function files

# If you saved them in a "functions" subfolder:
source(here("scripts", "blik_function.R"))
source(here("scripts", "bloglik_function.R"))
source(here("scripts", "plik_function.R"))
source(here("scripts", "ploglik_function.R"))

###################################
# script

# likelihood ratio compares how likely a parameter value is
# relative to the maximum likelihood estimate (MLE)
# likelihood ratio of 1.0 means that parameter value is exactly as likely as the MLE
# (i.e. is it the MLE)
# lower values indicate less plausible parameter values
# By default, in these functions, the cutoff is set to 0.1465 (just over 1/8), 
# which corresponds to parameter values within roughly 2 log-likelihood units of the maximum.

# In the plots given by these functions:
# - The shaded blue region shows the supported range
# - The vertical black line marks the MLE ($\hat{\pi}$)
# - The horizontal dashed red line shows the cutoff value
# - The black dots mark the MLE and the confidence limits 
# (where the curve crosses the cutoff)

## Question 1 ##

# plot likelihood for pi when 4events (d = disease) and 6 non-events (h = healthy)
# AND for 400 = d and 600 = h
# samex function to keep scales the same

# Create three plots with same x-axis scale
# Add scale_x_continuous() to ensure all plots use the same x-axis range (0 to 1)
q1a <- blik(d = 4, h = 6) + 
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1))
q1a

q1b <- blik(d = 40, h = 60) + 
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1))
q1b

q1c <- blik(d = 400, h = 600) + 
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1))
q1c

# Combine plots using patchwork
combined_plot <- q1a / q1b / q1c
print(combined_plot)

# Save combined plot
ggsave("q1_combined.png", combined_plot, width = 8, height = 12)

# If you have an 'outputs' folder
ggsave(here("outputs", "q1_combined.png"), combined_plot, width = 8, height = 12)

# notice how the 95% confidence interval gets smaller as the total number of subjects increases
# this means a more e=precise estiamte of the parameter

## Question 2 ##

# Try Question 1 with the cut point 0.2585 instead of 0.1465 (option `cut`). 
# Hint: what is the cutpoint demarcating? 
# Think analogously to confidence intervals...
# premonition: the cutoff point is like the cutoff for a 95% confidence interval?

q2 <- blik(d = 4, h = 6, cut = 0.2585) + 
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) # specified the breaks (range of x-axis - bac this is risk
q2

combined_plot <- q1a / q2
print(combined_plot)
ggsave(here("outputs", "q1_cutoff_0.2585.png"), combined_plot, width = 8, height = 12)

# What happens to the supported range / confidence interval?
# the supported range gets more narrow due to the more stringent cutoff

## Question 3 ##

# use blik() with more extreme splits (1:9 ot 1:99)
# curve is no longer symmetrical and bell-shaped
# we did not specify x-axis - specify what they are!
# ggplot decide what the range of the x-axis is based on that data

q3a <- blik(d = 1, h = 9)
q3a

q3b <- blik(d = 1, h = 99)
q3b

## Question 4 ##

# What happens when the number of events is zero?
# the MLE is zero
# and the supported range only goes into positive values to the right 

q4 <- blik(d = 0, h = 100)
q4

## Question 5 ##

# Use `plik()` to plot the likelihood 
# for the rate λ when 7 events are observed for 500 person years
# Do the same for 70 events and 5000 person years, 
# and for 700 events and 50000 person years. 
# 95% CI should get narrower as the total number of subjects increases

q5a <- plik(d = 7, y = 500)
q5a

q5b <- plik(d = 70, y = 5000)
q5b

q5c <- plik(d = 700, y = 50000)
q5c

## Question 6 ##

# Try using `plik()` with a very low number of events.

q6 <- plik(d = 1, y = 500)
q6

# quite a left skew, due to the MLE best matching a lower value

## Question 7 ##

# 25 participants are asked to choose between treatments A and B; 
# 15 prefer A and 10 prefer B. 
# Use `blik()` to plot the likelihood for π, that is, 
# the probability of preferring A. What is the most likely value of π?

# Use the `null = 0.5` option to obtain the likelihood ratio 
# for the null value π = 0.5 (i.e. equally likely to prefer A or B)
# and use the `pval` option to obtain an approximate P-value 
# for the test of the null hypothesis that π = 0.5. 

q7 <- blik(d = 15, h = 10, null = 0.5, pval = TRUE)
q7

# Note: the green line represents the value of the null hypothesis, 
# and the orange line where the null hypothesis meets the likelihood.

# GOOD EXPLANATION
# What are your conclusions?
# we posit the null hypothesis that ppl are equally likely
# to choose A or B
# the MLE estimate is closer to 0.6
# the null hypothesis intersecting with the likelihood ratio farther down the curve
# is indicative of a parameter not highly supported by the data
# but still within the supported range

## Question 8 ##

# Repeat Question 7 for 40 participants, 
# of whom 24 prefer A and 16 prefer B. 

q8 <- blik(d = 24, h = 16, null = 0.5, pval = TRUE)
q8

# What are your conclusions?

# The MLE is still 0.6
# but because the null hypothesis intersects with the likelihood ratio
# even lower down on the curve
# and the likelihood ratio is a measure of support for the model by the data
# there is even less reason to think that the null hypothesis
# is supported by the data 

## Question 9 ##

# Repeat Question 7 for 60 participants, 
# of whom 36 prefer A and 24 prefer B. 

q9 <- blik(d = 36, h = 24, null = 0.5, pval = TRUE)
q9

# What are your conclusions?
# even more similar to above, there the measure of the support for the null hypothesis
# (the likelihood ratio)
# intersects with the null value for the parameter even farther down the curve
# meaning there is poor evidence from the data to support 
# the null value being the parameter for the model

## Log Likelihoods ##

# Use `bloglik()` to plot the exact and approximate log likelihoods for π 
# (the probability of event) when there are 4 events and 6 non-events. 
# Make sure you understand the meaning of all of the output

q10 <- bloglik(d = 4, h = 6)
q10

# The meaning of all of the output is:
# following the binomial distribution, as this is for risk,
# the blue line is the approximate estimate of the log likelihood ratio
# of the model from the data provided, using the quadratic expression
# the 95% confidence interval, or "supported range"
# which lies between where the blue line (approximate) intersects with
# the cutoff on either side of the curve
# the red line is the exact curve for the log likelihood ratio for the risk parameter

## Question 11 ##

# Use `bloglik()` to plot the exact and approximate log likelihoods 
# for π when 40 events and 60 non-events are observed

q11 <- bloglik(d = 40, h = 60)
q11

# how does the 95% confidence interval (supported range) change?
# it goes from the parameter pi being 0.32 to 0.5
# this provides a better approximation due to having more events

## Question 12 ##

# repeat the above with 400 and 600 events

q12 <- bloglik(d = 400, h = 600)
q12

# the approximation is now almost perfect! 

## Question 13 ##

# Use `bloglik()` to plot the exact and approximate log likelihoods for π 
# when 2 events and 18 non-events are observed. 

q13a <- bloglik(d = 2, h = 18)
q13a

q13b <- bloglik(d = 20, h = 180)
q13b

# What are the 95% confidence limits from the approximate log likelihood?
# for the one with greater number of events,
# the range is 0.6 to 0.15 ish

# the approximation for the 2:18 split is quite poor bc the true log likelihood curve
# is for from quadratic in shape
# the lower limit is negative (not well shown in graph)
# which is an impossible value for the risk

# the approximation is better for more events

## Question 14 ##

# repeat Q13 using the log likelihood plotted against log odds parameter

q14 <- bloglik(d = 2, h = 18, logodds = TRUE)
q14

# what are the 95% confidence limits for the approximate log likelihood?
# look at the console for the parameter value estimates for the original risk scale
# between zero and one

## Question 15 ##

# Use `ploglik()` to plot the log likelihood for λ 
# when 7 events are observed for 500 person-years. 

q15 <- ploglik(d = 7, y = 500)
q15

# What are the 95% confidence limits from this approximation?
# 6 events per person-year to 27 events per person-year (?)

# Repeat Question 15 using the log likelihood plotted 
# against the log rate parameter.

q16 <- ploglik(d = 7, y = 500, lograte = TRUE)
q16

# What are the 95% limits from the approximation?

## Question 17 ##

# Use `ploglik()` to plot the log likelihood for λ 
# when there is 1 event for 1000 person years.

q17a <- ploglik(d = 1, y = 1000)
q17a

# How good is the approximation? 
# What is the 95% confidence interval for λ?

q17b <- ploglik(d = 1, y = 1000, lograte = TRUE)
q17b

# Repeat using the log rate scale. 
# What is the 95% confidence interval for λ?






