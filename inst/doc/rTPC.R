## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  #tidy.opts=list(width.cutoff=60),
  #tidy=TRUE,
  fig.align = 'center'
)

## ----setup, message=FALSE-----------------------------------------------------
# load packages
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)

## ----get_mod_names------------------------------------------------------------
get_model_names()

## ----first_plot, fig.width=6, fig.height = 4----------------------------------
# load in data
data("chlorella_tpc")

# keep just a single curve
d <- filter(chlorella_tpc, curve_id == 1)

# show the data
ggplot(d, aes(temp, rate)) +
  geom_point() +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures')


## -----------------------------------------------------------------------------
# choose model
mod = 'sharpschoolhigh_1981'

# get start vals
start_vals <- get_start_vals(d$temp, d$rate, model_name = 'sharpeschoolhigh_1981')

# get limits
low_lims <- get_lower_lims(d$temp, d$rate, model_name = 'sharpeschoolhigh_1981')
upper_lims <- get_upper_lims(d$temp, d$rate, model_name = 'sharpeschoolhigh_1981')

start_vals
low_lims
upper_lims

## ----fit_model----------------------------------------------------------------
# fit model
fit <- nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                                     data = d,
                                                     iter = 500,
                                                     start_lower = start_vals - 10,
                                                     start_upper = start_vals + 10,
                                                     lower = low_lims,
                                                     upper = upper_lims,
                                                     supp_errors = 'Y')

fit

## -----------------------------------------------------------------------------
# calculate additional traits
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)

## ----pred_and_plot, fig.width=6, fig.height = 4-------------------------------
# predict new data
new_data <- data.frame(temp = seq(min(d$temp), max(d$temp), 0.5))
preds <- augment(fit, newdata = new_data)

# plot data and model fit
ggplot(d, aes(temp, rate)) +
  geom_point() +
  geom_line(aes(temp, .fitted), preds, col = 'blue') +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures')


