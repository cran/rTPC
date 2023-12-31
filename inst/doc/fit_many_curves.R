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

## ----fit_multiple_models, warning=FALSE---------------------------------------
# load in data
data("chlorella_tpc")
d <- chlorella_tpc

# when scaling up our code to fit hundreds of models, its nice to have a progress bar
# edit nls_multstart to allow for a progress bar
nls_multstart_progress <- function(formula, data = parent.frame(), iter, start_lower, 
                                   start_upper, supp_errors = c("Y", "N"), convergence_count = 100, 
                                   control, modelweights, ...){
  if(!is.null(pb)){
    pb$tick()
  }
  nls_multstart(formula = formula, data = data, iter = iter, start_lower = start_lower, 
                start_upper = start_upper, supp_errors = supp_errors, convergence_count = convergence_count, 
                control = control, modelweights = modelweights, ...)
}

# start progress bar and estimate time it will take
number_of_models <- 2
number_of_curves <- length(unique(d$curve_id))

# setup progress bar
pb <- progress::progress_bar$new(total = number_of_curves*number_of_models,
                                 clear = FALSE,
                                 format ="[:bar] :percent :elapsedfull")

# fit two chosen model formulation in rTPC
d_fits <- nest(d, data = c(temp, rate)) %>%
  mutate(gaussian = map(data, ~nls_multstart_progress(rate~gaussian_1987(temp = temp, rmax,topt,a),
                        data = .x,
                        iter = c(3,3,3),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
         sharpeschoolhigh = map(data, ~nls_multstart_progress(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                        data = .x,
                        iter = c(3,3,3,3),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)))

## ----preds--------------------------------------------------------------------
# create new list column of for high resolution data
d_preds <- mutate(d_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100)))) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(gaussian,sharpeschoolhigh)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(curve_id, growth_temp, process, flux, model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

glimpse(d_preds)


## ----plot, fig.height = 12, fig.width=8---------------------------------------
# plot
ggplot(d_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate), d) +
  facet_wrap(~curve_id, scales = 'free_y', ncol = 6) +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'All fitted thermal performance curves',
       subtitle = 'gaussian in green; sharpeschoolhigh in orange')


## ----params-------------------------------------------------------------------
# stack models and calculate extra params
d_params <- pivot_longer(d_fits, names_to = 'model_name', values_to = 'fit', c(gaussian,sharpeschoolhigh)) %>%
  mutate(params = map(fit, calc_params)) %>%
  select(curve_id, growth_temp, process, flux, model_name, params) %>%
  unnest(params)

glimpse(d_params)

