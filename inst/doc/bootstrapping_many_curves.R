## ----setup, message=FALSE-----------------------------------------------------
# load packages
library(boot)
library(car)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)


## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  #tidy.opts=list(width.cutoff=60),
  #tidy=TRUE,
  fig.align = 'center',
  warning=FALSE
)

## ----load_data, fig.height=3.5, fig.width=9-----------------------------------
# load in data
data("chlorella_tpc")

# keep just a single curve
d <- filter(chlorella_tpc, curve_id <= 3)

# fit
d_fits <- nest(d, data = c(rate, temp)) %>%
  mutate(gaussian = map(data, ~nls_multstart(rate ~ gaussian_1987(temp, rmax, topt, a),
                                             data = .x,
                                             iter = c(3,3,3),
                                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 1,
                                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 1,
                                             lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)))

# create high resolution predictions
d_preds <- mutate(d_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100)))) %>%
  select(., -data) %>%
  mutate(preds = map2(gaussian, new_data, ~augment(.x, newdata = .y))) %>%
  select(curve_id, growth_temp, process, flux, preds) %>%
  unnest(preds)

# show the data
ggplot(d, aes(temp, rate)) +
  geom_point(size = 2) +
  geom_line(aes(temp, .fitted), d_preds) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Metabolic rate across temperatures') +
  facet_wrap(~curve_id)


## ----refit_nlsLM--------------------------------------------------------------
# get coefs
d_fits <- mutate(d_fits, coefs = map(gaussian, coef))

# fit with nlsLM instead
d_fits <- mutate(d_fits, nls_fit = map2(data, coefs, ~nlsLM(rate ~ gaussian_1987(temp, rmax, topt, a),
                                        data = .x,
                                        start = .y,
                                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'))))

head(d_fits)

d_fits$nls_fit[[1]]


## ----boot_error, error=TRUE---------------------------------------------------
# try and bootstrap # THIS BREAKS
d_fits <- mutate(d_fits, bootstrap = map(nls_fit, ~Boot(.x, method = 'residual')))


## ----work_around--------------------------------------------------------------
# create empty list column
d_fits <- mutate(d_fits, bootstrap = list(rep(NA, n())))

# run for loop to bootstrap each refitted model
for(i in 1:nrow(d_fits)){
  temp_data <- d_fits$data[[i]]
  temp_fit <- nlsLM(rate ~ gaussian_1987(temp, rmax, topt, a),
               data = temp_data,
               start = d_fits$coefs[[i]],
               lower = get_lower_lims(temp_data$temp, temp_data$rate, model_name = 'gaussian_1987'),
               upper = get_upper_lims(temp_data$temp, temp_data$rate, model_name = 'gaussian_1987'))
  boot <- Boot(temp_fit, method = 'residual')
  d_fits$bootstrap[[i]] <- boot
  rm(list = c('temp_fit', 'temp_data', 'boot'))
}

d_fits


## ----get_preds_and_plot, fig.height=3.5, fig.width=9--------------------------
# get the raw values of each bootstrap
d_fits <- mutate(d_fits, output_boot = map(bootstrap, function(x) x$t))

# calculate predictions with a gnarly written function
d_fits <- mutate(d_fits, preds = map2(output_boot, data, function(x, y){
  temp <- as.data.frame(x) %>%
    drop_na() %>%
    mutate(iter = 1:n()) %>%
    group_by_all() %>%
    do(data.frame(temp = seq(min(y$temp), max(y$temp), length.out = 100))) %>%
    ungroup() %>%
    mutate(pred = gaussian_1987(temp, rmax, topt, a))
  return(temp)
}))

# select, unnest and calculate 95% CIs of predictions
boot_conf_preds <- select(d_fits, curve_id, preds) %>%
  unnest(preds) %>%
  group_by(curve_id, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975),
            .groups = 'drop')

ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d, size = 2) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Rate') +
  facet_wrap(~curve_id)

## ----get_params_and_plot, fig.height=3.5, fig.width=9-------------------------
# get tidied parameters using broom::tidy
# get confidence intervals of parameters
d_fits <- mutate(d_fits, params = map(nls_fit, broom::tidy),
                 cis = map(bootstrap, function(x){
                   temp <- confint(x, method = 'bca') %>%
                     as.data.frame() %>%
                     rename(conf_lower = 1, conf_upper = 2) %>%
                     rownames_to_column(., var = 'term')
                   return(temp)
                   }))

# join parameter and confidence intervals in the same dataset 
left_join(select(d_fits, curve_id, growth_temp, flux, params) %>% unnest(params),
          select(d_fits, curve_id, growth_temp, flux, cis) %>% unnest(cis)) %>%
  ggplot(., aes(curve_id, estimate)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~term, scales = 'free')


## ----calc_params_and_plot, fig_width = 8, fig_height = 5----------------------
# create empty list column
d_fits <- mutate(d_fits, ci_extra_params = list(rep(NA, n())))

# run for loop to bootstrap extra params from each model
for(i in 1:nrow(d_fits)){
  temp_data <- d_fits$data[[i]]
  temp_fit <- nlsLM(rate ~ gaussian_1987(temp, rmax, topt, a),
               data = temp_data,
               start = d_fits$coefs[[i]],
               lower = get_lower_lims(temp_data$temp, temp_data$rate, model_name = 'gaussian_1987'),
               upper = get_upper_lims(temp_data$temp, temp_data$rate, model_name = 'gaussian_1987'))
  boot <- Boot(temp_fit, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(temp_fit)), R = 20, method = 'case') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param')
  d_fits$ci_extra_params[[i]] <- boot
  rm(list = c('temp_fit', 'temp_data', 'boot'))
}

# calculate extra params for each model and put in long format to begin with
d_fits <- mutate(d_fits, extra_params = map(nls_fit, function(x){calc_params(x) %>% pivot_longer(everything(), names_to =  'param', values_to = 'estimate')}))

left_join(select(d_fits, curve_id, growth_temp, flux, extra_params) %>% unnest(extra_params),
          select(d_fits, curve_id, growth_temp, flux, ci_extra_params) %>% unnest(ci_extra_params)) %>%
  ggplot(., aes(as.character(curve_id), estimate)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  labs(y = 'estimate', x = "curve id") +
  facet_wrap(~param, scales = 'free') +
  labs(title = 'Calculation of confidence intervals for extra parameters')


