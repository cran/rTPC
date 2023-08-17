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
library(ggrepel)

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


## ----fit_multiple_models, warning=FALSE---------------------------------------
# fit five chosen model formulations in rTPC
d_fits <- nest(d, data = c(temp, rate)) %>%
  mutate(boatman = map(data, ~nls_multstart(rate~boatman_2017(temp = temp, rmax, tmin, tmax, a,b),
                        data = .x,
                        iter = c(4,4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
         gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                        data = .x,
                        iter = c(4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
         oneill = map(data, ~nls_multstart(rate~oneill_1972(temp = temp, rmax, ctmax, topt, q10),
                        data = .x,
                        iter = c(4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
         quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                        data = .x,
                        iter = c(4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 0.5,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 0.5,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
         rezende = map(data, ~nls_multstart(rate~rezende_2019(temp = temp, q10, a,b,c),
                        data = .x,
                        iter = c(4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
         sharpeschoolhigh = map(data, ~nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                        data = .x,
                        iter = c(4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') - 10,
                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') + 10,
                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)))

## ----preds_and_plot, fig.width = 6, fig.height = 4----------------------------
# stack models
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', boatman:sharpeschoolhigh)

# get predictions using augment
newdata <- tibble(temp = seq(min(d$temp), max(d$temp), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# take a random point from each model for labelling
d_labs <- filter(d_preds, temp < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
ggplot(d_preds, aes(temp, .fitted)) +
  geom_line(aes(col = model_name)) +
  geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
  geom_point(aes(temp, rate), d) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2)


## ----calculate_fit_measures---------------------------------------------------
d_ic <- d_stack %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, sigma, AIC, AICc, BIC, df.residual)

d_ic


## ----model_selection, fig.width=6, fig.height = 4-----------------------------
# filter for best model
best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)
best_model

# get colour code
col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]

# plot
ggplot(d_preds, aes(temp, .fitted)) +
  geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
  geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
  geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
  geom_point(aes(temp, rate), d) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures',
       subtitle= 'The Sharpe-Schoolfield model is the best model') +
  geom_hline(aes(yintercept = 0), linetype = 2) 


## ----model_weights, fig.width=6, fig.height = 4-------------------------------
# get model weights
# filtering on AIC score is hashtagged out in this example
d_ic <- d_ic %>%
  # filter(d_ic, aic - min(aic) <= 2) %>%
  mutate(., weight = MuMIn::Weights(AICc))

select(d_ic, model_name, weight) %>%
  arrange(., desc(weight))

# calculate average prediction
ave_preds <- left_join(d_preds, select(d_ic, model_name, weight)) %>%
  group_by(temp) %>%
  summarise(., .fitted = sum(.fitted*weight)) %>%
  ungroup() %>%
  mutate(model_name = 'model average')

# create label for averaged predictions
d_labs <- filter(ave_preds, temp < 30) %>% sample_n(., 1)

# plot these
ggplot(d_preds, aes(temp, .fitted)) +
  geom_line(aes(col = model_name), alpha = 0.3) +
  geom_line(data = ave_preds, col = 'blue') +
  geom_label_repel(aes(label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = d_labs, col = 'blue') +
  geom_point(aes(temp, rate), d) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures',
       subtitle= 'Model averaged predictions') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2)



## ----ave_calc_params----------------------------------------------------------
# calculate estimated parameters
params <- d_stack %>%
  mutate(., params = map(fit, calc_params)) %>%
  select(-fit) %>%
  unnest(params)

# get averaged parameters based on model weights
ave_params <- left_join(params, select(d_ic, model_name, weight)) %>%
  summarise(., across(rmax:skewness, function(x){sum(x*.$weight)})) %>%
  mutate(model_name = 'model average')

# show them
bind_rows(select(params, model_name, rmax:skewness), ave_params) %>%
  mutate_if(is.numeric, round, 2)


