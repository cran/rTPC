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

## ----get_data, fig.height=4, fig.width = 6------------------------------------
# get curve data
data("chlorella_tpc")

d_ave <- filter(chlorella_tpc, process == 'adaptation', growth_temp == 33, flux == 'photosynthesis') %>%
  group_by(temp, flux) %>%
  summarise(., sd = sd(rate),
            ave_rate = mean(rate)) %>%
  ungroup()

# plot
ggplot() +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), d_ave) +
  geom_point(aes(temp, ave_rate), d_ave, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2)

## ----weighted_regression------------------------------------------------------
# fit kamykowski model using rTPC with and without weights
d_fits <- nest(d_ave, data = c(temp, ave_rate, sd)) %>%
  mutate(standard = map(data, ~nls_multstart(ave_rate~kamykowski_1985(temp = temp, tmin, tmax, a,b,c),
                        data = .x,
                        iter = c(4,4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'kamykowski_1985') - 10,
                        start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'kamykowski_1985') + 10,
                        lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'kamykowski_1985'),
                        upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'kamykowski_1985'),
                        supp_errors = 'Y',
                        convergence_count = FALSE)),
         weighted = map(data, ~nls_multstart(ave_rate~kamykowski_1985(temp = temp, tmin, tmax, a,b,c),
                        data = .x,
                        iter = c(4,4,4,4,4),
                        start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'kamykowski_1985') - 10,
                        start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'kamykowski_1985') + 10,
                        lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'kamykowski_1985'),
                        upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'kamykowski_1985'),
                        supp_errors = 'Y',
                        convergence_count = FALSE,
                        # include weights here!
                        modelweights = 1/sd)))

d_fits

## ----preds_and_plot, fig.width = 6, fig.height = 4----------------------------
# stack models
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', standard:weighted)

# get predictions using augment
newdata <- tibble(temp = seq(min(d_ave$temp), max(d_ave$temp), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# take a random point from each model for labelling
d_labs <- filter(d_preds, temp < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# get topt for each model
d_topt <- mutate(d_stack, topt = map_dbl(fit, get_topt), rmax = map_dbl(fit, get_rmax),
                 topt_text = paste(topt, 'ºC', sep = ' '))

# plot
ggplot(d_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), d_ave) +
  geom_point(aes(temp, ave_rate), d_ave, size = 2, shape = 21, fill = 'green4') +
  geom_label_repel(aes(topt, rmax, label = topt_text, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_topt) +
  geom_point(aes(topt, rmax, col = model_name), size = 4, d_topt) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2) +
  ylim(c(0, 3))

