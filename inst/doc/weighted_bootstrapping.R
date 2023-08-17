## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  #tidy.opts=list(width.cutoff=60),
  #tidy=TRUE,
  fig.align = 'center',
  warning=FALSE
)

## ----setup, message=FALSE-----------------------------------------------------
# load packages
library(boot)
library(car) # to install development version of car install.packages("car", repos = c("https://r-forge.r-project.org"), dep = FALSE)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)


## ----get_data, fig.height=4, fig.width = 6------------------------------------
# get curve data
data("chlorella_tpc")

d_ave <- filter(chlorella_tpc, process == 'adaptation', growth_temp == 33, flux == 'photosynthesis') %>%
  group_by(temp, flux) %>%
  summarise(., sd = sd(rate),
            ave_rate = mean(rate),
            groups = 'drop')

d_fit <- nest(d_ave, data = c(temp, ave_rate, sd)) %>%
  mutate(weighted = map(data, ~nls_multstart(ave_rate~lactin2_1995(temp = temp, a, b, tmax, delta_t),
                        data = .x,
                        iter = c(3,3,3,3),
                        start_lower = get_start_vals(.x$temp, .x$ave_rate, model_name = 'lactin2_1995') - 10,
                        start_upper = get_start_vals(.x$temp, .x$ave_rate, model_name = 'lactin2_1995') + 10,
                        lower = get_lower_lims(.x$temp, .x$ave_rate, model_name = 'lactin2_1995'),
                        upper = get_upper_lims(.x$temp, .x$ave_rate, model_name = 'lactin2_1995'),
                        supp_errors = 'Y',
                        convergence_count = FALSE,
                        # include weights here!
                        modelweights = 1/sd)))

# get predictions using augment
newdata <- tibble(temp = seq(min(d_ave$temp), max(d_ave$temp), length.out = 100))
d_preds <- d_fit %>%
  mutate(., preds = map(weighted, augment, newdata = newdata)) %>%
  select(-weighted) %>%
  unnest(preds)

# plot
ggplot() +
  geom_line(aes(temp, .fitted), d_preds) +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), d_ave) +
  geom_point(aes(temp, ave_rate), d_ave, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  ylim(c(-0.25, 3.5))

## ----case_boot----------------------------------------------------------------
# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(ave_rate~lactin2_1995(temp = temp, a, b, tmax, delta_t),
                        data = d_ave,
                        start = coef(d_fit$weighted[[1]]),
                        lower = get_lower_lims(d_ave$temp, d_ave$ave_rate, model_name = 'lactin2_1995'),
                        upper = get_upper_lims(d_ave$temp, d_ave$ave_rate, model_name = 'lactin2_1995'),
                        weights = 1/sd)

# perform case bootstrap
boot1 <- Boot(fit_nlsLM, method = 'case')

## ----bootstrap1, fig.height=3,fig.width=8-------------------------------------
# predict over new data
boot1_preds <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d_ave$temp), max(d_ave$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = lactin2_1995(temp, a, b, tmax, delta_t))

# calculate bootstrapped confidence intervals
boot1_conf_preds <- group_by(boot1_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975),
            .groups = 'drop')

# plot bootstrapped CIs
p1 <- ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'black') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds, fill = 'black', alpha = 0.3) +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), d_ave) +
  geom_point(aes(temp, ave_rate), d_ave, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  ylim(c(-3, 3.5))

# plot bootstrapped predictions
p2 <- ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'black') +
  geom_line(aes(temp, pred, group = iter), boot1_preds, col = 'black', alpha = 0.007) +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), d_ave) +
  geom_point(aes(temp, ave_rate), d_ave, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  ylim(c(-3, 3.5))

p1 + p2


## ----bootstrap2, fig.height=3,fig.width=8-------------------------------------
# perform residual bootstrap
boot2 <- Boot(fit_nlsLM, method = 'residual')

# predict over new data
boot2_preds <- boot2$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d_ave$temp), max(d_ave$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = lactin2_1995(temp, a, b, tmax, delta_t))

# calculate bootstrapped confidence intervals
boot2_conf_preds <- group_by(boot2_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975),
            .groups = 'drop')

# plot bootstrapped CIs
p1 <- ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'black') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds, fill = 'black', alpha = 0.3) +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), d_ave) +
  geom_point(aes(temp, ave_rate), d_ave, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  ylim(c(-3, 3.5))

# plot bootstrapped predictions
p2 <- ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'black') +
  geom_line(aes(temp, pred, group = iter), boot2_preds, col = 'black', alpha = 0.007) +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), d_ave) +
  geom_point(aes(temp, ave_rate), d_ave, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  ylim(c(-3, 3.5))

p1 + p2


## ----conf_int, error=TRUE, fig.height = 5, fig.width = 7----------------------
# get parameters of fitted model
param <- broom::tidy(fit_nlsLM) %>%
  select(param = term, estimate)

# calculate confidence intervals of models
ci1 <- nlstools::confint2(fit_nlsLM, method = 'asymptotic') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'asymptotic')
ci2 <- nlstools::confint2(fit_nlsLM, method = 'profile')
# profiling method fails
ci2 <- mutate(ci1, method = 'profile',
                    conf_lower = NA,
                    conf_upper = NA)

# CIs from case resampling
ci3 <- confint(boot1, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

# CIs from residual resampling
ci4 <- confint(boot2, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci <- bind_rows(ci1, ci2, ci3, ci4) %>%
  full_join(., param)

ggplot(ci, aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
  labs(title = 'Calculation of confidence intervals for model parameters',
       subtitle = 'For the chlorella TPC; profile method failes')


