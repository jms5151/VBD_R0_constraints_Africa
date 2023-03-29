# Use rTPC package to compare compare multiple metabolic curves for Zika pMI 

# load packages
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(ggrepel)
library(AICcmodavg)

# write function to label ggplot2 panels
label_facets_num <- function(string){
  len <- length(string)
  string = paste('(', 1:len, ') ', string, sep = '')
  return(string)
}

# load data
traits.df <- read.csv('../VBD-Data/aegypti_traits_temp_formatted.csv')

# for each trait, order by temperature from low to high 
traits.df <- traits.df %>% arrange(trait_name_new, Temperature)

d <- data.frame('temp' =  traits.df$Temperature[traits.df$trait_name_new == 'pMI_zikv']
                , 'rate' =traits.df$trait_value_new[traits.df$trait_name_new == 'pMI_zikv'])

# show the data
ggplot(d, aes(temp, rate)) +
  geom_point() +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures')

# fit every model formulation in rTPC
d_fits <- nest(d, data = c(temp, rate)) %>%
  mutate(beta = map(data, ~nls_multstart(rate~beta_2012(temp = temp, a, b, c, d, e),
                                         data = .x,
                                         iter = c(6,6,6,6,6),
                                         start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') - 10,
                                         start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') + 10,
                                         lower = get_lower_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                         upper = get_upper_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                         supp_errors = 'Y',
                                         convergence_count = FALSE)),
         # boatman = map(data, ~nls_multstart(rate~boatman_2017(temp = temp, rmax, tmin, tmax, a,b),
         #                                    data = .x,
         #                                    iter = c(4,4,4,4,4),
         #                                    start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') - 10,
         #                                    start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') + 10,
         #                                    lower = get_lower_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
         #                                    upper = get_upper_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
         #                                    supp_errors = 'Y',
         #                                    convergence_count = FALSE)),
         # briere2 = map(data, ~nls_multstart(rate~briere2_1999(temp = temp, tmin, tmax, a,b),
         #                                    data = .x,
         #                                    iter = c(4,4,4,4),
         #                                    start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') - 10,
         #                                    start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') + 10,
         #                                    lower = get_lower_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
         #                                    upper = get_upper_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
         #                                    supp_errors = 'Y',
         #                                    convergence_count = FALSE)),
         # delong = map(data, ~nls_multstart(rate~delong_2017(temp = temp, c, eb, ef, tm, ehc),
         #                                   data = .x,
         #                                   iter = c(4,4,4,4,4),
         #                                   start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'delong_2017') - 10,
         #                                   start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'delong_2017') + 10,
         #                                   lower = get_lower_lims(.x$temp, .x$rate, model_name = 'delong_2017'),
         #                                   upper = get_upper_lims(.x$temp, .x$rate, model_name = 'delong_2017'),
         #                                   supp_errors = 'Y',
         #                                   convergence_count = FALSE)),
         # flinn = map(data, ~nls_multstart(rate~flinn_1991(temp = temp, a, b, c),
         #                                  data = .x,
         #                                  iter = c(5,5,5),
         #                                  start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') - 10,
         #                                  start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') + 10,
         #                                  lower = get_lower_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
         #                                  upper = get_upper_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
         #                                  supp_errors = 'Y',
         #                                  convergence_count = FALSE)),
         gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                             data = .x,
                                             iter = c(4,4,4),
                                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                             lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
         # hinshelwood = map(data, ~nls_multstart(rate~hinshelwood_1947(temp = temp, a, e, b, eh),
         #                                        data = .x,
         #                                        iter = c(5,5,5,5),
         #                                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'hinshelwood_1947') - 1,
         #                                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'hinshelwood_1947') + 1,
         #                                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'hinshelwood_1947'),
         #                                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'hinshelwood_1947'),
         #                                        supp_errors = 'Y',
         #                                        convergence_count = FALSE)),
         # joehnk = map(data, ~nls_multstart(rate~joehnk_2008(temp = temp, rmax, topt, a, b, c),
         #                                   data = .x,
         #                                   iter = c(4,4,4,4, 4),
         #                                   start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') - 10,
         #                                   start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') + 10,
         #                                   lower = get_lower_lims(.x$temp, .x$rate, model_name = 'joehnk_2008'),
         #                                   upper = get_upper_lims(.x$temp, .x$rate, model_name = 'joehnk_2008'),
         #                                   supp_errors = 'Y',
         #                                   convergence_count = FALSE)),
         # johnson_lewin = map(data, ~suppressWarnings(nls_multstart(rate~ johnsonlewin_1946(temp = temp, r0, e, eh, topt),
         #                                                           data = .x,
         #                                                           iter = c(4,4,4,4),
         #                                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') - 1,
         #                                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') + 1,
         #                                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'johnsonlewin_1946'),
         #                                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'johnsonlewin_1946'),
         #                                                           supp_errors = 'Y',
         #                                                           convergence_count = FALSE))),
         # kamykowski = map(data, ~nls_multstart(rate~kamykowski_1985(temp = temp, tmin, tmax, a, b, c),
         #                                       data = .x,
         #                                       iter = c(4,4,4,4,4),
         #                                       start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') - 10,
         #                                       start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') + 10,
         #                                       lower = get_lower_lims(.x$temp, .x$rate, model_name = 'kamykowski_1985'),
         #                                       upper = get_upper_lims(.x$temp, .x$rate, model_name = 'kamykowski_1985'),
         #                                       supp_errors = 'Y',
         #                                       convergence_count = FALSE)),
         # lactin2 = map(data, ~nls_multstart(rate~lactin2_1995(temp = temp, a, b, tmax, delta_t),
         #                                    data = .x,
         #                                    iter = c(4,4,4,4),
         #                                    start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') - 10,
         #                                    start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') + 10,
         #                                    lower = get_lower_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
         #                                    upper = get_upper_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
         #                                    supp_errors = 'Y',
         #                                    convergence_count = FALSE)),
         # lrf = map(data, ~nls_multstart(rate~lrf_1991(temp = temp, rmax, topt, tmin, tmax),
         #                                data = d,
         #                                iter = c(3,3,3,3),
         #                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lrf_1991') - 10,
         #                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lrf_1991') + 10,
         #                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'lrf_1991'),
         #                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'lrf_1991'),
         #                                supp_errors = 'Y',
         #                                convergence_count = FALSE)),
         # modifiedgaussian = map(data, ~nls_multstart(rate~modifiedgaussian_2006(temp = temp, rmax, topt, a, b),
         #                                             data = .x,
         #                                             iter = c(4,4,4,4),
         #                                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') - 10,
         #                                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') + 10,
         #                                             lower = get_lower_lims(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006'),
         #                                             upper = get_upper_lims(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006'),
         #                                             supp_errors = 'Y',
         #                                             convergence_count = FALSE)),
         oneill = map(data, ~nls_multstart(rate~oneill_1972(temp = temp, rmax, ctmax, topt, q10),
                                           data = .x,
                                           iter = c(4,4,4,4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref, e, eh, topt, tref = 15),
                                          data = .x,
                                          iter = c(4,4,4,4),
                                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'pawar_2018') - 10,
                                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'pawar_2018') + 10,
                                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'pawar_2018'),
                                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'pawar_2018'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
         # quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
         #                                      data = .x,
         #                                      iter = c(4,4,4),
         #                                      start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 0.5,
         #                                      start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 0.5,
         #                                      lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
         #                                      upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
         #                                      supp_errors = 'Y',
         #                                      convergence_count = FALSE)),
         # ratkowsky = map(data, ~nls_multstart(rate~ratkowsky_1983(temp = temp, tmin, tmax, a, b),
         #                                      data = .x,
         #                                      iter = c(4,4,4,4),
         #                                      start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') - 10,
         #                                      start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') + 10,
         #                                      lower = get_lower_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
         #                                      upper = get_upper_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
         #                                      supp_errors = 'Y',
         #                                      convergence_count = FALSE)),
         # rezende = map(data, ~nls_multstart(rate~rezende_2019(temp = temp, q10, a,b,c),
         #                                    data = .x,
         #                                    iter = c(4,4,4,4),
         #                                    start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') - 10,
         #                                    start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') + 10,
         #                                    lower = get_lower_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
         #                                    upper = get_upper_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
         #                                    supp_errors = 'Y',
         #                                    convergence_count = FALSE)),
         # sharpeschoolfull = map(data, ~nls_multstart(rate~sharpeschoolfull_1981(temp = temp, r_tref,e,el,tl,eh,th, tref = 15),
         #                                             data = .x,
         #                                             iter = c(4,4,4,4,4,4),
         #                                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') - 10,
         #                                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') + 10,
         #                                             lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981'),
         #                                             upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981'),
         #                                             supp_errors = 'Y',
         #                                             convergence_count = FALSE)),
         sharpeschoolhigh = map(data, ~nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') - 10,
                                                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') + 10,
                                                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE))#,
         # sharpeschoollow = map(data, ~nls_multstart(rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15),
         #                                            data = .x,
         #                                            iter = c(4,4,4,4),
         #                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') - 10,
         #                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') + 10,
         #                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
         #                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
         #                                            supp_errors = 'Y',
         #                                            convergence_count = FALSE)),
         # spain = map(data, ~nls_multstart(rate~spain_1982(temp = temp, a,b,c,r0),
         #                                  data = .x,
         #                                  iter = c(4,4,4,4),
         #                                  start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') - 1,
         #                                  start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') + 1,
         #                                  lower = get_lower_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
         #                                  upper = get_upper_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
         #                                  supp_errors = 'Y',
         #                                  convergence_count = FALSE)),
         # thomas1 = map(data, ~nls_multstart(rate~thomas_2012(temp = temp, a,b,c,topt),
         #                                    data = .x,
         #                                    iter = c(4,4,4,4),
         #                                    start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2012') - 1,
         #                                    start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2012') + 2,
         #                                    lower = get_lower_lims(.x$temp, .x$rate, model_name = 'thomas_2012'),
         #                                    upper = get_upper_lims(.x$temp, .x$rate, model_name = 'thomas_2012'),
         #                                    supp_errors = 'Y',
         #                                    convergence_count = FALSE)),
         # thomas2 = map(data, ~nls_multstart(rate~thomas_2017(temp = temp, a,b,c,d,e),
         #                                    data = .x,
         #                                    iter = c(3,3,3,3,3),
         #                                    start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') - 10,
         #                                    start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') + 10,
         #                                    lower = get_lower_lims(.x$temp, .x$rate, model_name = 'thomas_2017'),
         #                                    upper = get_upper_lims(.x$temp, .x$rate, model_name = 'thomas_2017'),
         #                                    supp_errors = 'Y',
         #                                    convergence_count = FALSE)),
         # weibull = map(data, ~nls_multstart(rate~weibull_1995(temp = temp, a,topt,b,c),
         #                                    data = .x,
         #                                    iter = c(4,4,4,4),
         #                                    start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') - 10,
         #                                    start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') + 10,
         #                                    lower = get_lower_lims(.x$temp, .x$rate, model_name = 'weibull_1995'),
         #                                    upper = get_upper_lims(.x$temp, .x$rate, model_name = 'weibull_1995'),
         #                                    supp_errors = 'Y',
         #                                    convergence_count = FALSE))
         )
# stack models
# d_stack <- select(d_fits, -data) %>%
#   pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:weibull)

d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:sharpeschoolhigh)


# get parameters using tidy
params <- d_stack %>%
  mutate(., est = map(fit, tidy)) %>%
  select(-fit) %>%
  unnest(est)

# get predictions using augment
newdata <- tibble(temp = seq(min(d$temp), max(d$temp), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# plot all curves individually
ggplot(d_preds, aes(temp, rate)) +
  geom_point(aes(temp, rate), d) +
  geom_line(aes(temp, .fitted), col = 'blue') +
  facet_wrap(~model_name, labeller = labeller(model_name = label_facets_num), scales = 'free', ncol = 5) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Fits of every model available in rTPC') +
  geom_hline(aes(yintercept = 0), linetype = 2)


# take a random point from each model for labelling
d_labs <- filter(d_preds, temp < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot all curves overlaid (too busy with lots of curves, fine with subset)
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

# model selection
d_ic <- d_stack %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, AICcmodavg::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, sigma, AIC, AICc, BIC, df.residual)

# look at results
d_ic

# filter for best model
best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)
best_model

# get colour code
col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]

# plot overlaid curves, color best fit with orange
ggplot(d_preds, aes(temp, .fitted)) +
  geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
  geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
  geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
  geom_point(aes(temp, rate), d) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'pMI',
       title = 'Zika pMI across temperatures',
       subtitle= paste0('The ', best_model, ' is the best model')) +
  geom_hline(aes(yintercept = 0), linetype = 2) 

# get parameters for best_model
xx <- as.data.frame(params)
xx <- xx[xx$model_name == best_model,]

# get model equation
gaussian_1987
