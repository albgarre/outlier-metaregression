
library(tidyverse)
library(readxl)
library(RColorBrewer)
library(broom)
library(rstan)


# Importing the data

full_data <- read_excel("data/e.coli_ina.xlsx", sheet = "coli1_2", na = "NA") %>%
    rename(D_val = "D (min)", temp = "T") %>%
    mutate(logD = log10(D_val)) %>%
    unite("full_ref", c("Ref", "Year")) %>%
    mutate(T0 = temp - 60) %>%
    mutate(omit = ifelse(omit == "*Included", "Included", omit))

## Datasets

d_refsout <- full_data %>% filter(!is.na(omit))  # References with some points removed
d_final <- full_data %>% filter(omit %in% c(NA, "Included"))

d_out <- full_data %>%  # points that were removed
    filter(!is.na(omit)) %>%
    filter(omit != "*Included")

d_1 <- d_out %>%  # Data after the 1st round
    filter(!(full_ref %in% c("Herceg et al._2011", "Lang et al._2017", "Read et al._1961"
    ))) %>%
    # ggplot() + geom_point(aes(temp, logD, colour = full_ref, shape = omit))
    bind_rows(., d_final)

d_2 <- d_out %>%  # Data after the 2nd round
    filter(full_ref == "Evans et al._1970") %>%
    bind_rows(., d_final)

ggplot(d_final) +
    geom_point(aes(temp, logD))

## Figure 1

point_cols <- brewer.pal(9, "Set1")
line_cols <- brewer.pal(4, "Blues")

d_refsout_plot <- d_refsout %>%
    mutate(full_ref = gsub("_", ", ", full_ref)) %>%
    mutate(
        full_ref = ifelse(
            grepl("E.Malin", full_ref),
            "Malinowska-Panczyketal et al., 2019",
            full_ref
            )
        )

p <- full_data %>% 
    filter(is.na(omit)) %>%
    ggplot() +
    geom_point(aes(x = temp, y = logD),
               shape = 1, size = 4) +
    geom_point(aes(x = temp, y = logD, colour = full_ref, shape = omit),
               data = d_refsout_plot, size = 4) +
    scale_shape_manual(values = c(16, 15, 17, 18, 4)) +
    scale_color_manual(values = point_cols)

p <- p + 
    geom_smooth(aes(x = temp, y = logD), method = "lm",
                se = FALSE,
                data = full_data,
                colour = line_cols[1],
                # color = "black",
                linetype = 1,
                size = 1,
                ) +
    geom_smooth(aes(x = temp, y = logD), method = "lm",
                se = FALSE,
                data = d_1,
                colour = line_cols[2],
                # color = "black",
                linetype = 1,
                size = 1.5,
                ) +
    geom_smooth(aes(x = temp, y = logD), method = "lm",
                se = FALSE,
                data = d_2,
                colour = line_cols[3],
                # color = "black",
                linetype = 3,
                size = 2
                ) +
    geom_smooth(aes(x = temp, y = logD), method = "lm",
                se = FALSE,
                data = d_final,
                colour = line_cols[4],
                # color = "black",
                linetype = 2,
                size = 2.5
                ) +
    xlab("Temperature (ºC)") +
    ylab("Logarithm of the D-value (log min)") +
    theme_gray(base_size = 14) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 14),
          panel.grid = element_line(colour = "gray80"),
          panel.background = element_rect(fill = "gray85")
          # panel.border = element_rect(colour = "green")
          ) 

p

ggsave(p, filename = "Figure1.png",
       width = 12, height = 6)

## Table 1

list(
    full_data,
    d_1, 
    d_2,
    d_final
) %>%
    map(., 
        ~ nls(logD ~ logDref - (temp - 60)/z,
              start = list(logDref = 0, z = 5),
              data = .)
        ) %>%
    map(tidy)
 
## Rapid experimental bias - Figure 5

d_final$temp %>% range()

nls(logD ~ logDref - (temp - 60)/z,
    start = list(logDref = 0, z = 5), 
    data = d_final)

n_steps <- 5
step_size <- 2

p <- (0:n_steps) %>%
    set_names(., (0:n_steps)*step_size) %>%
    map(.,
        ~ list(model_left = filter(d_final, temp > (50 + .*step_size)),
               model_right = filter(d_final, temp < (70 - .*step_size))
               )
        ) %>%
    # map(., 
    #     ~ list(model_left = nls(logD ~ logDref - (temp - 60)/z,
    #                             start = list(logDref = 0, z = 5), 
    #                             data = .$model_left) %>% tidy,
    #            model_right = nls(logD ~ logDref - (temp - 60)/z,
    #                              start = list(logDref = 0, z = 5), 
    #                              data = .$model_right) %>% tidy
    #            )
    #     ) %>%
    map(.,
        ~ list(model_left = lm(logD ~ temp, data = .$model_left) %>% tidy,
               model_right = lm(logD ~ temp, data = .$model_right) %>% tidy)
        ) %>%
    map(.,
        ~ imap_dfr(., ~ mutate(.x, model = .y))
        ) %>%
    imap_dfr(., ~ mutate(.x, step = as.numeric(.y))) %>%
    # filter(term == "z") %>%
    filter(term == "temp") %>%
    ggplot(aes(x = step, y = estimate, colour = model)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error),
                  width = .5, size = 1) +
    xlab("Size of the shift (ºC)") + 
    ylab("Estimated slope (1/ºC)") +
    theme_bw(base_size = 14) +
    theme(legend.position = "none")

ggsave(p, filename = "Figure5.png",
       width = 9, height = 6)
    
# ggplot(d_final) +
#     geom_point(aes(temp, logD))

## Application of the heuristic algorithm to define U and L --------------------

#-  Step 1: Fixing L to the minimum value in the dataset

set.seed(91791)
x_ref <- 60

this_L <- d_final %>%
    pull(logD) %>%
    min()

U_models <- seq(.6, 2, length = 8) %>%
    set_names(., .) %>%
    map(., 
        ~ filter(d_final, 
                 logD > this_L,
                 logD < .)
    ) %>%
    imap(.,
         ~ stan(file = "truncated_model.stan",
                data = list(temperature = .x$temp,
                            logD = .x$logD,
                            refTemp = x_ref,
                            N = nrow(.x),
                            L = this_L, U = as.numeric(.y)),
                iter = 2000, chains = 1, cores = 1
         )
         
    )

U_models %>%
    map(as.data.frame) %>%
    map(.,
        ~ select(., -lp__)
    ) %>%
    map(.,
        ~ pivot_longer(., everything(), names_to = "var")
    ) %>%
    map(.,
        ~ group_by(., var)
    ) %>%
    map(.,
        ~ summarize(.,
                    m_value = median(value),
                    q10 = quantile(value, probs = .1),
                    q90 = quantile(value, probs = .9)
        )
    ) %>%
    imap_dfr(., ~ mutate(.x, U = as.numeric(.y))) %>%
    ggplot(aes(x = U, y = m_value)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = q10, ymax = q90)) +
    facet_wrap("var", scales = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.title = element_blank(),
          legend.position = "top") +
    xlab("Value of U (log min)") + ylab("Parameter estimate") 

#-  Step 2: Fixing U to the value identified in the previous step --------------

set.seed(91791)

my_U <- 1.5

L_models <- seq(-1.3, -0, length = 8) %>%
    set_names(., .) %>%
    map(., 
        ~ filter(d_final, 
                 logD > .,
                 logD < my_U)
    ) %>%
    imap(.,
         ~ stan(file = "truncated_model.stan",
                data = list(temperature = .x$temp,
                            logD = .x$logD,
                            refTemp = x_ref,
                            N = nrow(.x),
                            L = as.numeric(.y), U = my_U),
                iter = 2000, chains = 1, cores = 1
         )
         
    )

L_models %>%
    map(as.data.frame) %>%
    map(.,
        ~ select(., -lp__)
    ) %>%
    map(.,
        ~ pivot_longer(., everything(), names_to = "var")
    ) %>%
    map(.,
        ~ group_by(., var)
    ) %>%
    map(.,
        ~ summarize(.,
                    m_value = median(value),
                    q10 = quantile(value, probs = .1),
                    q90 = quantile(value, probs = .9)
        )
    ) %>%
    imap_dfr(., ~ mutate(.x, L = as.numeric(.y))) %>%
    ggplot(aes(x = L, y = m_value)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = q10, ymax = q90)) +
    facet_wrap("var", scales = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.title = element_blank(),
          legend.position = "top") +
    xlab("Value of L (log min)") + ylab("Parameter estimate") 

#-  Step 3 (last): Fixing L to the value identified above ----------------------

set.seed(91791)

my_L <- -1.3

U_models_2 <- seq(.8, 1.5, length = 8) %>%
    set_names(., .) %>%
    map(., 
        ~ filter(d_final, 
                 logD > my_L,
                 logD < .)
    ) %>%
    imap(.,
         ~ stan(file = "truncated_model.stan",
                data = list(temperature = .x$temp,
                            logD = .x$logD,
                            refTemp = x_ref,
                            N = nrow(.x),
                            L = my_L, U = as.numeric(.y)),
                iter = 2000, chains = 1, cores = 1
         )
         
    )

U_models_2 %>%
    map(as.data.frame) %>%
    map(.,
        ~ select(., -lp__)
    ) %>%
    map(.,
        ~ pivot_longer(., everything(), names_to = "var")
    ) %>%
    map(.,
        ~ group_by(., var)
    ) %>%
    map(.,
        ~ summarize(.,
                    m_value = median(value),
                    q10 = quantile(value, probs = .1),
                    q90 = quantile(value, probs = .9)
        )
    ) %>%
    imap_dfr(., ~ mutate(.x, U = as.numeric(.y))) %>%
    ggplot(aes(x = U, y = m_value)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = q10, ymax = q90)) +
    facet_wrap("var", scales = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.title = element_blank(),
          legend.position = "top") +
    xlab("Value of U (log min)") + ylab("Parameter estimate") 


#- The algorithm seems to have converged, so we move on

## Analysis of the fitted model ------------------------------------------------

my_U <- 1

trunc_model <- U_models_2[[3]]  # The 3rd model is the best one according to our criteria
trunc_model

set.seed(12412)

lin_model <- stan(file = "lineal_inactivation.stan",
                  data = list(temperature = d_final$temp,
                              logD = d_final$logD,
                              refTemp = x_ref,
                              N = nrow(d_final)
                  ),
                  iter = 2000
)

## Comparison of the posterior distributions

list(Truncated = trunc_model,
     Classical = lin_model) %>%
    map(as.data.frame) %>%
    imap_dfr(., ~ mutate(.x, model = .y)) %>%
    select(-lp__) %>%
    gather(par, value, -model) %>%
    ggplot() +
    geom_density(aes(value, colour = model)) +
    facet_wrap(~par, scales = "free")


list(Truncated = trunc_model,
     Classical = lin_model) %>%
    map(as.data.frame) %>%
    imap_dfr(., ~ mutate(.x, model = .y)) %>%
    select(-lp__) %>%
    gather(par, value, -model) %>%
    ggplot() +
    geom_boxplot(aes(x = par, y = value, colour = model)) +
    facet_wrap("par", scales = "free") +
    xlab("") + ylab("") +
    theme(legend.title = element_blank())


## Table 2

set.seed(12412)

list(Truncated = trunc_model,
     Classical = lin_model) %>%
    map(as.data.frame) %>%
    imap_dfr(., ~ mutate(.x, model = .y)) %>%
    select(-lp__) %>%
    # mutate(logD130 = logDref - (130-x_ref)/z) %>%
    gather(par, value, -model) %>%
    group_by(par, model) %>%
    summarize(mean(value), sd(value))

## Figure 6

set.seed(12412)

p <- list(`Truncated model` = trunc_model,
     `Classical model` = lin_model) %>%
    map(as.data.frame) %>%
    imap_dfr(., ~ mutate(.x, model = .y)) %>%
    select(-lp__) %>%
    mutate(sim = row_number()) %>%
    split(.$sim) %>%
    map_dfr(.,
            ~ tibble(temp = seq(48, 70, length = 101),
                     logD = .$logDref - (temp - x_ref)/.$z,
                     sigma = .$sigma,
                     model = .$model)
    ) %>%
    group_by(model, temp) %>%
    summarize(m_logD = median(logD),
              q10 = quantile(logD, .1),
              q90 = quantile(logD, .9),
              sigma = mean(sigma)) %>% 
    ggplot(aes(x = temp, colour = model, fill = model)) +
    # geom_point(aes(x = temp, y = logD), inherit.aes = FALSE,
    #            data = my_data, shape = 1) +
    geom_line(aes(y = q10), linetype = 2, size = 1) +
    geom_line(aes(y = q90), linetype = 2, size = 1) +
    geom_line(aes(y = m_logD-1.96*sigma), linetype = 3, size = 1) +
    geom_line(aes(y = m_logD+1.96*sigma), linetype = 3, size = 1) +
    geom_line(aes(y = m_logD), linetype = 1, size = 1) +
    geom_hline(yintercept = c(my_L, my_U), linetype = 2) +
    geom_point(aes(x = temp, y = logD), shape=1, 
               data = mutate(d_final, included = between(logD, my_L, my_U)), 
               inherit.aes = FALSE) +
    # geom_point(aes(x = temp, y = logD, colour = included), 
    #            data = mutate(my_data, included = between(logD, my_L, my_U)), 
    #            inherit.aes = FALSE) +
    theme_bw() +
    xlab("Temperature (ºC)") + ylab("Logarithm of the D-value (log min)") +
    theme(legend.position = "none",
          legend.title = element_blank(),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14)) 


ggsave(p, filename = "Figure6.png",
       width = 9, height = 6)

###

# d_final %>%
#     group_by(full_ref) %>%
#     mutate(n = n()) %>%
#     filter(n > 3) %>%
#     split(.$full_ref) %>%
#     map(.,
#         ~ lm(logD ~ temp, data = .) 
#         ) %>%
#     map(.,
#         ~ coef(.)[[2]]
#         ) %>%
#     imap_dfr(., ~ tibble(full_ref = .y, z = -1/.x)) 
    










