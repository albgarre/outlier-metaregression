
## Load libraries

library(tidyverse)
library(readxl)
library(broom)
library(plotly)
library(ellipse)

library(ggrepel)
library(RColorBrewer)

source("functions-meta.R")

# # Importing the data
# 
# full_data <- read_excel("data/e.coli_ina.xlsx", sheet = "coli1_2") %>%
#     rename(D_val = "D (min)", temp = "T") %>%
#     mutate(logD = log10(D_val)) %>%
#     unite("full_ref", c("Ref", "Year")) %>%
#     mutate(T0 = temp - 60)
# 
# my_data <- full_data
# 
# p <- ggplot(my_data) +
#     geom_point(aes(x = T0, y = logD, colour = full_ref))
# 
# ggplotly(p)

# Importing the data

full_data <- read_excel("data/e.coli_ina.xlsx", sheet = "coli1_2", na = "NA") %>%
    rename(D_val = "D (min)", temp = "T") %>%
    mutate(logD = log10(D_val)) %>%
    unite("full_ref", c("Ref", "Year")) %>%
    mutate(T0 = temp - 60)

## Datasets

d_refsout <- full_data %>% filter(!is.na(omit))  # References with some points removed
d_final <- full_data %>% filter(omit %in% c(NA, "*Included"))

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

## Fit the big model

my_data <- full_data

big_model <- fit_big(my_data)

## Fit leaving one ref out

models_ref_out <- fit_one_out(my_data, "full_ref")

## Calculation of DFBETAs

get_DFBETAs(big_model, models_ref_out)
p1 <- get_DFBETAs(big_model, models_ref_out) %>%
    ggplot(aes(x = DFBETA_int, y = DFBETA_T0)) +
    geom_point() +
    geom_label_repel(aes(label = ref), max.overlaps = 100) +
    xlab("Intercept's DFBETA (log min)") + ylab("Slope's DFBETA (1/ºC)") +
    theme_gray(base_size = 14)

## Calculation of RESRATIO

p2 <- get_RESRATIO(big_model, models_ref_out) %>%
    arrange(RESRATIO) %>%
    mutate(ref = factor(ref, levels = ref)) %>%
    ggplot(aes(x = ref, y = RESRATIO)) +
    geom_col(aes(fill = test)) +
    geom_point(aes(y = F_crit), colour = "gold", shape = 18, size = 5) +
    coord_flip() + 
    xlab("") +
    scale_fill_manual(values = c("gray30", "red4")) +
    theme_gray(base_size = 14) +
    theme(legend.position = "none",
          axis.text = element_text(size = 14))

## Figure 2

p <- cowplot::plot_grid(p1, p2, labels = "AUTO")

ggsave(p, filename = "Figure2.png",
       width = 15, height = 6)

## Removing Lang, Read and Herceg ----------------------------------------------

my_data <- d_1
# my_data <- my_data %>%
#     filter(!(full_ref %in% c("Lang et al._2017", "Herceg et al._2011", "Read et al._1961")))

p <- ggplot(my_data) +
    geom_point(aes(x = temp, y = logD, colour = full_ref))

ggplotly(p)

## Fit the big model

big_model <- fit_big(my_data)
# all_models$`2` <- big_model
summary(big_model)

## Fit leaving one ref out

models_ref_out <- fit_one_out(my_data, "full_ref")

## Calculation of DFBETAs

get_DFBETAs(big_model, models_ref_out)
p1 <- get_DFBETAs(big_model, models_ref_out) %>%
    ggplot(aes(x = DFBETA_int, y = DFBETA_T0)) +
    geom_point() +
    geom_label_repel(aes(label = ref), max.overlaps = 100) +
    xlab("Intercept's DFBETA (log min)") + ylab("Slope's DFBETA (1/ºC)") +
    theme_gray(base_size = 14)

## Calculation of RESRATIO

p2 <- get_RESRATIO(big_model, models_ref_out) %>%
    arrange(RESRATIO) %>%
    mutate(ref = factor(ref, levels = ref)) %>%
    ggplot(aes(x = ref, y = RESRATIO)) +
    geom_col(aes(fill = test)) +
    geom_point(aes(y = F_crit), colour = "gold", shape = 18, size = 5) +
    coord_flip() + 
    xlab("") +
    scale_fill_manual(values = c("gray30", "red4")) +
    theme_gray(base_size = 14) +
    theme(legend.position = "none",
          axis.text = element_text(size = 14))

## Figure 3

p <- cowplot::plot_grid(p1, p2, labels = "AUTO")

ggsave(p, filename = "Figure3.png",
       width = 15, height = 6)

## Removing Trevisani, one point by Malinowska & Pereira --------------------------------

my_data <- d_2
# my_data <- my_data %>%
#     filter(!(full_ref %in% c("Trevisani et al._2014")),
#            logD > -2.5  # A lazy way of removing Malinowska's point
#            )

p <- ggplot(my_data) +
    geom_point(aes(x = temp, y = logD, colour = full_ref))

ggplotly(p)

## Fit the big model

big_model <- fit_big(my_data)
summary(big_model)

## Fit leaving one ref out

models_ref_out <- fit_one_out(my_data, "full_ref")

## Calculation of DFBETAs

get_DFBETAs(big_model, models_ref_out)
p1 <- get_DFBETAs(big_model, models_ref_out) %>%
    ggplot(aes(x = DFBETA_int, y = DFBETA_T0)) +
    geom_point() +
    geom_label_repel(aes(label = ref), max.overlaps = 100) +
    xlab("Intercept's DFBETA (log min)") + ylab("Slope's DFBETA (1/ºC)") +
    theme_gray(base_size = 14)

## Calculation of RESRATIO

p2 <- get_RESRATIO(big_model, models_ref_out) %>%
    arrange(RESRATIO) %>%
    mutate(ref = factor(ref, levels = ref)) %>%
    ggplot(aes(x = ref, y = RESRATIO)) +
    geom_col(aes(fill = test)) +
    geom_point(aes(y = F_crit), colour = "gold", shape = 18, size = 5) +
    coord_flip() + 
    xlab("") +
    scale_fill_manual(values = c("gray30", "red4")) +
    theme_gray(base_size = 14) +
    theme(legend.position = "none",
          axis.text = element_text(size = 14))

## Figure 4

p <- cowplot::plot_grid(p1, p2, labels = "AUTO")

ggsave(p, filename = "Figure4.png",
       width = 15, height = 6)

## Removing Pereira  ----------------------------------------------

my_data <- my_data %>%
    filter(!(full_ref %in% c("Pereira et al._2007")))

p <- ggplot(my_data) +
    geom_point(aes(x = temp, y = logD, colour = full_ref))

ggplotly(p)

## Fit the big model

big_model <- fit_big(my_data)
# all_models$`4` <- big_model
summary(big_model)

## Fit leaving one ref out

models_ref_out <- fit_one_out(my_data, "full_ref")

## Calculation of DFBETAs

get_DFBETAs(big_model, models_ref_out)
get_DFBETAs(big_model, models_ref_out) %>%
    ggplot(aes(x = DFBETA_int, y = DFBETA_T0)) +
    geom_point() +
    geom_label_repel(aes(label = ref), max.overlaps = 100) +
    xlab("DFBETA on the intercept") + ylab("DFBETA on the slope") +
    theme_gray(base_size = 14)

## Calculation of RESRATIO

get_RESRATIO(big_model, models_ref_out) %>%
    arrange(RESRATIO) %>%
    mutate(ref = factor(ref, levels = ref)) %>%
    ggplot(aes(x = ref, y = RESRATIO)) +
    geom_col() +
    geom_point(aes(y = F_crit), colour = "red", shape = 18) +
    coord_flip()

##

ggplot(my_data) +
    geom_point(aes(x = temp, y = logD, colour = full_ref)) 

## Removing Evans  ----------------------------------------------

my_data <- my_data %>%
    filter(logD > -2)

p <- ggplot(my_data) +
    geom_point(aes(x = temp, y = logD, colour = full_ref))

ggplotly(p)

## Fit the big model

big_model <- fit_big(my_data)
# all_models$`5` <- big_model
summary(big_model)

## Fit leaving one ref out

models_ref_out <- fit_one_out(my_data, "full_ref")

## Calculation of DFBETAs

get_DFBETAs(big_model, models_ref_out)
get_DFBETAs(big_model, models_ref_out) %>%
    ggplot(aes(x = DFBETA_int, y = DFBETA_T0)) +
    geom_point() +
    geom_label_repel(aes(label = ref), max.overlaps = 100) +
    xlab("DFBETA on the intercept") + ylab("DFBETA on the slope") +
    theme_gray(base_size = 14)

## Calculation of RESRATIO

get_RESRATIO(big_model, models_ref_out) %>%
    arrange(RESRATIO) %>%
    mutate(ref = factor(ref, levels = ref)) %>%
    ggplot(aes(x = ref, y = RESRATIO)) +
    geom_col() +
    geom_point(aes(y = F_crit), colour = "red", shape = 18) +
    coord_flip()

##

ggplot(my_data) +
    geom_point(aes(x = temp, y = logD, colour = full_ref)) +
    geom_point(aes(x = temp, y = logD),
               shape = 1, size = 6,
               data = tibble(temp = c(55),
                             logD = log10(c(4.71, 5.38, 5.76, 5.90)))
    )







