
library(tidyverse)
library(readxl)
library(plotly)
library(minpack.lm)
library(wesanderson)
library(ggthemes)

## Load data

d <- read_excel("./data/data_Ingrid/final_data.xlsx") %>%
    rename(t = `Time (min`, logN = Log)

##

p <- d %>%
    # group_by(exp) %>%
    # mutate(logN0 = mean(ifelse(t == 0, logN, NA), na.rm = TRUE),
    #        logS = logN - logN0) %>%
    # filter(strain == "Ec16") %>%
    ggplot(aes(x = t, y = logN)) +
    geom_point(aes(colour = factor(exp))) +
    geom_line(aes(colour = factor(exp))) +
    geom_smooth(method = "lm") +
    facet_wrap(media ~ strain)

p

# ggplotly(p)

# ## Figure 1
# 
# d %>%
#     ggplot(aes(x = t, y = logN, colour = strain)) +
#     geom_point() +
#     geom_smooth(method = "lm", se = FALSE) +
#     facet_wrap("media")

## Model fitting - Metselaar

Delta <- 4

get_residuals <- function(p, this_data) {
    
    p <- as.list(p)
    
    this_data %>%
        mutate(pred = p$logN0 - Delta*(t/Delta/p$D)^p$p,
               res = pred - logN) %>%
        pull(res)
    
}

my_models <- d %>%
    mutate(strain.media = paste(strain, media, sep = "|")) %>%
    split(.$strain.media) %>%
    map(.,
        ~ nls.lm(c(logN0 = 7, D = 6, p = 1),
                 lower = NULL,
                 upper = NULL,
                 get_residuals,
                 this_data = .)
        )

my_models %>% map(summary)    

## Model fitting - Bigelow

get_residuals <- function(p, this_data) {
    
    p <- as.list(p)
    
    this_data %>%
        mutate(pred = p$logN0 - (t/p$D),
               res = pred - logN) %>%
        pull(res)
    
}

my_models <- d %>%
    mutate(strain.media = paste(strain, media, sep = "|")) %>%
    split(.$strain.media) %>%
    map(.,
        ~ nls.lm(c(logN0 = 7, D = 6),
                 lower = NULL,
                 upper = NULL,
                 get_residuals,
                 this_data = .)
    )

my_models %>% map(summary)    


## Figure 1

d %>%
    ggplot(aes(x = t, y = logN, colour = strain)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap("media")


preds <- my_models %>% 
    map(coef) %>%
    map(.,
        ~ tibble(t = c(0, 20),
                 logN = .[["logN0"]] - t/.[["D"]])
        ) %>%
    imap_dfr(., ~ mutate(.x, strain.media = .y)) %>%
    separate(strain.media, into = c("strain", "media")) 

my_palette <- wes_palette("Darjeeling1", 2)

preds %>%
    mutate(strain = ifelse(strain == "Ec16", "DSM498 (K12)", "DSM1103 (ATCC 25922)")) %>%
    ggplot(aes(x = t, y = logN, colour = media)) +
    geom_line() +
    geom_point(data = d %>% 
                   mutate(strain = ifelse(strain == "Ec16", "DSM498 (K12)", "DSM1103 (ATCC 25922)"))
               ) +
    facet_wrap("strain") +
    theme_clean(base_size = 14) +
    scale_colour_manual(values = my_palette) +
    theme(legend.position = "none") +
    xlab("Treatment time (min)") +
    ylab("Microbial concentration (log CFU/ml)")









