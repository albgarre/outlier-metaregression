
## High influence points detection based on Regression diagnostics: identifying influential data and sources of collinearity


#' Meta-regression model
#' 
fit_big <- function(d) {
    lm(logD ~ T0, data = d)
}

#' Meta-regression model leaving one out according to column var_out
#' 
fit_one_out <- function(d, var_out) {
    
    d[[var_out]] %>%
        unique() %>%
        set_names(., .) %>%
        map(.,
            ~ filter(d, .data[[var_out]] != .)
        ) %>%
        map(.,
            ~ lm(logD ~ T0, data = .)
        )
    
}

#' Scatter plot of the intercept versus the slope
#' 
scatter_estimates <- function(model_list, model0) {
    
    model_list %>%
        map(tidy) %>%
        map(.,
            ~ select(., term, estimate)
        ) %>%
        map(.,
            ~ pivot_wider(., names_from = term, values_from = estimate)
        ) %>%
        map2(., model_list,
             ~ bind_cols(.x,
                         tibble(sigma = summary(.y)$sigma)
             )
        ) %>%
        imap_dfr(.,
                 ~ mutate(.x, ref = .y)
        ) %>%
        ggplot(aes(x = `(Intercept)`, y = T0)) +
        # geom_point(aes(size = 1/sigma), shape = 1) +
        geom_label(aes(label = ref, fill = sigma)) +
        geom_point(aes(x = coef(model0)[1], y = coef(model0)[2]),
                   colour = "red", shape = 16, size = 6)
    
}

#' Plot of the confidence ellipsoids
#' 
ellipsoid_plot <- function(model_list, model0) {
    
    model_list %>%
        map(ellipse) %>%
        map(as_tibble) %>%
        imap_dfr(.,
                 ~ mutate(.x, ref = .y)
        ) %>%
        ggplot() +
        geom_path(aes(x = `(Intercept)`, y = T0, colour = ref)) +
        geom_path(aes(x = `(Intercept)`, y = T0),
                  data = as_tibble(ellipse(model0)),
                  size = 2, linetype = 2)
    
}

#' Calculation of DFBETAs
#' 
get_DFBETAs <- function(model0, model_list) {

model_list %>%
    map(coef) %>%
    imap_dfr(.,
             ~ tibble(par = names(.x),
                      value = .x,
                      ref = .y)
    ) %>%
    pivot_wider(names_from = par, values_from = value) %>%
    mutate(DFBETA_int = coef(model0)[1] - `(Intercept)`,
           DFBETA_T0 = coef(model0)[2] - T0) 

}

#' Calculation of scaled DFBETAs
#' 
get_scaled_DFBETAs <- function(model0, model_list) {
    
    model_list %>%
        map(coef) %>%
        imap(.,
             ~ tibble(par = names(.x),
                      value = .x)
        ) %>%
        map(.,
            ~ mutate(., 
                     pars_0 = coef(model0),
                     DFBETA = pars_0 - value
            )
        ) %>%
        map2_dbl(., model_list,
                ~ matrix(.x$DFBETA, nrow = 1) %*% t(model.matrix(.y)) %*% model.matrix(.y) %*% matrix(.x$DFBETA, ncol = 1)
        ) %>%
        sort()
    
}

#' Calculation of DFFIT
#' 
get_DFFIT <- function(model0, model_list, target_T0) {
    
    pred_0 <- predict(model0, newdata = tibble(T0 = target_T0))
    
    model_list %>%
        map(.,
            ~ tibble(pred = predict(., newdata = tibble(T0 = target_T0)))
        ) %>%
        map(.,
            ~ mutate(., DFFIT = pred_0 - pred)
        ) %>%
        imap_dfr(., 
                 ~ mutate(.x, ref = .y)
        ) 
}

#' Calculation of DFFITS
#' 
get_scaled_DFFIT <- function(model0, model_list, target_T0) {
    
    pred_0 <- predict(model0, newdata = tibble(T0 = target_T0))
    
    model_list %>%
        map(.,
            ~ tibble(pred = predict(., newdata = tibble(T0 = target_T0)),
                     sigma = summary(.)$sigma)
        ) %>%
        map(.,
            ~ mutate(., DFFITS = (pred_0 - pred)/sigma/1)
        ) %>%
        imap_dfr(., 
                 ~ mutate(.x, ref = .y)
        )
    
}

#' #' Calculation of RESRATIO
#' #' 
#' get_RESRATIO <- function(model0, model_list, only_top = TRUE) {
#'     
#'     SSR_big <- sum(model0$residuals^2)
#'     n <- length(model0$residuals)
#'     p <- length(model0$coefficients)
#'     
#'     d <- model_list %>%
#'         map(.,
#'             ~ tibble(SSR_old = SSR_big,
#'                      SSR_new = sum(.$residuals^2),
#'                      m = n - length(.$residuals),
#'                      RESRATIO = ((SSR_old - SSR_new)/m)/(SSR_new/(n-p-m)),
#'                      F_crit = qf(.95, m, n-p-m)
#'             )
#'         ) %>%
#'         imap_dfr(.,
#'                  ~ mutate(.x, ref = .y)
#'         ) %>%
#'         arrange(desc(RESRATIO)) %>%
#'         mutate(test = RESRATIO > F_crit)
#'     
#'     if (only_top) {  # Apply the statistical test
#'         
#'         flag <- FALSE
#'         
#'         for (i in 1:nrow(d)) {
#'             
#'             if (!d$test[[i]]) {
#'                 flag <- TRUE
#'             }
#'             
#'             if (flag) {
#'                 d$F_crit[[i]] <- NA
#'                 d$test[[i]] <- FALSE
#'             }
#'             
#'         }
#'         
#'     }
#'     
#'     ## Return 
#'     
#'     d
#'     
#' }

#' Calculation of RESRATIO
#' 
get_RESRATIO <- function(model0, model_list, only_top = TRUE) {
    
    SSR_big <- sum(model0$residuals^2)
    n <- length(model0$residuals)
    p <- length(model0$coefficients)
    
    d <- model_list %>%
        map(.,
            ~ tibble(SSR_old = SSR_big,
                     SSR_new = sum(.$residuals^2),
                     m = n - length(.$residuals),
                     RESRATIO = ((SSR_old - SSR_new)/m)/(SSR_new/(n-p-m)),
                     F_crit = qf(.95, m, n-p-m)
            )
        ) %>%
        imap_dfr(.,
                 ~ mutate(.x, ref = .y)
        ) %>%
        arrange(desc(RESRATIO)) %>%
        mutate(test = RESRATIO > F_crit)
    
    if (only_top) {  # Apply the statistical test
        
        flag <- FALSE
        
        for (i in 1:nrow(d)) {

            if (flag) {
                d$F_crit[[i]] <- NA
                d$test[[i]] <- FALSE
            }
            
            if (!d$test[[i]]) {
                flag <- TRUE
            }
            
        }
        
    }
    
    ## Return 
    
    d
    
}




