# meta-analysis helper functions

# load required packages
if (!require("pacman")) install.packages("pacman");
library(pacman) 
p_load(broom, dplyr, metafor, purrr, tibble, tidyr)



#collate sub-analyses into a table of results
collate_models <- function(model, digits = 3, p_digits = 3, row = 1, mv = T){
  fmt <- paste0("%.", digits, "f")
  
  out <- tidy(model, conf.int = TRUE)[row, ]
  
  if(mv){
    out <- out |>
      mutate(
        sigma_b = model$sigma2[1],
        sigma_w = model$sigma2[2])
  } else {
    out <- out |>
      mutate(
        tau_2 = model$tau2)
  }
  
  out <- out |>
    mutate(
      q = model$QE,
      q_p = metafor::fmtp(model$QEp, p_digits),
      across(where(is.numeric), \(x) sprintf(fmt = fmt, x)),
      ci = paste0("[", conf.low, ", ", conf.high, "]"),
      k = model$k,
      q_df = ifelse(mv, model$QEdf, NA),
    )
  
  if(mv){
    out |> 
      select(k, sigma_b, sigma_w, q, q_df, q_p, estimate, ci, statistic, p.value)
  } else {
    out |> 
      select(k, tau_2, q, q_p, estimate, ci, statistic, p.value)
  }
}

## run subanalyses for each level of a categorical predictor and collate into a nice tabulated form

get_mas <- function(tibble, predictor, es_id, study_id, digits = 2, p_digits = 3){
  mas <- tibble  |>
    dplyr::arrange({{predictor}}) |> 
    dplyr::group_by({{predictor}})  |>  
    tidyr::nest()  |> 
    dplyr::mutate(
      model = purrr::map(.x = data, 
                         .f = \(es_tib) rma.mv(yi = g, V = v_g, random = ~1|{{study_id}}/{{es_id}}, data = es_tib)),
      coefs = purrr::map(.x = model, .f = \(m) collate_models(m, digits = digits, p_digits = p_digits))
    )
  
  mas
}

mas_coefs <- function(mas, predictor){
  mas  |>
    dplyr::select(c({{predictor}}, coefs)) |> 
    tidyr::unnest(coefs)
}


## run moderation subanalyses for each level of a categorical predictor and collate into a nice tabulated form
# note row = 2 in collate models because we want to extract the moderation effect

get_mod_mas <- function(tibble, predictor, moderator, es_id, study_id, digits = 2, p_digits = 3){
  mas <- tibble  |>
    dplyr::arrange({{predictor}}) |> 
    dplyr::group_by({{predictor}})  |>  
    tidyr::nest()  |> 
    dplyr::mutate(
      model = purrr::map(.x = data, 
                         .f = \(es_tib) rma.mv(yi = g, V = v_g, mods = ~{{moderator}}, random = ~1|{{study_id}}/{{es_id}}, data = es_tib)),
      coefs = purrr::map(.x = model, .f = \(m) collate_models(m, row = 2, digits = digits, p_digits = p_digits))
    )
  
  mas
}



# report residual heterogeneity stats

report_het <- function(metaobject, digits = 2, p_digits = 3){
  paste0("$Q_E$(", metaobject$QEdf, ") = ", metafor::fmtx(metaobject$QE, digits), ", *p* ", metafor::fmtp(metaobject$QEp, p_digits, sep = T, equal = T))
}


# get row of tidy output from metafor

get_row <- function(tidyobject, row = 2, digits = 2, p_digits = 3){
  n <- nrow(tidyobject)
  
  tidyobject |>
    dplyr::mutate(
      row_number = 1:n,
      p.value = metafor::fmtp(p.value, p_digits),
      ci = paste0("[", metafor::fmtx(conf.low, digits), ", ", metafor::fmtx(conf.high, digits), "]"),
      across(.cols = where(is.double), \(x) metafor::fmtx(x, digits))
    ) |>
    dplyr::filter(row_number == row)
}


#report moderator omnibus stat
report_mod <- function(metaobject, digits = 2, p_digits = 3){
  
  if(is.null(metaobject$robumethod)){
    paste0("$Q_E$(", metaobject$QMdf[1], ") = ", metafor::fmtx(metaobject$QM, digits), ", *p* ", metafor::fmtp(metaobject$QMp, p_digits, sep = T, equal = T))
  } else {
    paste0("$F$(", metaobject$QMdf[1], ", ", metaobject$QMdf[2], ") = ", metafor::fmtx(metaobject$QM, digits), ", *p* ", metafor::fmtp(metaobject$QMp, p_digits, sep = T, equal = T))
  }
}


report_pars <- function(metaobject, row = 1, digits = 2, p_digits = 3, robust = T, mod = F){
  tidyobject <- tidy(metaobject, conf.int = T)
  row <- get_row(tidyobject, row, digits, p_digits)
  
  symbol <- ifelse(mod == TRUE, 
                 paste0("$\\hat{\\beta}$ = "),
                 paste0("$\\hat{\\theta}$ = "))
  
  stat <- ifelse(robust == TRUE, 
                 paste0(", *t*(", metaobject$dfs, ") = "),
                 paste0(", *z* = "))

  
  paste0(symbol, row$estimate, " ", row$ci, stat, row$statistic, ", *p* ", metafor::fmtp(row$p.value, p_digits, sep = T, equal = T))
}

report_par_tbl <- function(metaobject, digits = 2, p_digits = 3, pred_labels){
  fmt <- paste0("%.", digits, "f")
  
  broom::tidy(metaobject, conf.int = T) |>
    mutate(
      p = metafor::fmtp(p.value, 3),
      across(where(is.numeric), \(x) sprintf(fmt = fmt, x)),
      ci = paste0("[", conf.low, ", ", conf.high, "]"),
      predictor = pred_labels
    ) |> 
    dplyr::select(predictor, estimate, ci, statistic, p)
}


# collate regtest info

tidy_regtest <- function(regtest, pb = F){
  tibble(
    estimate = ifelse(pb, regtest$b, regtest$est),
    conf.low = regtest$ci.lb,
    conf.high = regtest$ci.ub,
    p.value = regtest$pval,
    statistic = regtest$zval
  )
}

collate_regtests <- function(model, digits = 3, p_digits = 3){
  fmt <- paste0("%.", digits, "f")
  
  model |>
    tidy_regtest() |> 
    mutate(
      across(where(is.numeric), \(x) sprintf(fmt = fmt, x)),
      ci = paste0("[", conf.low, ", ", conf.high, "]")
    ) |> 
      select(estimate, ci, statistic, p.value)
}

# function to compute full sample SD using O'Neill (2014) result 1 (variance decomposition)
pooled_var <- function(nx, ny, sdx, sdy, ux, uy, sd = F){
  n <- nx + ny
  uxy <- (nx*ux + ny*uy)/n
  pooled_var <- ((nx-1)*sdx^2 + (ny-1)*sdy^2 + (nx*ny*(ux-uy)^2)/n)/(n-1)
  if(sd){
    sqrt(pooled_var)
  } else {
    pooled_var
  }
}

# function to apply correction to d to make it unbiased (g*)
d_to_g <- function(n, d){
  (1-3/(4*(n)-9))*d
}

# function to estimate d based on r. When biserial = T the standard formula for 
# converting a point biserial correlation is used, when biserial = F it uses Mathur & VanderWeele (2019) 
# to treat the measures as continuous and delta (unit change required - which should be set to the mean difference between groups
# of the key continuous measure, X) and sx (standard deviation of the key continuous measure, X). X is the continuous
# variable version of the variable on which d is computed for groups. For example, if d compares anxious vs controls, X is a continuous measure of anxiety.

d_from_r <- function(delta, sx, r, biserial = T){
  if(biserial){
    2*r/(sqrt(1-r^2))
  } else {
    delta*r/(sx*sqrt(1-r^2))
  }
}

# function to estimate the sampling variance of d based on r. When biserial = T the standard formula for 
# converting a point biserial correlation is used, when biserial = F it uses Mathur & VanderWeele (2019) 
# to treat the measures as continuous and d must be specified. n is the total sample size. When var = T the sampling variance is returned, when FALSE the stancdard error is returned.
var_d_from_r <- function(n, r, d, biserial = T, var = T){
  if(biserial){
    se <- 2/sqrt((n-1)*(1-r^2))
  } else {
    se <- abs(d)*sqrt(1/(r^2*(n-3)) + 1/(2*(n-1)))
  }
  
  if(var){
    se^2
  } else {
    se
  }
}


### create bubble plots from a regtest object

get_bubble_data <- function(model){
  data <- model |> metafor::regplot()
  predict <- predict(model)
  tibble::tibble(
    xi = data$xi,
    yi = data$yi,
    psize = data$psize,
    pred = predict$pred,
    pil = predict$pi.lb,
    piu = predict$pi.ub,
    cil = predict$ci.lb,
    ciu = predict$ci.ub
  )
}

plot_bubble <- function(model, xlim = NULL, ylim = NULL, xbreaks = NULL, ybreaks = NULL, xlab = "Predictor", ylab = "Effect size", title = NULL, pi = TRUE){
  tibble <- get_bubble_data(model)
  gg <- ggplot2::ggplot(tibble, aes(x = xi, y = yi)) +
    geom_point(aes(size = psize), fill = "grey80", alpha = 0.2) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_line(aes(x = xi, y = pred), colour = "grey10") +
    labs(x = xlab, y = ylab, title = title) +
    theme_minimal() +
    theme(legend.position = "none")
  
  if(pi){
    gg <- gg +
      geom_line(aes(x = xi, y = pil), colour = "grey10", linetype = "longdash", linewidth = 0.25) +
      geom_line(aes(x = xi, y = piu), colour = "grey10", linetype = "longdash", linewidth = 0.25)
  } else {
    gg <- gg +
      geom_line(aes(x = xi, y = cil), colour = "grey10", linetype = "longdash", linewidth = 0.25) +
      geom_line(aes(x = xi, y = ciu), colour = "grey10", linetype = "longdash", linewidth = 0.25)
  }
  
  if(!is.null(xlim)){
    if(!is.null(ylim)){
      gg <- gg + coord_cartesian(xlim = xlim, ylim = ylim)
    } else {
      gg <- gg + coord_cartesian(xlim = xlim)
    }
  } else {
    if(!is.null(ylim)){
      gg <- gg + coord_cartesian(ylim = ylim)
    }
  }
  
  try(
    if(!is.null(xbreaks)){
      gg <- gg + scale_x_continuous(breaks = seq(xlim[1], xlim[2], xbreaks))
    },
    stop("You must set xlim to use xbreaks")
  )
  
  try(
    if(!is.null(ybreaks)){
      gg <- gg + scale_y_continuous(breaks = seq(ylim[1], ylim[2], ybreaks))
    },
    stop("You must set ylim to use ybreaks")
  )
  
  gg
}

## helper functions for publication bias
# Collates models of moderate nad severe publication bias across levels of a factor and tabulates them

collate_pub_bias <- function(model, pb_mod, pb_sev, digits = 3, p_digits = 3){
  fmt <- paste0("%.", digits, "f")
  
  pb_mod <- tidy(pb_mod, conf.int = T) |> 
    select(estimate, conf.low, conf.high, p.value)
  
  pb_sev <- tidy(pb_sev, conf.int = T) |> 
    select(estimate, conf.low, conf.high, p.value)
  
  out <- pb_mod |> 
    bind_rows(pb_sev) |> 
    mutate(theta = model$b[1],
           pb = c("Moderate", "Severe")) |> 
    mutate(
      p.value = metafor::fmtp(p.value, p_digits),
      across(where(is.numeric), \(x) sprintf(fmt = fmt, x)),
      ci = paste0("[", conf.low, ", ", conf.high, "]")
    ) |> 
    select(pb, theta, estimate, ci, p.value)
  
  out
}

# not use of pmap because collate_pub_bias() has multiple inputs

get_pbm <- function(tibble, predictor, digits = 2, p_digits = 3, a, moderate, severe){
  mas <- tibble  |>
    dplyr::arrange({{predictor}}) |> 
    dplyr::group_by({{predictor}})  |>  
    tidyr::nest()  |> 
    dplyr::mutate(
      model = purrr::map(.x = data, 
                         .f = \(es_tib) rma(yi = g, vi = v_g, data = es_tib)),
      pb_mod = purrr::map(.x = model, 
                          .f = \(m) selmodel(m, type = "stepfun", steps = a, delta = moderate)),
      pb_sev = purrr::map(.x = model, 
                          .f = \(m) selmodel(m, type = "stepfun", steps = a, delta = severe)),
      coefs = purrr::pmap(.l = list(model, pb_mod, pb_sev),
                          .f = \(model, pb_mod, pb_sev) collate_pub_bias(model, pb_mod, pb_sev, digits = digits, p_digits = p_digits)))
  
  mas
}
