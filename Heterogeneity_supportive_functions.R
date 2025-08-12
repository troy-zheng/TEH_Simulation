#Hete Supportive files

########################################
####Baseline dependency##################
########################################

get_cor <- function(x, y, type_x, type_y) {
  ## function to calculate the correlation between variables of mixed
  ## type. In case both x,y are continuous, the standard correlation
  ## is returned. Otherwise an alternative that would reduce to
  ## Pearson's correlation for continuous data
  Pearson_C <- function(X2, N) {
    X2 <- max(X2, 0)
    sqrt(X2 / (N + X2))
  }
  N <- length(x)
  type_x_cont <- type_x == "continuous"
  type_y_cont <- type_y == "continuous"
  if (type_x_cont & type_y_cont) {
    fit <- lm(y ~ x)
    fit0 <- lm(y ~ 1)
    X2 <- as.numeric(2 * (logLik(fit) - logLik(fit0)))
    return(Pearson_C(X2, N)) # should roughly equal cor(x,y)
  }
  if (!type_x_cont & !type_y_cont) { # both categorical
    fit <- nnet::multinom(y ~ x, trace = FALSE)
    fit0 <- nnet::multinom(y ~ 1, trace = FALSE)
    X2 <- as.numeric(-2 * (logLik(fit0) - logLik(fit)))
    return(Pearson_C(X2, N))
  }
  if (type_x_cont & !type_y_cont) {
    fit <- nnet::multinom(y ~ x, trace = FALSE) ## multinomial regression
    fit0 <- nnet::multinom(y ~ 1, trace = FALSE)
    X2 <- as.numeric(-2 * (logLik(fit0) - logLik(fit)))
    P1 <- Pearson_C(X2, N)
    fit <- lm(x ~ y)
    fit0 <- lm(x ~ 1)
    X2 <- as.numeric(-2 * (logLik(fit0) - logLik(fit)))
    P2 <- Pearson_C(X2, N)
    ## cat(sprintf("P1: %.2f, P2: %.2f\n", P1, P2))
    return((P1 + P2) * 0.5)
  }
  if (!type_x_cont & type_y_cont) {
    fit <- nnet::multinom(x ~ y, trace = FALSE) ## multinomial regression
    fit0 <- nnet::multinom(x ~ 1, trace = FALSE)
    X2 <- as.numeric(-2 * (logLik(fit0) - logLik(fit)))
    P1 <- Pearson_C(X2, N)
    fit <- lm(y ~ x)
    fit0 <- lm(y ~ 1)
    X2 <- as.numeric(-2 * (logLik(fit0) - logLik(fit)))
    P2 <- Pearson_C(X2, N)
    ## cat(sprintf("P1: %.2f, P2: %.2f\n", P1, P2))
    return((P1 + P2) * 0.5)
  }
}

get_dep <- function(dat) {
  var_type <- sapply(dat, function(x) {
    cl <- class(x)
    if (cl %in% c("factor", "ordered", "logical", "character")) {
      return("categorical")
    }
    if (cl %in% c("numeric", "integer", "double")) {
      return("continuous")
    }
  })
  n_vars <- ncol(dat)
  out <- matrix(0, nrow = n_vars, ncol = n_vars)
  res <- vector("list", 0.5 * n_vars * (n_vars - 1))
  nams <- names(dat)
  colnames(out) <- rownames(out) <- names(dat)
  z <- 1
  for (i in 1:(n_vars - 1)) {
    for (j in (i + 1):n_vars) {
      tmp <- get_cor(dat[, i, drop = TRUE], dat[, j, drop = TRUE], var_type[i], var_type[j])
      if (inherits(tmp, "try-error")) {
        txt <- sprintf("Problem for variables: %s and %s\n", nams[i], nams[j])
        cat(txt, attr(tmp, "condition")$message, "\n")
        cat("Contingency table for variables looks like this (missing values removed)\n")
        print(table(dat[, i], dat[, j]))
        break
      } else {
        out[i, j] <- out[j, i] <- tmp
        if (var_type[i] == var_type[j]) {
          comparison <- var_type[i]
        } else {
          comparison <- "mixed"
        }
        res[[z]] <- c(nams[i], nams[j], tmp, comparison)
        z <- z + 1
        get_cor <- function(x, y, type_x, type_y) {
          ## function to calculate the correlation between variables of mixed
          ## type. In case both x,y are continuous, the standard correlation
          ## is returned. Otherwise an alternative that would reduce to
          ## Pearson's correlation for continuous data
          Pearson_C <- function(X2, N) {
            X2 <- max(X2, 0)
            sqrt(X2 / (N + X2))
          }
          N <- length(x)
          type_x_cont <- type_x == "continuous"
          type_y_cont <- type_y == "continuous"
          if (type_x_cont & type_y_cont) {
            fit <- lm(y ~ x)
            fit0 <- lm(y ~ 1)
            X2 <- as.numeric(2 * (logLik(fit) - logLik(fit0)))
            return(Pearson_C(X2, N)) # should roughly equal cor(x,y)
          }
          if (!type_x_cont & !type_y_cont) { # both categorical
            fit <- nnet::multinom(y ~ x, trace = FALSE)
            fit0 <- nnet::multinom(y ~ 1, trace = FALSE)
            X2 <- as.numeric(-2 * (logLik(fit0) - logLik(fit)))
            return(Pearson_C(X2, N))
          }
          if (type_x_cont & !type_y_cont) {
            fit <- nnet::multinom(y ~ x, trace = FALSE) ## multinomial regression
            fit0 <- nnet::multinom(y ~ 1, trace = FALSE)
            X2 <- as.numeric(-2 * (logLik(fit0) - logLik(fit)))
            P1 <- Pearson_C(X2, N)
            fit <- lm(x ~ y)
            fit0 <- lm(x ~ 1)
            X2 <- as.numeric(-2 * (logLik(fit0) - logLik(fit)))
            P2 <- Pearson_C(X2, N)
            ## cat(sprintf("P1: %.2f, P2: %.2f\n", P1, P2))
            return((P1 + P2) * 0.5)
          }
          if (!type_x_cont & type_y_cont) {
            fit <- nnet::multinom(x ~ y, trace = FALSE) ## multinomial regression
            fit0 <- nnet::multinom(x ~ 1, trace = FALSE)
            X2 <- as.numeric(-2 * (logLik(fit0) - logLik(fit)))
            P1 <- Pearson_C(X2, N)
            fit <- lm(y ~ x)
            fit0 <- lm(y ~ 1)
            X2 <- as.numeric(-2 * (logLik(fit0) - logLik(fit)))
            P2 <- Pearson_C(X2, N)
            ## cat(sprintf("P1: %.2f, P2: %.2f\n", P1, P2))
            return((P1 + P2) * 0.5)
          }
        }
        
        get_dep <- function(dat) {
          var_type <- sapply(dat, function(x) {
            cl <- class(x)
            if (cl %in% c("factor", "ordered", "logical", "character")) {
              return("categorical")
            }
            if (cl %in% c("numeric", "integer", "double")) {
              return("continuous")
            }
          })
          n_vars <- ncol(dat)
          out <- matrix(0, nrow = n_vars, ncol = n_vars)
          res <- vector("list", 0.5 * n_vars * (n_vars - 1))
          nams <- names(dat)
          colnames(out) <- rownames(out) <- names(dat)
          z <- 1
          for (i in 1:(n_vars - 1)) {
            for (j in (i + 1):n_vars) {
              tmp <- get_cor(dat[, i, drop = TRUE], dat[, j, drop = TRUE], var_type[i], var_type[j])
              if (inherits(tmp, "try-error")) {
                txt <- sprintf("Problem for variables: %s and %s\n", nams[i], nams[j])
                cat(txt, attr(tmp, "condition")$message, "\n")
                cat("Contingency table for variables looks like this (missing values removed)\n")
                print(table(dat[, i], dat[, j]))
                break
              } else {
                out[i, j] <- out[j, i] <- tmp
                if (var_type[i] == var_type[j]) {
                  comparison <- var_type[i]
                } else {
                  comparison <- "mixed"
                }
                res[[z]] <- c(nams[i], nams[j], tmp, comparison)
                z <- z + 1
              }
            }
          }
          res <- as.data.frame(do.call("rbind", res))
          names(res) <- c("Variable 1", "Variable 2", "Correlation", "Comparison")
          res[, 3] <- as.numeric(res[, 3])
          list(cor_mat = out, results = res)
        }
        
        low_freq_categories <- function(dat, perc) {
          var_type <- sapply(dat, function(x) {
            cl <- class(x)
            if (cl %in% c("factor", "ordered", "logical", "character")) {
              return("categorical")
            }
            if (cl %in% c("numeric", "integer", "double")) {
              return("continuous")
            }
          })
          ind <- var_type == "categorical"
          vars <- sapply(dat[, ind], function(x) {
            tab <- table(x)
            freq_tab <- tab / sum(tab)
            min(freq_tab)
          })
          vars <- vars[order(vars, decreasing = FALSE)]
          data.frame(var = names(vars), size_sparsest_cat = as.numeric(vars))
        }
      }
    }
  }
  res <- as.data.frame(do.call("rbind", res))
  names(res) <- c("Variable 1", "Variable 2", "Correlation", "Comparison")
  res[, 3] <- as.numeric(res[, 3])
  list(cor_mat = out, results = res)
}

low_freq_categories <- function(dat, perc) {
  var_type <- sapply(dat, function(x) {
    cl <- class(x)
    if (cl %in% c("factor", "ordered", "logical", "character")) {
      return("categorical")
    }
    if (cl %in% c("numeric", "integer", "double")) {
      return("continuous")
    }
  })
  ind <- var_type == "categorical"
  vars <- sapply(dat[, ind], function(x) {
    tab <- table(x)
    freq_tab <- tab / sum(tab)
    min(freq_tab)
  })
  vars <- vars[order(vars, decreasing = FALSE)]
  data.frame(var = names(vars), size_sparsest_cat = as.numeric(vars))
}

########################################
####################Baseline Plots#######
########################################
plot_density_bar <- function(x.var, data) {
  n_miss <- sum(is.na(data[, x.var]))
  miss_perc <- round(n_miss / dim(data)[1] * 100)
  title <- sprintf("%s, missing: %s (%s%%)", x.var, n_miss, miss_perc)
  data_x <- unlist(data %>% dplyr::select(x = all_of(x.var)) %>% na.omit())
  if (is.numeric(data_x) & length(unique(data_x)) > 5) {
    pp <- data %>%
      dplyr::select(x.var = all_of(x.var)) %>%
      ggplot(mapping = aes(x = (x.var))) +
      geom_histogram(bins = 200) +
      geom_rug()
  } else {
    data_x <- data %>%
      dplyr::select(x = all_of(x.var)) %>%
      na.omit()
    n_miss <- sum(is.na(data[, x.var]))
    range_max <- max(table(data_x))
    pp <- ggplot(data_x, aes(x)) +
      geom_bar() +
      geom_text(stat = "count", aes(label = sprintf("%s (%s %%)", ..count.., round(..count.. / sum(..count..) * 100))), vjust = -0.1) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      ylim(c(0, range_max * 1.1))
  }
  pp +
    theme_bw() +
    xlab(x.var) +
    ylab("") +
    ggtitle(title)
}

### plot boxplot (continuous vars) or barplot (categorical vars) by stratification
plot_boxbar_by <- function(x.var, y.var, data) {
  data_x_y <- data %>%
    dplyr::select(x = all_of(x.var), y = all_of(y.var)) %>%
    na.omit()
  if (unlist(is.numeric(data_x_y$x)) & length(unlist(unique(data_x_y$x))) > 5) { ## continuous
    ## side by side boxplot for x.var using y.var as category
    nns <- table(data_x_y$y)
    pp <- ggplot(data_x_y, aes(x = y, y = x)) +
      stat_boxplot(geom = "errorbar") +
      geom_boxplot() +
      geom_jitter(shape = 16, position = position_jitter(0.2), color = "red", alpha = 0.2) +
      xlab("") +
      ylab(x.var) +
      ggtitle(x.var) +
      theme(strip.text.x = element_text(size = 3))
  } else { ## categorical
    pp <- data_x_y %>%
      count(y, x) %>%
      group_by(y) %>%
      mutate(pct = round(prop.table(n) * 100)) %>%
      mutate(n_p = sprintf("%s (%s%%)", n, pct)) %>%
      ggplot(aes(factor(y), n, fill = factor(x))) +
      geom_bar(position = "stack", stat = "identity", alpha = 0.8) +
      geom_text(aes(label = n_p), position = position_stack(vjust = 0.5)) +
      labs(title = x.var, fill = x.var) +
      xlab(x.var) +
      ylab("count") +
      theme(legend.position = "top")
  }
  pp +
    theme_bw()
}


########################################
####################Helper trt eff#######
########################################
## return a plot, with smoothing curve (y versus x) as well as histogram on bottom (based on x)
## if x categorical, return bar plot with count
#' Plot (smoothed) treatment difference against baseline covariate
#'
#' This function plots the (smoothed) treatment difference (in average
#' outcomes) against a baseline covariate. For numeric baseline
#' covariates a regression spline is used (with user-specified df). For
#' categorical baseline covariates just means are presented.
#'
#' @param data output from p_diffx_cal, include columns x, trt_eff, lower, upper 
#' (lower and upper bound for treatment effect), lev(optional), pred(optional): smoothed phi
#' @param dat_raw include x with lev variable
#' @param x.name Name of baseline variable in data
#' @param y.name Name of treatment variable in data (code assumes
#'   there are only two treatments)
#' @param lev.name facet variable name
#' @param data Data set
#' @param df Degrees of freedom for the regression spline
#' @param family Family to use for outcome variable (passed to glm)
#' @param span Span to use for smoothing the phis
#' @export
#' @examples

comb_plot <- function(data, dat_raw, x.name, y.name, lev.name, 
                      diff_lab = "probability diff"){
  ## defined commonly used theme
  common_theme <- theme_bw() + theme(legend.position = "top")
  
  dat <- data 
  
  if(is.null(lev.name)){
    title_txt <- paste0(diff_lab, " of ", y.name, " on ", x.name)
  }else{
    title_txt <- paste0(diff_lab, " of ", y.name, " on ", x.name, " by ", lev.name)
  }
  
  
  if(is.numeric(dat$x)){
    p_effect <- ggplot(aes(x, trt_eff), data = dat) +
      geom_line() +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
      geom_hline(yintercept = 0, linetype = "longdash") +
      common_theme +
      labs(
        # title = title_txt,
        # subtitle = "Average with 95% CI",
        x = x.name, y = diff_lab)
    if(!is.null(dat$phi_pred)){
      p_effect <- p_effect + geom_line(aes(x, phi_pred), lty = 2, alpha = 0.5)
    }
    
    p_dist <- ggplot(dat_raw, aes(x = x)) +
      geom_histogram() +
      theme_bw() +
      xlab(x.name)
    
    if(!is.null(lev.name)){
      dat$lev <- data[, "lev"]
      p_effect <- p_effect + facet_grid(~lev)
      p_dist <- p_dist + facet_grid(~x2)
    }
    p_plot <- cowplot::plot_grid(p_effect +
                                   theme(
                                     axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(),
                                     axis.title.x = element_blank()
                                   ),
                                 p_dist +
                                   theme(strip.background = element_blank(), strip.text.x = element_blank()),
                                 nrow = 2, rel_heights = c(3, 1), align = "v", axis = "lr"
    )
  }else{
    # browser()
    p_plot <- dat %>%
      mutate( x = factor(x, levels = unique(dat$x))) %>%
      ggplot(aes(x, trt_eff)) +
      geom_point() +
      theme_bw() +
      geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.5, width = 0.1) +
      geom_hline(yintercept = 0, linetype = "longdash") +
      labs(
        # title = title_txt,
        # subtitle = paste0("Average with 95% CI"),
        x = x.name, y = diff_lab)
    
    if(!is.null(dat$phi_pred)){
      p_plot <- p_plot + geom_point(aes(x, phi_pred), pch = 4, alpha = 0.5, data = dat)
    }
    
    if(!is.null(lev.name)){
      dat$lev <- data[, "lev"]
      p_plot <- p_plot + facet_grid(~lev, scale = "free")
    }
  }
  return(p_plot)
}


## functions to help calculate the dataframe used for p_diffx and p_diffxx
#'
#' This function calculate treatment difference 
#' (response difference or loglogs ratio difference) given x and trt
#'
#' @param dat a data frame with variable x, trt, and y
#' @param df Degrees of freedom for the regression spline
#' @param family Family to use for outcome variable (passed to glm)
#' @param phi_predictions Predicted phi's obtained from the DR learner (optional)
#' @param span Span to use for smoothing the phis
#' @export
#' @examples
p_diffx_cal <- function(dat,
                        df,
                        family,
                        span) {
  levels_trt <- levels(dat$trt)
  if (is.numeric(dat$x)) {
    ## code based on glm and contrast
    x_seq <-
      seq(min(dat$x, na.rm = T), max(dat$x, na.rm = T), length.out = 100)
    glm_fit1 <-
      glm(
        y ~ splines::ns(x, df = df),
        data = filter(dat, trt == levels_trt[1]),
        family = family
      )
    glm_fit2 <-
      glm(
        y ~ splines::ns(x, df = df),
        data = filter(dat, trt == levels_trt[2]),
        family = family
      )
    pred_dat <- data.frame(x = x_seq)
    pred1 <-
      predict(glm_fit1, pred_dat, type = "response", se.fit = TRUE)
    pred2 <-
      predict(glm_fit2, pred_dat, type = "response", se.fit = TRUE)
    plot_dat <-
      data.frame(
        x = x_seq,
        trt_eff = pred1$fit - pred2$fit,
        se = sqrt(pred1$se ^ 2 + pred2$se ^ 2)
      ) %>%
      mutate(lower = trt_eff - 2 * se, upper = trt_eff + 2 * se)
  } else {
    n_count <- dat %>%
      mutate( x = factor(x, levels = unique(dat$x))) %>%
      group_by(x) %>%
      summarise(n = n()) %>%
      mutate(x_label = paste0(x, "(n = ", n, ")")) %>%
      ungroup()
    
    dat <- dat %>%
      mutate(x = factor(x, levels = n_count$x, labels = n_count$x_label))
    
    x_lev <- n_count$x_label
    
    glm_fit <- glm(y ~ trt * x, data = dat, family = family)
    preds <- pred_se <- numeric(length(x_lev))
    for (i in 1:length(x_lev)) {
      glm_fit <-
        glm(y ~ trt,
            data = filter(dat, x == x_lev[i]),
            family = family)
      pred1 <- predict(
        glm_fit,
        newdata = data.frame(trt = levels_trt[1], x = x_lev[i]),
        se.fit = TRUE,
        type = "response"
      )
      pred2 <- predict(
        glm_fit,
        newdata = data.frame(trt = levels_trt[2], x = x_lev[i]),
        se.fit = TRUE,
        type = "response"
      )
      preds[i] <- pred1$fit - pred2$fit
      pred_se[i] <- sqrt(pred1$se.fit ^ 2 + pred2$se.fit ^ 2)
    }
    plot_dat <-
      data.frame(
        x = factor(x_lev, levels  = x_lev),
        trt_eff = preds,
        se = pred_se
      ) %>%
      mutate(lower = trt_eff - 2 * se, upper = trt_eff + 2 * se)
  }
  
  if (is.numeric(dat$x) & (!is.null(dat$phi))) {
    ## refit phi using smoothing curve
    fit <- loess(phi ~ x, data = dat, span = span)
    plot_dat$phi_pred <- predict(fit, newdata = plot_dat)
  }else if(!is.numeric(dat$x) & (!is.null(dat$phi))){
    mns <- tapply(dat$phi, dat$x, mean)
    pdat2 <- data.frame(x = names(mns), phi_pred = mns)
    plot_dat <- full_join(pdat2, plot_dat, by = "x")
  }
  
  return(plot_dat)
}

comb_plot_yx <- function(data, x.name, y.name, trt.name, lev.name, family, df = 2){
  ## defined commonly used theme
  common_theme <- theme_bw() + theme(legend.position = "top")
  
  dat <- data 
  
  if(is.null(lev.name)){
    title_txt <- paste0(y.name, " on ", x.name)
  }else{
    title_txt <- paste0(y.name, " on ", x.name, " by ", lev.name)
  }
  
  
  if(is.numeric(dat$x)){
    
    p_effect <- dat %>%
      ggplot(aes(x, y, col = trt)) +
      geom_point() +
      geom_smooth(
        method = "glm", method.args = list(family = family),
        formula = y ~ splines::ns(x, df = df)
      ) +
      common_theme +
      labs(
        # title = title_txt,
        # subtitle = "Average with 95% CI",
        x = x.name, y = y.name, color = trt.name
      )
    p_dist <- ggplot(dat, aes(x = x)) +
      geom_histogram() +
      theme_bw() +
      xlab(x.name) 
    
    if(!is.null(lev.name)){
      dat$lev <- dat[, "lev"]
      p_effect <- p_effect + facet_grid(~lev)
      p_dist <- p_dist + facet_grid(~lev)
    }
    p_plot <- cowplot::plot_grid(p_effect +
                                   theme(
                                     axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(),
                                     axis.title.x = element_blank()
                                   ),
                                 p_dist +
                                   theme(strip.background = element_blank(), strip.text.x = element_blank()),
                                 nrow = 2, rel_heights = c(3, 1), align = "v", axis = "lr"
    )
  }else{
    # browser()
    
    p_plot <- dat %>%
      ggplot(aes(x, y, col = trt, group = trt)) +
      stat_summary(fun = mean, geom = "point") +
      stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.1) +
      common_theme +
      labs(
        # title = title_txt,
        # subtitle = "Average with 95% CI",
        x = x.name, y = y.name, color = trt.name
      )
    
    if(!is.null(lev.name)){
      dat$lev <- data[, "lev"]
      p_plot <- p_plot + facet_grid(~lev, scale = "free")
    }
  }
  return(p_plot)
}



########################################
####################trt eff plots#######
########################################

## functions to generate plot for response by treatment/treatment effect v.s. covariates(one or two variables)
## to explore the relationship

#' Plot (smoothed) average outcome against baseline covariate by
#' treatment group
#'
#' This function plots the (smoothed) average outcome against baseline
#' covariate by treatment group. For numeric baseline covariates a
#' regression spline is used (with user-specified df). For categorical
#' baseline covariates just means are presented.
#'
#' @param y.name Name of outcome variable in data
#' @param x.name Name of baseline variable in data
#' @param trt.name Name of treatment variable in data
#' @param data Data set
#' @param df Degrees of freedom for the regression spline
#' @param family Family to use for outcome variable (passed to glm)
#' @export
#' @examples
p_yx <- function(y.name, x.name, trt.name, data, df = 2, family) {
  ## defined commonly used theme
  common_theme <- theme_bw() + theme(legend.position = "top")
  
  # browser()
  dat <- data %>% dplyr::select(all_of(c(y.name, trt.name, x.name))) %>% na.omit()
  colnames(dat) <- c("y", "trt", "x")
  dat$trt <- as.factor(dat$trt)
  levels_trt <- levels(dat$trt)
  if (length(levels_trt) > 2) {
    stop("code only works for two treatments for now")
  }
  
  if(!is.numeric(dat$x)){
    n_count <- dat %>%
      mutate( x = factor(x, levels = unique(dat$x))) %>%
      group_by(x) %>%
      summarise(n = n()) %>%
      mutate(x_label = paste0(x, "(n = ", n, ")")) %>%
      ungroup()
    
    dat <- dat %>%
      mutate(x = factor(x, levels = n_count$x, labels = n_count$x_label))
  }
  
  comb_plot_yx(data = dat, x.name, y.name, trt.name, lev.name = NULL, family = family, df = df)
}


#' Plot (smoothed) treatment difference against baseline covariate
#'
#' This function plots the (smoothed) treatment difference (in average
#' outcomes) against a baseline covariate. For numeric baseline
#' covariates a regression spline is used (with user-specified df). For
#' categorical baseline covariates just means are presented.
#'
#' @param y.name Name of outcome variable in data
#' @param x.name Name of baseline variable in data
#' @param trt.name Name of treatment variable in data (code assumes
#'   there are only two treatments)
#' @param data Data set
#' @param df Degrees of freedom for the regression spline
#' @param family Family to use for outcome variable (passed to glm)
#' @param facet_title Include a facet (useful for p_trt_effect_xx function)
#' @param label_y Include the label of the y axis (useful for p_trt_effect_xx function)
#' @export
#' @examples
p_diffx <- function(y.name, x.name, trt.name, data, family, df = 3,
                    phi_predictions = NULL, span = 0.75) {
  ## defined commonly used theme
  common_theme <- theme_bw() + theme(legend.position = "top")
  
  # browser()
  dat <- data %>% dplyr::select(all_of(c(y.name, trt.name, x.name)))
  colnames(dat) <- c("y", "trt", "x")
  dat$trt <- as.factor(dat$trt)
  levels_trt <- levels(dat$trt)
  if (length(levels_trt) > 2) {
    stop("code only works for two treatments for now")
  }
  
  if(!is.null(phi_predictions)){
    dat$phi <- phi_predictions
  }
  plot_dat <- p_diffx_cal(dat, df, family, span = span)
  
  
  comb_plot(data = plot_dat, dat_raw = dat %>% select(x),
            x.name = x.name, y.name = y.name, lev.name = NULL, 
            diff_lab = "Probability difference")           
}



###############################################################################
p_diffxx <- function(data, x.names, y.name, trt.name, family, df = 3,
                     phi_predictions = NULL, span = 0.75) {
  common_theme <- theme_bw() + theme(legend.position = "top")
  
  # browser()
  dat <- data %>% dplyr::select(all_of(c(y.name, trt.name, x.names))) 
  colnames(dat) <- c("y", "trt", "x1", "x2")
  
  dat$trt <- as.factor(dat$trt)
  levels_trt <- levels(dat$trt)
  if (length(levels_trt) > 2) {
    stop("code only works for two treatments for now")
  }
  
  if(!is.null(phi_predictions)){
    dat$phi <- phi_predictions
  }
  
  dat <- dat %>% drop_na()
  
  if(is.numeric(dat$x1)){
    ## x1 continuous, x2 continuous, categorize x2 and then facet by x2, facet_name is the original name of x2
    dat = as.data.frame(dat) %>% rename(x = x1)
    if(is.numeric(dat$x2)){
      dat[,'x2'] = cut( dat[,'x2'] , breaks=quantile( dat[,'x2'], seq(0.00,1,1/3), na.rm = TRUE),  include.lowest = TRUE)
    }
    
    n_count <- dat %>%
      group_by(x2) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      mutate(x2_label = paste0(x.names[2], ": ", x2, " (n = ", n, ")"))
    
    dat <- dat %>% mutate(x2 = factor(x2, levels = n_count$x2, labels = n_count$x2_label))
    
    dat_plot <- lapply(split(dat, dat[, "x2"]), function(xx_dat){
      p_diffx_cal(dat = xx_dat, df = df, family = family,
                  span = span) %>% mutate(lev = xx_dat[1, "x2"])
    }) %>% bind_rows()
    
    comb_plot(data = dat_plot, dat_raw = dat %>% select(x, x2),
              x.name = x.names[1], y.name = y.name, lev.name = x.names[2], 
              diff_lab = "Probability difference")
    
  }else if (!is.numeric(dat$x1) & is.numeric(dat$x2)){
    ## x1 categorical, x2 continuous, facet by x1, facet_name is the original name of x1
    dat <- as.data.frame(dat) %>% rename(x = x2, x2 = x1)
    ## x continuous, x2 categorical
    
    n_count <- dat %>%
      group_by(x2) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      mutate(x2_label = paste0(x.names[1], ": ", x2, " (n = ", n, ")"))
    dat <- dat %>% mutate(x2 = factor(x2, levels = n_count$x2, labels = n_count$x2_label))
    
    dat_plot <- lapply(split(dat, dat[, "x2"]), function(xx_dat){
      p_diffx_cal(dat = xx_dat, df = df, family = family,
                  span = span) %>% mutate(lev = xx_dat[1, "x2"])
    }) %>% bind_rows()
    
    comb_plot(data = dat_plot, dat_raw = dat %>% select(x, x2),
              x.name = x.names[2], y.name = y.name, lev.name = x.names[1], 
              diff_lab = "Probability difference")
    
  } else{
    dat <- as.data.frame(dat) %>% rename(x = x1)
    
    n_count <- dat %>%
      group_by(x2) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      mutate(x2_label = paste0(x.names[2], ": ", x2, " (n = ", n, ")"))
    dat <- dat %>% mutate(x2 = factor(x2, levels = n_count$x2, labels = n_count$x2_label))
    
    ## x categorical, x2 categorical
    dat_plot <- lapply(split(dat, dat[, "x2"]), function(xx_dat){
      p_diffx_cal(dat = xx_dat, df = df, family = family,
                  span = span) %>% mutate(lev = xx_dat[1, "x2"])
    }) %>% bind_rows() %>% mutate(x = as.character(x)) %>% arrange(lev, x)
    
    comb_plot(data = dat_plot, dat_raw = dat %>% select(x, x2),
              x.name = x.names[1], y.name = y.name, lev.name = x.names[2], 
              diff_lab = "Probability difference")
  } 
}


## bivariate outcome plots
## combination of two categorical variables or one continuous variable and one categorical variable

#' Plot (smoothed) treatment difference against two baseline 
#' covariate, to assess potential interactions
#'
#' This function plots the estimated outcome estimated from original outcome against two
#' baseline covariates. If both baseline variables are categorical, the phi's are smoothed using loess
#' and displayed using a contour plotit is simple average with CI. If one baseline variable is continuous and one categorical
#' a regression spline is used (with user-specified df) for each category of the second variable.
#' 
#'
#' @param y.name Name of outcome variable in data
#' @param x.names Names of baseline variable in data
#' @param trt.name Name of treatment variable in data (code assumes
#'   there are only two treatments)
#' @param data Data set
#' @param df Degrees of freedom for the regression spline
#' @param family Family to use for outcome variable (passed to glm)
#' @export
#' @examples
p_yxx <- function(y.name, x.names, trt.name, data, df = 2, family = "binomial") {
  common_theme <- theme_bw() + theme(legend.position = "top")
  
  dat <- data %>% dplyr::select(all_of(c(y.name, trt.name, x.names))) %>% drop_na()
  colnames(dat) <- c("y", "trt", "x1", "x2")
  
  dat$trt <- as.factor(dat$trt)
  levels_trt <- levels(dat$trt)
  if (length(levels_trt) > 2) {
    stop("code only works for two treatments for now")
  }
  
  if(is.numeric(dat$x1)){
    ## x1 continuous, x2 continuous, categorize x2 and then facet by x2, facet_name is the original name of x2
    dat_new <- as.data.frame(dat) %>% rename(x = x1)
    if(is.numeric(dat$x2)){
      dat_new[,'x2'] = cut( dat_new[,'x2'] , breaks=quantile( dat_new[,'x2'], seq(0.00,1,1/3), na.rm = TRUE),  include.lowest = TRUE)
    }
    
    n_count <- dat_new %>%
      group_by(x2) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      mutate(x2_label = paste0(x.names[2], ": ", x2, " (n = ", n, ")"))
    
    dat_new <- dat_new %>% mutate(x2 = factor(x2, levels = n_count$x2, labels = n_count$x2_label)) %>% rename(lev = x2)
    
    comb_plot_yx(data = dat_new, x.name = x.names[1], y.name, trt.name, 
                 lev.name = x.names[2], family = family, df = df)
    
  }else if (!is.numeric(dat$x1) & is.numeric(dat$x2)){
    ## x1 categorical, x2 continuous, facet by x1, facet_name is the original name of x1
    dat_new <- as.data.frame(dat) %>% rename(x = x2, x2 = x1)
    ## x continuous, x2 categorical
    
    n_count <- dat_new %>%
      group_by(x2) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      mutate(x2_label = paste0(x.names[1], ": ", x2, " (n = ", n, ")"))
    dat_new <- dat_new %>% mutate(x2 = factor(x2, levels = n_count$x2, labels = n_count$x2_label)) %>% rename(lev = x2)
    
    comb_plot_yx(data = dat_new, x.name = x.names[2], y.name, trt.name, 
                 lev.name = x.names[1], family = family, df = df)
    
  } else{
    dat_new <- as.data.frame(dat) %>% rename(x = x1)
    
    dat_new <- lapply(split(dat_new, dat_new$x2), function(xx){
      n_count <- xx %>%
        group_by(x) %>%
        summarise(n = n()) %>%
        ungroup() %>%
        mutate(x_label = paste0(x, " (n = ", n, ")"))
      xx %>% mutate(x = factor(x, levels = n_count$x, labels = n_count$x_label)) %>% 
        mutate(x2 = paste0(x.names[2], ": ", xx$x2[1], " (n = ", nrow(xx), ")"))
    }) %>% bind_rows() %>% rename(lev = x2)
    
    comb_plot_yx(data = dat_new, x.name = x.names[1], y.name, trt.name, 
                 lev.name = x.names[2], family = family, df = df)
  } 
}


########################################################################
## helper functions to create exhaustive subgroup plot (funnel plot)


#' function to guess variable type (numeric or categorical)
#'
#' @param x Data frame
#' @export
#' @examples
guess_type <- function(x) {
  N <- length(x)
  uN <- length(unique(x))
  ind <- !is.na(x) ## remove missings
  oldw <- getOption("warn")
  options(warn = -1)
  xx <- as.double(x[ind])
  options(warn = oldw)
  mn <- mean(is.na(xx))
  if (((uN > 8) | (N < 40 & uN / N > 0.4)) & mn < 0.01) { ## assume it's numeric
    return("numeric")
  } else {
    if (uN < 40) {
      return("categorical")
    } else {
      return("undigestable")
    }
  }
}

#' function to categorize a continuous variable in X categorizes
#' according to quantiles and returns in a list (i) the categorized
#' variable (ii) the category labels and (iii) character variable with
#' describing the link between labels and cut-offs
#'
#' @param x Data frame
#' @param nr Number of categories to categorize a numerical variable
#' @param var_nam Intended name for the variable
#' @export
#' @examples

categ_var <- function(x, nr, var_nam) {
  psq <- seq(0, 1, length = nr + 1)
  qq <- quantile(x, psq, na.rm = TRUE)
  delta <- diff(range(x, na.rm = TRUE))
  qq[1] <- min(x, na.rm = TRUE) - 0.02 * delta
  qq[nr + 1] <- max(x, na.rm = TRUE) + 0.02 * delta
  brks <- unique(qq)
  if(length(brks) == 3){
    labs <- c("low", "high")
    char <- sprintf("%s: low (<= %s) | high (> %s)",
                    var_nam, brks[2], brks[2])
  }
  if(length(brks) == 4){
    labs <- c("low", "mid", "high")
    char <- sprintf("%s: low (<= %s) | mid (%s, %s] | high (> %s)",
                    var_nam, brks[2], brks[2], brks[3], brks[3])
  }
  ct <- cut(x, brks, labels = labs)
  list(ct, labs, char)
}

#' function to produce categorical subgroup variables
#'
#' @param data Data frame
#' @export
#' @examples
get_subg_var <- function(data) {
  types <- sapply(data, guess_type)
  nams <- names(types)
  ind <- types == "undigestable"
  if (any(ind)) {
    types <- types[-ind]
    nams <- nams[-ind]
  }
  lst <- levs <- vector("list", length = length(nams))
  chars <- character(length(nams))
  for (i in 1:length(nams)) {
    if (types[i] == "numeric") 
      tmp <- categ_var(data[, get(nams[i])], 3, nams[i])
    if (types[i] == "categorical"){
      labs <- as.character(unique(data[, get(nams[i])]))
      tmp <- list(data[, get(nams[i])],
                  labs,
                  sprintf("%s: %s", nams[i], paste0(labs, collapse=" | ")))
    }
    lst[[i]] <- tmp[[1]]
    levs[[i]] <- tmp[[2]]
    chars[i] <- tmp[[3]]
  }
  names(lst) <- nams
  names(levs) <- nams
  list(subgr_vars = lst, subgr_levels = levs, chars = chars)
}


