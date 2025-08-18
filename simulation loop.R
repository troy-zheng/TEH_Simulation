# purer.R (modified)
library(tidyverse)
library(nnls)
library(SuperLearner)
library(permimp)
library(dplyr)
library(clustermq)
library(ggplot2)
library(ggpubr)
library(GGally)
library(vivid)
library(tableone)
library(quadprog)
library(naniar)
library(caret)
library(xgboost)
library(benchtm)
library(dplyr)

run_simulation <- function(seed) {
  set.seed(seed)
  
  ## Step 1
  scal <- 1
  X <- generate_X_syn(n=500* scal)
  trt <- generate_trt(n=500* scal, p_trt = 0.5)
  plot(density(X$X11))
  
  
  x <- X$X11
  m <- sum(!is.na(x))
  k <- ceiling(0.50 * m)   # number of 1's desired
  
  # rank from largest to smallest; break ties randomly to hit target proportion
  rk <- rank(-x, ties.method = "random")
  X$X31 <- as.integer(rk <= k)   # 1 for top 25%, else 0
  
  # check proportion and correlation
  mean(X$X11_bin, na.rm = TRUE)                  # ~0.25
  cor(X$X11, X$X31, use = "complete.obs")        # typically ~0.85–0.95
  
  plot(X$X11, X$X31)
  summary(X$X31)
  
  #X31 is the variable region. Here we decide to let 25% of the sample to be a same region (eg: North America) while the rest become the other group.
  
  #Define prog and pred function and calculate b1
  prog <- "2.32*(0.5*(X1=='Y')+X11)"
  pred <- "pnorm(20*(X11-0.5)) "
  coeff_all <- get_b(X, scal,prog, pred, trt,
                     type = "continuous",
                     power = c(0.5, 0.8), alpha = c(0.025, 0.1),
                     start = c(0, 0), sigma_error = 1
  )
  
  b1_star <- coeff_all[2]
  # we need to consider five situations: b1=0,b1=0.5*(b1*),b1=1*(b1*),b1=1.5*(b1*),b1=2*(b1*).
  
  #We simulate b1=2*(b1*).
  b0_star <- get_b0(X, scal, prog, pred, trt,
                    b1 = b1_star*2, type = "continuous",
                    power = 0.5, alpha = 0.025, interval = c(-2, 2), sigma_error = 1
  )
  #Calculate response variable based on b0 and b1
  dat <- generate_y(X, trt, prog ,
                    pred , b0 = b0_star, b1 =b1_star*2,
                    type = "continuous", sigma_error = 3)
  
  #Recode categorical variables as 0/1
  dat$X1 <- ifelse(dat$X1 == "Y", 1, 0)
  dat$X2 <- ifelse(dat$X2 == "Y", 1, 0)
  dat$X4 <- ifelse(dat$X4 == "Y", 1, 0)
  dat$X6 <- ifelse(dat$X6 == "Y", 1, 0)
  dat$X7 <- ifelse(dat$X7 == "Y", 1, 0)
  dat$X8 <- ifelse(dat$X8 == "Y", 1, 0)
  dat$X9 <- ifelse(dat$X9 == "F", 1, 0)
  levels <- c("a", "b", "c", "d", "e")
  dat$X3 <- match(dat$X3, levels) - 1
  
  
  X <- dat %>% dplyr::select(starts_with("X"))
  Y <- dat$Y
  trt <- dat$trt
  
  ## randomly generate missing value for baseline covariate with missing proportion 5%, 10%, 20% for X2, X5, X6 
  ## (for demonstration of data preprocessing)
  X2_miss_index <- rbinom(n = nrow(X), size = 1, p = 0.05)
  X5_miss_index <- rbinom(n = nrow(X), size = 1, p = 0.10)
  X6_miss_index <- rbinom(n = nrow(X), size = 1, p = 0.20)
  
  X[which(X2_miss_index == 1), 2] <- NA
  X[which(X5_miss_index == 1), 5] <- NA
  X[which(X6_miss_index == 1), 6] <- NA
  
  dat <- cbind(Y, trt, X) %>% mutate_at("trt", as.factor)
  
  saveRDS(dat, file = "/Users/sarah/Downloads/Figures/analysis_data.rds")
  
  
  ## Step 2
  
  
  adat <- dat
  readRDS(file = "/Users/sarah/Downloads/Figures/analysis_data.rds")
  vars <- c("X5", "X6")
  pp_plots <- lapply(vars, plot_density_bar, data = adat)
  res <- ggpubr::ggarrange(pp_plots[[1]], pp_plots[[2]], nrow = 1, ncol = 2)
  
  
  #dir.create(file.path("/Users/sarah/Downloads/Figures"))
  ggsave("/Users/sarah/Downloads/Figures/ida_a.pdf",
         res,
         width = 6, height = 4, units = "in"
  )
  
  ## text based summary
  tab <- tableone::CreateTableOne(strata = "trt", data = adat)
  print(tab)
  # tab <- tableone::CreateTableOne(strata = "STUDYID", data = adat)
  # print(tab)
  
  ## side-by-side boxplot or stacked bar-plots
  vars <- c("X5", "X6")
  pp_plots <- lapply(vars, plot_boxbar_by, y.var = "trt", data = adat)
  
  res <- ggpubr::ggarrange(pp_plots[[1]], pp_plots[[2]], nrow = 1, ncol = 2)
  
  ggsave("/Users/sarah/Downloads/Figures/ida_b.pdf",
         res,
         width = 6, height = 4, units = "in"
  )
  
  
  ## missing variables
  p1 <- naniar::gg_miss_var(adat %>% select(X1, X2, X3, X4, X5, X6, X7, X8), 
                            show_pct = TRUE)
  ## missing variable patterns
  p2 <- naniar::vis_miss(adat %>% select(X1, X2, X3, X4, X5, X6, X7, X8)) +
    coord_flip()
  
  res <- ggpubr::ggarrange(p1, p2, nrow = 1, ncol = 2)
  
  ggsave("/Users/sarah/Downloads/Figures/ida_c.pdf",
         res,
         width = 6, height = 4, units = "in"
  )
  
  
  
  ## identification of uninformative variables
  nzv <- caret::nearZeroVar(adat[, -1], saveMetrics = TRUE)
  head(nzv %>% arrange(desc(nzv)), n = 20)
  
  adat_dep <- adat[, 3:dim(adat)[2]]
  dependencies <- get_dep(adat_dep)
  
  ## assess dependencies by type of variable comparison
  dependencies$results %>%
    filter(Comparison == "continuous" & Correlation > 0.7) %>%
    arrange(desc(Correlation))
  ## for categorical use lower correlation threshold (as max achievable value can be <1)
  dependencies$results %>%
    filter(Comparison == "categorical" & Correlation > 0.5) %>%
    arrange(desc(Correlation))
  ## mixed comparisons
  dependencies$results %>%
    filter(Comparison == "mixed" & Correlation > 0.5) %>%
    arrange(desc(Correlation))
  
  ## use hierarchical clustering to assess similarity of variables
  hc <- hclust(as.dist(1 - dependencies$cor_mat), method = "average")
  pdf(file = "/Users/sarah/Downloads/Figures/ida_d.pdf", height = 4, width = 6)
  plot(hc, hang = -1, xlab = NA, sub = NA, cex = 0.67)
  dev.off()
  
  
  
  ## Step 3
  
  
  #adat <- readRDS(file = "/Users/sarah/Downloads/analysis_data.rds")
  
  #######################################################################
  ## Variable removal or transformations
  #######################################################################
  ## remove variables due to large missingness
  adat <- adat %>% select(-c())
  ## remove variables due to uninformativeness
  adat <- adat %>% select(-c())
  ## remove duplicated or highly correlated variables
  adat <- adat %>% select(-c())
  ## merge sparse categories
  adat <- adat 
  ## variable transformations
  adat <- adat %>%
    mutate(
      log_X5 = log(X5 + 0.01),
      # log_X13 = log(X13),
      log_X13 = log(X13 + 0.01),
      log_X25 = log(X25 + 0.01)
    ) %>% 
    select(-c(X5, X13, X25))
  
  
  ## save transformed data
  saveRDS(adat, file = "/Users/sarah/Downloads/Figures/analysis_data2.rds")
  
  adat <- adat %>%
    mutate_if(
      sapply(adat, class) %in% c("integer", "numeric"),
      as.numeric
    ) %>%
    mutate_if(
      sapply(adat, class) %in% c("factor", "character"),
      as.factor
    )
  
  
  ## Missing Value variables check again
  missing_value_count <- adat %>%
    summarise(., across(everything(), ~ sum(is.na(.)))) %>% # Use this under R >= 4.0.0
    # summarise_all(funs(sum(is.na(.)))) %>%  # Use this line under R < 4
    as.data.frame() %>%
    `rownames<-`("na_count") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    arrange(desc(na_count))
  
  head(missing_value_count)
  
  vars_impute <- c(missing_value_count %>%
                     filter(na_count > 0 & rowname != "Y") %>%
                     select(rowname)) # select covariates includes missing value
  
  
  ## optional step: Add indicator variable for missingness
  # adat <- adat %>%
  #  mutate_at(vars_impute$rowname, list(missing = ~ is.na(.) * 1.))
  
  ## Median/Mode imputation based on same variable
  df_univariate_impute <- adat
  
  # Median imputation for numeric variables
  for (rn in vars_impute$rowname) {
    if (class(df_univariate_impute[[rn]]) == "numeric") {
      df_univariate_impute[rn][is.na(df_univariate_impute[rn])] <- median(df_univariate_impute[[rn]], na.rm = TRUE)
    }
  }
  
  # Mode imputation for factor variables
  df_univariate_impute <- randomForest::na.roughfix(df_univariate_impute)
  exclude_vars <- c("Y", "trt")
  df_regress_impute <- adat
  ## you can add remove.collinear=FALSE if you have collinear variables in dataset to
  ## prevent imputed dataset contains NA
  imputation <- mice::mice(df_regress_impute %>% dplyr::select(-all_of(exclude_vars)),
                           m = 10, maxit = 5, method = "pmm", printFlag = TRUE, seed = 2020
  )
  # Find mean/mode of 10 imputations
  for (v in vars_impute$rowname) {
    if (class(df_regress_impute[[v]]) == "numeric") {
      imp_median <- apply(imputation$imp[[v]], 1, median)
      df_regress_impute[[v]][as.numeric(names(imp_median))] <- imp_median
    }
    if (class(df_regress_impute[[v]]) == "factor") {
      imp_most <- apply(
        imputation$imp[[v]], 1,
        function(x) {
          names(sort(table(x), decreasing = TRUE)[1])
        }
      )
      df_regress_impute[[v]][as.numeric(names(imp_most))] <- imp_most
    }
  }
  
  
  ## Save imputed data
  saveRDS(df_univariate_impute, file = "/Users/sarah/Downloads/Figures/analysis_univ_imputed.rds")
  saveRDS(df_regress_impute, file = "/Users/sarah/Downloads/Figures/analysis_regress_imputed.rds")
  
  
  ## Step 4
  
  adat <- readRDS(file = "/Users/sarah/Downloads/Figures/analysis_regress_imputed.rds")
  write.csv(adat, file = "/Users/sarah/Downloads/Figures/your_dataset.csv", row.names = FALSE)
  
  y <- adat$Y
  trt <- select(adat, trt)$trt
  trt <- as.numeric(trt == "1")
  
  
  covs_out <- setdiff(colnames(adat), c("Y"))
  X_out <- select(adat, all_of(covs_out))
  ## select covariates to use for modelling treatment (make sure to remove outcome and treatment)
  covs_pi <- setdiff(colnames(adat), c("Y", "trt"))
  X_pi <- select(adat, all_of(covs_pi))
  
  ## super-learners to include
  SL.library <- list(
    "SL.glmnet", "SL.xgboost", "SL.ranger"
  )
  
  
  k <- 10
  ## create k splits, stratified by treatment group
  ids <- caret::createFolds(y = adat$trt, k = k)
  N <- nrow(adat)
  fit_one_fold <- function(i) {
    library(nnls)
    library(SuperLearner)
    train_ids <- setdiff(1:N, ids[[i]])
    test_ids <- ids[[i]]
    ## SuperLearner for outcome
    SL_out <- SuperLearner::SuperLearner(
      Y = y[train_ids],
      X = X_out[train_ids, ],
      SL.library = SL.library, family = "gaussian",
      verbose = FALSE, method = c("method.CC_LS")
    )
    ## SuperLearner for treatment
    SL_pi <- SuperLearner::SuperLearner(
      Y = trt[train_ids],
      X = X_pi[train_ids, ],
      SL.library = SL.library, family = "binomial",
      verbose = FALSE, method = c("method.CC_LS")
    )
    
    ## predictions needed for outcome
    pred_X1 <- pred_X0 <- pred_trt <- X_out[test_ids, ]
    pred_X1$trt <- factor("1", levels = levels(adat$trt))
    pred_X0$trt <- factor("0", levels = levels(adat$trt))
    ## prediction under control
    m0 <- as.numeric(predict(SL_out, newdata = pred_X0)$pred)
    ## prediction under treatment
    m1 <- as.numeric(predict(SL_out, newdata = pred_X1)$pred)
    ## prediction under observed arm
    mtrt <- as.numeric(predict(SL_out, newdata = pred_trt)$pred)
    ## predictions needed for treatment
    pi <- as.numeric(predict(SL_pi, newdata = X_pi[test_ids, ])$pred)
    list(
      i = i, train_ids = train_ids, test_ids = test_ids,
      m0 = m0, m1 = m1, mtrt = mtrt, pi = pi
    )
  }
  ## perform cross-fitting on grid
  export <- list(
    adat = adat, X_out = X_out, X_pi = X_pi,
    ids = ids, N = N, SL.library = SL.library, trt = trt, y = y
  )
  res <- clustermq::Q(fit_one_fold, i = 1:k, n_jobs = k, export = export)
  
  ## Calculate pseudo-observations phi
  m1 <- m0 <- mtrt <- pi <- numeric(N)
  for (i in 1:k) {
    ids <- res[[i]]$test_ids
    m0[ids] <- res[[i]]$m0
    m1[ids] <- res[[i]]$m1
    mtrt[ids] <- res[[i]]$mtrt
    pi[ids] <- res[[i]]$pi
  }
  
  ## double robust pseudo observations for the treatment difference (see https://arxiv.org/abs/2004.14497)
  phi <- m1 - m0 + (trt - pi) / (pi * (1 - pi)) * (y - mtrt)
  saveRDS(phi, file = "/Users/sarah/Downloads/Figures/phi.rds")
  
  n_cov <- ncol(X_pi)
  control_cforest <- party::cforest_unbiased(mtry = 5, ntree = 500)
  fit <- party::cforest(y ~ ., data.frame(X_pi, y = phi), control = control_cforest)
  ## standard deviation of observed "individual" treatment effects
  sd_obs <- sd(predict(fit, OOB = TRUE))
  ## assess variable importance (to get more stable VI rankings increase nperm)
  cf_vi <- permimp::permimp(fit, conditional = FALSE, nperm = 10)
  ## Optional: The conditional  permutation importance is prohibitive in terms of computational cost.
  # cf_vi_cond <- permimp::permimp(fit, conditional = TRUE, nperm = 10)
  
  
  ## assess variability across trees
  plot(cf_vi)
  plot(cf_vi, type = "box")
  
  
  pval <- numeric(n_cov)
  lm_fit_null <- lm(y ~ ., data = data.frame(y=phi))
  LR_null <- logLik(lm_fit_null)
  for (j in 1:n_cov) {
    lm_fit <- lm(y ~ ., data = data.frame(x=X_pi[,j], y=phi))
    pval[j] <- anova(lm_fit, lm_fit_null)[["Pr(>F)"]][2]
  }
  ## present p-value on "surprise" scale
  temp<-data.frame(variable = colnames(X_pi), surprise = pval)
  
  p_value <- temp %>%
    filter(variable == "X11") %>%
    select(variable, surprise)
  
  importance_scores <- data.frame(
    variable = names(cf_vi$values),
    vi = cf_vi$values #, cond_vi = cf_vi_cond$values
  )
  saveRDS(importance_scores, file = "/Users/sarah/Downloads/Figures/importance.rds")
  
  # 
  # ### Optional: interaction variable importance based on partial dependence function
  # ### refer paper: https://arxiv.org/pdf/1805.04755.pdf
  # ## only based on top nx variables from importance score
  nx <- 5
  importance_scores <- as_tibble(importance_scores)
  topn <- importance_scores %>% 
    dplyr::arrange(desc(vi)) %>%  # Sort by importance (descending)
    dplyr::slice(1:5) %>%         # Take top 5 rows
    dplyr::select(variable, vi)    # Keep only 'variable' and 'vi' columns
  
  
  #topn <- importance_scores[1:5,2] %>% arrange(desc(vi)) %>% slice(1:nx) %>% pull(variable)
  comb_var <- combn(topn, 2)
  
  p_value1 <- temp %>%
    filter(variable == "X11") %>%
    select(variable, surprise)
  
  return(list(p_value1 = p_value1, topn = topn))
}


### Simulate for 100 time
set.seed(123)
seeds <- sample(1:10000, 100)
results <- lapply(seeds, run_simulation)
x11_pvalues <- sapply(results, function(res) res$p_value1$surprise)

# Combine topn variables into one dataframe
top_variables <- do.call(rbind, lapply(results, function(res) res$topn))

saveRDS(results, file = "/Users/sarah/Downloads/Figures/simulation_results.rds")
write.csv(data.frame(seed = seeds, x11_pvalue = x11_pvalues),
          "/Users/sarah/Downloads/Figures/x11_pvalues.csv", row.names = FALSE)

# Count how often each variable appeared in top 5
top_var_freq <- top_variables %>%
  count(variable, sort = TRUE)

write.csv(top_var_freq, "/Users/sarah/Downloads/Figures/top_variable_frequency.csv", row.names = FALSE)
summary(x11_pvalues)
hist(x11_pvalues, breaks = 20, main = "Distribution of X11 p-values", xlab = "p-value")


###Export to the result to csv file
flat_results <- do.call(rbind, lapply(seq_along(results), function(i) {
  res <- results[[i]]
  # Extract p-value
  pval <- res$p_value1$surprise
  # Extract top 5 variables and importance
  top_vars <- res$topn
  top_vars$rank <- 1:nrow(top_vars)
  top_vars$seed <- i
  top_vars$x11_pvalue <- pval
  return(top_vars)
}))

flat_results <- flat_results %>%
  select(seed, rank, variable, vi, x11_pvalue)

head(flat_results, 10)
write.csv(flat_results, file = "/Users/sarah/Downloads/Figures/flat_simulation_results.csv", row.names = FALSE)

