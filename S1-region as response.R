
#### generate data from benchtm

library(benchtm)
library(dplyr)

# define seeds
seeds <- 1:15

# empty list to collect results
all_results <- list()

for (s in seeds) {
  set.seed(s)

#### data part: data generated on different scenarios saved in "scen_param" from library(benchtm)
data(scen_param)

cases <- scen_param %>% filter(type == "continuous" &
                                 pred == "pnorm(20*(X11-0.5))" &
                                 b1_rel == 0)

## generate data using benchtm
dat <- generate_scen_data(scen = cases,  include_truth = F) %>%
  mutate_if(is.character, as.factor) %>%
  mutate_at('trt', as.factor)

X <- dat %>% dplyr::select(starts_with("X"))
rk <- rank(dat$X11, ties.method = "random") #may change to a better rank method
X$X31 <- ifelse(rk<=125, 1, 0)
#dat$X31<-X$X31  
dat$Y<-X$X31
# rank from largest to smallest; break ties randomly to hit target proportion
#mean(X$X11_bin, na.rm = TRUE)                  # ~0.25
#cor(X$X11, X$X31, use = "complete.obs")        # typically ~0.85–0.95
#summary(X$X31)   # 1 for top 25%, else 0

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

# Check is the data folder exists, and if not create it
#dir.create(file.path("data"))

saveRDS(dat, file = "/Users/sarah/Downloads/Figures/analysis_data.rds")

library(tidyverse)

## assumes the code in 01_read_data.R has already be used to create data/analysis_data.rds
adat <- readRDS(file = "/Users/sarah/Downloads/Figures/analysis_data.rds")

## histogram/barplots of all variables (save plots as pdf under reports/ for easier review)
## for paper, only draw plot with 2 variables
vars <- c("X5", "X6")
pp_plots <- lapply(vars, plot_density_bar, data = adat)
res <- ggpubr::ggarrange(pp_plots[[1]], pp_plots[[2]], nrow = 1, ncol = 2)

# Check is the reports folder exists, and if not create it
dir.create(file.path("reports"))
# Check is the reports/figures folder exists, and if not create it
dir.create(file.path("reports/figures"))
ggsave("/Users/sarah/Downloads/Figures/ida_a.pdf",
       res,
       width = 6, height = 4, units = "in"
)

#######################################################################
## Stratified univariate summary of baseline variables
#######################################################################
## Observe baseline summaries stratified by another categorical factor.
## Typical stratification factors are "study or "treatment group".
## In this case we want to compare placebo vs Cosentyx 300mg
## but Cosentyx 300mg is not available in all studies. So both
## treatment and study are of interest as stratification factors

## text based summary
tab <- tableone::CreateTableOne(strata = "trt", data = adat)
print(tab)
# tab <- tableone::CreateTableOne(strata = "STUDYID", data = adat)
# print(tab)

## side-by-side boxplot or stacked bar-plots
vars <- c("X5", "X6")
pp_plots <- lapply(vars, plot_boxbar_by, y.var = "trt", data = adat)

res <- ggpubr::ggarrange(pp_plots[[1]], pp_plots[[2]], nrow = 1, ncol = 2)

ggsave("reports/figures/ida_b.pdf",
       res,
       width = 6, height = 4, units = "in"
)

#######################################################################
## Evalulate missing values/non-informative baseline variables
#######################################################################
## Here missingness & missingness patterns are explored further, in
## addition to assessment of variables with low information (e.g. all
## observations of one variable equal to one value). Both may suggest
## removal of certain variables. In addition for categorical variables
## it is observed whether there are variables with sparsely populated
## categories (may merge those categories with other categories).

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

## identification of categorical variables with sparse categories
## (code in src/util/baseline_dependency.R)
low_freq_categories(adat)

#######################################################################
## Assess dependencies across baseline variables
#######################################################################
## Provides a better understanding of the joint distribution of baseline
## variables and helps identify duplicate (or close-to-duplicate) variables in
## the data. In addition it may help with interpretation of final results.
## To assess dependency between two variables X and Y we calculate
## sqrt(X2/(X2+N)), where X2 is the chi-squared statistic comparing a model for
## p(X|Y) versus p(X).
## If both X and Y are continuous linear models can be used (adjusting for Y linearly)
## and the result is very close to the Pearson correlation.
## If both X and Y are categorical multinomial regression models can be used
## (adjusting for Y as categorical variable). This will give results close to
## Pearson's contingency coefficient.
## For mixed data both multinomial and continuous regression can be used (the
## results assessing p(Y|X) vs p(Y) and p(X|Y) vs p(X) are usually very similar.
## Downside of this approach: For non-continuous data, the maximum achievable
## value can be <1, so assess the results by data-type.

## remove outcome and some problematic variables (e.g. just 1 value)
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


## Including Plots


library(tidyverse)

adat <- readRDS(file = "/Users/sarah/Downloads/Figures/analysis_data.rds")

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

#######################################################################
## Impute remaining missing covariates
#######################################################################
## For missing data in the *outcome* variable we propose (for consistency)
## to follow the analytical strategy used in the main pre-specified
## clinical trial analyses for this endpoint, following the decided
## intercurrent event strategies and missing data handling approaches
## (see ICH E9 addendum). In situations, where these approaches are
## complex (e.g. various multiple imputation strategies used) subsequent
## analyses may become time-consuming or infeasible. In this case we
## suggest to use simpler analyses that are still in the spirit of
## the main estimand targeted in the clinical trial analyses.
## This may for example mean using single imputation approaches.

## For missing data in *baseline* variables, we provide following
## considerations and suggestions.
## * Consider to drop covariates for missingness > 10%-20%
## * We don't recommend using multiple imputation for baseline variables
##   as it would result in multiple analysis data-sets and thus more
##   time-consuming subsequent analyses.
## * For single imputation we recommend two methods:
##   (i) Imputation of median or mode of the non-missing values from the
##       same baseline variable.
##   or better
##   (ii) Regression imputation of the missing baseline
##        covariates/biomarkers based on all other variables.
## * While there are different ways of doing (ii) above we
##   recommend to perform multiple imputation to get multiple data-sets
##   with complete baseline variables, but then, rather than using all of 
##   those data-sets in subsequent analyses, taking median/mode on the imputed 
##   values for each missing to get a single dataset. The [`mice`](https://CRAN.R-project.org/package=mice)
##   provides a lot of flexibility and ease of use for this purpose.
## * The imputation of baseline covariates should be independent of
##   outcome or any other post-baseline variables.
## * The missing value indicator method creates extra binary covariates for
##   each covariate with missing data, indicating where a missing value was
##   observed. Running subsequent analyses could then be done using these
##   missing indicators. This allows to assess if the missingness is
##   informative [Groenwold et al (2012)](https://doi.org/10.1503/cmaj.110977).

## We encourage user to run analysis with different imputation methods for a stability check.

## Change variable class to either numeric or factor
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

## Single imputation based on regression multiple imputation

# Here multiple imputation of the baseline covariates based on all
# other baseline covariates is performed using the mice package. The
# median (continuous variable) or mode (categorical variable) of the 10
# multiply imputed values are finally used for the single imputed data-set.

## remove outcome and trt (as recommended)
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



## Derive variable importance scores and perform global heterogeneity
library(nnls)
library(SuperLearner)
library(permimp)
library(dplyr)
library(clustermq)
library(ggplot2)

adat <- readRDS(file = "/Users/sarah/Downloads/Figures/analysis_regress_imputed.rds")
y <- adat$Y
trt <- select(adat, trt)$trt
trt <- as.numeric(trt == "1")

## First models are build to model the outcome (and for observational
## studies, or when pooling across multiple studies for the treatment
## assignement), based on the stacking procedure as implemented in the
## ['SuperLearner'](https://CRAN.R-project.org/package=SuperLearner) R
## package. Different base models can be utilized for stacking. These
## models are then used to provide pseudo-observations $\phi$ for the
## treatment effect for each patient as in the double-robust learning
## algorithm [Kennedy (2022)](https://arxiv.org/abs/2004.14497). To be
## able to perform statistical inference later the $\phi$ values for
## each patient are created using cross-fitting based on 10-fold
## cross-validation. The prediction for the pseudo-observation for
## every patient is only based on the models where this patient was
## not used for fitting the models (out-of-fold predictions). Note
## that for small data-sets (and or low number of covariates) the
## utilized base models for the super-learner can/should be changed
## (e.g. standard unpenalized regression can also be used).

##############################################################################
## first create pseudo-observations phi (see https://arxiv.org/abs/2004.14497)
##############################################################################

## select covariates to use for modelling the outcome
## make sure to remove outcome
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

## Based on the obtained $\phi$ values a global heterogeneity test
## is performed and variable importance is assessed. This is done by
## fitting a random forest based on conditional inference trees (as
## implemented in the ['party'](https://CRAN.R-project.org/package=party)
## R package). One advantage of this approach versus the approach
## implemented in the randomForest package is that conditional inference
## trees do not suffer from variable selection bias (towards variables
## with many split possibilities) when choosing the variables to split
## [Hothorn et al (2006)](https://doi.org/10.1198/106186006X133933).

## Based on the fitted forest then the standard deviation of the
## model-based predictions of the treatment effect for all patients
## are extracted (a large standard deviation would indicate heterogeneity).
## Then the $\phi$ values are permuted against the considered covariates
## and for each permutation the model above is re-fitted and the standard
## deviation of the model-based predictions of the treatment effect is
## calculated. Under the permutation distribution a low standard
## deviation is expected, as the covariates will not explain the $\phi$
## values. A p-value can be extracted based on the proportion of
## permutations standard deviations that are larger than the observed
## standard deviation for the unpermuted data.

## Variable importance is calculated based on the fitted conditional
## random forest using the approaches outlined in [Debeer and Strobl
## (2020)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03622-2),
## as implemented in the
## ['party'](https://CRAN.R-project.org/package=party) and
## ['permimp'](https://CRAN.R-project.org/package=permimp) R
## packages. The basic idea of standard variable importance is to permute
## a covariate and then to observe how much the out-of-bag mean squared
## error for the treatment effect increases (averaged over the trees in
## the forest). To better quantify the unique contribution of a variable
## also conditional permutation variable importance is available, where
## a variable is permuted within strata of variables that are correlated
## with the variable of interest.

########################################################################
## perform global heterogeneity test and obtain variable importance
########################################################################
## Fit the conditional forest model to pseudo observations
n_cov <- ncol(X_pi)
## control parameters for conditional forest
## can assess sensitivity of results to mtry (mtry=5 is the default in cforest,
## n_cov/3 is the default in the randomForest package)
## larger ntree will reduce variability in results (e.g. VI rankings)
control_cforest <- party::cforest_unbiased(mtry = 5, ntree = 500)
fit <- party::cforest(y ~ ., data.frame(X_pi, y = Y), control = control_cforest)
## standard deviation of observed "individual" treatment effects
sd_obs <- sd(predict(fit, OOB = TRUE))
## assess variable importance (to get more stable VI rankings increase nperm)
cf_vi <- permimp::permimp(fit, conditional = FALSE, nperm = 10)
## Optional: The conditional  permutation importance is prohibitive in terms of computational cost.
# cf_vi_cond <- permimp::permimp(fit, conditional = TRUE, nperm = 10)


## assess variability across trees
plot(cf_vi)
plot(cf_vi, type = "box")


## Global heterogeneity test using coin
test <- coin::independence_test(y ~ ., data.frame(X_pi, y = phi), teststat="quadratic")
p_value <- coin::pvalue(test)
p_value

saveRDS(p_value,
        file = "/Users/sarah/Downloads/Figures/p_value.rds"
)

#
importance_scores <- data.frame(
  variable = names(cf_vi$values),
  vi = cf_vi$values #, cond_vi = cf_vi_cond$values
)
saveRDS(importance_scores, file = "/Users/sarah/Downloads/Figures/importance.rds")


## Summary plots see docs/FURTHER_INFO.md for more information.
library(ggplot2)
library(dplyr)
library(ggpubr)
library(GGally)
library(vivid)

## load data without missing covariates imputed
adat <- readRDS(file = "/Users/sarah/Downloads/Figures/analysis_data2.rds")
importance_scores <- readRDS(file = "/Users/sarah/Downloads/Figures/importance.rds")
phi <- readRDS(file = "/Users/sarah/Downloads/Figures/phi.rds")
p_value <- readRDS(file = "/Users/sarah/Downloads/Figures/p_value.rds")

# create a template, users can change it according to their preferences
plot_theme <- theme(
  # Hide panel borders and remove grid lines
  axis.text=element_text(size=12),
  axis.title=element_text(size=14),
  legend.title = element_text(size = 14),
  legend.text = element_text( size = 14),
  strip.text.x = element_text( size = 14),
  plot.title = element_text(size=14))

#####################
## Evidence of TEH ##
#####################
# No plot for visualizing the vidence, just the p-value

######################
## Effect modifiers ##
######################
## variable importance plot
## standard permutation variable importance
p_importance <- ggplot(head(importance_scores %>% arrange(-vi) %>% filter(vi>0), 5),
                       aes(reorder(variable, vi), vi)) +
  geom_col(color="black", fill = "lightblue" ) + 
  xlab("Variables") + ylab("") +
  # ylab("Variable importance\n(average increase in OOB-MSE\n under permutation)")+
  theme_classic() + coord_flip() + plot_theme 

ggsave("/Users/sarah/Downloads/Figures/var_importance.pdf", p_importance, width = 6, height = 4, units = "in")
# saveRDS(p_importance, file = "/Users/tianyuzheng/Documents/PKU RA/Novartis/Figures/var_importance.rds")


## select top 10 variables of interest
vars_to_plot <- rownames(head(importance_scores %>% arrange(-vi), 5))

vars_to_plot_df <- data.frame(
  variable = vars_to_plot,
  seed = s
)

colnames(vars_to_plot_df) <- "Top variables"
all_results[[s]] <- vars_to_plot_df

}

final_results <- bind_rows(all_results)
saveRDS(final_results, file = "/Users/sarah/Downloads/Figures/Top_ranking_result_regionY.rds")

head(final_results, 10)
sum(final_results$`Top variables` == "X11")

#final_results <- readRDS("/Users/sarah/Downloads/Figures/Top_ranking_result_0.rds")

############TEH Plot############
library(ggplot2)
df <- data.frame(
  heterogeneity = c(0, 0.5, 1, 1.5, 2),
  prob_X2 = c(0.133333333, 0.266666667, 0.466666667, 0.866666667, 1)
)

# Plot
p <- ggplot(df, aes(x = heterogeneity, y = prob_X2)) +
  geom_line(color = "blue", size = 1.2) +
  geom_point(color = "blue", size = 3) +
  labs(
    x = "Degree of Treatment Effect Heterogeneity",
    y = "Probability That X2 is Truly Predictive",
    title = "Predictive Power of X2 Across Heterogeneity Levels"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2))

p 


############
library(dplyr)

# Rename for clarity
colnames(final_results) <- c("variable", "sim_id")

# Count how many simulations each variable appeared in
selection_counts <- final_results %>%
  group_by(variable, sim_id) %>%
  summarise(appeared = 1, .groups = "drop") %>%   # one row per variable per sim
  group_by(variable) %>%
  summarise(simulations = n(), .groups = "drop")  # total number of sims where variable appeared

# Add probability (divide by total number of unique simulations)
n_sims <- n_distinct(final_results$sim_id)

selection_probs <- selection_counts %>%
  mutate(probability = round(simulations / n_sims, 3)) %>%
  arrange(desc(probability))

selection_probs$probability

df <- data.frame(
  method = rep(c( "DR-Perm(cfforest)"), times = 16),
  estimate =  selection_probs$probability)

# Compute group means for red circles
df_mean <- df %>%
  group_by(method) %>%
  summarise(mean_est = mean(estimate), .groups = "drop")

max_val <- max(df$estimate)
min_val <- min(df$estimate)

ggplot(df, aes(x = estimate, y = method)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(shape = 3, color = "skyblue", size = 2, stroke = 0.7) +        
  geom_point(data = df_mean, aes(x = mean_est, y = method),
             color = "red", size = 3, shape = 21, stroke = 1.2) +  
  labs(
    title = "Scenario 1:\n2.32*(0.5*(X1 == 'Y') + X11) + A*β0\n \nβ2*/β2=2",
    x = "Probability being selected", y = NULL
  ) +
  scale_x_continuous(
    limits = c(0, 1),  
    breaks = seq(0.05, 1, by = 0.1),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 9, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 9, angle = 45, hjust = 1, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks = element_line(color = "black")
  )
