# load data and source code
library(bootstrap) # for jackknife standard error
dt <- read.csv("data_example_252pts.csv")
source("fcn_median_estimation.R")
names(dt)

# Specify baseline covariate data frame
x_df <- subset(dt, select = -c(r_vec, y_vec, e_vec))

# Estimate working model of expected Y given X
# Here we use linear regression, using all the patients with Y observed
terms <- paste(colnames(x_df), collapse = "+")
formula <- as.formula(paste("y_vec~", terms))  # just an easier way to throw 
# all the X into the regression model
print(formula)
fit_glm <- lm(formula, data = subset(dt, r_vec == 1))

# Specify the initial working model coefficients
beta_init <- coefficients(fit_glm)

print(beta_init)

# Specify the form of distribution of Y given X
# Here we assume Y|X follows a normal distribution with mean X*beta and standard error 1.
py_vec_fcn = function(t, x, beta){
  pnorm(t, mean = x%*%beta[-1] + beta[1], sd = 1) # beta[1] is intercept
}


est <- deduct_est_median(x_df, dt$r_vec, dt$y_vec,
                         dt$e_vec, py_vec_fcn,
                         beta_init)
print(est)

deduct_est_median(x_df, dt$r_vec, dt$y_vec,
                  dt$e_vec, py_vec_fcn,
                  beta_init, free_idx = "sf36_phy_scr") # sf36_phy_scr is continuous

deduct_est_median(x_df, dt$r_vec, dt$y_vec,
                  dt$e_vec, py_vec_fcn,
                  beta_init, free_idx = "female") # female is binary

est_w_se <- deduct_est_median(x_df, dt$r_vec, dt$y_vec,
                              dt$e_vec, py_vec_fcn,
                              beta_init, jack_se = TRUE)
print(est_w_se)