# Clear environment
rm(list = ls())

# Load required libraries
library(MASS)
library(Matrix)
library(writexl)
library(corrplot)

# Load the dataset
data <- read.csv("C:\\Users\\hnayem\\OneDrive - Florida International University\\Documents\\Private\\PhD work\\Folder for Simulation work\\Collab_Paper\\example_data2.csv")

# Extract response and predictor variables
y <- data$Y
X <- as.matrix(data[, c("X1", "X2", "X3")])

# Save correlation plot
png("plot2.png", width = 6, height = 6, units = "in", res = 300)
cor_matrix <- cor(X)
corrplot(cor_matrix, method = "number")
dev.off()

# Compute condition number of X
condition_number <- kappa(X)
print(condition_number)

# Set parameters
p <- ncol(X)  # Number of predictors
n <- nrow(X)  # Number of observations
k_values <- seq(0.1, 0.9, by = 0.1)  # Ridge parameter
sigma_squared <- var(y - mean(y))  # Estimate error variance

# Define true beta
beta <- rep(1, p)  # Hypothetical true beta values

# Function to compute MSE
calculate_mse <- function(beta_hat, true_beta) {
  mean((beta_hat - true_beta)^2)
}

# Identity matrix
I <- diag(p)
S <- t(X) %*% X

# Initialize results data frame
results_df <- data.frame()

# Loop over ridge parameter values
for (k in k_values) {

  cat("Running estimators for k =", k, "\n")

  # Fit the OLS model
  ols_model <- lm(y ~ ., data = as.data.frame(X))
  beta_ols <- coef(ols_model)[-1]  # Exclude intercept

  # Compute MSE
  ols_mse <- calculate_mse(beta_ols, beta)

  # Ordinary Ridge Regression (ORR)
  beta_orr <- solve(S + k * I) %*% t(X) %*% y
  orr_mse <- calculate_mse(beta_orr, beta)

  # Unbiased Ridge Regression (URR)
  var_j_urr <- (sigma_squared / k) * I
  var_j_urr_pd <- nearPD(var_j_urr)$mat
  chol_var_j_urr <- chol(var_j_urr_pd)
  J_urr <- mvrnorm(1, mu = beta, Sigma = chol_var_j_urr)
  beta_urr <- solve(S + k * I) %*% (t(X) %*% y + k * J_urr)
  urr_mse <- calculate_mse(beta_urr, beta)

  # Restricted Least Squares (RLS) - Adjust R Matrix
  R <- matrix(runif(2 * p, -1, 1), nrow = 2, ncol = p)  # Random constraints
  r <- matrix(rep(0, 2), nrow = 2)

  R <- matrix(c(-1, -2,0,0,-1,2), nrow = 2, ncol = 3)  # Adjust restriction based on p = 4
  r <- matrix(c(0,0), nrow = 2,ncol=1)
  beta_rls <- beta_ols + solve(S) %*% t(R) %*% solve(R %*% solve(S) %*% t(R)) %*% (r - R %*% beta_ols)
  rls_mse <- calculate_mse(beta_rls, beta)

  # Convex Restricted Estimator (CR)
  M <- I - solve(S) %*% t(R) %*% solve(R %*% solve(S) %*% t(R)) %*% R
  W <- solve(I + k * solve(S))
  var_j <- sigma_squared * W %*% M %*% solve(S) %*% solve(I - W)
  var_j_pd <- nearPD(var_j)$mat
  chol_var_j <- chol(var_j_pd)
  J_cr <- mvrnorm(1, mu = beta, Sigma = chol_var_j)
  beta_cr <- W %*% beta_rls + (I - W) %*% J_cr
  cr_mse <- calculate_mse(beta_cr, beta)

  # Store results
  temp_df <- data.frame(
    k = k,
    URR_MSE = urr_mse,
    RLS_MSE = rls_mse,
    CR_MSE = cr_mse,
    ORR_MSE = orr_mse,
    OLS_MSE = ols_mse
  )

  results_df <- rbind(results_df, temp_df)
}

# Save results to Excel
#write_xlsx(results_df, "C:\\Users\\Student\\OneDrive - Florida International University\\Documents\\Private\\Folder for Simulation work\\Collab_Paper\\estimator_results_data_10_26.xlsx")

# Display results
print(results_df)
