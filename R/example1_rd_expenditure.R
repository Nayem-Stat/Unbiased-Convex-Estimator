# Clear environment
rm(list = ls())

# Load required libraries
library(MASS)
library(Matrix)
library(writexl)

# Load the dataset
data <- read.csv("C:\\Users\\hnayem\\OneDrive - Florida International University\\Documents\\Private\\Folder for Simulation work\\Collab_Paper\\example_data.csv")

# Extract response and predictor variables from the dataset
y <- data$Y
X <- as.matrix(data[, c("X1", "X2", "X3", "X4")])

# Load the corrplot package
library(corrplot)
# Save plot as PNG file with 300 DPI
png("plot2.png", width = 6, height = 6, units = "in", res = 300)

# Calculate the correlation matrix for X
cor_matrix <- cor(X)

# Generate the correlation plot
corrplot(cor_matrix, method = "number")

# Close the PNG device
dev.off()


# Load the necessary package
library(Matrix)

# Compute the condition number of X
condition_number <- kappa(X)

# Print the condition number
print(condition_number)



# Set parameters
p <- ncol(X)  # Number of predictors
n <- nrow(X)  # Number of observations
k_values <- seq(0.1, 0.9, by = 0.1)  # Ridge parameter
sigma_squared <- var(y - mean(y))  # Estimate error variance from data

# Adjust true beta based on the number of predictors (p)
beta <- rep(1, p)  # Set hypothetical true beta values: 1, 1, ..., 1

# Initialize results data frame
results_df <- data.frame()

# Function to compute MSE
calculate_mse <- function(beta_hat, true_beta) {
  mean((beta_hat - true_beta)^2)
}

# Ensure no negative values in X (replace negative with absolute value)
X[X < 0] <- abs(X[X < 0])

# Identity matrix
I <- diag(p)
S <- t(X) %*% X

# Loop over ridge parameter values
for (k in k_values) {

  cat("Running estimators for k =", k, "\n")

  # 1. Ordinary Least Squares (OLS)
  beta_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  ols_mse <- calculate_mse(beta_ols, beta)

  # 2. Ordinary Ridge Regression (ORR)
  beta_orr <- solve(t(X) %*% X + k * I) %*% t(X) %*% y
  orr_mse <- calculate_mse(beta_orr, beta)

  # 3. Unbiased Ridge Regression (URR)
  var_j_urr <- (sigma_squared / k) * I
  var_j_urr_pd <- nearPD(var_j_urr)$mat
  chol_var_j_urr <- chol(var_j_urr_pd)
  J_urr <- mvrnorm(1, mu = beta, Sigma = chol_var_j_urr)
  beta_urr <- solve(t(X) %*% X + k * I) %*% (t(X) %*% y + k * J_urr)
  urr_mse <- calculate_mse(beta_urr, beta)

  # 4. Restricted Least Squares (RLS)
  R <- matrix(c(-1, -2, 1, 0,0,0,-1,2), nrow = 2, ncol = 4)  # Adjust restriction based on p = 4
  r <- matrix(c(0,0), nrow = 2,ncol=1)
  beta_rls <- beta_ols + solve(S) %*% t(R) %*% solve(R %*% solve(S) %*% t(R)) %*% (r - R %*% beta_ols)
  rls_mse <- calculate_mse(beta_rls, beta)

  # 5. Convex Restricted Estimator (CR)
  M <- I - solve(S) %*% t(R) %*% solve(R %*% solve(S) %*% t(R)) %*% R
  W <- solve(I + k * solve(S))
  var_j <- sigma_squared * W %*% M %*% solve(S) %*% solve(I - W)
  var_j_pd <- nearPD(var_j)$mat
  chol_var_j <- chol(var_j_pd)
  J_cr <- mvrnorm(1, mu = beta, Sigma = chol_var_j)
  beta_cr <- W %*% beta_rls + (I - W) %*% J_cr
  cr_mse <- calculate_mse(beta_cr, beta)

  # Store results for each k
  temp_df <- data.frame(
    k = k,
    URR_MSE = urr_mse,
    RLS_MSE = rls_mse,
    CR_MSE = cr_mse,
    ORR_MSE = orr_mse,
    OLS_MSE = ols_mse
  )

  # Append to the results data frame
  results_df <- rbind(results_df, temp_df)
}

results_df
# Save results to Excel
write_xlsx(results_df, "C:\\Users\\Student\\OneDrive - Florida International University\\Documents\\Private\\Folder for Simulation work\\Collab_Paper\\estimator_results_data_10_26.xlsx")
