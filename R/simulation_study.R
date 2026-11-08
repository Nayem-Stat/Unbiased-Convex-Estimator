# Clear environment
rm(list = ls())

# Load required libraries
library(MASS)
library(Matrix)
library(writexl)

# Function to generate the design matrix X
generate_X = function(n, p, rho) {
  Z = matrix(rnorm(n * p), nrow = n, ncol = p)
  X = sqrt(1 - rho^2) * Z + rho * matrix(rnorm(n), nrow = n, ncol = p)
  return(X)
}

# Function to compute MSE
calculate_mse = function(beta_hat, true_beta) {
  mean((beta_hat - true_beta)^2)
}

# Function to run the simulation for a single combination of n, p, rho, k, and sigma_squared
run_simulation = function(n, p, rho, k, sigma_squared, beta, num_sim) {

  mse_results = matrix(NA, nrow = num_sim, ncol = 5)
  colnames(mse_results) = c("URR", "RLS", "CR", "ORR", "OLS")

  for (sim in 1:num_sim) {

    # Generate the design matrix X
    X = generate_X(n, p, rho)

    # Generate response variable y
    y = X %*% beta + rnorm(n, mean = 0, sd = sqrt(sigma_squared))

    # Ensure no negative values in X (not y)
    X[X < 0] = abs(X[X < 0])

    I = diag(p)
    S = t(X) %*% X

    # 1. Ordinary Least Regression (OLR)
    beta_ols = solve(t(X) %*% X) %*% t(X) %*% y
    ols_mse = calculate_mse(beta_ols, beta)

    # 2. Ordinary Ridge Regression (ORR)
    beta_orr = solve(t(X) %*% X + k * I) %*% t(X) %*% y
    orr_mse = calculate_mse(beta_orr, beta)

    # 3. Unbiased Ridge Regression (URR)
    var_j_urr = (sigma_squared / k) * I
    var_j_urr_pd = nearPD(var_j_urr)$mat
    chol_var_j_urr = chol(var_j_urr_pd)
    J_urr = mvrnorm(1, mu = beta, Sigma = chol_var_j_urr)
    beta_urr = solve(t(X) %*% X + k * I) %*% (t(X) %*% y + k * J_urr)
    urr_mse = calculate_mse(beta_urr, beta)

    # 4. Restricted Least Squares (RLS)
    if (p == 3) {
      R = matrix(c(1, 1, 0), nrow = 1, ncol = 3)
    } else if (p == 5) {
      R = matrix(c(1, 1, -2, 0, 0), nrow = 1, ncol = 5)
    } else if (p == 10) {
      R = matrix(c(1, 1, 0, 0, 1, 0, -2, -1, 1, 0), nrow = 1, ncol = 10)
    }

    r = matrix(c(0), nrow = 1, ncol = 1)
    beta_rls = beta_ols + solve(S) %*% t(R) %*% solve(R %*% solve(S) %*% t(R)) %*% (r - R %*% beta_ols)
    rls_mse = calculate_mse(beta_rls, beta)

    # 5. Convex Restricted Estimator (CR)
    M = I - solve(S) %*% t(R) %*% solve(R %*% solve(S) %*% t(R)) %*% R
    W = solve(I + k * solve(S))
    var_j = sigma_squared * W %*% M %*% solve(S) %*% solve(I - W)
    var_j_pd = nearPD(var_j)$mat
    chol_var_j = chol(var_j_pd)
    J_cr = mvrnorm(1, mu = beta, Sigma = chol_var_j)
    beta_cr = W %*% beta_rls + (I - W) %*% J_cr
    cr_mse = calculate_mse(beta_cr, beta)

    # Store MSE results for this simulation
    mse_results[sim, ] = c(urr_mse, rls_mse, cr_mse, orr_mse, ols_mse)
  }

  # Return the average MSEs across all simulations
  colMeans(mse_results)
}

# Set parameters
n_values = c(20, 30, 50, 100)  # Sample sizes
p_values = c(3, 5, 10)  # Number of predictors
k_values = seq(0.1, 0.9, by = 0.1)  # Ridge parameter
num_sim = 10000  # Number of simulations
rho_values = c(0.5, 0.7, 0.8)  # Correlation level
sigma_squared_values = c(1, 5)  # Error variance

# To store results
results_df = data.frame()

# Loop over sample sizes, number of predictors, rho values, and k values
for (n in n_values) {
  for (p in p_values) {

    # Adjust true beta based on the number of predictors (p)
    beta_p = rep(1, p)  # True beta values: 1, 1, ..., 1 (length p)

    for (rho in rho_values) {
      for (sigma_squared in sigma_squared_values) {
        for (k in k_values) {

          cat("Running simulation for n =", n, ", p =", p, ", rho =", rho, ", sigma_squared =", sigma_squared, ", k =", k, "\n")

          # Run the simulation for the given parameters
          mse_averages = run_simulation(n, p, rho, k, sigma_squared, beta_p, num_sim)

          # Create a temporary data frame to hold results for this combination
          temp_df = data.frame(
            n = n,
            p = p,
            rho = rho,
            sigma_squared = sigma_squared,
            k = k,
            URR_MSE = mse_averages[1],
            RLS_MSE = mse_averages[2],
            CR_MSE = mse_averages[3],
            ORR_MSE = mse_averages[4],
            OLS_MSE = mse_averages[5]
          )

          # Append to the results data frame
          results_df = rbind(results_df, temp_df)
        }
      }
    }
  }
}

# Save results to Excel
write_xlsx(results_df, "C:\\Users\\hnayem\\OneDrive - Florida International University\\Documents\\simulation_result_1_28_26.xlsx")
