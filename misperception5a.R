library(MASS) 
library(mvtnorm)

set.seed(123)
n <- 1000
n_simulations <- 500

# Initialize vectors to store coefficients
beta_X1e_values <- numeric(n_simulations)
beta_X2_values <- numeric(n_simulations)
beta_X2_value_model2 <- numeric(n_simulations)
beta_X1e_value_model3 <- numeric(n_simulations)

for (i in 1:n_simulations) {
  corr <- 0.8
  
  Sigma <- matrix(c(1, corr, corr, 1), ncol = 2) ; Sigma

  x <- mvtnorm::rmvnorm(n = n, mean = c(0, 0), sigma = Sigma)
  x1 <-  (x[,1])
  x1e <-  x1 + rnorm(n=n,sd=5)
  
  x2 <- x[,2]
  # Simulate response
  beta0 <- 1. 
  beta1 <- 2.5 # x1 has non-zero coefficient
  beta2 <- 0. # x2 has zero coefficient
  y <- beta0 + beta1*x1 + beta2*x2 + rnorm(n=n,mean=0,sd=1)
  
  # Fit regression model
  model1 <- lm(y ~  x1e + x2)
  model2 <- lm(y ~ x2)
  model3 <- lm(y ~ x1e)
  # Check model coefficients
  summary(model1)
  summary(model2)
  summary(model3)
  
  # coefficients for beta_X1e and beta_X2
  coeff1 <- model1$coefficients
  coeff2 <- model2$coefficients
  coeff3 <- model3$coefficients
  beta_X1e_values[i] <- coeff1['x1e']
  beta_X2_values[i] <- coeff1['x2']
  beta_X2_value_model2[i] <- coeff2['x2']
  beta_X1e_value_model3[i] <- coeff3['x1e']
  
}

# Calculate the mean values for beta_X1e and beta_X2 in each model
mean_beta_X1e <- mean(beta_X1e_values)
mean_beta_X2 <- mean(beta_X2_values)
mean_beta_X2_model2 <- mean(beta_X2_value_model2)
mean_beta_X1e_model3 <- mean(beta_X1e_value_model3)

# Print the mean values
cat("Mean value for beta_X1e in model 1:", mean_beta_X1e, "\n")
cat("Mean value for beta_X2 in model 1:", mean_beta_X2, "\n")
cat("Mean value for beta_X2 in model 2:", mean_beta_X2_model2, "\n")
cat("Mean value for beta_X1e in model 3:", mean_beta_X1e_model3, "\n")

mean_values_df <- data.frame(
  Variable = c("beta_X1e in model 1", "beta_X2 in model 1", "beta_X2 in model 2", "beta_X1e in model 3"),
  Mean_Value = c(mean_beta_X1e, mean_beta_X2, mean_beta_X2_model2, mean_beta_X1e_model3)
)

mean_values_df

