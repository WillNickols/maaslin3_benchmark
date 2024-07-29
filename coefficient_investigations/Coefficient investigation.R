n <- 100
features <- 50
nsims <- 1000

outputs <- vector(length = nsims)
for (i in 1:nsims) {
  x <- matrix(rnorm(n, 0, 1), nrow = n)
  beta_abs <- rnorm(features, 0, 1)
  a <- exp(beta_abs %*% t(x) + matrix(rnorm(n * features, 0, 0.1), nrow = features)) # feature x samples
  drop_out_prob <- 0.7
  a <- t(apply(a, 1, function(x){ifelse(runif(length(x)) < 1 - drop_out_prob, x, 0)}))
  d <- colSums(a) # total depth
  y <- t(t(a) / d) # relative abundances (features x samples)
  
  betas_fit <- vector(length = features)
  for (feature_id in 1:features) {
    knock_outs <- sample(1:n, 10, replace = F)
    current_y <- y[feature_id,]
    current_a <- a[feature_id,]
    betas_fit[feature_id] <- 1/sum(x[current_y != 0]^2) * sum(x[current_y != 0] * log(current_y[current_y != 0]))# + 
       #1/sum(x[current_y != 0]^2) * sum(x[current_y != 0] * log(d[current_y != 0]))
    beta_abs[feature_id] <- 1/sum(x[current_a != 0]^2) * sum(x[current_a != 0] * log(current_a[current_a != 0]))
  }
  
  outputs[i] <- coef(lm(betas_fit ~ beta_abs))['beta_abs']
}
hist(outputs)
