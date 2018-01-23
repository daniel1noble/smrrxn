// generated with brms 2.0.1
functions { 
} 
data { 
  int<lower=1> N;  // total number of observations 
  vector[N] Y_t_22;  // response variable 
  int<lower=1> K_t_22;  // number of population-level effects 
  matrix[N, K_t_22] X_t_22;  // population-level design matrix 
  vector[N] Y_t_24;  // response variable 
  int<lower=1> K_t_24;  // number of population-level effects 
  matrix[N, K_t_24] X_t_24;  // population-level design matrix 
  vector[N] Y_t_26;  // response variable 
  int<lower=1> K_t_26;  // number of population-level effects 
  matrix[N, K_t_26] X_t_26;  // population-level design matrix 
  vector[N] Y_t_28;  // response variable 
  int<lower=1> K_t_28;  // number of population-level effects 
  matrix[N, K_t_28] X_t_28;  // population-level design matrix 
  vector[N] Y_t_30;  // response variable 
  int<lower=1> K_t_30;  // number of population-level effects 
  matrix[N, K_t_30] X_t_30;  // population-level design matrix 
  vector[N] Y_t_32;  // response variable 
  int<lower=1> K_t_32;  // number of population-level effects 
  matrix[N, K_t_32] X_t_32;  // population-level design matrix 
  int<lower=1> nresp;  // number of responses
  int nrescor;  // number of residual correlations
  // data for group-level effects of ID 1 
  int<lower=1> J_1[N]; 
  int<lower=1> N_1; 
  int<lower=1> M_1; 
  vector[N] Z_1_t_22_1; 
  vector[N] Z_1_t_24_2; 
  vector[N] Z_1_t_26_3; 
  vector[N] Z_1_t_28_4; 
  vector[N] Z_1_t_30_5; 
  vector[N] Z_1_t_32_6; 
  int<lower=1> NC_1; 
  int prior_only;  // should the likelihood be ignored? 
} 
transformed data { 
  int Kc_t_22 = K_t_22 - 1; 
  matrix[N, K_t_22 - 1] Xc_t_22;  // centered version of X_t_22 
  vector[K_t_22 - 1] means_X_t_22;  // column means of X_t_22 before centering 
  int Kc_t_24 = K_t_24 - 1; 
  matrix[N, K_t_24 - 1] Xc_t_24;  // centered version of X_t_24 
  vector[K_t_24 - 1] means_X_t_24;  // column means of X_t_24 before centering 
  int Kc_t_26 = K_t_26 - 1; 
  matrix[N, K_t_26 - 1] Xc_t_26;  // centered version of X_t_26 
  vector[K_t_26 - 1] means_X_t_26;  // column means of X_t_26 before centering 
  int Kc_t_28 = K_t_28 - 1; 
  matrix[N, K_t_28 - 1] Xc_t_28;  // centered version of X_t_28 
  vector[K_t_28 - 1] means_X_t_28;  // column means of X_t_28 before centering 
  int Kc_t_30 = K_t_30 - 1; 
  matrix[N, K_t_30 - 1] Xc_t_30;  // centered version of X_t_30 
  vector[K_t_30 - 1] means_X_t_30;  // column means of X_t_30 before centering 
  int Kc_t_32 = K_t_32 - 1; 
  matrix[N, K_t_32 - 1] Xc_t_32;  // centered version of X_t_32 
  vector[K_t_32 - 1] means_X_t_32;  // column means of X_t_32 before centering 
  vector[nresp] Y[N];  // response matrix
  for (i in 2:K_t_22) { 
    means_X_t_22[i - 1] = mean(X_t_22[, i]); 
    Xc_t_22[, i - 1] = X_t_22[, i] - means_X_t_22[i - 1]; 
  } 
  for (i in 2:K_t_24) { 
    means_X_t_24[i - 1] = mean(X_t_24[, i]); 
    Xc_t_24[, i - 1] = X_t_24[, i] - means_X_t_24[i - 1]; 
  } 
  for (i in 2:K_t_26) { 
    means_X_t_26[i - 1] = mean(X_t_26[, i]); 
    Xc_t_26[, i - 1] = X_t_26[, i] - means_X_t_26[i - 1]; 
  } 
  for (i in 2:K_t_28) { 
    means_X_t_28[i - 1] = mean(X_t_28[, i]); 
    Xc_t_28[, i - 1] = X_t_28[, i] - means_X_t_28[i - 1]; 
  } 
  for (i in 2:K_t_30) { 
    means_X_t_30[i - 1] = mean(X_t_30[, i]); 
    Xc_t_30[, i - 1] = X_t_30[, i] - means_X_t_30[i - 1]; 
  } 
  for (i in 2:K_t_32) { 
    means_X_t_32[i - 1] = mean(X_t_32[, i]); 
    Xc_t_32[, i - 1] = X_t_32[, i] - means_X_t_32[i - 1]; 
  } 
  for (n in 1:N) {
    Y[n] = [Y_t_22[n], Y_t_24[n], Y_t_26[n], Y_t_28[n], Y_t_30[n], Y_t_32[n]]';
  }
} 
parameters { 
  vector[Kc_t_22] b_t_22;  // population-level effects 
  real temp_t_22_Intercept;  // temporary intercept 
  real<lower=0> sigma_t_22;  // residual SD 
  vector[Kc_t_24] b_t_24;  // population-level effects 
  real temp_t_24_Intercept;  // temporary intercept 
  real<lower=0> sigma_t_24;  // residual SD 
  vector[Kc_t_26] b_t_26;  // population-level effects 
  real temp_t_26_Intercept;  // temporary intercept 
  real<lower=0> sigma_t_26;  // residual SD 
  vector[Kc_t_28] b_t_28;  // population-level effects 
  real temp_t_28_Intercept;  // temporary intercept 
  real<lower=0> sigma_t_28;  // residual SD 
  vector[Kc_t_30] b_t_30;  // population-level effects 
  real temp_t_30_Intercept;  // temporary intercept 
  real<lower=0> sigma_t_30;  // residual SD 
  vector[Kc_t_32] b_t_32;  // population-level effects 
  real temp_t_32_Intercept;  // temporary intercept 
  real<lower=0> sigma_t_32;  // residual SD 
  // parameters for multivariate linear models 
  cholesky_factor_corr[nresp] Lrescor; 
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations 
  matrix[M_1, N_1] z_1;  // unscaled group-level effects 
  // cholesky factor of correlation matrix 
  cholesky_factor_corr[M_1] L_1; 
} 
transformed parameters { 
  // group-level effects 
  matrix[N_1, M_1] r_1 = (diag_pre_multiply(sd_1, L_1) * z_1)'; 
  vector[N_1] r_1_t_22_1 = r_1[, 1]; 
  vector[N_1] r_1_t_24_2 = r_1[, 2]; 
  vector[N_1] r_1_t_26_3 = r_1[, 3]; 
  vector[N_1] r_1_t_28_4 = r_1[, 4]; 
  vector[N_1] r_1_t_30_5 = r_1[, 5]; 
  vector[N_1] r_1_t_32_6 = r_1[, 6]; 
} 
model { 
  vector[N] mu_t_22 = Xc_t_22 * b_t_22 + temp_t_22_Intercept; 
  vector[N] mu_t_24 = Xc_t_24 * b_t_24 + temp_t_24_Intercept; 
  vector[N] mu_t_26 = Xc_t_26 * b_t_26 + temp_t_26_Intercept; 
  vector[N] mu_t_28 = Xc_t_28 * b_t_28 + temp_t_28_Intercept; 
  vector[N] mu_t_30 = Xc_t_30 * b_t_30 + temp_t_30_Intercept; 
  vector[N] mu_t_32 = Xc_t_32 * b_t_32 + temp_t_32_Intercept; 
  // multivariate linear predictor matrix 
  vector[nresp] Mu[N]; 
  vector[nresp] sigma = [sigma_t_22, sigma_t_24, sigma_t_26, sigma_t_28, sigma_t_30, sigma_t_32]';
  // cholesky factor of residual covariance matrix 
  matrix[nresp, nresp] LSigma = diag_pre_multiply(sigma, Lrescor); 
  for (n in 1:N) { 
    mu_t_22[n] = mu_t_22[n] + (r_1_t_22_1[J_1[n]]) * Z_1_t_22_1[n]; 
    mu_t_24[n] = mu_t_24[n] + (r_1_t_24_2[J_1[n]]) * Z_1_t_24_2[n]; 
    mu_t_26[n] = mu_t_26[n] + (r_1_t_26_3[J_1[n]]) * Z_1_t_26_3[n]; 
    mu_t_28[n] = mu_t_28[n] + (r_1_t_28_4[J_1[n]]) * Z_1_t_28_4[n]; 
    mu_t_30[n] = mu_t_30[n] + (r_1_t_30_5[J_1[n]]) * Z_1_t_30_5[n]; 
    mu_t_32[n] = mu_t_32[n] + (r_1_t_32_6[J_1[n]]) * Z_1_t_32_6[n]; 
    Mu[n] = [mu_t_22[n], mu_t_24[n], mu_t_26[n], mu_t_28[n], mu_t_30[n], mu_t_32[n]]';
  } 
  // priors including all constants 
  target += student_t_lpdf(temp_t_22_Intercept | 3, -5.71, 10); 
  target += student_t_lpdf(sigma_t_22 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10); 
  target += student_t_lpdf(temp_t_24_Intercept | 3, -5.53, 10); 
  target += student_t_lpdf(sigma_t_24 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10); 
  target += student_t_lpdf(temp_t_26_Intercept | 3, -5.33, 10); 
  target += student_t_lpdf(sigma_t_26 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10); 
  target += student_t_lpdf(temp_t_28_Intercept | 3, -5.12, 10); 
  target += student_t_lpdf(sigma_t_28 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10); 
  target += student_t_lpdf(temp_t_30_Intercept | 3, -4.89, 10); 
  target += student_t_lpdf(sigma_t_30 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10); 
  target += student_t_lpdf(temp_t_32_Intercept | 3, -4.79, 10); 
  target += student_t_lpdf(sigma_t_32 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10); 
  target += lkj_corr_cholesky_lpdf(Lrescor | 1); 
  target += student_t_lpdf(sd_1 | 3, 0, 10)
    - 6 * student_t_lccdf(0 | 3, 0, 10); 
  target += lkj_corr_cholesky_lpdf(L_1 | 1); 
  target += normal_lpdf(to_vector(z_1) | 0, 1); 
  // likelihood including all constants 
  if (!prior_only) { 
    target += multi_normal_cholesky_lpdf(Y | Mu, LSigma); 
  } 
} 
generated quantities { 
  // actual population-level intercept 
  real b_t_22_Intercept = temp_t_22_Intercept - dot_product(means_X_t_22, b_t_22); 
  // actual population-level intercept 
  real b_t_24_Intercept = temp_t_24_Intercept - dot_product(means_X_t_24, b_t_24); 
  // actual population-level intercept 
  real b_t_26_Intercept = temp_t_26_Intercept - dot_product(means_X_t_26, b_t_26); 
  // actual population-level intercept 
  real b_t_28_Intercept = temp_t_28_Intercept - dot_product(means_X_t_28, b_t_28); 
  // actual population-level intercept 
  real b_t_30_Intercept = temp_t_30_Intercept - dot_product(means_X_t_30, b_t_30); 
  // actual population-level intercept 
  real b_t_32_Intercept = temp_t_32_Intercept - dot_product(means_X_t_32, b_t_32); 
  matrix[nresp, nresp] Rescor = multiply_lower_tri_self_transpose(Lrescor); 
  vector<lower=-1,upper=1>[nrescor] rescor; 
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1); 
  vector<lower=-1,upper=1>[NC_1] cor_1; 
  // take only relevant parts of residual correlation matrix 
  rescor[1] = Rescor[1, 2]; 
  rescor[2] = Rescor[1, 3]; 
  rescor[3] = Rescor[2, 3]; 
  rescor[4] = Rescor[1, 4]; 
  rescor[5] = Rescor[2, 4]; 
  rescor[6] = Rescor[3, 4]; 
  rescor[7] = Rescor[1, 5]; 
  rescor[8] = Rescor[2, 5]; 
  rescor[9] = Rescor[3, 5]; 
  rescor[10] = Rescor[4, 5]; 
  rescor[11] = Rescor[1, 6]; 
  rescor[12] = Rescor[2, 6]; 
  rescor[13] = Rescor[3, 6]; 
  rescor[14] = Rescor[4, 6]; 
  rescor[15] = Rescor[5, 6]; 
  // take only relevant parts of correlation matrix 
  cor_1[1] = Cor_1[1,2]; 
  cor_1[2] = Cor_1[1,3]; 
  cor_1[3] = Cor_1[2,3]; 
  cor_1[4] = Cor_1[1,4]; 
  cor_1[5] = Cor_1[2,4]; 
  cor_1[6] = Cor_1[3,4]; 
  cor_1[7] = Cor_1[1,5]; 
  cor_1[8] = Cor_1[2,5]; 
  cor_1[9] = Cor_1[3,5]; 
  cor_1[10] = Cor_1[4,5]; 
  cor_1[11] = Cor_1[1,6]; 
  cor_1[12] = Cor_1[2,6]; 
  cor_1[13] = Cor_1[3,6]; 
  cor_1[14] = Cor_1[4,6]; 
  cor_1[15] = Cor_1[5,6]; 
} 