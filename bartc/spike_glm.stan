data {
    int<lower = 0> N; // number of observations
    int<lower = 0> P; // number of regressors
    int y[N]; // response variable
    matrix[N, P] X; // predictor variable
}

parameters {
    real b0; // intercept
    vector[P] beta; // slope
}

transformed parameters {
    vector[N] log_lambda; // Poisson rate per bin
    log_lambda <- b0 + X * beta;
}
 
model {
    b0 ~ normal(0, 2.5);
    for (idx in 1:P) {
        beta[idx] ~ student_t(1, 0, 0.01);
    }
    y ~ poisson_log(log_lambda);
}
