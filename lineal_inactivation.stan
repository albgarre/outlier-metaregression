//
//

data {
  int<lower=0> N;
  vector[N] temperature;
  vector[N] logD;
  real refTemp;
}

parameters {
  real logDref;
  real<lower = 0> z;
  real<lower=0> sigma;
}


model {
  
  z ~ normal(6, 2);
  logDref ~ normal(0, 2);
  sigma ~ exponential(2);
  
  
  for (n in 1:N)
    logD[n] ~ normal(logDref - (temperature[n] - refTemp)/z, sigma);
}

