//
//

data {
  int<lower=0> N;
  real U;
  real L;
  vector[N] temperature;
  vector<lower=L,upper=U>[N] logD;
  real refTemp;
}

parameters {
  real logDref;
  real<lower = 0> z;
  real<lower=0> sigma;
}


model {
  
  z ~ normal(6, 1);
  logDref ~ normal(0, 2);
  sigma ~ exponential(2);
  
  
  for (n in 1:N)
    logD[n] ~ normal(logDref - (temperature[n] - refTemp)/z, sigma) T[L,U];
}
