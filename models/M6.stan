functions {
  // ODE for country-level transmission model --------------------------------------------------------------
  real[] ode(real t,
              real[] y, //  [5] S I T J U + 4 dummies
              real[] theta, // [9] beta, tau, nu, xi, eta, delta, omega, iota, kappa
              real[] x_r, // [S+4] life_expectancy
              int[] x_i // dummy
  ) {
    real dydt[11];
    real tau_t = theta[2] / ( 1 + exp(-theta[4]*(t+1-theta[3]))); // progressive introduction of ART using a logistic function
    real mu = 1.0/(x_r[x_i[1]+4]-15); // approximate life expectancy for adults aged 15+ by life expectancy at birth - 15
    dydt[1] = theta[5] - theta[1]*y[1]*(y[2]+y[4]+y[5]) / sum(y[1:6]) - mu*y[1];
    dydt[2] = theta[1]*y[1]*y[2] / sum(y[1:6]) - tau_t*y[2] - theta[6]*y[2] - mu*y[2];
    dydt[3] = tau_t*y[2] - theta[7]*y[3] - mu*y[3];
    dydt[4] = theta[1]*y[1]*(y[4]+y[5]) / sum(y[1:6]) - tau_t*y[4] - theta[6]*y[4] - mu*y[4];
    dydt[5] = tau_t*y[4] + theta[7]*y[3] - theta[6]*y[5] - mu*y[5] - theta[9]*y[5];
    dydt[6] = theta[9]*y[5] - mu*y[6];
    // dummy: incidence 
    dydt[7] = theta[1] * y[1] * (y[2]+y[4]+y[5]) / sum(y[1:6]);
    // dummy: AIDS deaths 
    dydt[8] = theta[6] * (y[2]+y[4]+y[5]);
    // dummy: art initiators
    dydt[9] = tau_t * (y[2] + y[4]);
    // dummy: resistance among art initiators
    dydt[10] = tau_t * y[4];
    // dummy: resistance among incidence
    dydt[11] = theta[1] * y[1] * (y[4]+y[5]) / sum(y[1:6]);
    return(dydt);
  }
  // parallelized ODE solver  ----------------------------------------------------------------------
  vector country_ode(vector hypertheta, vector theta, real[] xrs, int[] xis) {
    int S = xis[1]; // simulation length
    int C = xis[2]; // number of compartments
    real y[S,C];
    vector[S*C] y_rect;
    real init[C] = {(xrs[S+2]-xrs[S+3])/xrs[S+2],xrs[S+3]*(1-theta[8])/xrs[S+2],0.0,xrs[S+3]*theta[8]/xrs[S+2],0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    y = integrate_ode_rk45(
        ode, 
        init, // initial states
        xrs[1], // t0
        xrs[2:(S+1)], // evaluation dates (ts)  
        to_array_1d(theta), // parameters
        xrs, // real data: dummy
        xis, // integer data: dummy
        1.0E-8, 1.0E-8, 1.0E3); // controls from https://discourse.mc-stan.org/t/running-hierarchical-model-consisting-of-ode/2013/2
    for(i in 1:C) y_rect[((i-1)*S+1):(i*S)] = to_vector(y[,i]) + 1.0E-7; // removes any negative value due to numerical approximation of small values 
    return y_rect;
  }
}
data {
  // controls
  int inference;
  // data structure
  int t0; // initial time
  int S; // duration of simulation (in years)
  real ts[S]; // evaluation dates (yearly)
  int C; // number of compartments
  int D; // number of parameters
  int E; // number of fitted compartments
  int G; // number of countries
  int L; // number of country-level data points from start_year (L<=S)
  int M; // number of surveys
  int N; // number of suppression data points
  // country data
  real init_pop[G]; // initial population size
  real init_prev[G]; // initial prevalence
  int rollout[G]; // year of ART rollout
  real life_expectancy[G]; // life_expectancy at birth
  // indicator data (A-E)
  vector[L] A_est[G]; // prevalence
  vector[L] B_est[G]; // treatment
  vector[L] C_est[G]; // mortality
  vector[L] D_est[G]; // population
  vector[L] E_est[G]; // pmtct
  // survey data (F)
  int F_country[M]; // country indicator
  int F_n[M]; // sample size
  int F_k[M]; // number with NNRTI resistance
  int F_t[M]; // median date of sampling
  // priors
  real p_beta;
  real p_tau;
  real p_nu;
  real p_xi;
  real p_eta;
  real p_delta;
  real p_kappa;
  real p_sigma;
  real p_mu_omega;
  real p_sigma_omega;
  real p_mu_iota[2];
  real p_sigma_iota;
}
transformed data {
  int xis[G,2]; 
  real xrs[G,S+4];
  
  for(i in 1:G) {
    // integer data
    xis[i] = {S,C};
    // real data
    xrs[i,1] = t0;
    xrs[i,2:(S+1)] = ts;
    xrs[i,(S+2):(S+3)] = {init_pop[i],init_prev[i]};
    xrs[i,S+4] = life_expectancy[i];
  } 
}
parameters {
  // parameters for the ODE system
  real<lower=0> beta_raw[G]; // transmission rate
  real<lower=0> tau_raw[G]; // treatment rate
  real<lower=0> nu_raw[G]; // time to 50% ART roll-out (after 2000)
  real<lower=0,upper=1> xi_raw[G]; // slope of logistic
  real<lower=0> eta_raw[G];
  real<lower=0> delta_raw[G];
  real<lower=0> kappa_raw[G];
  real<lower=0> mu_omega;
  real<lower=0> sigma_omega;
  real<lower=0> omega[G];
  real<lower=0,upper=1> mu_iota;
  real<lower=0> sigma_iota;
  real iota_raw[G];
  // variance and covariance parameters
  vector<lower=0>[E] sigma_raw[G];
}
transformed parameters {
  // transformed parameters
  real beta[G];
  real tau[G];
  real xi[G];
  real nu[G];
  real eta[G];
  real delta[G];
  real kappa[G];
  real<lower=0,upper=1> iota[G];
  vector[E] sigma[G];
  // formatted objects for map_rect
  vector[0] hypertheta;
  vector[D] theta[G];
  vector[G*S*C] y_rect;
  real y[G,S,C];
  // model outputs
  vector[L] A_output[G]; 
  vector[L] B_output[G];
  vector[L] C_output[G];
  vector[L] D_output[G];
  vector[L] F_output[G];
  vector[M] F_study_logit; // proportion of resistance among ART initiators by study (logit-transformed)
  matrix[L,E] output[G];
  // sqrt-transform indicator data
  matrix[L,E] indicator_data[G];
  
  // transform parameters
  for(i in 1:G) {
    beta[i] = beta_raw[i] / p_beta;
    tau[i] = tau_raw[i] / p_tau;
    xi[i] = 0.5 + xi_raw[i];
    nu[i] = nu_raw[i] / p_nu;
    eta[i] = eta_raw[i] / p_eta;
    delta[i] = delta_raw[i] / p_delta;
    kappa[i] = kappa_raw[i] / p_kappa;
    iota[i] = inv_logit(iota_raw[i]);
    sigma[i] = sigma_raw[i] / p_sigma;
  }
  // concatenate parameters
  for(i in 1:G) theta[i] = to_vector({beta[i],tau[i],nu[i],xi[i],eta[i],delta[i],omega[i],iota[i],kappa[i]});
  // solve ODE
  y_rect = map_rect(country_ode,hypertheta,theta,xrs,xis);
  // unpack outputs from rectangular shape
  for(i in 1:G) {
    for(j in 1:C) y[i,,j] = to_array_1d(y_rect[(((j-1)*S+1)+S*C*(i-1)):((j*S)+S*C*(i-1))]);
    // compute model outputs with sqrt transformation
    for(j in 1:L) {
      A_output[i,j] = sum(y[i,j,2:6]) * init_pop[i]; // prevalence (sum comps 2:6)
      B_output[i,j] = (y[i,j,3]+y[i,j,5]+y[i,j,6]) * init_pop[i]; // treatment (sum comps 3+5+6)
      C_output[i,j] = (y[i,j,8] - (j==1 ? 0 : y[i,j-1,8])) <= 0 ? 0.0 : (y[i,j,8] - (j==1 ? 0 : y[i,j-1,8])) * init_pop[i]; // yearly AIDS-related mortality (lag comp 8, insure no negative values)
      D_output[i,j] = sum(y[i,j,1:6]) * init_pop[i]; // population (sum comps 1:6)
      F_output[i,j] = y[i,j,4] / (y[i,j,2]+y[i,j,4]); // proportion of resistance among untreated (comp 4 / comps 2+4)
    }
    // output to be fitted into matrix format
    output[i,,1] = sqrt(A_output[i,1:L]);
    output[i,,2] = sqrt(B_output[i,1:L]);    
    output[i,,3] = sqrt(C_output[i,1:L]);   
    output[i,,4] = sqrt(D_output[i,1:L]);   
    // sqrt-transform indicator data
    indicator_data[i,,1] = sqrt(A_est[i]);
    indicator_data[i,,2] = sqrt(B_est[i]);
    indicator_data[i,,3] = sqrt(C_est[i]);
    indicator_data[i,,4] = sqrt(D_est[i]);
  }
  for(i in 1:M) F_study_logit[i] = logit(F_output[F_country[i],F_t[i]]);
}
model {
  // priors
  for(i in 1:G) {  
    beta_raw[i] ~ exponential(1);
    tau_raw[i] ~ exponential(1);
    nu_raw[i] ~ exponential(1);
    xi_raw[i] ~ beta(p_xi,p_xi);
    eta_raw[i] ~ exponential(1);
    delta_raw[i] ~ exponential(1);
    kappa_raw[i] ~ exponential(1);
    omega[i] ~ gamma(mu_omega^2/sigma_omega^2,mu_omega/sigma_omega^2);
    iota_raw[i] ~ normal(logit(mu_iota),sigma_iota);
    sigma_raw[i] ~ exponential(1);
  }
  mu_omega ~ exponential(p_mu_omega);
  sigma_omega ~ exponential(p_sigma_omega);
  mu_iota ~ beta(p_mu_iota[1],p_mu_iota[2]);
  sigma_iota ~ exponential(p_sigma_iota);
  // inference
  if(inference==1) {
    for(i in 1:G) for(j in 1:E) target += normal_lpdf( indicator_data[i,,j] | output[i,,j], sigma[i,j]) ;
    target += binomial_logit_lpmf(F_k | F_n, F_study_logit) ;
  }
}
generated quantities {
  vector[L] A_pred[G];
  vector[L] B_pred[G];
  vector[L] C_pred[G];
  vector[L] D_pred[G];
  vector[L] tau_t[G];
  matrix[L,E] log_lik[G];
  vector[M] log_lik2;
  
  // predicted data
  for(i in 1:G) {
    for(j in 1:L) {
      A_pred[i,j] = square(normal_rng(output[i,j,1], sigma[i,1]));
      B_pred[i,j] = square(normal_rng(output[i,j,2], sigma[i,2]));
      C_pred[i,j] = square(normal_rng(output[i,j,3], sigma[i,3]));
      D_pred[i,j] = square(normal_rng(output[i,j,4], sigma[i,4]));
      tau_t[i,j] =  tau[i] / (1+exp(-xi[i]*(j-nu[i])));
      for(e in 1:E) log_lik[i,j,e] =  normal_lpdf( indicator_data[i,j,e] | output[i,j,e], sigma[i,e]) ;
    }
  }
  for(i in 1:M) log_lik2[i] = binomial_logit_lpmf(F_k[i] | F_n[i], F_study_logit[i]);
}
