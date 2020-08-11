functions {
  // ODE for country-level transmission model --------------------------------------------------------------
    real[] ode(real t,
               real[] y, // [1] S I T
               real[] theta, // [4] beta, tau, nu, xi
               real[] x_r, // dummy
               int[] x_i // dummy
    ) {
      real dydt[5];
      real tau_t = theta[2] / ( 1 + exp(-theta[4]*(t+1-theta[3]))); // progressive introduction of ART using a logistic function
      dydt[1] = - theta[1] * y[1] * y[2] ;
      dydt[2] = theta[1] * y[1] * y[2] / sum(y[1:3]) - tau_t * y[2];
      dydt[3] = tau_t * y[2];
      // dummy: incidence 
      dydt[4] = theta[1] * y[1] * y[2];
      // dummy: AIDS deaths 
      dydt[5] = 0.0; 
      return(dydt);
    }
  // parallelized ODE solver  ----------------------------------------------------------------------
    vector country_ode(vector hypertheta, vector theta, real[] xrs, int[] xis) {
      int S = xis[1]; // simulation length
      int C = xis[2]; // number of compartments
      real y[S,C];
      vector[S*C] y_rect;
      y = integrate_ode_rk45(
        ode, 
        xrs[(S+2):(S+C+1)], // initial states
        xrs[1], // t0
        xrs[2:(S+1)], // evaluation dates (ts)  
        to_array_1d(theta), // parameters
        xrs, // real data: dummy
        xis, // integer data: dummy
        1.0E-8, 1.0E-8, 1.0E3); // controls from https://discourse.mc-stan.org/t/running-hierarchical-model-consisting-of-ode/2013/2
      for(i in 1:C) y_rect[((i-1)*S+1):(i*S)] = to_vector(y[,i]) + 1.0E-7; // add more than tolerance to avoid approximation leading to negative values
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
  // initial values
  real init_pop[G]; // initial population size
  real init_prev[G]; // initial prevalence
  int rollout[G];
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
  real p_sigma;
}
transformed data {
  int xis[G,3]; 
  real xrs[G,S+C+1];
  
  for(i in 1:G) {
    // integer data
    xis[i] = {S,C,rollout[i]};
    // real data
    xrs[i,1] = t0;
    xrs[i,2:(S+1)] = ts;
    xrs[i,(S+2):(S+C+1)] = {(init_pop[i]-init_prev[i])/init_pop[i],init_prev[i]/init_pop[i],0.0,0.0,0.0};
  } 
}
parameters {
  // parameters for the ODE system
  real<lower=0> beta_raw[G]; // transmission rate
  real<lower=0> tau_raw[G]; // treatment rate
  real<lower=0> nu_raw[G]; // time to 50% ART roll-out (after 2000)
  real<lower=0,upper=1> xi_raw[G]; // slope of logistic
  // variance and covariance parameters
  vector<lower=0>[E] sigma_raw[G];
}
transformed parameters {
  // transformed parameters
  real beta[G];
  real tau[G];
  real xi[G];
  real nu[G];
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
  matrix[L,E] output[G];
  // sqrt-transform indicator data and add uncertainty 
  matrix[L,E] indicator_data[G];
  
  // transform parameters
  for(i in 1:G) {
    beta[i] = beta_raw[i] / p_beta;
    tau[i] = tau_raw[i] / p_tau;
    xi[i] = 0.5 + xi_raw[i];
    nu[i] = nu_raw[i] / p_nu;
    sigma[i] = sigma_raw[i] / p_sigma;
 }
  // concatenate parameters
  for(i in 1:G) theta[i] = to_vector({beta[i],tau[i],nu[i],xi[i]});
  // solve ODE
  for(i in 1:G) y_rect[((i-1)*S*C+1):(S*C*i)] = country_ode(hypertheta, theta[i], xrs[i], xis[i]);
  // unpack outputs from rectangular shape
  for(i in 1:G) {
    for(j in 1:C) y[i,,j] = to_array_1d(y_rect[(((j-1)*S+1)+S*C*(i-1)):((j*S)+S*C*(i-1))]);
    // compute model outputs with sqrt transformation
    // (need to check for 0 values to prevent numerical failure in chain rule derivation)
    for(j in 1:L) {
      A_output[i,j] = sum(y[i,j,2:3]) * init_pop[i]; // prevalence (sum comps 2:3)
      B_output[i,j] = y[i,j,3] * init_pop[i]; // treatment (comp 3)
      C_output[i,j] = 0.0; // yearly AIDS-related mortality
      D_output[i,j] = sum(y[i,j,1:3]) * init_pop[i]; // population (sum comps 1:3)
    }
    // output to be fitted into matrix format
    output[i,,1] = sqrt(A_output[i,1:L]);
    output[i,,2] = sqrt(B_output[i,1:L]);    
    // sqrt-transform indicator data
    indicator_data[i,,1] = sqrt(A_est[i]);
    indicator_data[i,,2] = sqrt(B_est[i]);
  }
}
model {
  // priors
  for(i in 1:G) {  
    beta_raw[i] ~ exponential(1);
    tau_raw[i] ~ exponential(1);
    nu_raw[i] ~ exponential(1);
    xi_raw[i] ~ beta(p_xi,p_xi);
    sigma_raw[i] ~ exponential(1);
  }
  // inference
  if(inference==1) {
    for(i in 1:G) for(j in 1:E) target += normal_lpdf( indicator_data[i,,j] | output[i,,j], sigma[i,j]) ;
  }
}
generated quantities {
  vector[L] A_pred[G];
  vector[L] B_pred[G];
  vector[L] C_pred[G];
  vector[L] D_pred[G];
  vector[L] tau_t[G];
  matrix[L,E] log_lik[G];
  
  // predicted data
  for(i in 1:G) {
    for(j in 1:L) {
      A_pred[i,j] = square(normal_rng(output[i,j,1], sigma[i,1]));
      B_pred[i,j] = square(normal_rng(output[i,j,2], sigma[i,2]));
      C_pred[i,j] = C_output[i,j];    
      D_pred[i,j] = D_output[i,j];    
      tau_t[i,j] = j+1 < rollout[i] ? 0 : tau[i];
      for(e in 1:E) log_lik[i,j,e] =  normal_lpdf( indicator_data[i,j,e] | output[i,j,e], sigma[i,e]) ;
    }
  }
}
