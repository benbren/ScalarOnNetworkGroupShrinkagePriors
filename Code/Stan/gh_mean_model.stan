data{
  int Q; // number of possible connections 
  int G; // number of groups 
  int group_sizes[G]; // group_sizes 
  int sum_gs[G];
  vector[Q] y; // outcome 
}

parameters{
  vector [Q] beta; // regression coefficients ordered by group 
  real logsigma2; 
  real <lower=0> tau2; 
  real <lower=0> at; 
  vector <lower=0> [Q] al; 
  vector <lower=0> [Q]lam;
  vector <lower=0> [G] ag; 
  vector <lower=0> [G] gam;
}

transformed parameters{
  real sigma2;  
  vector[Q] beta_vars; 
  sigma2 = exp(logsigma2);
  
  for (g in 1:G){
    for (gs in 1:group_sizes[g]){
      beta_vars[sum_gs[g] + gs] = tau2*sigma2*lam[sum_gs[g] + gs]*gam[g] ; 
    }
  }
  
 }

model{ 
  
 // likelihood  
 y ~ normal(beta, sigma2);
 
 // priors 
 at ~ inv_gamma(0.5, 1) ;
 tau2 ~ inv_gamma(0.5, 1/at) ; 
  
 for (g in 1:G){
    ag[g] ~ inv_gamma(0.5, 1) ;
    gam[g] ~ inv_gamma(0.5, 1/ag[g]) ; 
  }

  
 for (q in 1:Q){
    al[q] ~ inv_gamma(0.5,1) ;
    lam[q] ~ inv_gamma(0.5, 1/al[q]);
  }

 beta ~ normal(0,beta_vars);
  
}

generated quantities{
  
  vector[Q] w;
  
   for (g in 1:G){
    for (gs in 1:group_sizes[g]){
      w[sum_gs[g] + gs] = 1 - 1/(1+tau2*lam[sum_gs[g] + gs]*gam[g]) ; 
    }
  }
  
 }
  

