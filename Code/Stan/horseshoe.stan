data{
  int Q; // number of possible connections 
  vector[Q] y; // outcome 
}

parameters{
  vector [Q] beta; // regression coefficients ordered by group 
  real logsigma2; 
  real <lower=0> tau2; 
  real <lower=0> at; 
  vector <lower=0> [Q] al; 
  vector <lower=0> [Q]lam;

}

transformed parameters{
  real sigma2;  
  sigma2 = exp(logsigma2);
  
} 

model{ 
  
 // likelihood  
 y ~ normal(beta, sigma2);
 
 // priors 
 at ~ inv_gamma(0.5, 1) ;
 tau2 ~ inv_gamma(0.5, 1/at) ; 

 for (q in 1:Q){
    al[q] ~ inv_gamma(0.5,1) ;
    lam[q] ~ inv_gamma(0.5, 1/al[q]);
  }

 for (q in 1:Q){
   beta[q] ~ normal(0,tau2*sigma2*lam[q]);
 } 
  
}

generated quantities{
  vector[Q] w;
   for (q in 1:Q){
      w[q] = 1 - 1/(1+tau2*lam[q]) ; 
  }
 }
  