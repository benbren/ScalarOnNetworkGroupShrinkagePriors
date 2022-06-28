write("data{
  int N; // number of subjects 
  int Q; // number of possible connections 
  int G; // number of groups 
  int group_sizes[G]; // group_sizes 
  int sum_gs[G];
  vector[N] y; // outcome 
  matrix[N,Q] X; //matrix of lower triangular networks 
}

parameters{
  vector [Q] beta; // regression coefficients ordered by group 
  real <lower=0> logsigma2; 
  real <lower=0> tau2; 
  real <lower=0> at; 
  vector <lower=0> [Q] al; 
  vector<lower=0> [Q]lam;
  vector <lower=0> [G] ag; 
  vector<lower=0> [G] gam;
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

  at ~ inv_gamma(0.5, 1) ;
  tau2 ~ inv_gamma(0.5, 1/at) ; 
  
  //for (g in 1:G){
    //ag[g] ~ inv_gamma(0.5, 1) ;
    //gam[g] ~ inv_gamma(0.5, 1/ag[g]) ; 
  //}

  ag ~ inv_gamma(0.5, 1)
  gam ~ inv_gamma(0.5, 1/ag) 
  
  //for (n in 1:N){
    //al[n] ~ inv_gamma(0.5,1) ;
    //lam[n] ~ inv_gamma(0.5, 1/al[n]) ;
  //}

  al ~ inv_gamma(0.5, 1)
  lam ~ inv_gamma(0.5, 1/al)
  
  //for(q in 1:Q){
    //beta[q] ~ normal(0,beta_vars[q]) ;
  //}

  beta ~ normal(0,beta_vars)

  y ~ normal(beta*x, sigma2)
  
  //for (n in 1:N){
    //y[n] ~ normal(X[n,]*beta, sigma2) ;
  //}
  
}

", 
"group_horseshoe.stan")
