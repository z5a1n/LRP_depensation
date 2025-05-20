functions {
   real c_prior_log(real c) { //used in Perala and Kuparinen 2017. Modified for c over [0,5]
      return log(1+0.62269*((1-2001^(2*c-2))/(1+2001^(2*c-2))));
   }
}
data
{
   int Y;               // Number of years
   vector[Y] R;         // Number of recruits
   vector[Y] S;         // Spawning stock biomass
   real kLow;           // Constraints of the parameters
   real kUp;
   real SkLow;
   real SkUp;
   real cLow;
   real cUp;
   real sigmaLow;
   real sigmaUp;
}
parameters
{
   real<lower=kLow,upper=kUp> k;
   real<lower=SkLow,upper=SkUp> Sk;
   real<lower=cLow,upper=cUp> c;
   real<lower=sigmaLow,upper=sigmaUp> sigma;
}
transformed parameters
{
   vector[Y] mu;
   vector[Y] ER;

   for (y in 1:Y) {
      ER[y] = k * pow(S[y]/Sk,c) * exp(c*(1-S[y]/Sk));
   }
   mu = log(ER) - 0.5*sigma^2;
}
model
{
   // Priors
   k ~ uniform(kLow,kUp);
   Sk ~ uniform(SkLow,SkUp);
   c ~ c_prior();
   sigma ~ uniform(sigmaLow,sigmaUp);

   // Likelihood
   R ~ lognormal(mu, sigma);
}
generated quantities
{
   vector[Y] log_lik;

   for (y in 1:Y) {
      log_lik[y] = lognormal_lpdf(R[y]| mu[y], sigma);
   }
}
