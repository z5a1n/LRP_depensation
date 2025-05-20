data
{
   int Y;               // Number of years
   vector[Y] R;         // Number of recruits
   vector[Y] S;         // Spawning stock biomass
   real kLow;           // Constraints of the parameters
   real kUp;
   real SkLow;
   real SkUp;
   real sigmaLow;
   real sigmaUp;
}
parameters
{
   real<lower=kLow,upper=kUp> k;
   real<lower=SkLow,upper=SkUp> Sk;
   real<lower=sigmaLow,upper=sigmaUp> sigma;
}
transformed parameters
{
   vector[Y] mu;
   vector[Y] ER;

   for (y in 1:Y) {
      ER[y] = k * S[y]/Sk * exp(1-S[y]/Sk);
   }
   mu = log(ER) - 0.5*sigma^2;
}
model
{
   // Priors
   k ~ uniform(kLow,kUp);
   Sk ~ uniform(SkLow,SkUp);
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
