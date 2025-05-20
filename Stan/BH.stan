data
{
   int Y;               // Number of years
   vector[Y] R;         // Number of recruits
   vector[Y] S;         // Spawning stock biomass
   real RinfLow;        // Constraints of the parameters
   real RinfUp;
   real SqLow;
   real SqUp;
   real sigmaLow;
   real sigmaUp;
}
parameters
{
   real<lower=RinfLow, upper=RinfUp> Rinf;
   real<lower=SqLow, upper=SqUp> Sq;
   real<lower=sigmaLow, upper=sigmaUp> sigma;
}
transformed parameters
{
   vector[Y] mu;
   vector[Y] ER;

   for (y in 1:Y) {
      ER[y] = Rinf / ((Sq / S[y]) + 1);
   }
   mu = log(ER) - 0.5*sigma^2;
}
model
{
   // Priors
   Rinf ~ uniform(RinfLow, RinfUp);
   Sq ~ uniform(SqLow, SqUp);
   sigma ~ uniform(sigmaLow, sigmaUp);

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

