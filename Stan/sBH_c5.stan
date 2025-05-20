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
   real RinfLow;        // Constraints of the parameters
   real RinfUp;
   real SqLow;
   real SqUp;
   real cLow;
   real cUp;
   real sigmaLow;
   real sigmaUp;
}
parameters
{
   real<lower=RinfLow, upper=RinfUp> Rinf;
   real<lower=SqLow, upper=SqUp> Sq;
   real<lower=cLow, upper=cUp> c;
   real<lower=sigmaLow, upper=sigmaUp> sigma;
}
transformed parameters
{
   vector[Y] mu;
   vector[Y] ER;

   for (y in 1:Y) {
      ER[y] = Rinf / (pow(Sq / S[y] , c) + 1);
   }
   mu = log(ER) - 0.5*sigma^2;
}
model
{
   // Priors
   Rinf ~ uniform(RinfLow, RinfUp);
   Sq ~ uniform(SqLow, SqUp);
   c ~ c_prior();
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

