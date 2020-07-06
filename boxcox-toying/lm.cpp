#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);
  DATA_VECTOR(Y);
  PARAMETER(alpha);
  PARAMETER(beta);
  PARAMETER(logSigma);
  Type sigma = exp(logSigma);
  PARAMETER(logLambda);
  Type lambda = exp(logLambda);
  vector<Type> Ytrans=(pow(Y,lambda)-1)/lambda;
  Type nll = -sum(dnorm(Ytrans, alpha*x+beta, sigma, true));
  nll += -sum(vector<Type>(log(pow(Y,lambda-Type(1)))));
  ADREPORT(lambda);
  return nll;
}
