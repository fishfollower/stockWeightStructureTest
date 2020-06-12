#include <TMB.hpp>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(Wr)
  DATA_MATRIX(Wc)
  DATA_MATRIX(Wd)
  DATA_MATRIX(Y)

  matrix<Type> pred(Y.rows(),Y.cols());
    
  PARAMETER_VECTOR(logPhi)
  PARAMETER_VECTOR(mu)
  PARAMETER(logSdProc)   
  PARAMETER(logSdObs)   
  PARAMETER_MATRIX(z)
  vector<Type> phi=exp(logPhi);
    
  matrix<Type> I(Wr.rows(),Wr.cols());
  I.setIdentity();

  matrix<Type> Q=I-phi(0)*Wr-phi(1)*Wc-phi(2)*Wd;
    
  using namespace density;
  Type nll = SCALE(GMRF(asSparseMatrix(Q)),exp(logSdProc))(z.vec());

  for(int i=0; i<Y.rows(); ++i){
    for(int j=0; j<Y.cols(); ++j){
      pred(i,j)=z(i,j)+mu(j);
      if(!isNA(Y(i,j))){
        nll += -dnorm(Y(i,j),pred(i,j),exp(logSdObs),true);
      }
    }
  }
  ADREPORT(pred);
  return(nll);
}
