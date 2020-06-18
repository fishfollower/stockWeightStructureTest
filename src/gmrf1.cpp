#include <TMB.hpp>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(mode)
  Type nll = 0;
    
  if(mode==1){ // GMRF with neighbors in year, age, and cohort direction  
    DATA_MATRIX(Wr)
    DATA_MATRIX(Wc)
    DATA_MATRIX(Wd)
    DATA_ARRAY(Y)
    DATA_ARRAY_INDICATOR(keep,Y);
    matrix<Type> pred(Y.dim[0],Y.dim[1]);
    
    PARAMETER_VECTOR(logPhi)
    PARAMETER_VECTOR(mu)
    PARAMETER(logSdProc)   
    PARAMETER_VECTOR(logSdObs)   
    PARAMETER_MATRIX(z)
    vector<Type> phi=exp(logPhi);
    
    matrix<Type> I(Wr.rows(),Wr.cols());
    I.setIdentity();

    matrix<Type> Q=I-phi(0)*Wr-phi(1)*Wc-phi(2)*Wd;
    
    using namespace density;
    nll += SCALE(GMRF(asSparseMatrix(Q)),exp(logSdProc))(z.vec());

    for(int i=0; i<Y.dim[0]; ++i){
      for(int j=0; j<Y.dim[1]; ++j){
        pred(i,j)=z(i,j)+mu(j);
        if(!isNA(Y(i,j))){
          nll += -dnorm(Y(i,j),pred(i,j),exp(logSdObs(j)),true)*keep(i,j);
        }
      }
    }
    ADREPORT(pred);
  }

  
  if(mode==2){ // AR(time)xAR(age)+added cohort effect  
    DATA_MATRIX(Wr)
    DATA_MATRIX(Wc)
    DATA_MATRIX(Wd)
    DATA_ARRAY(Y)  
    DATA_ARRAY_INDICATOR(keep,Y);
    
    matrix<Type> pred(Y.dim[0],Y.dim[1]);
    
    PARAMETER_VECTOR(logitRho)
    PARAMETER_VECTOR(mu)
    PARAMETER_VECTOR(logSdProc)
      
    PARAMETER_VECTOR(logSdObs)   
    PARAMETER_ARRAY(omega)
    PARAMETER_VECTOR(z)      
    vector<Type> rho=invlogit(logitRho);
    using namespace density;
    nll += SCALE(SEPARABLE(AR1(rho(1)),AR1(rho(0))),exp(logSdProc(0)))(omega);
    nll += SCALE(AR1(rho(2)),exp(logSdProc(1)))(z);
    for(int i=0; i<Y.dim[0]; ++i){
      for(int j=0; j<Y.dim[1]; ++j){
        pred(i,j)=mu(j)+omega(i,j)+z(i-j+(Y.dim[1]-1));
        if(!isNA(Y(i,j))){
          nll += -dnorm(Y(i,j),pred(i,j),exp(logSdObs(j)),true)*keep(i,j);
        }
      }
    }
    ADREPORT(pred);
  }
  return(nll);
}
