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


  if(mode==3){ // as mode 2, but with correlated obs errors

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
    PARAMETER(logitRhoObs)
 
    vector<Type> rho=invlogit(logitRho);
    Type rhoObs = invlogit(logitRhoObs);
    vector<Type> sdObs = exp(logSdObs);
    using namespace density;
    int nA = Y.dim[1];
    
    matrix<Type> covObs(nA,nA);
    for(int i=0;i<nA;i++)
      for(int j=0;j<nA;j++){
	covObs(i,j)=pow(rhoObs,Type(abs(i-j)))*sdObs[i]*sdObs[j];
      }
    
    nll += SCALE(SEPARABLE(AR1(rho(1)),AR1(rho(0))),exp(logSdProc(0)))(omega);
    nll += SCALE(AR1(rho(2)),exp(logSdProc(1)))(z);

    MVNORM_t<Type> neg_log_density(covObs);
    
    vector<Type> obsvec(nA); 
    vector<Type> predvec(nA);
    data_indicator<vector<Type>,Type> keepvec(nA);
    
    for(int i=0; i<Y.dim[0]; ++i){
      for(int j=0; j<nA; ++j){
        pred(i,j)=mu(j)+omega(i,j)+z(i-j+(nA-1));
	obsvec(j) = Y(i,j);
	predvec(j) = pred(i,j);
	keepvec(j) = keep(i,j);
      }
      if(!isNA(Y(i,0))){ 
	nll += neg_log_density( obsvec - predvec, keepvec);
      }
    }
    ADREPORT(pred);
  }

  if(mode==4){ // as mode 3, but where cohorts and only increase in weight

    DATA_MATRIX(Wr)
    DATA_MATRIX(Wc)
    DATA_MATRIX(Wd)
    DATA_ARRAY(Y)  
    DATA_ARRAY_INDICATOR(keep,Y);
    DATA_SCALAR(trans); // 0 = log, >0 = x^(trans)
    matrix<Type> pred(Y.dim[0],Y.dim[1]);
    
    PARAMETER_VECTOR(logitRho)
    PARAMETER_VECTOR(mu)
    PARAMETER_VECTOR(logSdProc)
      
    PARAMETER_VECTOR(logSdObs)   
    PARAMETER_ARRAY(omega)
    PARAMETER_VECTOR(z)
    PARAMETER(logitRhoObs)
 
    vector<Type> rho=invlogit(logitRho);
    Type rhoObs = invlogit(logitRhoObs);
    vector<Type> sdObs = exp(logSdObs);
    using namespace density;
    int nY = Y.dim[0];
    int nA = Y.dim[1];
    
    matrix<Type> covObs(nA,nA);
    for(int i=0;i<nA;i++)
      for(int j=0;j<nA;j++){
	covObs(i,j)=pow(rhoObs,Type(abs(i-j)))*sdObs[i]*sdObs[j];
      }

    SCALE_t< SEPARABLE_t<AR1_t<N01<Type> > , AR1_t<N01<Type> > > > d1=SCALE(SEPARABLE(AR1(rho(1)),AR1(rho(0))),exp(logSdProc(0)));
    nll += d1(omega);
    SIMULATE{
      SEPARABLE(AR1(rho(1)),AR1(rho(0))).simulate(omega);
      omega*=exp(logSdProc(0));
    }
    SCALE_t<AR1_t<N01<Type> > >  d2=SCALE(AR1(rho(2)),exp(logSdProc(1)));
    nll += d2(z);
    SIMULATE{
      d2.simulate(z);
    }
    
    matrix<Type> dW(nY,nA);
    dW.setZero();
        
    for(int i=0; i<nY; i++){
      for(int j=0; j<nA; j++){
	dW(i,j) = exp( omega(i,j) + mu(j) + z(i-j+(nA-1)) );
      }
    }
    
    for(int i=0; i<nY; i++) pred(i,0) = exp(mu(0)); 
    for(int j=1; j<nA; j++) pred(0,j) = pred(0,j-1) + exp(mu(j));
    
    for(int i=0; i<nY; i++) pred(i,0) = pred(i,0)*exp(omega(i,0) + z(i+nA-1));
    for(int j=1; j<nA; j++) pred(0,j) = pred(0,j)*exp(omega(0,j) + z(-j+nA-1));
    
    for(int i=1; i<nY; i++){
      for(int j=1; j<nA; j++){ 
	pred(i,j) = pred(i-1,j-1) + dW(i,j);
      }
    }
      
    MVNORM_t<Type> neg_log_density(covObs);
    
    vector<Type> obsvec(nA); 
    vector<Type> predvec(nA);
    data_indicator<vector<Type>,Type> keepvec(nA);
    
    for(int i=0; i<Y.dim[0]; ++i){
      for(int j=0; j<nA; ++j){
	obsvec(j) = Y(i,j);
	if(trans<=0){
	  pred(i,j) = log(pred(i,j));
	} else {
	  pred(i,j) = pow(pred(i,j),trans);
	}
	predvec(j) = pred(i,j); 
	keepvec(j) = keep(i,j);
      }
      if(!isNA(Y(i,0))){ 
	nll += neg_log_density( obsvec - predvec, keepvec);
	SIMULATE{
          obsvec=neg_log_density.simulate()+predvec;
	  for(int j=0; j<nA; ++j){
            Y(i,j)=obsvec(j);
	  }  
	}
      }
    }
    ADREPORT(pred);
    
    SIMULATE{
      REPORT(Y)
      REPORT(omega)
      REPORT(z)
    }  
  }

  
  
  return(nll);
}
