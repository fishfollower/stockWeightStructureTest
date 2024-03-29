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

  if(mode==0){  
    DATA_MATRIX(Wr)
    DATA_MATRIX(Wc)
    DATA_MATRIX(Wd)
    DATA_ARRAY(Y)
    DATA_ARRAY_INDICATOR(keep,Y);
    DATA_INTEGER(aveYears);
    matrix<Type> pred(Y.dim[0],Y.dim[1]);
    
    PARAMETER_VECTOR(logSdObs)
    Type s;
    int c;
    for(int i=0; i<Y.dim[0]; ++i){
      for(int j=0; j<Y.dim[1]; ++j){
	s=0;
	c=0;
        if(i<aveYears){
          for(int ii=0; ii<aveYears; ++ii){
	    if(!isNA(Y(ii,j))){
  	      s+=Y(ii,j);
	      ++c;
	    }
	  }
	}else{
          for(int ii=(i-aveYears); ii<i; ++ii){
	    if(!isNA(Y(ii,j))){
	      s+=Y(ii,j);
              ++c;
	    }
	  }
	}
	pred(i,j)=s/c;
        if(!isNA(Y(i,j))){
          nll += -dnorm(Y(i,j),pred(i,j),exp(logSdObs(j)),true)*keep(i,j);
        }
      }
    }
    ADREPORT(pred);
  }
  
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

  if(mode==4){ // as mode 3, but where cohorts can only increase in weight

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

  if(mode==5){ // as model 11, but with correlated errors

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
    PARAMETER(logitRhoObs)  
    vector<Type> phi=exp(logPhi);
    
    matrix<Type> I(Wr.rows(),Wr.cols());
    I.setIdentity();

    matrix<Type> Q=I-phi(0)*Wc-phi(1)*Wd;
    
    using namespace density;
    nll += SCALE(GMRF(asSparseMatrix(Q)),exp(logSdProc))(z.vec());

    Type rhoObs = invlogit(logitRhoObs);
    vector<Type> sdObs = exp(logSdObs) +1.0e-5;
    int nA = Y.dim[1];
    matrix<Type> covObs(nA,nA);
    for(int i=0;i<nA;i++)
      for(int j=0;j<nA;j++){
	covObs(i,j)=pow(rhoObs,Type(abs(i-j)))*sdObs[i]*sdObs[j];
      }
    MVNORM_t<Type> neg_log_density(covObs);
    
    vector<Type> obsvec(nA); 
    vector<Type> predvec(nA);
    data_indicator<vector<Type>,Type> keepvec(nA);
    
    for(int i=0; i<Y.dim[0]; ++i){
      for(int j=0; j<Y.dim[1]; ++j){
        pred(i,j)=z(i,j)+mu(j);
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

  if(mode==6){ // mode 1 and mode 4 mix
    DATA_MATRIX(Wr)
    DATA_MATRIX(Wc)
    DATA_MATRIX(Wd)
    DATA_ARRAY(Y)
    DATA_ARRAY_INDICATOR(keep,Y);
    DATA_SCALAR(trans); // 0 = log, >0 = x^(trans)
    matrix<Type> pred(Y.dim[0],Y.dim[1]);
    
    PARAMETER_VECTOR(logPhi)
    PARAMETER_VECTOR(mu)
    PARAMETER(logSdProc)   
    PARAMETER_VECTOR(logSdObs)   
    PARAMETER_MATRIX(z)
    PARAMETER(logitRhoObs)  
    vector<Type> phi=exp(logPhi);
    
    matrix<Type> I(Wr.rows(),Wr.cols());
    I.setIdentity();

    matrix<Type> Q=I-phi(0)*Wr-phi(1)*Wc-phi(2)*Wd;
    
    using namespace density;
    nll += SCALE(GMRF(asSparseMatrix(Q)),exp(logSdProc))(z.vec());

    Type rhoObs = invlogit(logitRhoObs);
    vector<Type> sdObs = exp(logSdObs);
    int nY = Y.dim[0];
    int nA = Y.dim[1];
    matrix<Type> covObs(nA,nA);
    for(int i=0;i<nA;i++)
      for(int j=0;j<nA;j++){
	covObs(i,j)=pow(rhoObs,Type(abs(i-j)))*sdObs[i]*sdObs[j];
      }
    MVNORM_t<Type> neg_log_density(covObs);

    matrix<Type> dW(nY,nA);
    dW.setZero();
        
    for(int i=0; i<nY; i++){
      for(int j=0; j<nA; j++){
	dW(i,j) = exp( z(i,j)+mu(j));
      }
    }
    
    for(int i=0; i<nY; i++) pred(i,0) = exp(mu(0)); 
    for(int j=1; j<nA; j++) pred(0,j) = pred(0,j-1) + exp(mu(j));
    
    for(int i=0; i<nY; i++) pred(i,0) = pred(i,0)*exp(z(i,0));
    for(int j=1; j<nA; j++) pred(0,j) = pred(0,j)*exp(z(0,j));
    
    for(int i=1; i<nY; i++){
      for(int j=1; j<nA; j++){ 
	pred(i,j) = pred(i-1,j-1) + dW(i,j);
      }
    }
    
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
      }
    }
    ADREPORT(pred);
  }

  if(mode==10){ // GMRF with neighbors in year, age, and cohort direction with box-cox  
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
    PARAMETER(logLamBC)
    vector<Type> phi=exp(logPhi);
    Type lamBC = exp(logLamBC);
    array<Type> Ytrans=Y*0;
    for(int i=0; i<Y.dim[0]; ++i){
      for(int j=0; j<Y.dim[1]; ++j){
        Ytrans(i,j)=(pow(Y(i,j),lamBC)-1)/lamBC;
	if(!isNA(Ytrans(i,j))){
          nll += -log(pow(Y(i,j),lamBC-Type(1)));
	}
      }
    }
    matrix<Type> I(Wr.rows(),Wr.cols());
    I.setIdentity();

    matrix<Type> Q=I-phi(0)*Wr-phi(1)*Wc-phi(2)*Wd;
    
    using namespace density;
    nll += SCALE(GMRF(asSparseMatrix(Q)),exp(logSdProc))(z.vec());

    for(int i=0; i<Ytrans.dim[0]; ++i){
      for(int j=0; j<Ytrans.dim[1]; ++j){
        pred(i,j)=z(i,j)+mu(j);
        if(!isNA(Ytrans(i,j))){
          nll += -dnorm(Ytrans(i,j),pred(i,j),exp(logSdObs(j)),true)*keep(i,j);
        }
      }
    }

    for(int i=0; i<Ytrans.dim[0]; ++i){
      for(int j=0; j<Ytrans.dim[1]; ++j){
        pred(i,j)=pow(lamBC*pred(i,j)+Type(1),Type(1)/lamBC);
      }
    }
    
    ADREPORT(pred);
  }

  if(mode==11){ // GMRF with neighbors in year, age, and cohort direction  
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

    matrix<Type> Q=I-phi(0)*Wc-phi(1)*Wd;
    
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

  if(mode==12){ // as mode 4, but no cohort RE, no obs. corr., and option to switch on/off exp on dW
    DATA_MATRIX(Wr)
    DATA_MATRIX(Wc)
    DATA_MATRIX(Wd)
    DATA_ARRAY(Y)  
    DATA_ARRAY_INDICATOR(keep,Y);
    DATA_INTEGER(expIt);
    int nY = Y.dim[0];
    int nA = Y.dim[1];
    
    matrix<Type> pred(Y.dim[0],Y.dim[1]);
    
    PARAMETER_VECTOR(logitRho)
    PARAMETER_VECTOR(mu)
    PARAMETER_VECTOR(logSdProc)
      
    PARAMETER_VECTOR(logSdObs)   
    PARAMETER_ARRAY(omega)
    vector<Type> rho=invlogit(logitRho);
    using namespace density;
    nll += SCALE(SEPARABLE(AR1(rho(1)),AR1(rho(0))),exp(logSdProc(0)))(omega);

    int od = nY + nA;
    matrix<Type> dW(od,nA);
    dW.setZero();
        
    for(int i=0; i<od; i++){
      for(int j=0; j<nA; j++){
	if(expIt==1){
	  dW(i,j) = exp( omega(i,j) + mu(j) );
	} else {
	  dW(i,j) = omega(i,j) + mu(j);
	}
      }
    }
    matrix<Type> predExt(od,nA);

    // initial conditions
    if(expIt==1){
      for(int i=0; i<od; i++) predExt(i,0) = exp(mu(0));
      for(int j=1; j<nA; j++) predExt(0,j) = predExt(0,j-1) + exp(mu(j));
      for(int i=0; i<od; i++) predExt(i,0) = predExt(i,0)*exp(omega(i,0));
      for(int j=1; j<nA; j++) predExt(0,j) = predExt(0,j)*exp(omega(0,j));
    } else {
      for(int i=0; i<od; i++) predExt(i,0) = mu(0);
      for(int j=1; j<nA; j++) predExt(0,j) = predExt(0,j-1) + mu(j);
      for(int i=0; i<od; i++) predExt(i,0) = predExt(i,0) + omega(i,0);
      for(int j=1; j<nA; j++) predExt(0,j) = predExt(0,j) + omega(0,j);
    }
    
    for(int i=1; i<od; i++){
      for(int j=1; j<nA; j++){ 
	predExt(i,j) = predExt(i-1,j-1) + dW(i,j);
      }
    }
    
    for(int i=0; i<Y.dim[0]; ++i){
      for(int j=0; j<Y.dim[1]; ++j){
	if(expIt){
	  pred(i,j) = log(predExt(i+nA,j)); 
	} else {
	  pred(i,j) = predExt(i+nA,j); 
	}
        if(!isNA(Y(i,j))){
          nll += -dnorm(Y(i,j),pred(i,j),exp(logSdObs(j)),true)*keep(i,j);
        }
      }
    }
    ADREPORT(pred);
  }

  if(mode==13){ // GMRF with neighbors in age, and cohort direction no obs error 
    DATA_MATRIX(Wr)
    DATA_MATRIX(Wc)
    DATA_MATRIX(Wd)
    DATA_ARRAY(Y)
    DATA_ARRAY_INDICATOR(keep,Y)
    matrix<Type> pred(Y.dim[0],Y.dim[1]);
    
    PARAMETER_VECTOR(logPhi)
    PARAMETER_VECTOR(mu)
    PARAMETER(logSdProc)   
    PARAMETER_VECTOR(missing)
    matrix<Type> z(Y.dim[0],Y.dim[1]); z.setZero();
    vector<Type> phi=exp(logPhi);
    matrix<Type> I(Wr.rows(),Wr.cols());
    I.setIdentity();
    matrix<Type> Q=I-phi(0)*Wc-phi(1)*Wd;

    int misCount=0;
    for(int i=0; i<Y.dim[0]; ++i){
      for(int j=0; j<Y.dim[1]; ++j){
        if(!isNA(Y(i,j))){
          z(i,j)=Y(i,j)-mu(j); 
        }else{
          z(i,j)=missing(misCount++)-mu(j); 
	}
      }
    }
    
    using namespace density;
    nll += SCALE(GMRF(asSparseMatrix(Q)),exp(logSdProc))(z.vec());

    for(int i=0; i<Y.dim[0]; ++i){
      for(int j=0; j<Y.dim[1]; ++j){
        pred(i,j)=z(i,j)+mu(j);
      }
    }
    ADREPORT(pred);
  }

  
  if(mode==14){ // As 12 but using N's
    DATA_MATRIX(Wr)
    DATA_MATRIX(Wc)
    DATA_MATRIX(Wd)
    DATA_ARRAY(Y)  
    DATA_ARRAY_INDICATOR(keep,Y);
    DATA_INTEGER(expIt);
    DATA_ARRAY(N);
    int nY = Y.dim[0];
    int nA = Y.dim[1];
    
    matrix<Type> pred(Y.dim[0],Y.dim[1]);
    
    PARAMETER_VECTOR(logitRho)
    PARAMETER_VECTOR(mu)
    PARAMETER_VECTOR(logSdProc)
      
    PARAMETER_VECTOR(logSdObs)   
    PARAMETER_ARRAY(omega)
    PARAMETER_VECTOR(lalpha);

    vector<Type> alpha=exp(lalpha);
    vector<Type> rho=invlogit(logitRho);
    using namespace density;
    nll += SCALE(SEPARABLE(AR1(rho(1)),AR1(rho(0))),exp(logSdProc(0)))(omega);

    int od = nY + nA;
    matrix<Type> dW(od,nA);
    dW.setZero();
        
    for(int i=0; i<od; i++){
      for(int j=0; j<nA; j++){
	if(expIt==1){
	  dW(i,j) = exp( omega(i,j) + mu(j) );
	} else {
	  dW(i,j) = omega(i,j) + mu(j);
	}
      }
    }
    matrix<Type> predExt(od,nA);

    // initial conditions
    if(expIt==1){
      for(int i=0; i<od; i++) predExt(i,0) = exp(mu(0));
      for(int j=1; j<nA; j++) predExt(0,j) = predExt(0,j-1) + exp(mu(j));
      for(int i=0; i<od; i++) predExt(i,0) = predExt(i,0)*exp(omega(i,0));
      for(int j=1; j<nA; j++) predExt(0,j) = predExt(0,j)*exp(omega(0,j));
    } else {
      for(int i=0; i<od; i++) predExt(i,0) = mu(0);
      for(int j=1; j<nA; j++) predExt(0,j) = predExt(0,j-1) + mu(j);
      for(int i=0; i<od; i++) predExt(i,0) = predExt(i,0) + omega(i,0);
      for(int j=1; j<nA; j++) predExt(0,j) = predExt(0,j) + omega(0,j);
    }
    
    for(int i=1; i<od; i++){
      for(int j=1; j<nA; j++){
	if(i<nA){
	  predExt(i,j) = predExt(i-1,j-1) + dW(i,j);
	} else {
	  predExt(i,j) = predExt(i-1,j-1) + dW(i,j)*exp(-alpha(j-1)*log(N(i-nA,j)+1));
	}
      }
    }
    
    for(int i=0; i<Y.dim[0]; ++i){
      for(int j=0; j<Y.dim[1]; ++j){
	if(expIt==1){
	  pred(i,j) = log(predExt(i+nA,j)); // + lalpha(j)*log(N(i,j)); 
	} else {
	  pred(i,j) = predExt(i+nA,j); // + lalpha(j)*log(N(i,j));  
	}
        if(!isNA(Y(i,j))){
          nll += -dnorm(Y(i,j),pred(i,j),exp(logSdObs(j)),true)*keep(i,j);
        }
      }
    }
    ADREPORT(pred);
  }
  
  if(mode==15){ // GMRF with neighbors in age, and cohort direction no obs error 
    DATA_MATRIX(Wr)
    DATA_MATRIX(Wc)
    DATA_MATRIX(Wd)
    DATA_ARRAY(Y)
    DATA_ARRAY_INDICATOR(keep,Y)
    matrix<Type> pred(Y.dim[0],Y.dim[1]);
    
    PARAMETER_VECTOR(logPhi)
    PARAMETER_VECTOR(mu)
    PARAMETER_VECTOR(logSdProc)   
    PARAMETER_VECTOR(missing)
    matrix<Type> z(Y.dim[0],Y.dim[1]); z.setZero();
    vector<Type> phi=exp(logPhi);
    matrix<Type> I(Wr.rows(),Wr.cols());
    I.setIdentity();
    matrix<Type> Q=I-phi(0)*Wc-phi(1)*Wd;

    int misCount=0;
    for(int i=0; i<Y.dim[0]; ++i){
      for(int j=0; j<Y.dim[1]; ++j){
        if(!isNA(Y(i,j))){
          z(i,j)=Y(i,j)-mu(j); 
        }else{
          z(i,j)=missing(misCount++)-mu(j); 
	}
      }
    }
    
    using namespace density;
    nll += VECSCALE(GMRF(asSparseMatrix(Q)),exp(logSdProc))(z.vec());

    for(int i=0; i<Y.dim[0]; ++i){
      for(int j=0; j<Y.dim[1]; ++j){
        pred(i,j)=z(i,j)+mu(j);
      }
    }
    ADREPORT(pred);
  }

  //if(mode==16){ // GMRF with neighbors in age, and cohort direction no obs error 
  //  DATA_MATRIX(Wr)
  //  DATA_MATRIX(Wc)
  //  DATA_MATRIX(Wd)
  //  DATA_ARRAY(Y)
  //  DATA_ARRAY(logN)
  //  DATA_ARRAY(Z)
  //  DATA_ARRAY_INDICATOR(keep,Y)
  //  matrix<Type> pred(Y.dim[0],Y.dim[1]);
  //  
  //  PARAMETER_VECTOR(logPhi)
  //  PARAMETER_VECTOR(mu)
  //  PARAMETER(logSdProc)
  //  PARAMETER_VECTOR(logSdLogN) 
  //  PARAMETER_VECTOR(missing)
  //  matrix<Type> z(Y.dim[0],Y.dim[1]); z.setZero();
  //  vector<Type> phi=exp(logPhi);
  //  matrix<Type> I(Wr.rows(),Wr.cols());
  //  I.setIdentity();
  //  matrix<Type> Q=I-phi(0)*Wc-phi(1)*Wd;
  //
  //  int nrow=logN.dim[0];
  //  int ncol=logN.dim[1];
  //  Type predN;
  //  for(int y=1; y<nrow; ++y){
  //    predN=logN(y-1,0);
  //    nll += -dnorm(logN(y,0),predN,exp(logSdLogN(0)),true)
  //    for(int a=1; a<ncol; ++a){
  //      predN=logN(y-1,a-1)-Z(y-1,a-1);
  //      if(a==(ncol-1)){
  //        predN=log(exp(predN)+exp(logN(y-1,a)-N(y-1,a)));
  //      }
  //      jnll += -dnorm(logN(y,a),predN,exp(logSdLogN(1)),true);
  //    }
  //  }
  //
  //
  //  
  //  int misCount=0;
  //  for(int i=0; i<Y.dim[0]; ++i){
  //    for(int j=0; j<Y.dim[1]; ++j){
  //      if(!isNA(Y(i,j))){
  //        z(i,j)=Y(i,j)-mu(j); 
  //      }else{
  //        z(i,j)=missing(misCount++)-mu(j); 
  //	}
  //	pred(i,j)=z(i,j)+mu(j);
  //    }
  //  }
  //  
  //  using namespace density;
  //  nll += SCALE(GMRF(asSparseMatrix(Q)),exp(logSdProc))(z.vec());
  //
  //  ADREPORT(pred);
  //}
  if(mode==17){  // VonB type model
    DATA_MATRIX(Wr);
    DATA_MATRIX(Wc);
    DATA_MATRIX(Wd);
    DATA_ARRAY(Y);
    DATA_VECTOR(SSB);
    DATA_ARRAY_INDICATOR(keep,Y);
    matrix<Type> pred(Y.dim[0],Y.dim[1]);
    
    PARAMETER(logWinf);
    PARAMETER(logk);
    PARAMETER(logitRho0);
    PARAMETER_VECTOR(logSdObs);
    PARAMETER(logalpha); // density dependence
    PARAMETER_VECTOR(logitRho1); //
    PARAMETER_VECTOR(logSdOmega);
    PARAMETER_VECTOR(omegaw); // Winf deviations
    PARAMETER_VECTOR(omegak); // k deviations

    // Priors
    nll -= dnorm(logSdOmega(0),Type(-3),Type(3),true);
    nll -= dnorm(logSdOmega(1),Type(-3),Type(3),true);
    
    Type rho0 = invlogit(logitRho0);
    Type rho1 = invlogit(logitRho1(0));
    Type rho2 = invlogit(logitRho1(1));
    
    Type alpha = exp(logalpha);
    
    
    using namespace density;
    nll += SCALE(AR1(rho1),exp(logSdOmega(0)))(omegaw);
    nll += SCALE(AR1(rho2),exp(logSdOmega(1)))(omegak);
    
    
    for(int i=0; i<Y.dim[0]; ++i){
      for(int j=0; j<Y.dim[1]; ++j){
	int cohort = i - j + Y.dim[1] - 1;
	Type Winfc  = exp(logWinf - alpha*SSB(i) + omegaw(cohort));
	Type kc = exp(logk + omegak(cohort));
        pred(i,j) = log(Winfc) + Type(3)*log(Type(1) - (Type(1) - rho0)*exp(-kc*Type(j)));  //log(Winf*pow(Type(1) - exp(-k*(j-t0)),3));
	  
        if(!isNA(Y(i,j))){
          nll += -dnorm(Y(i,j),pred(i,j),exp(logSdObs(j)),true)*keep(i,j);
        }
      }
    }
    ADREPORT(pred);
    
  }
  
  return(nll);
}
