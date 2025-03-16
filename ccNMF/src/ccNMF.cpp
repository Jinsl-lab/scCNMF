#ifndef UTIL_H
#define UTIL_H
#include <Rcpp.h>
#include<RcppEigen.h>
#include<math.h>
#include<float.h>
#include<stdlib.h>
#include<vector>
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::plugins(openmp)]]
using namespace Rcpp;
#endif


//[[Rcpp::export]]
double eps(double a){
  int i=1;
  int b =(int)abs(a);
  double c=b;
  double epsilon=DBL_EPSILON;
  while((c/2)>1){
    i++;
    c/=2;
    epsilon*=2;
  }
  return epsilon;
}

Eigen::MatrixXd CppoperationMA_demo(Eigen::MatrixXd M, NumericVector A,int type){
  Eigen::MatrixXd ans(M.rows(),M.cols());
  switch (type)
  {
  case 0:
    for (int i = 0; i < M.rows(); i++){
      for (int j = 0; j < M.cols(); j++){
        ans(i,j)=A[i]+M(i,j);
      }
    }
    break;
  case 1:
    for (int i = 0; i < M.rows(); i++){
      for (int j = 0; j < M.cols(); j++){
        ans(i,j)=M(i,j)-A[i];
      }
    }
    break;
  case 2:
    for (int i = 0; i < M.rows(); i++){
      for (int j = 0; j < M.cols(); j++){
        ans(i,j)=A[i]*M(i,j);
      }
    }
    break;
  case 3:
    for (int i = 0; i < M.rows(); i++){
      for (int j = 0; j < M.cols(); j++){
        ans(i,j)=M(i,j)/A[i];
      }
    }
    break;
  }
  return ans;
}

Eigen::MatrixXd ifzeroMCPP(Eigen::MatrixXd M){
  Eigen::MatrixXd ans(M.rows(),M.cols());
  for (int i = 0; i < M.rows(); i++)
  {
    for (int j = 0; j < M.cols(); j++)
    {
      if(M(i,j)>0)    ans(i,j)=M(i,j);
      else    ans(i,j)=0;
    }
  }
  return ans;
}

Eigen::MatrixXd chooesVinMCPP(const Eigen::MatrixXd& M, const Rcpp::NumericVector& A, int type) {
  Eigen::VectorXd A_eigen = Rcpp::as<Eigen::VectorXd>(A);
  if (type == 0) {
    Eigen::MatrixXd ans(A_eigen.size(), M.cols());
    for (int i = 0; i < A_eigen.size(); i++) {
      for (int j = 0; j < ans.cols(); j++) {
        ans(i, j) = M(static_cast<int>(A_eigen(i)) - 1, j);
      }
    }
    return ans;
  }
  if (type == 1) {
    Eigen::MatrixXd ans(M.rows(), A_eigen.size());
    for (int i = 0; i < ans.rows(); i++) {
      for (int j = 0; j < A_eigen.size(); j++) {
        ans(i, j) = M(i, static_cast<int>(A_eigen(j)) - 1);
      }
    }
    return ans;
  }
  return Eigen::MatrixXd();
}

Eigen::MatrixXd epsMCpp(Eigen::MatrixXd M){
  Eigen::MatrixXd ans(M.rows(),M.cols());
  for (int i = 0; i < M.rows(); i++)
  {
    for (int j = 0; j < M.cols(); j++)
    {
      ans(i,j)=eps(M(i,j));
    }
  }
  return ans;
}


double Cvar(NumericVector AA, double Mean, int length){
  double Sum=sum(AA)/(length-1);
  Sum-=pow(Mean,2.0)*length/(length-1);
  return Sum;
}

double Cttest(NumericVector A, NumericVector B,NumericVector AA,NumericVector BB,int lengthA, int lengthB){
  double meanA=sum(A)/lengthA;
  double meanB=sum(B)/lengthB;
  double vartotal = Cvar(AA, meanA, lengthA) / lengthA + Cvar(BB, meanB, lengthB) / lengthB;
  return (meanA - meanB) / sqrt(vartotal);
}


// ccNMF_cpp_function
Rcpp::List ccNMF_cpp_function(Eigen::Map<Eigen::MatrixXd> X1, Eigen::Map<Eigen::MatrixXd> X2, Eigen::Map<Eigen::MatrixXd> Reg, int K, int maxIter,int stop_rule, double alpha, double beta, double gamma, Eigen::Map<Eigen::MatrixXd> W10, Eigen::Map<Eigen::MatrixXd> W20, Eigen::Map<Eigen::MatrixXd> W30, Eigen::Map<Eigen::MatrixXd> H0, Eigen::Map<Eigen::MatrixXd> θ0, Eigen::Map<Eigen::MatrixXd> O0, NumericVector c1, NumericVector c2, NumericVector Reg_w, int core);
RcppExport SEXP _RegNMF_ccNMF_cpp_function(SEXP X1SEXP, SEXP X2SEXP, SEXP RegSEXP, SEXP KSEXP, SEXP maxIterSEXP,SEXP stop_ruleSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP W10SEXP, SEXP W20SEXP, SEXP W30SEXP, SEXP H0SEXP, SEXP θ0SEXP, SEXP O0SEXP, SEXP c1SEXP, SEXP c2SEXP, SEXP Reg_wSEXP, SEXP coreSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Reg(RegSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< int >::type stop_rule(stop_ruleSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type W10(W10SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type W20(W20SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type W30(W30SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type H0(H0SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type θ0(θ0SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type O0(O0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type c1(c1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type c2(c2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Reg_w(Reg_wSEXP);
    Rcpp::traits::input_parameter< int >::type core(coreSEXP);
    rcpp_result_gen = Rcpp::wrap(ccNMF_cpp_function(X1, X2, Reg, K, maxIter,stop_rule, alpha, beta, gamma, W10, W20, W30, H0, θ0, O0, c1, c2, Reg_w, core));
    return rcpp_result_gen;
END_RCPP
}
// Fold_RE_TG_MultiAdjustCore
NumericMatrix Fold_RE_TG_MultiAdjustCore(NumericMatrix E2, NumericMatrix O2, NumericMatrix Symbol_location, NumericMatrix Peak_location);
RcppExport SEXP _RegNMF_Fold_RE_TG_MultiAdjustCore(SEXP E2SEXP, SEXP O2SEXP, SEXP Symbol_locationSEXP, SEXP Peak_locationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type E2(E2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type O2(O2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Symbol_location(Symbol_locationSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Peak_location(Peak_locationSEXP);
    rcpp_result_gen = Rcpp::wrap(Fold_RE_TG_MultiAdjustCore(E2, O2, Symbol_location, Peak_location));
    return rcpp_result_gen;
END_RCPP
}
// eps
double eps(double a);
RcppExport SEXP _RegNMF_eps(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(eps(a));
    return rcpp_result_gen;
END_RCPP
}


//[[Rcpp::export]]
double Cvar2(NumericVector AA, double Mean, int length){
  double Sum=sum(AA)/(length-1);
  Sum-=pow(Mean,2.0)*length/(length-1);
  return Sum;
}

// [[Rcpp::export]]
Rcpp::List ccNMF_cpp_function(Eigen::Map<Eigen::MatrixXd> X1,
                             Eigen::Map<Eigen::MatrixXd> X2,
                             Eigen::Map<Eigen::MatrixXd> Reg,
                             int K,
                             int maxIter,
                             int stop_rule,
                             double alpha,
                             double beta,
                             double gamma,
                             Eigen::Map<Eigen::MatrixXd> W10,
                             Eigen::Map<Eigen::MatrixXd> W20,
                             Eigen::Map<Eigen::MatrixXd> W30,
                             Eigen::Map<Eigen::MatrixXd> H0,
                             Eigen::Map<Eigen::MatrixXd> θ0,
                             Eigen::Map<Eigen::MatrixXd> O0,
                             NumericVector c1,
                             NumericVector c2,
                             NumericVector Reg_w,
                             int core){
  Eigen::setNbThreads(core);
  int n = Eigen::nbThreads();
  Rprintf("Core=%d\n",n);
  double tolx=1e-4;
  double tolfun=1e-6;
  double sqrteps = sqrt(DBL_EPSILON);
  double dn=0,dn_old=0;
  Eigen::MatrixXd numerR ,W1, W2, W3, H,θ, Tmp1,Tmp2;
  Eigen::MatrixXd Ht,θO,temp1,temp2,temp3,temp4,OX2tW2H,ORegtW3H,XtX2θOO,RegtRegθOO;
  Eigen::MatrixXd W1t,W2t,W3t,θ0t,X2t,Regt;
  dn_old=(X1-(W10*H0)).squaredNorm()+(alpha*(X2*(θ0*O0)-(W20*H0)).squaredNorm())+(beta*(Reg-(W30*H0)).squaredNorm())-(gamma*(θ0-(H0.transpose()*H0)).squaredNorm());
  for(int iter=0;iter<maxIter;iter++){
    Rcpp::checkUserInterrupt();
    Rprintf("iter=%d\n",iter);
     for (int i = 0; i < H0.rows(); i++) {
      double sumrow = 0;
       for(int j=0;j<H0.cols();j++){
        sumrow=sumrow+H0(i,j);
       }
       for(int k=0;k<H0.cols();k++){
        H0(i,k)=H0(i,k)/sumrow;
       }
    }

    /*update W1*/
    Ht=H0.transpose();
    W1=ifzeroMCPP((W10.array()*((X1*Ht).array()/(W10*H0*Ht+epsMCpp(X1*Ht)).array())).matrix());
    /*update W2*/
    θO = θ0;
    for (int i = 0; i < O0.rows(); i++) {
        for (int j = 0; j < O0.cols(); j++) {
            if (O0(i,j) == 0) {θO(i,j) = 0;}}
    }
    W2=ifzeroMCPP((W20.array()*(((X2*θO)*Ht).array()/((W20*H0*Ht)+epsMCpp((X2*θO)*Ht)).array())).matrix());
    /*update W3*/
    Tmp2=chooesVinMCPP(W2,c2,0);
    Tmp1=chooesVinMCPP(W1,c1,0);
    numerR=Tmp1+Tmp2;
    numerR=CppoperationMA_demo(numerR,Reg_w,2);
    numerR=Reg*Ht+numerR;
    W3=ifzeroMCPP((W30.array()*((numerR).array()/(W30*H0*Ht+W30+epsMCpp(numerR)).array())).matrix());
    /*update H*/
    W1t=W1.transpose();
    W2t=W2.transpose();
    W3t=W3.transpose();
    θ0t=θ0.transpose();
    temp1=(W1t*X1) + (alpha * (W2t*X2*θO)) + (beta * (W3t*Reg)) + (gamma * (H0*(θ0+θ0t)));
    temp2=((W1t*W1) + (alpha * (W2t*W2)) + (beta * (W3t*W3)) + (2 * gamma * (H0*Ht)))*H0;
    H=ifzeroMCPP((H0.array()*(temp1.array()/(temp2+epsMCpp(temp1)).array())).matrix());
    /*update θ*/
    X2t=X2.transpose();
    Regt=Reg.transpose();
    Ht=H.transpose();
    OX2tW2H = X2t*W2*H;
    ORegtW3H =Regt*W3*H;
    XtX2θOO = X2t*X2*θO;
    RegtRegθOO = Regt*Reg*θO;
    for (int i = 0; i < O0.rows(); i++) {
        for (int j = 0; j < O0.cols(); j++) {
            if (O0(i,j) == 0) {
                OX2tW2H(i,j) = 0;
                ORegtW3H(i,j) = 0;
                XtX2θOO(i,j) = 0;
                RegtRegθOO(i,j) = 0;}}
    }
    temp3=(alpha * OX2tW2H) + (gamma * Ht*H);
    temp4=(alpha * XtX2θOO) + (gamma * θ0); 
    θ=ifzeroMCPP((θ0.array()*(temp3.array()/(temp4+epsMCpp(temp3)).array())).matrix());
    /*Determine if the iteration is stopped*/
    if(stop_rule==2){
      dn=(X1-(W1*H)).squaredNorm()+(alpha*(X2*(θ*O0)-(W2*H)).squaredNorm())+(beta*(Reg-(W3*H)).squaredNorm())-(gamma*(θ-(H.transpose()*H)).squaredNorm());
     if(iter>500){
       if(iter>1 and (dn_old-dn)/dn_old<=tolfun){
         Rprintf("dnorm_old-dnorm %f is small", dn_old-dn);
         Rprintf("iter：", iter);
         break;
        }
       else if(iter==maxIter)
         break;
       dn_old=dn;
      }    
    } 
    W10 = W1;
    H0 = H;
    θ0 = θ;
    W20 = W2;
    W30 = W3;
}
return Rcpp::List::create(Named("W1") = wrap(W1),
                          Named("W2") = wrap(W2),
                          Named("W3") = wrap(W3),
                          Named("H") = wrap(H),
                          Named("θ") = wrap(θ));
}

// [[Rcpp::export]]
NumericMatrix Fold_RE_TG_MultiAdjustCore(NumericMatrix E2,
                                         NumericMatrix O2,
                                         NumericMatrix Symbol_location,
                                         NumericMatrix Peak_location){
  LogicalVector location1,location2,preid;
  NumericMatrix E2sqrt (E2.rows(),E2.cols());
  NumericMatrix P_1 (E2.rows(),O2.rows());
  IntegerVector id,rowid,colid;
  for(int i=0;i<(E2.length());i++){
    if(E2[i]!=0){E2sqrt[i]=pow(E2[i],2);}
    }
  Function w("which");
  for(int i=0;i<O2.rows();i++){
    Rcpp::checkUserInterrupt();
    location1=Symbol_location(_,0)==Peak_location(i,0);
    location2=abs(Symbol_location(_,1)-Peak_location(i,1))<1000000;
    preid=location1&location2;
    id=w(preid==TRUE);
    id=id-1;
    int tsum=0;
    for(int j=0;j<O2.cols();j++){
      if(O2(i,j)>0) tsum++;
    }
    for(int j=0;j<id.length();j++){
      int cnt0=0;
      int cnt1=0;
      for(int k=0;k<O2.cols();k++){
        if(E2(id(j),k)!=0){
          if(O2(i,k)>0){
            cnt1++;
          }
          else if (O2(i,k)==0) {
            cnt0++;
          }
        }
      }
      NumericVector set0 (cnt0);
      NumericVector set0sqrt (cnt0);
      NumericVector set1 (cnt1);
      NumericVector set1sqrt (cnt1);
      for(int k=0;k<O2.cols();k++){
      if(E2(id(j),k)!=0){
        if(O2(i,k)>0){
          cnt1--;
          set1(cnt1)=E2(id(j),k);
          set1sqrt(cnt1)=E2sqrt(id(j),k);
          }
        else if (O2(i,k)==0) {
          cnt0--;
          set0(cnt0)=E2(id(j),k);
          set0sqrt(cnt0)=E2sqrt(id(j),k);
          }
      }
      }
      P_1(id(j),i)=Cttest(set1,set0,set1sqrt,set0sqrt,tsum,O2.cols()-tsum);
    }
  }
  return wrap(P_1);
}

