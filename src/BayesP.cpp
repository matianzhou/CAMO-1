#include <RcppArmadillo.h>
#include <Rmath.h>
#include <cmath>
#include <math.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace std;
using namespace Rcpp;
using namespace arma;
//using namespace Rcpp::sugar;

//[[Rcpp::depends(RcppArmadillo)]]

const double mu0 = 0.0;
const double sigma0 = 10.0;
const double beta0 = 0.5;
const double alpha = 1.0;

int signC(int x) {
  if (x > 0) {
    return 1;
  } else if (x == 0) {
    return 0;
  } else {
    return -1;
  }
}

double integral_pos(double z,
                    double mu1, double sigma1){
  double factor1 = ( (R::pnorm(0,(z/(1*1) + mu1/(sigma1*sigma1))/
                      (1/(1*1) + 1/(sigma1*sigma1) + 1e-100),
                sqrt(1/(1/(1*1) + 1/(sigma1*sigma1))),0,0)) ) /
                     ((R::pnorm(0,mu1,sigma1,0,0)) + 1e-100);
  double factor2 = R::dnorm(z,  mu1, sqrt(1*1 + sigma1*sigma1),0);
  return(factor1*factor2);
}

double integral_neg(double z,
                    double mu1, double sigma1){
  double factor1 = ((R::pnorm(0,(z/(1*1) + mu1/(sigma1*sigma1))/
                                 (1/(1*1) + 1/(sigma1*sigma1)),
                        sqrt(1/(1/(1*1) + 1/(sigma1*sigma1))),1,0)))/
                     ((R::pnorm(0,mu1,sigma1,1,0)) + 1e-100 );
  double factor2 = R::dnorm(z,  mu1, sqrt(1*1 + sigma1*sigma1),0);
  return(factor1*factor2);
}

double h_zero(double z) {
  return(R::dnorm(z,0,1,0));
}

double h_old_pos(double z, double zc,
                 double zk, double nk, int G, double k) {

  if(zc==k) {
    nk = nk - 1;
    zk = zk - zc;
  }
  double factor1 = nk/(G-1 + alpha);
  double factor2 = integral_pos(z,
                             (zk/(1*1) + mu0/(sigma0*sigma0))/
                             (nk/(1*1) + 1/(sigma0*sigma0)) ,
                             1/(nk/(1*1) + 1/(sigma0*sigma0))+1  );
  return(factor1*factor2);
}

double h_new_pos(double z, int G){
  double factor1 = alpha/(G-1 + alpha);
  double factor2 = integral_pos(z,
                                mu0,
                                sigma0);
  return(factor1*factor2);
}

double h_old_neg(double z, double zc,
                 double zk, double nk, int G, int k) {

  if(zc==k) {
    nk = nk - 1;
    zk = zk - zc;
  }
  double factor1 = nk/(G-1 + alpha);
  double factor2 = integral_neg(z,
                                (zk/(1*1) + mu0/(sigma0*sigma0))/
                                  (nk/(1*1) + 1/(sigma0*sigma0)) ,
                                  1/(nk/(1*1) + 1/(sigma0*sigma0))+1  );
  return(factor1*factor2);
}

double h_new_neg(double z, int G){
  double factor1 = alpha/(G-1 + alpha);
  double factor2 = integral_neg(z,
                                mu0,
                                sigma0);
  return(factor1*factor2);
}


NumericVector updateC(NumericVector Z, NumericVector C,
                      const arma::vec pi, const arma::vec delta) {
  int G = Z.size();
  NumericVector Cunique = Rcpp::unique(C);
  double Kunique = Cunique.size();

  NumericVector nk(Kunique);
  NumericVector zk(Kunique);
  for(int g=0; g<G; g++){
    for(int k=0; k<Kunique; k++){
      if(C(g)==Cunique(k)){
        zk(k) = zk(k) + Z(g);
        nk(k)++;
      }
    }
  }

  for(int g=0; g<G; g++){

    double z=Z(g);
    double zc=C(g);
    double pipos=pi(g)*delta(g);
    double pineg=pi(g)*(1-delta(g));
    double pizero = 1-pi(g);

    NumericVector prob(Kunique);
    for(int k=0; k<Kunique; k++){
      if(Cunique(k)==0) {
        prob(k) = pizero*h_zero(z);
      }
      if(Cunique(k)>0) {
        prob(k) = pipos*h_old_pos(z, zc, zk(k), nk(k), G, Cunique(k));
      }
      if(Cunique(k)<0) {
        prob(k) = pineg*h_old_neg(z, zc, zk(k), nk(k), G, Cunique(k));
      }

     if(std::isnan(prob(k))) {
        //cout << prob <<endl;
        prob(k) = 0;
      }
    }

    double newk_pos = max(Cunique) + 1;
    double probnew_pos = pipos*h_new_pos(z, G);
    double newk_neg = min(Cunique) - 1;
    double probnew_neg = pineg*h_new_neg(z, G);

    NumericVector stdprob(Kunique+2);
    NumericVector ktotal(Kunique+2);
      double probtotal = sum(prob) +  probnew_pos + probnew_neg;
      for(int k=0; k< Kunique; k++){
        stdprob(k) = prob(k)/probtotal;
        ktotal(k) = Cunique(k);
      }

      stdprob(Kunique) = probnew_pos/probtotal;
      stdprob(Kunique+1) = probnew_neg/probtotal;
      ktotal(Kunique) = newk_pos;
      ktotal(Kunique+1) = newk_neg;

      C(g) = sample(ktotal,1,0,stdprob)[0];

       if(C(g) == newk_pos || C(g) == newk_neg) {
         NumericVector Cunique2(Kunique+1);
         NumericVector zk2(Kunique+1);
         NumericVector nk2(Kunique+1);
         for(int k=0; k<Kunique; k++){
           Cunique2(k) = Cunique(k);
             zk2(k) = zk(k);
             nk2(k) = nk(k);
           if(C(g)==k){
             zk(k) = zk(k) - Z(g);
           }
         }
         Cunique2(Kunique) = C(g);
         zk2(Kunique) = z;
         nk2(Kunique) = 1;
         Cunique = Cunique2;
         zk = zk2;
         nk = nk2;
         Kunique = Cunique.size();
      }  else {
            for(int k=0; k<Kunique; k++){
              if(Cunique(k)==zc){
                   nk(k) = nk(k) - 1;
                   zk(k) = zk(k) - z;
                  }
              if(Cunique(k)==C(g)){
                   nk(k) = nk(k) + 1;
                   zk(k) = zk(k) + z;
                  }
             }
        }

      if(min(nk)==0){
        int ind = which_min(nk);
        Cunique.erase(ind);
        zk.erase(ind);
        nk.erase(ind);
        Kunique = Cunique.size();
      }

  }

  return(C);
}

// [[Rcpp::export]]

Rcpp::List MCMC(NumericVector Z, int iteration, double const gamma) {

  int G = Z.size();
  arma::mat statePi(G, iteration);
  arma::mat stateDelta(G, iteration);
  arma::mat stateY(G, iteration);
  arma::mat stateC(G, iteration);

  arma::vec pi = rbeta(G, gamma/(G-gamma),1);
  statePi.col(0) = pi;
  arma::vec delta = rbeta(G, beta0, beta0);
  stateDelta.col(0) = delta;

  arma::vec Y(G);
  NumericVector C(G);

  IntegerVector startpool = seq(-1,1);
  //cout << startpool <<endl;
  for(int g=0; g<G; g++){
    double pi_pos = pi(g)*delta(g);
    double pi_neg = pi(g)*(1-delta(g));
    double pi_zero = 1 - pi_pos - pi_neg;
    if(pi_pos>pi_neg && pi_pos>pi_zero){
      Y(g) = 1;
    } else if(pi_neg>pi_pos && pi_neg>pi_zero){
      Y(g) = -1;
    } else {
      Y(g) = 0;
    }
    C(g) = sample(startpool,1)[0];
  }

  stateY.col(0) = Y;
  arma::vec C2  = C;
  stateC.col(0) = C2;

  for (int iter=1; iter<iteration; iter++){
    if(iter % 100==0) {
      cout << "Iteration: "<< iter+1 << endl;
    }
    //update pi
    for(int g=0; g<G; g++){
      if(Y(g) !=0) {
        pi(g) = rbeta(1,gamma/(G-gamma) + 1,  1)[0];
      } else if(Y(g)== 0) {
        pi(g) = rbeta(1,gamma/(G-gamma),  2)[0];
      }
    }
    statePi.col(iter) = pi;

    //update delta
    for(int g=0; g<G; g++){
      if(Y(g) ==1) {
        delta(g) = rbeta(1,beta0 + 1,  beta0 )[0];
      } else if(Y(g)== -1) {
        delta(g) = rbeta(1,beta0,  beta0 +1 )[0];
      } else {
        delta(g) = rbeta(1, beta0, beta0)[0];
      }
    }
    stateDelta.col(iter) = delta;

    //update C
    C = updateC(Z,C,pi,delta);
    arma::vec C2  = C;
    stateC.col(iter) = C2;

   // cout << Rcpp::unique(C) <<endl;

    //update Y
    for(int g=0; g<G; g++){
      Y(g) = signC(C(g));
    }
    stateY.col(iter) = Y;

  }

  Rcpp::List state;

  state["pi"] = statePi;
  state["delta"] = stateDelta;
  state["C"]=stateC;
  state["Y"]=stateY;

  cout << "MCMC sampling completed." <<endl;

  return(state);

}
