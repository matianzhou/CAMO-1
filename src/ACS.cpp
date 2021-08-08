#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

/* ***** Function 1 ******** */
//' @export
// [[Rcpp::export]]
double calcSensC(NumericMatrix h, NumericMatrix m) {
  int iter = h.ncol();
  int G = h.nrow();
  double num = 0;
  double denom = 0;
  for(int i = 0; i < iter; ++i) {
      for(int g = 0; g < G; ++g)  {
        num += abs(h(g,i)) *(m(g,i) == h(g,i));
        denom += abs(h(g,i));
      }
  }
  return num / denom;
}

/* ***** Function 2 ******** */
//' @export
// [[Rcpp::export]]
double calcSpecC(NumericMatrix h, NumericMatrix m) {
  int iter = h.ncol();
  int G = h.nrow();
  double num = 0;
  double denom = 0;
  for(int i = 0; i < iter; ++i) {
      for(int g = 0; g < G; ++g)  {
        num += (1-abs(h(g,i))) *(m(g,i) == 0);
        denom += (1-abs(h(g,i)));
      }
  }
  return num / denom;
}

/* ***** Function 3 ******** */
//' @export
// [[Rcpp::export]]
double calcPrecC(NumericMatrix h, NumericMatrix m) {
  int iter = h.ncol();
  int G = h.nrow();
  double num = 0;
  double denom = 0;
  for(int i = 0; i < iter; ++i) {
      for(int g = 0; g < G; ++g)  {
        num += abs(h(g,i))*(m(g,i) == h(g,i));
        denom += abs(m(g,i));
      }
  }
  return num / denom;
}

/* ***** Function 4 ******** */
//' @export
// [[Rcpp::export]]
double calcESensC(NumericMatrix h, NumericMatrix m) {
  int iter = h.ncol();
  int G = h.nrow();
  double hmp = 0;
  double hmn = 0;
  double mmp = 0;
  double mmn = 0;

  for(int i = 0; i < iter; ++i) {
    for(int g = 0; g < G; ++g)  {
      if(h(g,i)==1) {
        hmp += 1 ;
      }
      if(h(g,i)== -1) {
        hmn += 1 ;
      }
      if(m(g,i)==1) {
        mmp += 1 ;
      }
      if(m(g,i)== -1) {
        mmn += 1 ;
      }
    }
  }

  return (mmp*hmp+mmn*hmn)/(G*iter*(hmp+hmn)) ;
}

/* ***** Function 5 ******** */
//' @export
// [[Rcpp::export]]
double calcESpecC(NumericMatrix h, NumericMatrix m) {
  int iter = h.ncol();
  int G = h.nrow();
  double mm0 = 0;
  for(int i = 0; i < iter; ++i) {
    for(int g = 0; g < G; ++g)  {
      if(m(g,i)==0){
        mm0 += 1;
      }
    }
  }
  return mm0/(G*iter);
}

/* ***** Function 6 ******** */
//' @export
// [[Rcpp::export]]
double calcEPrecC(NumericMatrix h, NumericMatrix m) {
  int iter = h.ncol();
  int G = h.nrow();
  double hmp = 0;
  double hmn = 0;
  double mmp = 0;
  double mmn = 0;

  for(int i = 0; i < iter; ++i) {
    for(int g = 0; g < G; ++g)  {
      if(h(g,i)==1) {
        hmp += 1 ;
      }
      if(h(g,i)== -1) {
        hmn += 1 ;
      }
      if(m(g,i)==1) {
        mmp += 1 ;
      }
      if(m(g,i)== -1) {
        mmn += 1 ;
      }
    }
  }

  return (mmp*hmp+mmn*hmn)/(G*iter*(mmp+mmn)) ;
}

