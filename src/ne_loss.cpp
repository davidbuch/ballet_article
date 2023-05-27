#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;


// [[Rcpp::export]]
double ne_loss_cpp(IntegerVector ghat,
                   IntegerVector g,
                   IntegerVector qhat,
                   IntegerVector q,
                   IntegerVector colormatch_indices, // not including 0,0
                   double sep_loss = 1,
                   double join_loss = 1,
                   double miscolor_loss = 1/4) {
  int N = ghat.length();
  int N_matches = colormatch_indices.length();
  
  double loss = (N - 1)*miscolor_loss*Rcpp::sum(qhat != q);
  
  for(int i = 0; i < N_matches - 1; i++){
    int id1 = colormatch_indices(i) - 1;
    for(int j = i + 1; j < N_matches; j++){
      // adjust indices to C++ convention
      int id2 = colormatch_indices(j) - 1;
      loss += (((ghat[id1] != ghat[id2])*(g[id1] == g[id2])) ? sep_loss: 0.) +
         (((ghat[id1] == ghat[id2])*(g[id1] != g[id2])) ? join_loss: 0.);
    }
  }
  return loss;
}

