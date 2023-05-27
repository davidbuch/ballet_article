#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerMatrix encode_subpartition(IntegerVector x) {
  int N = x.size();
  IntegerMatrix A(N, N); // initialized to all 0s
  
  for(int i = 0; i < N; i++){
    if(x[i] != 0){ // leave row empty for inactive items
      for(int j = i; j < N; j++){
        A(i,j) = x[i] == x[j];
      }
    }
  }
  return A;
}

// [[Rcpp::export]]
IntegerVector decode_subpartition(IntegerMatrix A) {
  int N = A.nrow();
  // create a length N vector with default value -1
  // this will help us track which items have been 'visited'
  IntegerVector x(N, -1);
  
  // track the largest cluster index used so far
  int max_cluster = 1;
  for(int i = 0; i < N; i++){
    if(x[i] != -1){ 
      // the node has already been visited
      continue;
    }else if (A(i,i) == 0){
      // the node is 'inactive'
      x[i] = 0;
    }else{
      // the node is 'active'
      x[i] = max_cluster;
      for(int j = (i + 1); j < N; j++){
        // iterate over the other nodes
        // if j is connected to i, give it the same label
        if(A(i,j) == 1){
          x[j] = max_cluster;
        }
      }
      // update label to assign to the next 'active' item
      max_cluster++;
    }
  }
  return x;
}


// [[Rcpp::export]]
IntegerMatrix min_subpartition_pair(IntegerVector x, IntegerMatrix A) {
  int N = x.size();
  // Step 1: Determine which nodes are going to be active
  for(int i = 0; i < N; i++){
    A(i,i) = ((x[i] == 0) && (A(i,i) == 0))?0:1;
  }

  // Step 2: Graph min. Newly activated nodes should have all their
  // edges zeroes out by default since they won't be connected to old active
  // nodes
  for(int i = 0; i < N; i++){
    for(int j = i + 1; j < N; j++){
      // get the boolean encoding of x
      int x_code = x[i] == x[j];
      // compute the minimum between the x encoding and A
      A(i,j) = (x_code < A(i,j))?x_code:A(i,j);
    }
  }
  return A;
}

// [[Rcpp::export]]
IntegerMatrix max_subpartition_pair(IntegerVector x, IntegerMatrix A) {
  int N = x.size();
  
  // Step 1: Determine which nodes are going to be active
  for(int i = 0; i < N; i++){
    A(i,i) = ((x[i] == 0) || (A(i,i) == 0))?0:1;
  }
  
  // Step 2: Graph max
  for(int i = 0; i < N; i++){
    if(A(i,i) != 0){
      for(int j = i + 1; j < N; j++){
        A(i,j) = (x[i] == x[j])?1:A(i,j);
      }
    }
  }
  return A;
}

// [[Rcpp::export]]
IntegerMatrix floyd_warshall(IntegerMatrix adj, IntegerVector changed){
  int N = adj.nrow();
  
  
  // make the matrix symmetric (convenient for readibility)
  for(int i = 0; i < N - 1; i++){
    for(int j = i + 1; j < N; j++){
      adj(j,i) = adj(i,j);
    }
  }
  
  /* Add all vertices one by one to
   the set of intermediate vertices.
   ---> Before start of an iteration,
   we have shortest distances between all
   pairs of vertices such that the
   shortest distances consider only the
   vertices in set {0, 1, 2, .. k-1} as
   intermediate vertices.
   ----> After the end of an iteration,
   vertex no. k is added to the set of
   intermediate vertices and the set becomes {0, 1, 2, ..
   k} */
  for (int k = 0; k < N; k++) {
    Rcpp::checkUserInterrupt();
    if(changed[k] == 1){
      // Pick all vertices as source one by one
      for (int i = 0; i < N - 1; i++) {
        // Pick all vertices as destination for the
        // above picked source
        for (int j = i + 1; j < N; j++) {
          // If vertex k is on the shortest path from
          // i to j, then update the value of
          // dist[i][j]
          if ( adj(i,k) == 1 && adj(k,j) == 1 ){
            adj(i,j) = 1;
            adj(j,i) = 1;
          }
        }
      }
    }
  }
  return adj;
}

// [[Rcpp::export]]
IntegerVector edgewise_lowerbound(NumericMatrix psm, double alpha){
  int N = psm.nrow();
  IntegerMatrix A(N,N);
  
  for(int i = 0; i < N - 1; i++){
    A(i,i) = 1;
    for(int j = i + 1; j < N; j++){
      A(i,j) = (psm(i,j) < 1 - alpha)?0:1;
    }
  }
  
  return decode_subpartition(A);
}

// [[Rcpp::export]]
IntegerVector edgewise_upperbound(NumericMatrix psm, double alpha){
  int N = psm.nrow();
  IntegerVector changed(N,1);
  IntegerMatrix A(N,N);
  
  for(int i = 0; i < N - 1; i++){
    A(i,i) = 1;
    for(int j = i + 1; j < N; j++){
      A(i,j) = (psm(i,j) > alpha)?1:0;
    }
  }
  
  return decode_subpartition(floyd_warshall(A, changed));
}

// // [[Rcpp::export]]
// IntegerMatrix max_subpartition_pair(IntegerVector x, IntegerMatrix A) {
//   int N = x.size();
//   IntegerVector active(N, 0); // initialized to all 0s
//   IntegerVector changed(N, 0); // initialized to all 0s
// 
//   // Step 1: Determine which nodes are going to be active
//   int active_count = 0;
//   for(int i = 0; i < N; i++){
//     if((x[i] == 0) || (A(i,i) == 0)){
//       A(i,i) = 0;
//     }else{
//       A(i,i) = 1;
//       active[active_count] = 1;
//       active_count++;
//     }
//   }
// 
//   // Make sure all the inactive observations are separated from the actives
//   // for(int i = 0; i < N; i++){
//   //   if(A(i,i) == 0){
//   //     for(int j = i + 1; j < N; j++){
//   //       A(i,j) = (A(j,j) == 2)?1:0;
//   //     }
//   //   }
//   // }
// 
//   // Step 2: Among pairs for which both observations are active
//   // perform 'graph maximum'
//   for(int ii = 0; ii < active_count - 1; ii++){
//     for(int jj = ii + 1; jj < active_count; jj++){
//       int i = active[ii];
//       int j = active[jj];
// 
//       // get the boolean encoding of x
//       int x_code = x[i] == x[j];
//       // compute the minimum between the x encoding and A
//       if(x_code > A(i,j)){
//         A(i,j) = x_code;
//         if(A(i,i)){
//           changed[i] = 1;
//           changed[j] = 1;
//         }
//       }
//     }
//   }
//   // Transitive completion
//   A = floyd_warshall(A, changed);
//   return A;
// }

// [[Rcpp::export]]
IntegerVector min_subpartition(IntegerMatrix clustering_samps) {
  int S = clustering_samps.nrow();

  // keep A as a 'running minimum' across all samples
  IntegerVector x1 = clustering_samps.row(0);
  IntegerMatrix A = encode_subpartition(x1);
  for(int s = 1; s < S; s++){
    A = min_subpartition_pair(clustering_samps.row(s), A);
  }

  return decode_subpartition(A);
}

// [[Rcpp::export]]
IntegerVector max_subpartition(IntegerMatrix clustering_samps) {
  int S = clustering_samps.nrow();
  int N = clustering_samps.ncol();
  IntegerVector changed(N,1);
  
  // keep A as a 'running maximum' across all samples
  IntegerVector x1 = clustering_samps.row(0);
  IntegerMatrix A = encode_subpartition(x1);
  for(int s = 1; s < S; s++){
    A = max_subpartition_pair(clustering_samps.row(s), A);
  }

  return decode_subpartition(floyd_warshall(A, changed));
}


