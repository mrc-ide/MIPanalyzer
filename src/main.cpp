
#include "main.h"
#include "misc_v2.h"

using namespace std;

//------------------------------------------------
// Dummy function to test Rcpp working as expected
// [[Rcpp::export]]
Rcpp::List dummy1_cpp(Rcpp::List args) {
  
  // print message to console
  Rcpp::Rcout << "running C++ dummy1_cpp function\n";
  
  // get inputs from Rcpp format to base C++ format
  vector<double> x = Rcpp::as<vector<double>>(args("x"));
  
  // square values
  for (int i=0; i<int(x.size()); i++) {
    x[i] *= x[i];
  }
  
  // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("x_squared") = x);
  return ret;
}

//------------------------------------------------
// estimate f by maximum likelihood
// [[Rcpp::export]]
Rcpp::List inbreeding_mle(Rcpp::List args) {
  
  // extract inputs
  vector<vector<int>> x = rcpp_to_matrix_int(args["x"]);
  vector<double> f = rcpp_to_vector_double(args["f"]);
  vector<double> p = rcpp_to_vector_double(args["p"]);
  bool report_progress = rcpp_to_bool(args["report_progress"]);
  
  return Rcpp::List::create(Rcpp::Named("ret") = -9);
  
  /*
  // get dimensions
  int n = x.size();
  int L = x[0].size();
  
  // create q = 1-p
  vector<double> q(L);
  for (int j=0; j<L; ++j) {
    q[j] = 1.0 - p[j];
  }
  
  // create lookup tables
  vector<vector<double>> lookup_homo1(f_breaks, vector<double>(L));
  vector<vector<double>> lookup_homo2(f_breaks, vector<double>(L));
  vector<vector<double>> lookup_het(f_breaks, vector<double>(L));
  for (int k=0; k<f_breaks; ++k) {
    double f = k/double(f_breaks-1);
    for (int j=0; j<L; ++j) {
      lookup_homo1[k][j] = log((1-f)*p[j]*p[j] + f*p[j]);
      lookup_homo2[k][j] = log((1-f)*q[j]*q[j] + f*q[j]);
      lookup_het[k][j] = log((1-f)*2*p[j]*q[j]);
    }
  }
  
  // create objects for storing results
  vector<double> loglike_vec(f_breaks);
  vector<vector<double>> ret(n, vector<double>(n));
  
  // loop through all pairwise samples
  for (int i1=0; i1<(n-1); ++i1) {
    
    // report progress
    if (report_progress) {
      print("sample", i1, "of", n);
    }
    
    for (int i2=(i1+1); i2<n; ++i2) {
      
      // calculate loglike for every value of f
      for (int k=0; k<f_breaks; ++k) {
        double loglike = 0;
        for (int j=0; j<L; ++j) {
          if (x[i1][j] == -1 || x[i2][j] == -1) {
            continue;
          }
          if (x[i1][j] == 1 && x[i2][j] == 1) {
            loglike += lookup_homo1[k][j];
          } else if (x[i1][j] == 0 && x[i2][j] == 0) {
            loglike += lookup_homo2[k][j];
          } else {
            loglike += lookup_het[k][j];
          }
        }
        loglike_vec[k] = loglike;
      }
      
      // store maximum likelihood f
      double best_loglike = loglike_vec[0];
      for (int k=1; k<f_breaks; ++k) {
        if (loglike_vec[k] > best_loglike) {
          best_loglike = loglike_vec[k];
          ret[i1][i2] = k/double(f_breaks-1);
        }
      }
      
    }  // end i2 loop
  }  // end i1 loop
  
  // return as Rcpp object
  return Rcpp::List::create(Rcpp::Named("ret") = ret);
  */
}
