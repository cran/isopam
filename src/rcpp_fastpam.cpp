// Rcpp wrapper for FastPAM
// Minimal version for isopam

// [[Rcpp::plugins("cpp11")]]

#include <Rcpp.h>
using namespace Rcpp;

#include "pam.h"

//' FastPAM clustering
//' 
//' @description FastPAM: An improved version of PAM, that is usually O(k) times faster.
//' Based on Schubert and Rousseeuw (2019).
//' 
//' @param rdist The distance matrix (lower triangular matrix, column wise storage)
//' @param n The number of observations
//' @param k The number of clusters to produce.
//' @param maxiter The maximum number of iterations (default: 0)
//' @param initializer Initializer: either "BUILD" or "LAB" (default)
//' @param fasttol Tolerance for fast swapping behavior. Default: 1.0
//' @param seed Seed for random number generator. Default: 123456789
//' @return A list with cost, medoids, and assignment
//' @keywords internal
// [[Rcpp::export]]
List cpp_fastpam(NumericVector rdist, int n, int k, int maxiter = 0, 
                 std::string initializer = "LAB", double fasttol = 1.0, 
                 int seed = 123456789) {
    std::vector<double> dist = as<std::vector<double> >(rdist);
    
    RDistMatrix dm(n, dist);
    
    PAMInitializer* init;
    if (initializer == "BUILD") {
        init = new BUILD(&dm);
    } else {
        init = new LAB(&dm, seed);
    }
    
    FastPAM pam(n, &dm, init, k, maxiter, fasttol);
    
    double cost = pam.run();
    
    std::vector<int> medoids = pam.getMedoids();
    std::vector<int> results = pam.getResults();
    
    delete init;
    
    return List::create(
        Named("cost") = cost,
        Named("medoids") = medoids,
        Named("assignment") = results
    );
}
