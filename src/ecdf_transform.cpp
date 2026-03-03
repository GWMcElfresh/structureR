#include <Rcpp.h>
#include <algorithm>
#include <numeric>
using namespace Rcpp;

//' Fast Rank-Based ECDF Transformation to Pseudo-Observations
//'
//' Transforms a numeric vector to pseudo-observations in \eqn{(0,1)} using
//' a rank-based approach: \eqn{u_i = \mathrm{rank}(x_i) / (n+1)}.
//' Average ranks are used for ties. Implemented in C++ via \pkg{Rcpp} for
//' speed on large single-cell datasets.
//'
//' @param x A numeric vector of observations.
//' @return A numeric vector of pseudo-observations in \eqn{(0,1)}.
//' @export
// [[Rcpp::export]]
NumericVector rcpp_ecdf_transform(NumericVector x) {
  int n = x.size();
  if (n == 0) return NumericVector(0);

  // Create index array sorted by value in O(n log n)
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&x](int i, int j) {
    return x[i] < x[j];
  });

  // Assign average ranks to tied groups
  NumericVector ranks(n);
  int i = 0;
  while (i < n) {
    // Find the end of this tie group
    int j = i;
    while (j < n && x[idx[j]] == x[idx[i]]) ++j;
    // Average rank (1-indexed) for the group [i, j)
    double avg_rank = (i + 1.0 + j) / 2.0;
    for (int k = i; k < j; ++k) {
      ranks[idx[k]] = avg_rank;
    }
    i = j;
  }

  // Scale to (0, 1): rank / (n + 1)
  return ranks / (n + 1.0);
}
