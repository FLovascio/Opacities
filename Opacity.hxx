#pragma once

#include "Conductivity.hxx"
#include "roots.hxx"
#include <cmath>
#include <complex>
#include <functional>
#include <vector>

namespace dust {
template <class T> class dustDistribution {
public:
  T smin;
  T smax;
  T rhograin;
  int nbin;
  std::vector<T> dustSizeBins;
  std::vector<T> dustSizeDensity;

  dustDistribution(T &smin, T &smax, T &rhograin, T &nbin)
};
} // namespace dust

namespace opacity {
template <class T> class opacity {

  dust::dustDistribution<T> dustDist;
  conductivity::mixedGrain<T> grainCond;
  std::vector<T> opacity;

  T Kappa_j(int i, T &H) {
    T Kappa = 0.333333333 * dustDist.dustSizeDensity[i] * H *
              dustDist.dustSizeBins * dustDist.dustSizeBins[i] *
              dustDist.dustSizeBins[i] / dustDist.rhograin;
    return Kappa;
  }
};

template <class T> T H_j(T &xj, T &e1, T &e2) {
  return 1.0 + (xj * xj * ((e1 + 2.0) * (e1 + 2.0) + (e2 * e2)));
}

template <class T> T sigma_jk(T & lambdak, T & e1, T & e2, T& ljk) {
  return 2.0*M_PI*e2/((lambdak*ljk*ljk)*((e1+1.0/ljk -1.0)*(e1+1.0/ljk -1.0)+e2*e2));
}

template <class T> T e1(std::complex<T> & n) {
    return n.real*n.real-n.imag*n.imag;
}

template <class T> T e2(std::complex<T> & n) {
    return 2.0*n.real*n.imag;
}
} // namespace opacity
