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

  dustDistribution(T &sMin, T &sMax, T &rhoGrain, T &nBin,
                   std::function<T(T)> &dustSizeFunction)
      : smin(sMin), smax(sMax), rhograin(rhoGrain) {
    dustSizeBins = makeSizeBins(smin, smax, nbin);
    dustSizeDensity = makeSizeDensity(dustSizeFunction);
  };
  std::vector<T> makeSizeBins(T &smin, T &smax, T &nbin) {
    std::vector<T> retVec(nbin, 0.0);
    for (int i = 0; i < nbin; ++i) {
      retVec[i] = smin + (smax - smin) * double(i) / double(nbin);
    }
    return retVec;
  }
  std::vector<T> makeSizeDensity(std::function<T(T)> &dustSizeFunction) {
    std::vector<T> dustDensity = dustSizeBins;
    for (auto bin = std::begin(dustDensity); bin < std::end(dustDensity);
         ++bin) {
      bin = dustSizeFunction(bin);
    }
    return dustDensity;
  }
};
}; // namespace dust

namespace opacity {
template <class T> T H_j(T &xj, T &e1, T &e2) {
  return 1.0 + (xj * xj * ((e1 + 2.0) * (e1 + 2.0) + (e2 * e2)));
}
template <class T> T  Kappa_j(int i, T & H,dust::dustDistribution<T> & dustDist) {
  T Kappa = 0.333333333 * dustDist.dustSizeDensity[i] * H *
            dustDist.dustSizeBins * dustDist.dustSizeBins[i] *
            dustDist.dustSizeBins[i] / dustDist.rhograin;
  return Kappa;
}
template <class T> void KappaDust_fast(std::vector<T> & output, conductivity::mixedGrain<T> & grain,dust::dustDistribution<T> & dustDist){
#ifdef WARN
  std::cerr << "Warning: Fast methods do not boundscheck and may unexpectedly cause crashes and out of bound access\n
  Hint: output (passed as an arg to this function may have different length than lambda_k (also passed as an arg)"
#endif 
  T fillValue =(T)0.0;
  std::fill(std::begin(output), std::end(output), fillValue);
  T lambda=0.0;
  T e1Var=0.0;
  T e2Var=0.0;
  T sigma=0.0;
  T xVar=0.0;
  T HVar=0.0;
  for(int k= 0; k<grain.lambda_k.size();++k){
    lambda=grain.lambda_k[k];
    e1Var=e1(grain.sigma_eff_j[k]);
    e2Var=e2(grain.sigma_eff_j[k]);
    sigma=sigma_jk(lambda, e1Var, e2Var, 0.3333333333333333);
    xVar=xj(sigma,lambda);
    HVar=H_j(xVar,e1Var,e2Var);
    for(int idust=0; idust<dustDist.nbin; ++idust){
      output[k]+=Kappa_j(idust,HVar,dustDist);
    }
  }

}
template <class T> T xj(T lambda, T sigma){
  return 2.0*M_PI/lambda*(0.3333333333333333*sigma+0.6666666666666666*sigma);
}

template <class T> T sigma_jk(T & lambdak, T & e1, T & e2, T & ljk) {
  return 2.0 * M_PI * e2 /
         ((lambdak * ljk * ljk) *
          ((e1 + 1.0 / ljk - 1.0) * (e1 + 1.0 / ljk - 1.0) + e2 * e2));
}

template <class T> T e1(std::complex<T> & n) {
  return n.real * n.real - n.imag * n.imag;
}

template <class T> T e2(std::complex<T> & n) { return 2.0 * n.real * n.imag; }
}; // namespace opacity
