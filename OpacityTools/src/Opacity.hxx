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
  T epsilon;
  int nbin;
  std::vector<T> dustSizeBins;
  std::vector<T> dustSizeDensity;

  dustDistribution(const T &sMin, const T &sMax, const T &rhoGrain, const int &nBin, const T &Epsilon,
                   const std::function<T(T)> &dustSizeFunction)
      : smin(sMin), smax(sMax), nbin(nBin), rhograin(rhoGrain), epsilon(Epsilon) {
    dustSizeBins = makeSizeBins(smin, smax, nbin);
    dustSizeDensity = makeSizeDensity(dustSizeFunction);
  };
  std::vector<T> makeSizeBins(T &smin, T &smax, int &nbin) {
    std::vector<T> retVec(nbin, 0.0);
    for (int i = 0; i < nbin; ++i) {
      retVec[i] = smin + (smax - smin) * double(i) / double(nbin);
    }
    return retVec;
  }
  std::vector<T> makeSizeDensity(const std::function<T(T)> &dustSizeFunction) {
    std::vector<T> dustDensity = dustSizeBins;
    for (auto bin = std::begin(dustDensity); bin < std::end(dustDensity);
         ++bin) {
      *bin = epsilon*dustSizeFunction(*bin);
    }
    return dustDensity;
  }
};

template <class T>
std::function<T(T)> MRN_Pollack = [](T r) {
  return pow(r,-3)*(r >= 5e-4
             ? 0.0
             : (r >= 1e-4 ? (pow(1 / 0.005e-4, 2) * pow(0.005e-4 / r, 5.5))
                          : (r >= 0.005e-4 ? pow(0.005e-4 / r, 3.5) : 1)));
};

}; // namespace dust

namespace opacity {
template <class T> T H_j(T &xj, T &e1, T &e2) {
  return 1.0 + (xj * xj * ((e1 + 2.0) * (e1 + 2.0) + (e2 * e2)));
}

template <class T> T xj(T lambda, T sigma) {
  return 2.0 * M_PI / lambda *
         (0.3333333333333333 * sigma + 0.6666666666666666 * sigma);
}

template <class T> T sigma_jk(T &lambdak, T &e1, T &e2, const T &ljk) {
  return 2.0 * M_PI * e2 /
         ((lambdak * ljk * ljk) *
          ((e1 + 1.0 / ljk - 1.0) * (e1 + 1.0 / ljk - 1.0) + e2 * e2));
}

template <class T> T e1(std::complex<T> &n) {
  return (n.real() * n.real()) - (n.imag() * n.imag());
}

template <class T> T e2(std::complex<T> &n) { return 2.0 * n.real() * n.imag(); }

template <class T> T Kappa_j(int i, T &H, dust::dustDistribution<T> &dustDist) {
  T Kappa = 0.333333333 * dustDist.dustSizeDensity[i] * H *
            dustDist.dustSizeBins[i] * dustDist.dustSizeBins[i] *
            dustDist.dustSizeBins[i] / dustDist.rhograin;
  return Kappa;
}

template <class T>
void KappaDust_fast(std::vector<T> &output, conductivity::mixedGrain<T> &grain,
                    dust::dustDistribution<T> &dustDist) {
#ifdef WARN_FAST
  std::cerr
      << 'Warning: Fast methods do not boundscheck and may unexpectedly cause crashes and out of bound access\n'
      << 'Hint: output (passed as an arg to this function may have different length than lambda_k (also passed as an arg)\n';
#endif
  T fillValue = (T)0.0;
  std::fill(std::begin(output), std::end(output), fillValue);
  T lambda = 0.0;
  T e1Var = 0.0;
  T e2Var = 0.0;
  T sigma = 0.0;
  T xVar = 0.0;
  T HVar = 0.0;
  for (int k = 0; k < grain.lambda_k.size(); ++k) {
    lambda = grain.lambda_k[k];
    e1Var = e1(grain.sigma_eff_j[k]);
    e2Var = e2(grain.sigma_eff_j[k]);
    sigma = sigma_jk(lambda, e1Var, e2Var, 0.3333333333333333);
    xVar = xj(sigma, lambda);
    HVar = H_j(xVar, e1Var, e2Var);
    for (int idust = 0; idust < dustDist.nbin; ++idust) {
      output[k] += Kappa_j(idust, HVar, dustDist);
    }
  }
}

template <class T>
std::vector<T> KappaDust(conductivity::mixedGrain<T> &grain,
                    dust::dustDistribution<T> &dustDist) {
  T fillValue = (T)0.0;
  T lambda = 0.0;
  T e1Var = 0.0;
  T e2Var = 0.0;
  T sigma = 0.0;
  T xVar = 0.0;
  T HVar = 0.0;
  std::vector<T> output(grain.lambda_k.size(),fillValue);
  for (int k = 0; k < grain.lambda_k.size(); ++k) {
    lambda = grain.lambda_k[k];
    e1Var = e1<T>(grain.sigma_eff_j[k]);
    e2Var = e2<T>(grain.sigma_eff_j[k]);
    sigma = sigma_jk<T>(lambda, e1Var, e2Var, 0.3333333333333333);
    xVar = xj<T>(sigma, lambda);
    HVar = H_j<T>(xVar, e1Var, e2Var);
    for (int idust = 0; idust < dustDist.nbin; ++idust) {
      output[k] += Kappa_j(idust, HVar, dustDist);
    }
  }
  return output;
}
}; // namespace opacity
