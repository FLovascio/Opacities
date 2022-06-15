#pragma once

#include "roots.hxx"
#include <cmath>
#include <complex>
#include <functional>
#include <vector>

namespace conductivity {
template <class T> class mixedGrain {
public:
  std::vector<std::vector<std::complex<T>>> sigma_ij;
  std::vector<T> delta_i;
  std::vector<T> lambda_k;
  std::vector<std::complex<T>> sigma_eff_j;

  mixedGrain(std::vector<std::vector<std::complex<T>>> sigma_ij_input,
             std::vector<T> delta_i_input, std::vector<T> lambda_k_input) {
    sigma_ij = sigma_ij_input;
    delta_i = delta_i_input;
    lambda_k = lambda_k_input;
    sigma_eff_j = sigma_ij[0];
  }

  std::complex<T> BrugemannSum(int lambda_k_i, std::complex<T> sigma_eff) {
    std::complex<T> sum_value(0.0, 0.0);
    int n = delta_i.size();
    T n_minus_1 = n - 1.0;
    for (int i = 0; i < n; ++i) {
      sum_value += delta_i[i] * (sigma_ij[lambda_k_i][i] - sigma_eff) /
                   (sigma_ij[lambda_k_i][i] + (n_minus_1 * sigma_eff));
    }
    return sum_value;
  }

  std::complex<T> BrugemannSumDerivative(int lambda_k_i,
                                         std::complex<T> sigma_eff) {
    std::complex<T> sum_value(0.0, 0.0);
    int n = delta_i.size();
    T n_minus_1 = n - 1.0;
    for (int i = 0; i < n; ++i) {
      sum_value += delta_i[i] *
                   ((-sigma_ij[lambda_k_i][i] - (n_minus_1 * sigma_eff)) -
                    n_minus_1 * (sigma_ij[lambda_k_i][i] - sigma_eff)) /
                   ((sigma_ij[lambda_k_i][i] + (n_minus_1 * sigma_eff)) *
                    (sigma_ij[lambda_k_i][i] + (n_minus_1 * sigma_eff)));
    }
    return sum_value;
  }
};

template <class T> void solveSystem(mixedGrain<T> &grain) {
  std::complex<T> current_best_guess(1.0, 1.0);
  for (int i = 0; i < grain.lambda_k.size(); ++i) {
    grain.sigma_eff_j[i] = complexRootFind::solve(
        current_best_guess, grain.BrugemannSum, grain.BrugemannSumDerivative,
        grain.lambda_k[i], 1e-8);
    current_best_guess = grain.sigma_eff_j[i];
  }
}
}; // namespace conductivities
