#pragma once

#include <cmath>
#include <vector>
#include <complex>
#include <functional>
#include "roots.hxx"

namespace conductivities{
    template <class T>
    class mixedGrain{
        std::vector<std::vector<std::complex<T>>> sigma_ij;
        std::vector<T> delta_i;
        std::vector<T> lambda_k;
        std::vector<std::complex<T>> sigma_eff_j;

    public:
        mixedGrain(std::vector<std::vector<std::complex<T>>> sigma_ij_input, std::vector<T> delta_i_input, std::vector<T> lambda_k_input){
            sigma_ij = sigma_ij_input;
            delta_i = delta_i_input;
            lambda_k = lambda_k_input;
            sigma_eff_j = sigma_ij[0];
        }

        std::complex<T> BrugemannSum(int lambda_k_i, std::complex<T> sigma_eff){
            std::complex<T> sum_value(0.0, 0.0);
            int n = delta_ij[lambda_k_i].size();
            T n_minus_1 = n - 1.0;
            for (int i = 0; i < n; ++i){
                sum_value += delta_i[i] * (sigma_ij[lambda_k_i][i] - sigma_eff) / (sigma_ij[lambda_k_i][i] + (n_minus_1 * sigma_eff));
            }
            return sum_value;
        }

        std::complex<T> BrugemannSumDerivative(int lambda_k_i, std::complex<T> sigma_eff){
            std::complex<T> sum_value(0.0, 0.0);
            int n = delta_i.size();
            T n_minus_1 = n - 1.0;
            for (int i = 0; i < n; ++i){
                sum_value += delta_i[i] * ((-sigma_ij[lambda_k_i][i] - (n_minus_1 * sigma_eff)) - n_minus_1 * (sigma_ij[lambda_k_i][i] - sigma_eff)) / ((sigma_ij[lambda_k_i][i] + (n_minus_1 * sigma_eff)) * (sigma_ij[lambda_k_i][i] + (n_minus_1 * sigma_eff)));
            }
            return sum_value;
        }

        void solveSystem(){
            std::complex<T> current_best_guess(1.0, 1.0);
            for (int i = 0; i < lambda_k.size(); ++i){
                sigma_eff_j[i] = complexRootFind::solve(current_best_guess, BrugemannSum, BrugemannSumDerivative, lambda_k[i], 1e-8);
                current_best_guess = sigma_eff_j[i];
            }
        }
    };
};
