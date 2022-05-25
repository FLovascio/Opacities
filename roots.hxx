#pragma once

#include <cmath>
#include <complex>
#include <vector>
#include <functional>
#include <thread>
#include <tuple>
#include <complex>
#include <iostream>

namespace complexRootFind{
    template <class T, class... argType>
    std::complex<T> newtonStep(std::complex<T> & x_0, std::function<std::complex<T>(std::tuple<argType...>,std::complex<T>)> & f,  std::function<std::complex<T>(std::tuple<argType...>,std::complex<T>)> & f_prime, std::tuple<argType...> & args){
        T x_1=x_0-f(args,x_0)/f_prime(args,x_0);
        return x_1;
    }
    template <class T, class... argType>
    std::complex<T> solve(std::complex<T> & x_0, std::function<std::complex<T>(std::tuple<argType...>,std::complex<T>)> & f,  std::function<std::complex<T>(std::tuple<argType...>,std::complex<T>)> & f_prime, std::tuple<argType...> & args, T & threshold){
        std::complex<T> x_n = x_0; 
        for(int i=0;i < 1000){
            x_n=newtonStep(x_n,f,f_prime,args);
            if(std::abs(f(args,x_n))<threshold){
                return x_n;
            }
        }
        std::cerr << "ROOTFINDER WARNING: The root finder did not converge in 1000 steps: returned x at last step" << std::endl;
    } 
}