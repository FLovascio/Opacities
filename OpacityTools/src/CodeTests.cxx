////////////////////////////////////////////////
// This file implements tests for the c++    //
// implementation of the opacity calculator  //
///////////////////////////////////////////////
#include "Conductivity.hxx"
#include "FileIO.hxx"
#include "Opacity.hxx"
#include "roots.hxx"

#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <vector>

double quadraticRoot(double x, double Solution) { return x * x - Solution; }
double quadraticDerivative(double x, double Solution) { return 2.0 * x; }

std::complex<double> quadraticComplexRoot(std::complex<double> x,
                                          std::complex<double> offset) {
  return x * x + offset;
}
std::complex<double> quadraticComplexDerivative(std::complex<double> x,
                                                std::complex<double> offset) {
  return 2.0 * x;
}

int main() {
  // simple code integrity test//
  std::cout << "compiled ok\n";
  // root-finder test//
  // real roots//
  std::function<double(double, double)> f = quadraticRoot;
  std::function<double(double, double)> fprime = quadraticDerivative;
  double solutionVal = 1.0;
  double Thresh = 1e-4;
  auto solved =
      rootFind::solve<double, double>(0.1, Thresh, f, fprime, solutionVal);
  std::cout << "real root found: " << solved << "  expected: " << solutionVal
            << "\n";
  if (solved > solutionVal - Thresh && solved < solutionVal + Thresh) {
    std::cout << "real solver ok :) \n";
  } else {
    std::cout << "real solver bugged :( \n";
  }
  // complex roots
  std::function<std::complex<double>(std::complex<double>,
                                     std::complex<double>)>
      fC = quadraticComplexRoot;
  std::function<std::complex<double>(std::complex<double>,
                                     std::complex<double>)>
      fCprime = quadraticComplexDerivative;
  std::complex<double> offset(1.0, 0.0);
  auto solvedC = complexRootFind::solve<double, std::complex<double>>(
      std::complex<double>(0.1, 0.1), Thresh, fC, fCprime, solutionVal);
  std::cout << "complex root found: " << solvedC
            << "  expected: " << std::complex<double>(0.0, 1.0) << "\n";
  if (abs(fC(solvedC, offset)) < Thresh) {
    std::cout << "complex solver ok :) \n";
  } else {
    std::cout << "abs(f(x))=" << abs(fC(solvedC, offset))
              << ": complex solver bugged :( \n";
  }
  // conductivity solver tests//
  conductivity::mixedGrain<double> testGrain(conductivity::readGrain<double>(
      "/Users/fra/Code/Opacity/new_cons/Normal_silicates/"));
  std::cout << "read grain ok!\n";
  std::cout << "testGrain delta_i[0]=" << testGrain.delta_i[0] << "\n";
  std::cout << "lambda[0],sigma_ij[0][0]=" << testGrain.lambda_k[0] << ","
            << testGrain.sigma_ij[0][0]
            << "should be:  0.10000E+00,(0.14453E+01,0.89017E+00)\n";
  conductivity::solveSystem<double>(testGrain);
  std::cout << "lambda[0],sigma_eff_j[0]=" << testGrain.lambda_k[0] << ","
            << testGrain.sigma_eff_j[0] << "\n";
  return 0;
}
