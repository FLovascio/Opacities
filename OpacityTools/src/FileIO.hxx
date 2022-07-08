#include <complex>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>

namespace Parser { // required for faster file io  with delimited files, will
                   // implement
template <class T, int realVal, int imagVal>
std::complex<T> parseLineToComplex(char *) {
  std::complex<T> value(0, 0);
  return value;
}
}; // namespace Parser

namespace delimitedFiles {
template <class T>
bool readDatToComplexVector(std::vector<std::complex<T>> &outputVector,
                            std::string filename) {
  std::ifstream input_file(filename);
  if (!input_file.is_open()) {
    std::cerr << "Could not open the file - '"
           << filename << "'\n";
    return false;
  }
  T numbers[3];
  int i =0;
  while (input_file >> numbers[0] >> numbers[1] >> numbers[2]) {
    outputVector[i]=std::complex<T>(numbers[1],numbers[2]);
    i++;
  }
  input_file.close();
  return true;
}
template <class T>
bool readDatToVector(std::vector<T> &outputVector, std::string filename) {
    std::ifstream input_file(filename);
  if (!input_file.is_open()) {
    std::cerr << "Could not open the file - '"
           << filename << "'\n";
    return false;
  }
  T numbers[3];
  int i =0;
  while (input_file >> numbers[0] >> numbers[1] >> numbers[2]) {
    outputVector[i]=numbers[0];
    i++;
  }
  input_file.close();
  return true;
}
template <class T>
bool writeComplexVectorToFile(std::vector<std::complex<T>> &inputVector,
                            std::string filename) {
  std::ofstream output_file(filename);
  if (!output_file.is_open()) {
    std::cerr << "Could not open the file - '"
           << filename << "'\n";
    return false;
  }
  int i = 0;
  while (output_file >> string(inputVector[i].real())+" "+string(inputVector[i].imag())+"\n") {
    i++;
  }
  output_file.close();
  return true;
}
template <class T>
bool writeVectorToFile(std::vector<T> &inputVector,
                            std::string filename) {
  std::ofstream output_file(filename);
  if (!output_file.is_open()) {
    std::cerr << "Could not open the file - '"
           << filename << "'\n";
    return false;
  }
  int i = 0;
  while (output_file >> string(inputVector[i])+"\n") {
    i++;
  }
  output_file.close();
  return true;
}
}; // namespace ReadFiles

namespace binaryFiles {
template <class T> bool writeComplexVectorToBinary() { return true; }
}; // namespace WriteFiles
