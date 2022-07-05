#include <complex>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace Parser {
template<class T, int realVal, int imagVal>
std::complex<T> parseLineToComplex(char*){
    std::complex<T> value(0,0);
    return value;
}
};

namespace ReadFiles {
template <class T, int realColumn, int imagColumn>
bool readDatToComplexVector(std::vector<std::complex<T>> &outputVector,
                            std::string filename) {
  return true;
}
template <class T, int Column>
bool readDatToVector(std::vector<T> &outputVector, std::string filename) {
  return true;
}
}; // namespace ReadFiles

namespace WriteFiles {
template <class T> bool writeComplexVectorToBinary() { return true; }
}; // namespace WriteFiles
