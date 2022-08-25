#pragma once
#include "Conductivity.hxx"
#include "complex.h"

template <typename T>
T * calculateConductivity(T** ReCond, T**ImCond, unsigned int nMaterial, unsigned int nLambda);
template <typename T>
T * calculateConductivityFromDir(char* dir);

