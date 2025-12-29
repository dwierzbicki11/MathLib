#ifndef COMPLEX_H
#define COMPLEX_H
#include "Common.h"
DLLEXPORT void Complex_Add(double r1, double i1, double r2, double i2, double* rOut, double* iOut);
DLLEXPORT void Complex_Subtract(double r1, double i1, double r2, double i2, double* rOut, double* iOut);
DLLEXPORT void Complex_Multiply(double r1, double i1, double r2, double i2, double* rOut, double* iOut);
DLLEXPORT int Complex_Divide(double r1, double i1, double r2, double i2, double* rOut, double* iOut);
DLLEXPORT double Complex_Modulus(double r, double i);
DLLEXPORT double Complex_Argument(double r, double i);
#endif
