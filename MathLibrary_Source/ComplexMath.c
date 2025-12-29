#include "ComplexMath.h"
#include <math.h>
DLLEXPORT void Complex_Add(double r1, double i1, double r2, double i2, double* rOut, double* iOut) { *rOut = r1 + r2; *iOut = i1 + i2; }
DLLEXPORT void Complex_Subtract(double r1, double i1, double r2, double i2, double* rOut, double* iOut) { *rOut = r1 - r2; *iOut = i1 - i2; }
DLLEXPORT void Complex_Multiply(double r1, double i1, double r2, double i2, double* rOut, double* iOut) {*rOut = r1*r2 - i1*i2; *iOut = r1*i2 + i1*r2;}
DLLEXPORT int Complex_Divide(double r1, double i1, double r2, double i2, double* rOut, double* iOut) {double den = r2*r2 + i2*i2; if(den == 0) return 0; *rOut = (r1*r2 + i1*i2) / den; *iOut = (i1*r2 - r1*i2) / den; return 1;}
DLLEXPORT double Complex_Modulus(double r, double i) { return sqrt(r*r + i*i); }
DLLEXPORT double Complex_Argument(double r, double i) { return atan2(i, r) * (180.0/M_PI); }
