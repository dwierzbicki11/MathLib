#include "Algebra.h"
#include <math.h>
#include <stdlib.h>
DLLEXPORT int SolveLinearEquation(double a, double b, double* x) { if (a == 0) return (b == 0) ? 2 : 0; *x = -b / a; return 1; }
DLLEXPORT int SolveLinearSystem2x2(double a1, double b1, double c1, double a2, double b2, double c2, double* x, double* y) {
    double W = a1*b2 - b1*a2, Wx = c1*b2 - b1*c2, Wy = a1*c2 - c1*a2;
    if (W == 0) return (Wx == 0 && Wy == 0) ? 2 : 0; *x = Wx/W; *y = Wy/W; return 1;
}
DLLEXPORT int SolveQuadraticEquation(double a, double b, double c, double* r1, double* r2) {
    if (a == 0) return -1; double d = b*b - 4*a*c;
    if (d < 0) return 0; if (d == 0) { *r1 = -b/(2*a); return 1; }
    *r1 = (-b - sqrt(d))/(2*a); *r2 = (-b + sqrt(d))/(2*a); return 2;
}
DLLEXPORT void GetQuadraticVertex(double a, double b, double c, double* p, double* q) { if(a!=0){*p = -b/(2*a); *q = -(b*b - 4*a*c)/(4*a);} }
DLLEXPORT double VieteSum(double a, double b) { return (a==0)?0:-b/a; }
DLLEXPORT double VieteProduct(double a, double c) { return (a==0)?0:c/a; }
DLLEXPORT double ArithmeticSequenceTerm(double a1, double r, int n) { return a1 + (n-1)*r; }
DLLEXPORT double ArithmeticSequenceSum(double a1, double an, int n) { return ((a1+an)/2.0)*n; }
DLLEXPORT double GeometricSequenceTerm(double a1, double q, int n) { return a1 * pow(q, n-1); }
DLLEXPORT double GeometricSequenceSum(double a1, double q, int n) { if(q==1) return a1*n; return a1*((1-pow(q,n))/(1-q)); }
DLLEXPORT double InfiniteGeometricSum(double a1, double q) { if(fabs(q)>=1) return 0; return a1/(1-q); }
DLLEXPORT double EvaluatePolynomialHorner(int d, double c[], double x) { double r = c[0]; for(int i=1; i<=d; i++) r = r*x + c[i]; return r; }
DLLEXPORT double LogarithmArbitraryBase(double b, double v) { if(b<=0||b==1||v<=0) return -1; return log(v)/log(b); }
DLLEXPORT double CompoundInterest(double K, double p, int n) { return K * pow(1.0 + p/100.0, n); }
DLLEXPORT void MatrixAdd3x3(double A[], double B[], double R[]){for(int i=0;i<9;i++)R[i]=A[i]+B[i];}
DLLEXPORT void MatrixMultiply3x3(double A[], double B[], double R[]){for(int r=0;r<3;r++) for(int c=0;c<3;c++){R[r*3+c] = A[r*3+0]*B[0*3+c] + A[r*3+1]*B[1*3+c] + A[r*3+2]*B[2*3+c];}}
DLLEXPORT double MatrixDeterminant3x3(double A[]){return A[0]*A[4]*A[8]+A[1]*A[5]*A[6]+A[2]*A[3]*A[7]-A[2]*A[4]*A[6]-A[0]*A[5]*A[7]-A[1]*A[3]*A[8];}
