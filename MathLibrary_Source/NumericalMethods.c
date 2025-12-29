#include "NumericalMethods.h"
#include <math.h>
DLLEXPORT double Integrate_Trapezoidal(int n, double y[], double h) {if(n < 2) return 0; double sum = (y[0] + y[n-1]) / 2.0;for(int i=1; i < n-1; i++) sum += y[i]; return sum * h;}
DLLEXPORT int LinearRegression(int n, double x[], double y[], double* a, double* b, double* r_squared) {
    if (n < 2) return 0; double sX=0, sY=0, sXY=0, sX2=0;
    for(int i=0; i<n; i++) { sX+=x[i]; sY+=y[i]; sXY+=x[i]*y[i]; sX2+=x[i]*x[i]; }
    double den = (n*sX2 - sX*sX); if(den == 0) return 0;
    *a = (n*sXY - sX*sY) / den; *b = (sY - *a*sX) / n;
    double mY = sY/n, ssT=0, ssR=0; for(int i=0; i<n; i++) { ssT+=pow(y[i]-mY,2); ssR+=pow(y[i]-(*a*x[i]+*b),2); }
    *r_squared = (ssT==0)?1:(1 - (ssR/ssT)); return 1;
}
double D3(double A[]){return A[0]*A[4]*A[8]+A[1]*A[5]*A[6]+A[2]*A[3]*A[7]-A[2]*A[4]*A[6]-A[0]*A[5]*A[7]-A[1]*A[3]*A[8];}
DLLEXPORT int MatrixInverse3x3(double A[], double R[]) {
    double det=D3(A); if(fabs(det)<1e-9) return 0; double iD=1.0/det;
    R[0]=(A[4]*A[8]-A[5]*A[7])*iD; R[3]=-(A[3]*A[8]-A[5]*A[6])*iD; R[6]=(A[3]*A[7]-A[4]*A[6])*iD;
    R[1]=-(A[1]*A[8]-A[2]*A[7])*iD; R[4]=(A[0]*A[8]-A[2]*A[6])*iD; R[7]=-(A[0]*A[7]-A[1]*A[6])*iD;
    R[2]=(A[1]*A[5]-A[2]*A[4])*iD; R[5]=-(A[0]*A[5]-A[2]*A[3])*iD; R[8]=(A[0]*A[4]-A[1]*A[3])*iD; return 1;
}
