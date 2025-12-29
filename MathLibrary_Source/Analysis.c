#include "Analysis.h"
#include <math.h>
double PV(int d, double c[], double x){double r=c[0];for(int i=1;i<=d;i++)r=r*x+c[i];return r;}
double PD(int d, double c[], double x){double r=0;for(int i=0;i<d;i++)r+=c[i]*(d-i)*pow(x,d-i-1);return r;}
DLLEXPORT double PolynomialDerivativeAtPoint(int d, double c[], double x){return PD(d,c,x);}
DLLEXPORT int OptimizeQuadraticFunction(double a, double b, double c, double* x, double* y){if(a==0)return 0; *x=-b/(2*a); *y=a*(*x)*(*x)+b*(*x)+c; return (a>0)?1:2;}
DLLEXPORT void GetTangentLine(int d, double c[], double x0, double* a, double* b){*a=PD(d,c,x0); *b=PV(d,c,x0)-(*a*x0);}
double EvalFunc(int type, double x, double p[]) {switch(type) {case 0:return sin(x);case 1:return cos(x);case 2:return tan(x);case 3:return exp(x);case 4:if(p[2]*x+p[3]==0)return 0;return (p[0]*x+p[1])/(p[2]*x+p[3]);default:return 0;}}
DLLEXPORT double CalculateNumericalDerivative(int type, double x, double p[]) {double h=0.000001; return (EvalFunc(type,x+h,p)-EvalFunc(type,x-h,p))/(2.0*h);}
