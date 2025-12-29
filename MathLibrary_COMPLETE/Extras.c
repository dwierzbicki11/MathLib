#include "Extras.h"
#include <math.h>
DLLEXPORT int SolveQuadraticInequality(double a, double b, double c, double* x1, double* x2){
    if(a==0)return-1; double d=b*b-4*a*c;
    if(d<0)return(a>0)?1:0; 
    *x1=(-b-sqrt(d))/(2*a); *x2=(-b+sqrt(d))/(2*a); if(*x1>*x2){double t=*x1;*x1=*x2;*x2=t;}
    return(d==0)?4:((a>0)?3:2);
}
DLLEXPORT double CalculateMean(int n, double d[]){double s=0;for(int i=0;i<n;i++)s+=d[i];return s/n;}
DLLEXPORT double CalculateStandardDeviation(int n, double d[]){double m=CalculateMean(n,d); double v=0; for(int i=0;i<n;i++)v+=pow(d[i]-m,2); return sqrt(v/n);}
DLLEXPORT int IsPrime(int n){if(n<2)return 0;for(int i=2;i<=sqrt(n);i++)if(n%i==0)return 0;return 1;}
