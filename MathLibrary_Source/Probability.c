#include "Probability.h"
#include <math.h>
DLLEXPORT unsigned long long Factorial(int n){unsigned long long r=1;for(int i=2;i<=n;i++)r*=i;return r;}
DLLEXPORT unsigned long long BinomialCoefficient(int n, int k){if(k<0||k>n)return 0; if(k>n/2)k=n-k; unsigned long long r=1; for(int i=1;i<=k;i++){r=r*(n-i+1);r/=i;} return r;}
DLLEXPORT double BernoulliTrial(int n, int k, double p){return BinomialCoefficient(n,k)*pow(p,k)*pow(1-p,n-k);}
DLLEXPORT unsigned long long VariationsNoRepetition(int n, int k){unsigned long long r=1;for(int i=0;i<k;i++)r*=(n-i);return r;}
DLLEXPORT double VariationsWithRepetition(int n, int k){return pow((double)n,(double)k);}
DLLEXPORT double ClassicalProbability(double f, double a){return (a==0)?0:f/a;}
