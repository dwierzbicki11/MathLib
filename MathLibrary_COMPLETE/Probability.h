#ifndef PROBABILITY_H
#define PROBABILITY_H
#include "Common.h"
DLLEXPORT unsigned long long Factorial(int n);
DLLEXPORT unsigned long long BinomialCoefficient(int n, int k);
DLLEXPORT double BernoulliTrial(int n, int k, double p);
DLLEXPORT unsigned long long VariationsNoRepetition(int n, int k);
DLLEXPORT double VariationsWithRepetition(int n, int k);
DLLEXPORT double ClassicalProbability(double f, double a);
#endif
