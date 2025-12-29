<# :
@echo off
setlocal
cd /d %~dp0
powershell -NoProfile -ExecutionPolicy Bypass -Command "Invoke-Expression (Get-Content '%~f0' | Out-String)"
goto :eof
#>

# --- Poczatek Czesci PowerShell ---

$folderName = "MathLibrary_Source"
if (Test-Path $folderName) { Remove-Item $folderName -Recurse -Force }
New-Item -ItemType Directory -Force -Path $folderName | Out-Null
Set-Location $folderName

Write-Host "Generowanie BIBLIOTEKI TOTALNEJ (Wszystkie wzory z tablic CKE + Studia)..." -ForegroundColor Cyan

# --- 1. Common.h ---
$codeCommonH = @"
#ifndef COMMON_H
#define COMMON_H
#ifdef _WIN32
    #define DLLEXPORT __declspec(dllexport)
#else
    #define DLLEXPORT
#endif
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#endif
"@
Set-Content -Path "Common.h" -Value $codeCommonH

# --- 2. Solid Geometry (TOTALNIE ROZBUDOWANA) ---
$codeSolidH = @"
#ifndef SOLID_GEOMETRY_H
#define SOLID_GEOMETRY_H
#include "Common.h"

// --- GRANIASTOSLUPY I WALCE ---
DLLEXPORT double GetCubeVolume(double a);
DLLEXPORT double GetCubeArea(double a);
DLLEXPORT double GetCuboidVolume(double a, double b, double c); // Prostopadloscian
DLLEXPORT double GetCuboidArea(double a, double b, double c);
DLLEXPORT double GetCylinderVolume(double r, double h);
DLLEXPORT double GetCylinderArea(double r, double h);
DLLEXPORT double GetRegularPrismVolume(double baseArea, double h); // Graniastoslup ogolny

// --- OSTROSLUPY I STOZKI ---
DLLEXPORT double GetConeVolume(double r, double h);
DLLEXPORT double GetConeArea(double r, double h); // Pole calkowite (z pitagorasa)
DLLEXPORT double GetPyramidVolume(double baseArea, double h); // Ostroslup dowolny
DLLEXPORT double GetTetrahedronVolume(double a); // Czworoscian foremny
DLLEXPORT double GetTetrahedronArea(double a);
DLLEXPORT double GetOctahedronVolume(double a); // Osmioscian foremny

// --- BRYLY SCIETE (Tablice + Studia) ---
DLLEXPORT double GetFrustumConeVolume(double r1, double r2, double h); // Stozek sciety
DLLEXPORT double GetFrustumConeArea(double r1, double r2, double h);   // Pole boczne + podstawy
DLLEXPORT double GetFrustumPyramidVolume(double areaBase1, double areaBase2, double h); // Ostroslup sciety

// --- KULA I SFERA ---
DLLEXPORT double GetSphereVolume(double r);
DLLEXPORT double GetSphereArea(double r);
DLLEXPORT double GetSphericalCapVolume(double r, double h); // Czasza kuli (odcinek kuli)
DLLEXPORT double GetSphericalCapArea(double r, double h);   // Czasza kuli (pole sfery)

// --- INNE ---
DLLEXPORT double GetEllipsoidVolume(double a, double b, double c); // Studia
#endif
"@
Set-Content -Path "SolidGeometry.h" -Value $codeSolidH

$codeSolidC = @"
#include "SolidGeometry.h"
#include <math.h>

// Graniastoslupy
DLLEXPORT double GetCubeVolume(double a){return pow(a,3);}
DLLEXPORT double GetCubeArea(double a){return 6*pow(a,2);}
DLLEXPORT double GetCuboidVolume(double a, double b, double c){return a*b*c;}
DLLEXPORT double GetCuboidArea(double a, double b, double c){return 2*(a*b + a*c + b*c);}
DLLEXPORT double GetCylinderVolume(double r, double h){return M_PI*r*r*h;}
DLLEXPORT double GetCylinderArea(double r, double h){return 2*M_PI*r*(r+h);}
DLLEXPORT double GetRegularPrismVolume(double baseArea, double h){return baseArea*h;}

// Ostroslupy
DLLEXPORT double GetConeVolume(double r, double h){return (M_PI*r*r*h)/3.0;}
DLLEXPORT double GetConeArea(double r, double h){double l=sqrt(r*r+h*h); return M_PI*r*(r+l);}
DLLEXPORT double GetPyramidVolume(double baseArea, double h){return (baseArea*h)/3.0;}

// Wielościany Foremne (Wzory z tablic)
DLLEXPORT double GetTetrahedronVolume(double a){return (pow(a,3)*sqrt(2))/12.0;}
DLLEXPORT double GetTetrahedronArea(double a){return pow(a,2)*sqrt(3);}
DLLEXPORT double GetOctahedronVolume(double a){return (pow(a,3)*sqrt(2))/3.0;}

// Bryly Sciete
DLLEXPORT double GetFrustumConeVolume(double R, double r, double h){
    // V = (PI * h / 3) * (R^2 + R*r + r^2)
    return (M_PI * h / 3.0) * (R*R + R*r + r*r);
}
DLLEXPORT double GetFrustumConeArea(double R, double r, double h){
    double l = sqrt(pow(R-r, 2) + h*h); // tworzaca stozka scietego
    return M_PI * (R*R + r*r + (R+r)*l);
}
DLLEXPORT double GetFrustumPyramidVolume(double A1, double A2, double h){
    // V = (h/3) * (A1 + A2 + sqrt(A1*A2))
    return (h/3.0) * (A1 + A2 + sqrt(A1*A2));
}

// Kula i Elementy
DLLEXPORT double GetSphereVolume(double r){return (4.0/3.0)*M_PI*pow(r,3);}
DLLEXPORT double GetSphereArea(double r){return 4*M_PI*r*r;}
DLLEXPORT double GetSphericalCapVolume(double r, double h){
    // V = (PI * h^2 / 3) * (3r - h)
    return (M_PI * h * h / 3.0) * (3*r - h);
}
DLLEXPORT double GetSphericalCapArea(double r, double h){
    // Pole samej czaszy (bez podstawy): 2 * PI * r * h
    return 2 * M_PI * r * h;
}
DLLEXPORT double GetEllipsoidVolume(double a, double b, double c){
    return (4.0/3.0) * M_PI * a * b * c;
}
"@
Set-Content -Path "SolidGeometry.c" -Value $codeSolidC

# --- 3. Geometry 2D (UZUPELNIONA O TABLICE) ---
$codeGeometryH = @"
#ifndef GEOMETRY_H
#define GEOMETRY_H
#include "Common.h"
// Figury
DLLEXPORT double GetTriangleAreaHeron(double a, double b, double c);
DLLEXPORT double GetCircleArea(double r);
DLLEXPORT double GetCircleCircumference(double r);
DLLEXPORT double GetTrapezoidArea(double a, double b, double h); // NOWOSC
DLLEXPORT double GetRhombusArea(double d1, double d2); // NOWOSC
DLLEXPORT double GetRegularPolygonArea(int n, double side); // NOWOSC
DLLEXPORT double GetEllipseArea(double a, double b); // NOWOSC
DLLEXPORT double GetPolygonAreaFromPoints(int n, double x[], double y[]); // Wzor Sznurowadlowy

// Analityczna
DLLEXPORT double GetDistance2D(double x1, double y1, double x2, double y2);
DLLEXPORT void GetMidPoint(double x1, double y1, double x2, double y2, double* mx, double* my);
DLLEXPORT int ArePointsCollinear(double x1, double y1, double x2, double y2, double x3, double y3);
DLLEXPORT double DistancePointToLine(double A, double B, double C, double x0, double y0);
DLLEXPORT int AreLinesParallel(double a1, double a2);
DLLEXPORT int AreLinesPerpendicular(double a1, double a2);
DLLEXPORT double GetLineSlope(double x1, double y1, double x2, double y2);
#endif
"@
Set-Content -Path "Geometry.h" -Value $codeGeometryH

$codeGeometryC = @"
#include "Geometry.h"
#include <math.h>
#include <stdlib.h>

DLLEXPORT double GetTriangleAreaHeron(double a, double b, double c) { if(a+b<=c||a+c<=b||b+c<=a) return -1; double p=(a+b+c)/2.0; return sqrt(p*(p-a)*(p-b)*(p-c)); }
DLLEXPORT double GetCircleArea(double r){return (r<0)?-1:M_PI*r*r;}
DLLEXPORT double GetCircleCircumference(double r){return (r<0)?-1:2*M_PI*r;}

// Nowe Figury
DLLEXPORT double GetTrapezoidArea(double a, double b, double h){return ((a+b)*h)/2.0;}
DLLEXPORT double GetRhombusArea(double d1, double d2){return (d1*d2)/2.0;}
DLLEXPORT double GetRegularPolygonArea(int n, double side){
    if (n < 3) return 0;
    // Area = (n * s^2) / (4 * tan(PI/n))
    return (n * side * side) / (4.0 * tan(M_PI / n));
}
DLLEXPORT double GetEllipseArea(double a, double b){return M_PI * a * b;}
DLLEXPORT double GetPolygonAreaFromPoints(int n, double x[], double y[]) {
    // Shoelace Formula
    double area = 0.0;
    int j = n - 1; 
    for (int i = 0; i < n; i++) { 
        area += (x[j] + x[i]) * (y[j] - y[i]); 
        j = i; 
    } 
    return fabs(area / 2.0);
}

// Analityczna
DLLEXPORT double GetDistance2D(double x1, double y1, double x2, double y2){return sqrt(pow(x2-x1,2)+pow(y2-y1,2));}
DLLEXPORT void GetMidPoint(double x1, double y1, double x2, double y2, double* mx, double* my){*mx=(x1+x2)/2.0;*my=(y1+y2)/2.0;}
DLLEXPORT int ArePointsCollinear(double x1, double y1, double x2, double y2, double x3, double y3){return (fabs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)) < 1e-9);}
DLLEXPORT double DistancePointToLine(double A, double B, double C, double x0, double y0){double d = sqrt(A*A+B*B); return (d==0)?-1:fabs(A*x0+B*y0+C)/d;}
DLLEXPORT int AreLinesParallel(double a1, double a2){return (fabs(a1-a2)<1e-9);}
DLLEXPORT int AreLinesPerpendicular(double a1, double a2){return (fabs(a1*a2+1)<1e-9);}
DLLEXPORT double GetLineSlope(double x1, double y1, double x2, double y2){return (x1==x2)?0:(y2-y1)/(x2-x1);}
"@
Set-Content -Path "Geometry.c" -Value $codeGeometryC

# --- 4. Algebra (BEZ ZMIAN - JUZ DOBRA) ---
$codeAlgebraH = @"
#ifndef ALGEBRA_H
#define ALGEBRA_H
#include "Common.h"
DLLEXPORT int SolveLinearEquation(double a, double b, double* x);
DLLEXPORT int SolveLinearSystem2x2(double a1, double b1, double c1, double a2, double b2, double c2, double* x, double* y);
DLLEXPORT int SolveQuadraticEquation(double a, double b, double c, double* r1, double* r2);
DLLEXPORT void GetQuadraticVertex(double a, double b, double c, double* p, double* q);
DLLEXPORT double VieteSum(double a, double b);
DLLEXPORT double VieteProduct(double a, double c);
DLLEXPORT double ArithmeticSequenceTerm(double a1, double r, int n);
DLLEXPORT double ArithmeticSequenceSum(double a1, double an, int n);
DLLEXPORT double GeometricSequenceTerm(double a1, double q, int n);
DLLEXPORT double GeometricSequenceSum(double a1, double q, int n);
DLLEXPORT double InfiniteGeometricSum(double a1, double q);
DLLEXPORT double EvaluatePolynomialHorner(int degree, double coefficients[], double x);
DLLEXPORT double LogarithmArbitraryBase(double base, double value);
DLLEXPORT double CompoundInterest(double capital, double rate, int years);
DLLEXPORT void MatrixAdd3x3(double A[], double B[], double R[]);
DLLEXPORT void MatrixMultiply3x3(double A[], double B[], double R[]);
DLLEXPORT double MatrixDeterminant3x3(double A[]);
#endif
"@
Set-Content -Path "Algebra.h" -Value $codeAlgebraH

$codeAlgebraC = @"
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
"@
Set-Content -Path "Algebra.c" -Value $codeAlgebraC

# --- 5. Trigonometry ---
$codeTrigH = @"
#ifndef TRIGONOMETRY_H
#define TRIGONOMETRY_H
#include "Common.h"
DLLEXPORT double LawOfCosines_GetSide(double a, double b, double gammaDeg);
DLLEXPORT double LawOfSines_GetSide(double side, double ang1, double ang2);
DLLEXPORT double GetCircumradius(double a, double alphaDeg);
DLLEXPORT double Cotangent(double deg);
#endif
"@
Set-Content -Path "Trigonometry.h" -Value $codeTrigH

$codeTrigC = @"
#include "Trigonometry.h"
#include <math.h>
double ToRad(double d){return d*(M_PI/180.0);}
DLLEXPORT double LawOfCosines_GetSide(double a, double b, double g){double c2=a*a+b*b-2*a*b*cos(ToRad(g)); return (c2<0)?0:sqrt(c2);}
DLLEXPORT double LawOfSines_GetSide(double s, double a1, double a2){return (sin(ToRad(a1))==0)?0:(s*sin(ToRad(a2)))/sin(ToRad(a1));}
DLLEXPORT double GetCircumradius(double a, double al){return (sin(ToRad(al))==0)?0:a/(2*sin(ToRad(al)));}
DLLEXPORT double Cotangent(double d){double r=ToRad(d); return (sin(r)==0)?0:cos(r)/sin(r);}
"@
Set-Content -Path "Trigonometry.c" -Value $codeTrigC

# --- 6. Analysis ---
$codeAnalysisH = @"
#ifndef ANALYSIS_H
#define ANALYSIS_H
#include "Common.h"
DLLEXPORT double PolynomialDerivativeAtPoint(int degree, double coefficients[], double x);
DLLEXPORT int OptimizeQuadraticFunction(double a, double b, double c, double* xOpt, double* yOpt);
DLLEXPORT void GetTangentLine(int degree, double coefficients[], double x0, double* a, double* b);
DLLEXPORT double CalculateNumericalDerivative(int type, double x, double params[]);
#endif
"@
Set-Content -Path "Analysis.h" -Value $codeAnalysisH

$codeAnalysisC = @"
#include "Analysis.h"
#include <math.h>
double PV(int d, double c[], double x){double r=c[0];for(int i=1;i<=d;i++)r=r*x+c[i];return r;}
double PD(int d, double c[], double x){double r=0;for(int i=0;i<d;i++)r+=c[i]*(d-i)*pow(x,d-i-1);return r;}
DLLEXPORT double PolynomialDerivativeAtPoint(int d, double c[], double x){return PD(d,c,x);}
DLLEXPORT int OptimizeQuadraticFunction(double a, double b, double c, double* x, double* y){if(a==0)return 0; *x=-b/(2*a); *y=a*(*x)*(*x)+b*(*x)+c; return (a>0)?1:2;}
DLLEXPORT void GetTangentLine(int d, double c[], double x0, double* a, double* b){*a=PD(d,c,x0); *b=PV(d,c,x0)-(*a*x0);}
double EvalFunc(int type, double x, double p[]) {switch(type) {case 0:return sin(x);case 1:return cos(x);case 2:return tan(x);case 3:return exp(x);case 4:if(p[2]*x+p[3]==0)return 0;return (p[0]*x+p[1])/(p[2]*x+p[3]);default:return 0;}}
DLLEXPORT double CalculateNumericalDerivative(int type, double x, double p[]) {double h=0.000001; return (EvalFunc(type,x+h,p)-EvalFunc(type,x-h,p))/(2.0*h);}
"@
Set-Content -Path "Analysis.c" -Value $codeAnalysisC

# --- 7. Probability ---
$codeProbH = @"
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
"@
Set-Content -Path "Probability.h" -Value $codeProbH

$codeProbC = @"
#include "Probability.h"
#include <math.h>
DLLEXPORT unsigned long long Factorial(int n){unsigned long long r=1;for(int i=2;i<=n;i++)r*=i;return r;}
DLLEXPORT unsigned long long BinomialCoefficient(int n, int k){if(k<0||k>n)return 0; if(k>n/2)k=n-k; unsigned long long r=1; for(int i=1;i<=k;i++){r=r*(n-i+1);r/=i;} return r;}
DLLEXPORT double BernoulliTrial(int n, int k, double p){return BinomialCoefficient(n,k)*pow(p,k)*pow(1-p,n-k);}
DLLEXPORT unsigned long long VariationsNoRepetition(int n, int k){unsigned long long r=1;for(int i=0;i<k;i++)r*=(n-i);return r;}
DLLEXPORT double VariationsWithRepetition(int n, int k){return pow((double)n,(double)k);}
DLLEXPORT double ClassicalProbability(double f, double a){return (a==0)?0:f/a;}
"@
Set-Content -Path "Probability.c" -Value $codeProbC

# --- 8. ComplexMath ---
$codeComplexH = @"
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
"@
Set-Content -Path "ComplexMath.h" -Value $codeComplexH

$codeComplexC = @"
#include "ComplexMath.h"
#include <math.h>
DLLEXPORT void Complex_Add(double r1, double i1, double r2, double i2, double* rOut, double* iOut) { *rOut = r1 + r2; *iOut = i1 + i2; }
DLLEXPORT void Complex_Subtract(double r1, double i1, double r2, double i2, double* rOut, double* iOut) { *rOut = r1 - r2; *iOut = i1 - i2; }
DLLEXPORT void Complex_Multiply(double r1, double i1, double r2, double i2, double* rOut, double* iOut) {*rOut = r1*r2 - i1*i2; *iOut = r1*i2 + i1*r2;}
DLLEXPORT int Complex_Divide(double r1, double i1, double r2, double i2, double* rOut, double* iOut) {double den = r2*r2 + i2*i2; if(den == 0) return 0; *rOut = (r1*r2 + i1*i2) / den; *iOut = (i1*r2 - r1*i2) / den; return 1;}
DLLEXPORT double Complex_Modulus(double r, double i) { return sqrt(r*r + i*i); }
DLLEXPORT double Complex_Argument(double r, double i) { return atan2(i, r) * (180.0/M_PI); }
"@
Set-Content -Path "ComplexMath.c" -Value $codeComplexC

# --- 9. NumericalMethods ---
$codeNumH = @"
#ifndef NUMERICAL_H
#define NUMERICAL_H
#include "Common.h"
DLLEXPORT double Integrate_Trapezoidal(int n, double y[], double h);
DLLEXPORT int LinearRegression(int n, double x[], double y[], double* a, double* b, double* r_squared);
DLLEXPORT int MatrixInverse3x3(double A[], double Result[]);
#endif
"@
Set-Content -Path "NumericalMethods.h" -Value $codeNumH

$codeNumC = @"
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
"@
Set-Content -Path "NumericalMethods.c" -Value $codeNumC

# --- 10. Vectors ---
$codeVectorsH = @"
#ifndef VECTORS_H
#define VECTORS_H
#include "Common.h"
DLLEXPORT void Vector2D_FromPoints(double xA, double yA, double xB, double yB, double* vx, double* vy);
DLLEXPORT double Vector2D_Length(double vx, double vy);
DLLEXPORT void Vector2D_Add(double ux, double uy, double vx, double vy, double* resX, double* resY);
DLLEXPORT double Vector2D_DotProduct(double ux, double uy, double vx, double vy);
DLLEXPORT double Vector2D_Angle(double ux, double uy, double vx, double vy);
DLLEXPORT double Vector3D_Length(double vx, double vy, double vz);
DLLEXPORT double Vector3D_DotProduct(double ux, double uy, double uz, double vx, double vy, double vz);
DLLEXPORT void Vector3D_CrossProduct(double ux, double uy, double uz, double vx, double vy, double vz, double* rx, double* ry, double* rz);
#endif
"@
Set-Content -Path "Vectors.h" -Value $codeVectorsH

$codeVectorsC = @"
#include "Vectors.h"
#include <math.h>
DLLEXPORT void Vector2D_FromPoints(double xA, double yA, double xB, double yB, double* vx, double* vy) { *vx = xB - xA; *vy = yB - yA; }
DLLEXPORT double Vector2D_Length(double vx, double vy) { return sqrt(vx*vx + vy*vy); }
DLLEXPORT void Vector2D_Add(double ux, double uy, double vx, double vy, double* rx, double* ry) { *rx = ux + vx; *ry = uy + vy; }
DLLEXPORT double Vector2D_DotProduct(double ux, double uy, double vx, double vy) { return ux*vx + uy*vy; }
DLLEXPORT double Vector2D_Angle(double ux, double uy, double vx, double vy) {double L1=Vector2D_Length(ux,uy), L2=Vector2D_Length(vx,vy); if(L1==0||L2==0)return 0;double c=Vector2D_DotProduct(ux,uy,vx,vy)/(L1*L2); if(c>1)c=1; if(c<-1)c=-1; return acos(c)*(180.0/M_PI);}
DLLEXPORT double Vector3D_Length(double x, double y, double z){return sqrt(x*x+y*y+z*z);}
DLLEXPORT double Vector3D_DotProduct(double ux, double uy, double uz, double vx, double vy, double vz){return ux*vx+uy*vy+uz*vz;}
DLLEXPORT void Vector3D_CrossProduct(double ux, double uy, double uz, double vx, double vy, double vz, double* rx, double* ry, double* rz){*rx=uy*vz-uz*vy; *ry=uz*vx-ux*vz; *rz=ux*vy-uy*vx;}
"@
Set-Content -Path "Vectors.c" -Value $codeVectorsC

# --- 11. Extras ---
$codeExtraH = @"
#ifndef EXTRAS_H
#define EXTRAS_H
#include "Common.h"
DLLEXPORT int SolveQuadraticInequality(double a, double b, double c, double* x1, double* x2);
DLLEXPORT double CalculateMean(int n, double d[]);
DLLEXPORT double CalculateStandardDeviation(int n, double d[]);
DLLEXPORT int IsPrime(int n);
#endif
"@
Set-Content -Path "Extras.h" -Value $codeExtraH

$codeExtraC = @"
#include "Extras.h"
#include <math.h>
DLLEXPORT int SolveQuadraticInequality(double a, double b, double c, double* x1, double* x2){if(a==0)return-1; double d=b*b-4*a*c;if(d<0)return(a>0)?1:0; *x1=(-b-sqrt(d))/(2*a); *x2=(-b+sqrt(d))/(2*a); if(*x1>*x2){double t=*x1;*x1=*x2;*x2=t;}return(d==0)?4:((a>0)?3:2);}
DLLEXPORT double CalculateMean(int n, double d[]){double s=0;for(int i=0;i<n;i++)s+=d[i];return s/n;}
DLLEXPORT double CalculateStandardDeviation(int n, double d[]){double m=CalculateMean(n,d); double v=0; for(int i=0;i<n;i++)v+=pow(d[i]-m,2); return sqrt(v/n);}
DLLEXPORT int IsPrime(int n){if(n<2)return 0;for(int i=2;i<=sqrt(n);i++)if(n%i==0)return 0;return 1;}
"@
Set-Content -Path "Extras.c" -Value $codeExtraC

# --- 12. C# Wrapper (PELNY) ---
$codeCSharp = @"
using System;
using System.Runtime.InteropServices;

namespace MathApp
{
    internal static class NativeMethods
    {
        private const string DllName = "MathLibrary.dll";
        private const CallingConvention Convention = CallingConvention.Cdecl;

        // --- SOLID GEOMETRY (NEW EXTENDED) ---
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetCubeVolume(double a);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetCubeArea(double a);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetCuboidVolume(double a, double b, double c);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetCuboidArea(double a, double b, double c);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetCylinderVolume(double r, double h);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetCylinderArea(double r, double h);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetRegularPrismVolume(double baseArea, double h);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetConeVolume(double r, double h);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetConeArea(double r, double h);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetPyramidVolume(double baseArea, double h);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetTetrahedronVolume(double a);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetTetrahedronArea(double a);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetOctahedronVolume(double a);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetFrustumConeVolume(double r1, double r2, double h);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetFrustumConeArea(double r1, double r2, double h);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetFrustumPyramidVolume(double a1, double a2, double h);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetSphereVolume(double r);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetSphereArea(double r);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetSphericalCapVolume(double r, double h);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetSphericalCapArea(double r, double h);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetEllipsoidVolume(double a, double b, double c);

        // --- GEOMETRY 2D (NEW EXTENDED) ---
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetTriangleAreaHeron(double a, double b, double c);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetCircleArea(double r);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetCircleCircumference(double r);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetTrapezoidArea(double a, double b, double h);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetRhombusArea(double d1, double d2);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetRegularPolygonArea(int n, double side);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetEllipseArea(double a, double b);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetPolygonAreaFromPoints(int n, double[] x, double[] y);

        // --- ALGEBRA & ANALYSIS & REST ---
        [DllImport(DllName, CallingConvention = Convention)] public static extern int SolveLinearEquation(double a, double b, out double x);
        [DllImport(DllName, CallingConvention = Convention)] public static extern int SolveLinearSystem2x2(double a1, double b1, double c1, double a2, double b2, double c2, out double x, out double y);
        [DllImport(DllName, CallingConvention = Convention)] public static extern int SolveQuadraticEquation(double a, double b, double c, out double r1, out double r2);
        [DllImport(DllName, CallingConvention = Convention)] public static extern void GetQuadraticVertex(double a, double b, double c, out double p, out double q);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double VieteSum(double a, double b);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double VieteProduct(double a, double c);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double ArithmeticSequenceTerm(double a1, double r, int n);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double ArithmeticSequenceSum(double a1, double an, int n);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GeometricSequenceTerm(double a1, double q, int n);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GeometricSequenceSum(double a1, double q, int n);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double InfiniteGeometricSum(double a1, double q);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double CompoundInterest(double K, double p, int n);
        [DllImport(DllName, CallingConvention = Convention)] public static extern void MatrixAdd3x3(double[] A, double[] B, double[] R);
        [DllImport(DllName, CallingConvention = Convention)] public static extern void MatrixMultiply3x3(double[] A, double[] B, double[] R);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double MatrixDeterminant3x3(double[] A);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double LawOfCosines_GetSide(double a, double b, double g);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double LawOfSines_GetSide(double s, double a1, double a2);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetCircumradius(double a, double al);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double Cotangent(double d);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double PolynomialDerivativeAtPoint(int d, double[] c, double x);
        [DllImport(DllName, CallingConvention = Convention)] public static extern int OptimizeQuadraticFunction(double a, double b, double c, out double x, out double y);
        [DllImport(DllName, CallingConvention = Convention)] public static extern void GetTangentLine(int d, double[] c, double x0, out double a, out double b);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double CalculateNumericalDerivative(int type, double x, double[] parameters);
        [DllImport(DllName, CallingConvention = Convention)] public static extern ulong Factorial(int n);
        [DllImport(DllName, CallingConvention = Convention)] public static extern ulong BinomialCoefficient(int n, int k);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double BernoulliTrial(int n, int k, double p);
        [DllImport(DllName, CallingConvention = Convention)] public static extern ulong VariationsNoRepetition(int n, int k);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double VariationsWithRepetition(int n, int k);
        [DllImport(DllName, CallingConvention = Convention)] public static extern int SolveQuadraticInequality(double a, double b, double c, out double x1, out double x2);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double CalculateMean(int n, double[] d);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double CalculateStandardDeviation(int n, double[] d);
        [DllImport(DllName, CallingConvention = Convention)] public static extern int IsPrime(int n);
        [DllImport(DllName, CallingConvention = Convention)] public static extern void Complex_Add(double r1, double i1, double r2, double i2, out double rOut, out double iOut);
        [DllImport(DllName, CallingConvention = Convention)] public static extern void Complex_Multiply(double r1, double i1, double r2, double i2, out double rOut, out double iOut);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double Integrate_Trapezoidal(int n, double[] y, double h);
        [DllImport(DllName, CallingConvention = Convention)] public static extern int LinearRegression(int n, double[] x, double[] y, out double a, out double b, out double r_sq);
        [DllImport(DllName, CallingConvention = Convention)] public static extern int MatrixInverse3x3(double[] A, double[] R);
        [DllImport(DllName, CallingConvention = Convention)] public static extern void Vector2D_FromPoints(double xA, double yA, double xB, double yB, out double vx, out double vy);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double Vector2D_Length(double vx, double vy);
        [DllImport(DllName, CallingConvention = Convention)] public static extern void Vector3D_CrossProduct(double ux, double uy, double uz, double vx, double vy, double vz, out double rx, out double ry, out double rz);
        
        // ANALITYCZNA
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetDistance2D(double x1, double y1, double x2, double y2);
        [DllImport(DllName, CallingConvention = Convention)] public static extern void GetMidPoint(double x1, double y1, double x2, double y2, out double mx, out double my);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double DistancePointToLine(double A, double B, double C, double x0, double y0);
        [DllImport(DllName, CallingConvention = Convention)] public static extern int AreLinesParallel(double a1, double a2);
        [DllImport(DllName, CallingConvention = Convention)] public static extern int AreLinesPerpendicular(double a1, double a2);
    }
}
"@
Set-Content -Path "NativeMethods.cs" -Value $codeCSharp

Write-Host "ZAKONCZONO POMYSLNIE!" -ForegroundColor Green
Write-Host "Dodano brakujace bryly: Ostroslupy, Sciete, Czworosciany, Czasze..."
Write-Host "Uzupelniono geometrie 2D: Trapez, Romb, Wielokat Foremny, Elipsa..."
Start-Sleep -Seconds 3