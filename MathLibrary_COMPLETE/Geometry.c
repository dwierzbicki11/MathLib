#include "Geometry.h"
#include <math.h>
DLLEXPORT double GetTriangleAreaHeron(double a, double b, double c) { if(a+b<=c||a+c<=b||b+c<=a) return -1; double p=(a+b+c)/2.0; return sqrt(p*(p-a)*(p-b)*(p-c)); }
DLLEXPORT double GetCircleArea(double r){return (r<0)?-1:M_PI*r*r;}
DLLEXPORT double GetCircleCircumference(double r){return (r<0)?-1:2*M_PI*r;}
DLLEXPORT double GetDistance2D(double x1, double y1, double x2, double y2){return sqrt(pow(x2-x1,2)+pow(y2-y1,2));}
DLLEXPORT void GetMidPoint(double x1, double y1, double x2, double y2, double* mx, double* my){*mx=(x1+x2)/2.0;*my=(y1+y2)/2.0;}
DLLEXPORT int ArePointsCollinear(double x1, double y1, double x2, double y2, double x3, double y3){return (fabs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)) < 1e-9);}
DLLEXPORT double DistancePointToLine(double A, double B, double C, double x0, double y0){double d = sqrt(A*A+B*B); return (d==0)?-1:fabs(A*x0+B*y0+C)/d;}
DLLEXPORT int AreLinesParallel(double a1, double a2){return (fabs(a1-a2)<1e-9);}
DLLEXPORT int AreLinesPerpendicular(double a1, double a2){return (fabs(a1*a2+1)<1e-9);}
DLLEXPORT double GetLineSlope(double x1, double y1, double x2, double y2){return (x1==x2)?0:(y2-y1)/(x2-x1);}
