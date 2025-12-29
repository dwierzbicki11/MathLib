#include "Trigonometry.h"
#include <math.h>
double ToRad(double d){return d*(M_PI/180.0);}
DLLEXPORT double LawOfCosines_GetSide(double a, double b, double g){double c2=a*a+b*b-2*a*b*cos(ToRad(g)); return (c2<0)?0:sqrt(c2);}
DLLEXPORT double LawOfSines_GetSide(double s, double a1, double a2){return (sin(ToRad(a1))==0)?0:(s*sin(ToRad(a2)))/sin(ToRad(a1));}
DLLEXPORT double GetCircumradius(double a, double al){return (sin(ToRad(al))==0)?0:a/(2*sin(ToRad(al)));}
DLLEXPORT double Cotangent(double d){double r=ToRad(d); return (sin(r)==0)?0:cos(r)/sin(r);}
