using System;
using System.Runtime.InteropServices;

namespace MathApp
{
    internal static class NativeMethods
    {
        private const string DllName = "MathLibrary.dll";
        private const CallingConvention Convention = CallingConvention.Cdecl;

        // --- COMPLEX NUMBERS (STUDIA) ---
        [DllImport(DllName, CallingConvention = Convention)] public static extern void Complex_Add(double r1, double i1, double r2, double i2, out double rOut, out double iOut);
        [DllImport(DllName, CallingConvention = Convention)] public static extern void Complex_Subtract(double r1, double i1, double r2, double i2, out double rOut, out double iOut);
        [DllImport(DllName, CallingConvention = Convention)] public static extern void Complex_Multiply(double r1, double i1, double r2, double i2, out double rOut, out double iOut);
        [DllImport(DllName, CallingConvention = Convention)] public static extern int Complex_Divide(double r1, double i1, double r2, double i2, out double rOut, out double iOut);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double Complex_Modulus(double r, double i);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double Complex_Argument(double r, double i);

        // --- NUMERICAL METHODS (STUDIA) ---
        [DllImport(DllName, CallingConvention = Convention)] public static extern double Integrate_Trapezoidal(int n, double[] y, double h);
        [DllImport(DllName, CallingConvention = Convention)] public static extern int LinearRegression(int n, double[] x, double[] y, out double a, out double b, out double r_sq);
        [DllImport(DllName, CallingConvention = Convention)] public static extern int MatrixInverse3x3(double[] A, double[] R);

        // --- VECTORS (LICEUM + STUDIA) ---
        [DllImport(DllName, CallingConvention = Convention)] public static extern void Vector2D_FromPoints(double xA, double yA, double xB, double yB, out double vx, out double vy);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double Vector2D_Length(double vx, double vy);
        [DllImport(DllName, CallingConvention = Convention)] public static extern void Vector2D_Add(double ux, double uy, double vx, double vy, out double rx, out double ry);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double Vector2D_DotProduct(double ux, double uy, double vx, double vy);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double Vector2D_Angle(double ux, double uy, double vx, double vy);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double Vector3D_Length(double vx, double vy, double vz);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double Vector3D_DotProduct(double ux, double uy, double uz, double vx, double vy, double vz);
        [DllImport(DllName, CallingConvention = Convention)] public static extern void Vector3D_CrossProduct(double ux, double uy, double uz, double vx, double vy, double vz, out double rx, out double ry, out double rz);

        // --- ALGEBRA (LICEUM + STUDIA) ---
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
        // Macierze
        [DllImport(DllName, CallingConvention = Convention)] public static extern void MatrixAdd3x3(double[] A, double[] B, double[] R);
        [DllImport(DllName, CallingConvention = Convention)] public static extern void MatrixMultiply3x3(double[] A, double[] B, double[] R);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double MatrixDeterminant3x3(double[] A);
        
        // --- GEOMETRY (LICEUM) ---
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetTriangleAreaHeron(double a, double b, double c);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetCircleArea(double r);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetDistance2D(double x1, double y1, double x2, double y2);
        [DllImport(DllName, CallingConvention = Convention)] public static extern void GetMidPoint(double x1, double y1, double x2, double y2, out double mx, out double my);
        [DllImport(DllName, CallingConvention = Convention)] public static extern int ArePointsCollinear(double x1, double y1, double x2, double y2, double x3, double y3);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double DistancePointToLine(double A, double B, double C, double x0, double y0);
        [DllImport(DllName, CallingConvention = Convention)] public static extern int AreLinesParallel(double a1, double a2);
        [DllImport(DllName, CallingConvention = Convention)] public static extern int AreLinesPerpendicular(double a1, double a2);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetLineSlope(double x1, double y1, double x2, double y2);

        // --- TRIGONOMETRY (LICEUM) ---
        [DllImport(DllName, CallingConvention = Convention)] public static extern double LawOfCosines_GetSide(double a, double b, double g);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double LawOfSines_GetSide(double s, double a1, double a2);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetCircumradius(double a, double al);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double Cotangent(double d);

        // --- ANALYSIS (LICEUM + STUDIA) ---
        [DllImport(DllName, CallingConvention = Convention)] public static extern double PolynomialDerivativeAtPoint(int d, double[] c, double x);
        [DllImport(DllName, CallingConvention = Convention)] public static extern int OptimizeQuadraticFunction(double a, double b, double c, out double x, out double y);
        [DllImport(DllName, CallingConvention = Convention)] public static extern void GetTangentLine(int d, double[] c, double x0, out double a, out double b);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double CalculateNumericalDerivative(int type, double x, double[] parameters);

        // --- PROBABILITY (LICEUM) ---
        [DllImport(DllName, CallingConvention = Convention)] public static extern ulong Factorial(int n);
        [DllImport(DllName, CallingConvention = Convention)] public static extern ulong BinomialCoefficient(int n, int k);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double BernoulliTrial(int n, int k, double p);
        [DllImport(DllName, CallingConvention = Convention)] public static extern ulong VariationsNoRepetition(int n, int k);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double VariationsWithRepetition(int n, int k);

        // --- SOLID GEOMETRY (LICEUM) ---
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetSphereVolume(double r);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetConeArea(double r, double h);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetCubeVolume(double a);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetCubeArea(double a);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetSphereArea(double r);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetCylinderVolume(double r, double h);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetCylinderArea(double r, double h);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double GetConeVolume(double r, double h);

        // --- EXTRAS ---
        [DllImport(DllName, CallingConvention = Convention)] public static extern int SolveQuadraticInequality(double a, double b, double c, out double x1, out double x2);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double CalculateMean(int n, double[] d);
        [DllImport(DllName, CallingConvention = Convention)] public static extern double CalculateStandardDeviation(int n, double[] d);
        [DllImport(DllName, CallingConvention = Convention)] public static extern int IsPrime(int n);
    }
}
