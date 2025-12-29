using System;

namespace MathApp
{
    class Program
    {
        // --- GŁÓWNA PĘTLA PROGRAMU ---
        static void Main(string[] args)
        {
            Console.Title = "MATH ENGINE ULTIMATE (C + C#)";

            while (true)
            {
                Console.Clear();
                DrawHeader("GŁÓWNE MENU");

                Console.ForegroundColor = ConsoleColor.White;
                Console.WriteLine(" 1. Algebra, Równania i Ciągi");
                Console.WriteLine(" 2. Geometria Płaska (2D) i Trygonometria");
                Console.WriteLine(" 3. Stereometria (Bryły 3D) [PEŁNA]");
                Console.WriteLine(" 4. Analiza Matematyczna (Pochodne)");
                Console.WriteLine(" 5. Wektory i Macierze");
                Console.WriteLine(" 6. Prawdopodobieństwo i Statystyka");
                Console.WriteLine(" 7. POZIOM STUDENCKI (Numeryczne, Zespolone)");
                Console.WriteLine(" 0. Wyjście");
                Console.ResetColor();
                Console.WriteLine("---------------------------------------------");
                Console.Write(" Wybierz opcję: ");

                var key = Console.ReadLine();

                switch (key)
                {
                    case "1": MenuAlgebra(); break;
                    case "2": MenuGeometry2D(); break;
                    case "3": MenuStereometry(); break;
                    case "4": MenuAnalysis(); break;
                    case "5": MenuVectors(); break;
                    case "6": MenuProbStat(); break;
                    case "7": MenuUniversity(); break;
                    case "0": return;
                    default: Console.WriteLine("Nieznana opcja!"); System.Threading.Thread.Sleep(500); break;
                }
            }
        }

        // =============================================================
        //                      DZIAŁ 1: ALGEBRA
        // =============================================================
        static void MenuAlgebra()
        {
            Console.Clear();
            DrawHeader("ALGEBRA I RÓWNANIA", ConsoleColor.Yellow);

            Console.WriteLine("1. Równanie Kwadratowe (Delta, Pierwiastki)");
            Console.WriteLine("2. Wzory Viete'a (Suma i Iloczyn)");
            Console.WriteLine("3. Układ Równań Liniowych 2x2");
            Console.WriteLine("4. Ciąg Arytmetyczny (Wyraz n i Suma)");
            Console.WriteLine("5. Ciąg Geometryczny (Wyraz n i Suma)");
            Console.WriteLine("6. Szereg Geometryczny Nieskończony");
            Console.WriteLine("7. Procent Składany (Lokaty)");
            Console.WriteLine("8. Nierówność Kwadratowa");

            Console.Write("\nWybór: ");
            switch (Console.ReadLine())
            {
                case "1":
                    double a = ReadDouble("a"), b = ReadDouble("b"), c = ReadDouble("c");
                    int count = NativeMethods.SolveQuadraticEquation(a, b, c, out double x1, out double x2);
                    if (count == 2) Console.WriteLine($"Dwa rozwiązania: x1={x1}, x2={x2}");
                    else if (count == 1) Console.WriteLine($"Jedno rozwiązanie: x0={x1}");
                    else Console.WriteLine("Brak rozwiązań rzeczywistych (Delta < 0)");

                    // Przy okazji wierzchołek
                    NativeMethods.GetQuadraticVertex(a, b, c, out double p, out double q);
                    Console.WriteLine($"Wierzchołek paraboli W(p,q): ({p}, {q})");
                    break;
                case "2":
                    double va = ReadDouble("a"), vb = ReadDouble("b"), vc = ReadDouble("c");
                    Console.WriteLine($"x1 + x2 = {NativeMethods.VieteSum(va, vb)}");
                    Console.WriteLine($"x1 * x2 = {NativeMethods.VieteProduct(va, vc)}");
                    break;
                case "3":
                    Console.WriteLine("Układ: { a1x + b1y = c1 }");
                    Console.WriteLine("       { a2x + b2y = c2 }");
                    double a1 = ReadDouble("a1"), b1 = ReadDouble("b1"), c1 = ReadDouble("c1");
                    double a2 = ReadDouble("a2"), b2 = ReadDouble("b2"), c2 = ReadDouble("c2");
                    int res = NativeMethods.SolveLinearSystem2x2(a1, b1, c1, a2, b2, c2, out double x, out double y);
                    if (res == 1) Console.WriteLine($"Wynik: x={x}, y={y}");
                    else Console.WriteLine("Układ sprzeczny lub nieoznaczony.");
                    break;
                case "4":
                    double a_ar = ReadDouble("a1"), r_ar = ReadDouble("r");
                    int n_ar = (int)ReadDouble("n");
                    Console.WriteLine($"Wyraz {n_ar}: {NativeMethods.ArithmeticSequenceTerm(a_ar, r_ar, n_ar)}");
                    Console.WriteLine($"Suma {n_ar}:  {NativeMethods.ArithmeticSequenceSum(a_ar, NativeMethods.ArithmeticSequenceTerm(a_ar, r_ar, n_ar), n_ar)}");
                    break;
                case "5":
                    double a_geo = ReadDouble("a1"), q_geo = ReadDouble("q");
                    int n_geo = (int)ReadDouble("n");
                    Console.WriteLine($"Wyraz {n_geo}: {NativeMethods.GeometricSequenceTerm(a_geo, q_geo, n_geo)}");
                    Console.WriteLine($"Suma {n_geo}:  {NativeMethods.GeometricSequenceSum(a_geo, q_geo, n_geo)}");
                    break;
                case "6":
                    double a_inf = ReadDouble("a1"), q_inf = ReadDouble("q (|q|<1)");
                    double s_inf = NativeMethods.InfiniteGeometricSum(a_inf, q_inf);
                    Console.WriteLine(s_inf == 0 && Math.Abs(q_inf) >= 1 ? "Szereg rozbieżny!" : $"Suma: {s_inf}");
                    break;
                case "7":
                    double K = ReadDouble("Kapitał"), proc = ReadDouble("Procent %"), lata = ReadDouble("Lata");
                    Console.WriteLine($"Kapitał końcowy: {NativeMethods.CompoundInterest(K, proc, (int)lata):F2}");
                    break;
                case "8":
                    double na = ReadDouble("a"), nb = ReadDouble("b"), nc = ReadDouble("c");
                    int ineqRes = NativeMethods.SolveQuadraticInequality(na, nb, nc, out double ix1, out double ix2);
                    Console.WriteLine($"Kod wyniku: {ineqRes} (opis w dokumentacji biblioteki)");
                    Console.WriteLine($"Miejsca zerowe do przedziału: {ix1}, {ix2}");
                    break;
            }
            Pause();
        }

        // =============================================================
        //                      DZIAŁ 2: GEOMETRIA 2D
        // =============================================================
        static void MenuGeometry2D()
        {
            Console.Clear();
            DrawHeader("GEOMETRIA PŁASKA I TRYGONOMETRIA", ConsoleColor.Green);

            Console.WriteLine("1. Pole Trójkąta (Heron)");
            Console.WriteLine("2. Pole Wielokąta Foremnego (np. 5-kąt)");
            Console.WriteLine("3. Pole Elipsy, Rombu lub Trapezu");
            Console.WriteLine("4. Pole dowolnego wielokąta (z wierzchołków)");
            Console.WriteLine("5. Odległość punktu od prostej");
            Console.WriteLine("6. Twierdzenie Cosinusów (Bok c)");
            Console.WriteLine("7. Promień Okręgu Opisanego (R)");

            Console.Write("\nWybór: ");
            switch (Console.ReadLine())
            {
                case "1":
                    Console.WriteLine($"Pole: {NativeMethods.GetTriangleAreaHeron(ReadDouble("a"), ReadDouble("b"), ReadDouble("c"))}");
                    break;
                case "2":
                    Console.WriteLine($"Pole: {NativeMethods.GetRegularPolygonArea((int)ReadDouble("Liczba boków n"), ReadDouble("Długość boku a"))}");
                    break;
                case "3":
                    Console.WriteLine("a) Romb  b) Trapez  c) Elipsa");
                    string sub = Console.ReadLine();
                    if (sub == "a") Console.WriteLine($"Pole: {NativeMethods.GetRhombusArea(ReadDouble("d1"), ReadDouble("d2"))}");
                    if (sub == "b") Console.WriteLine($"Pole: {NativeMethods.GetTrapezoidArea(ReadDouble("a"), ReadDouble("b"), ReadDouble("h"))}");
                    if (sub == "c") Console.WriteLine($"Pole: {NativeMethods.GetEllipseArea(ReadDouble("półoś a"), ReadDouble("półoś b"))}");
                    break;
                case "4":
                    Console.WriteLine("Ile wierzchołków?");
                    int n = (int)ReadDouble("n");
                    double[] px = new double[n]; double[] py = new double[n];
                    for (int i = 0; i < n; i++) { px[i] = ReadDouble($"x{i + 1}"); py[i] = ReadDouble($"y{i + 1}"); }
                    Console.WriteLine($"Pole wielokąta: {NativeMethods.GetPolygonAreaFromPoints(n, px, py)}");
                    break;
                case "5":
                    Console.WriteLine("Prosta Ax + By + C = 0");
                    Console.WriteLine($"Odległość: {NativeMethods.DistancePointToLine(ReadDouble("A"), ReadDouble("B"), ReadDouble("C"), ReadDouble("x0"), ReadDouble("y0"))}");
                    break;
                case "6":
                    Console.WriteLine($"Bok c: {NativeMethods.LawOfCosines_GetSide(ReadDouble("a"), ReadDouble("b"), ReadDouble("kąt gamma"))}");
                    break;
                case "7":
                    Console.WriteLine($"Promień R: {NativeMethods.GetCircumradius(ReadDouble("bok a"), ReadDouble("kąt naprzeciwko a"))}");
                    break;
            }
            Pause();
        }

        // =============================================================
        //                      DZIAŁ 3: STEREOMETRIA (Pełna)
        // =============================================================
        static void MenuStereometry()
        {
            Console.Clear();
            DrawHeader("STEREOMETRIA (BRYŁY 3D)", ConsoleColor.Cyan);

            Console.WriteLine("1. Sześcian i Prostopadłościan");
            Console.WriteLine("2. Kula i Czasza Kuli");
            Console.WriteLine("3. Walec i Stożek");
            Console.WriteLine("4. Stożek Ścięty (Frustum)");
            Console.WriteLine("5. Ostrosłup (Dowolny)");
            Console.WriteLine("6. Ostrosłup Ścięty");
            Console.WriteLine("7. Czworościan Foremny");
            Console.WriteLine("8. Ośmiościan Foremny");
            Console.WriteLine("9. Elipsoida");

            Console.Write("\nWybór: ");
            switch (Console.ReadLine())
            {
                case "1":
                    Console.WriteLine("a) Sześcian  b) Prostopadłościan");
                    if (Console.ReadLine() == "a")
                    {
                        double a = ReadDouble("a");
                        Console.WriteLine($"V={NativeMethods.GetCubeVolume(a)}, Pc={NativeMethods.GetCubeArea(a)}");
                    }
                    else
                    {
                        double a = ReadDouble("a"), b = ReadDouble("b"), c = ReadDouble("c");
                        Console.WriteLine($"V={NativeMethods.GetCuboidVolume(a, b, c)}, Pc={NativeMethods.GetCuboidArea(a, b, c)}");
                    }
                    break;
                case "2":
                    Console.WriteLine("a) Pełna Kula  b) Czasza Kuli");
                    double r = ReadDouble("Promień kuli R");
                    if (Console.ReadLine() == "a")
                    {
                        Console.WriteLine($"V={NativeMethods.GetSphereVolume(r)}, Pc={NativeMethods.GetSphereArea(r)}");
                    }
                    else
                    {
                        double h = ReadDouble("Wysokość czaszy h");
                        Console.WriteLine($"V={NativeMethods.GetSphericalCapVolume(r, h)}, Pole czaszy={NativeMethods.GetSphericalCapArea(r, h)}");
                    }
                    break;
                case "3":
                    Console.WriteLine("a) Walec  b) Stożek");
                    double r_cyl = ReadDouble("r"), h_cyl = ReadDouble("h");
                    if (Console.ReadLine() == "a")
                        Console.WriteLine($"V={NativeMethods.GetCylinderVolume(r_cyl, h_cyl)}, Pc={NativeMethods.GetCylinderArea(r_cyl, h_cyl)}");
                    else
                        Console.WriteLine($"V={NativeMethods.GetConeVolume(r_cyl, h_cyl)}, Pc={NativeMethods.GetConeArea(r_cyl, h_cyl)}");
                    break;
                case "4":
                    double R_fr = ReadDouble("R (dół)"), r_fr = ReadDouble("r (góra)"), h_fr = ReadDouble("h");
                    Console.WriteLine($"Objętość: {NativeMethods.GetFrustumConeVolume(R_fr, r_fr, h_fr):F2}");
                    Console.WriteLine($"Pole całk: {NativeMethods.GetFrustumConeArea(R_fr, r_fr, h_fr):F2}");
                    break;
                case "5":
                    Console.WriteLine($"V={NativeMethods.GetPyramidVolume(ReadDouble("Pole podstawy"), ReadDouble("Wysokość h"))}");
                    break;
                case "6":
                    Console.WriteLine($"V={NativeMethods.GetFrustumPyramidVolume(ReadDouble("Pole dolne"), ReadDouble("Pole górne"), ReadDouble("Wysokość"))}");
                    break;
                case "7":
                    double tet_a = ReadDouble("Krawędź a");
                    Console.WriteLine($"V={NativeMethods.GetTetrahedronVolume(tet_a)}, Pc={NativeMethods.GetTetrahedronArea(tet_a)}");
                    break;
                case "8":
                    Console.WriteLine($"V={NativeMethods.GetOctahedronVolume(ReadDouble("Krawędź a"))}");
                    break;
                case "9":
                    Console.WriteLine($"V={NativeMethods.GetEllipsoidVolume(ReadDouble("a"), ReadDouble("b"), ReadDouble("c"))}");
                    break;
            }
            Pause();
        }

        // =============================================================
        //                      DZIAŁ 4: ANALIZA
        // =============================================================
        static void MenuAnalysis()
        {
            Console.Clear();
            DrawHeader("ANALIZA MATEMATYCZNA", ConsoleColor.Magenta);

            Console.WriteLine("1. Styczna do wykresu wielomianu");
            Console.WriteLine("2. Optymalizacja (Max/Min funkcji kwadratowej)");
            Console.WriteLine("3. Pochodna NUMERYCZNA (sin, cos, exp...)");

            Console.Write("\nWybór: ");
            switch (Console.ReadLine())
            {
                case "1":
                    double[] c = ReadArray("Współczynniki wielomianu (np. 2 0 -1 dla 2x^2-1)");
                    double x0 = ReadDouble("Punkt x0");
                    NativeMethods.GetTangentLine(c.Length - 1, c, x0, out double ta, out double tb);
                    Console.WriteLine($"Styczna: y = {ta}x + {tb}");
                    break;
                case "2":
                    int type = NativeMethods.OptimizeQuadraticFunction(ReadDouble("a"), ReadDouble("b"), ReadDouble("c"), out double xOpt, out double yOpt);
                    Console.WriteLine($"{(type == 1 ? "MIN" : "MAX")} w x={xOpt}, y={yOpt}");
                    break;
                case "3":
                    Console.WriteLine("Typ: 0=sin, 1=cos, 2=tan, 3=exp");
                    int ftype = (int)ReadDouble("Typ");
                    double xv = ReadDouble("Punkt x");
                    Console.WriteLine($"Pochodna f'(x) ≈ {NativeMethods.CalculateNumericalDerivative(ftype, xv, new double[4])}");
                    break;
            }
            Pause();
        }

        // =============================================================
        //                      DZIAŁ 5: WEKTORY
        // =============================================================
        static void MenuVectors()
        {
            Console.Clear();
            DrawHeader("WEKTORY I MACIERZE", ConsoleColor.Blue);

            Console.WriteLine("1. Wektory 2D (Długość, Kąt, Iloczyn skalarny)");
            Console.WriteLine("2. Wektory 3D (Iloczyn Wektorowy / Cross Product)");
            Console.WriteLine("3. Operacje na Macierzach 3x3 (Mnożenie, Wyznacznik, Odwracanie)");

            Console.Write("\nWybór: ");
            switch (Console.ReadLine())
            {
                case "1":
                    double vx = ReadDouble("vx"), vy = ReadDouble("vy");
                    double ux = ReadDouble("ux"), uy = ReadDouble("uy");
                    Console.WriteLine($"Długość V: {NativeMethods.Vector2D_Length(vx, vy)}");
                    Console.WriteLine($"Iloczyn skalarny: {NativeMethods.Vector2D_DotProduct(vx, vy, ux, uy)}");
                    Console.WriteLine($"Kąt między wektorami: {NativeMethods.Vector2D_Angle(vx, vy, ux, uy)} stopni");
                    break;
                case "2":
                    Console.WriteLine("Wektor U [x y z]:");
                    double u3x = ReadDouble("x"), u3y = ReadDouble("y"), u3z = ReadDouble("z");
                    Console.WriteLine("Wektor V [x y z]:");
                    double v3x = ReadDouble("x"), v3y = ReadDouble("y"), v3z = ReadDouble("z");
                    NativeMethods.Vector3D_CrossProduct(u3x, u3y, u3z, v3x, v3y, v3z, out double rx, out double ry, out double rz);
                    Console.WriteLine($"Wynik U x V = [{rx}, {ry}, {rz}]");
                    break;
                case "3":
                    Console.WriteLine("a) Mnożenie A*B  b) Wyznacznik A  c) Odwracanie A^-1");
                    string sub = Console.ReadLine();
                    double[] mA = ReadArray("Podaj macierz A (9 liczb)");
                    if (sub == "a")
                    {
                        double[] mB = ReadArray("Podaj macierz B (9 liczb)");
                        double[] res = new double[9];
                        NativeMethods.MatrixMultiply3x3(mA, mB, res);
                        PrintMatrix(res);
                    }
                    if (sub == "b") Console.WriteLine($"Det = {NativeMethods.MatrixDeterminant3x3(mA)}");
                    if (sub == "c")
                    {
                        double[] inv = new double[9];
                        if (NativeMethods.MatrixInverse3x3(mA, inv) == 1) PrintMatrix(inv);
                        else Console.WriteLine("Macierz osobliwa (Det=0).");
                    }
                    break;
            }
            Pause();
        }

        // =============================================================
        //                      DZIAŁ 6: PROB & STAT
        // =============================================================
        static void MenuProbStat()
        {
            Console.Clear();
            DrawHeader("PRAWDOPODOBIEŃSTWO I STATYSTYKA", ConsoleColor.DarkYellow);

            Console.WriteLine("1. Kombinatoryka (Silnia, Newton, Wariacje)");
            Console.WriteLine("2. Schemat Bernoulliego");
            Console.WriteLine("3. Statystyka (Średnia, Odchylenie Std)");

            Console.Write("\nWybór: ");
            switch (Console.ReadLine())
            {
                case "1":
                    int n = (int)ReadDouble("n"), k = (int)ReadDouble("k");
                    Console.WriteLine($"n! = {NativeMethods.Factorial(n)}");
                    Console.WriteLine($"Newton (n po k) = {NativeMethods.BinomialCoefficient(n, k)}");
                    Console.WriteLine($"Wariacje bez powt. = {NativeMethods.VariationsNoRepetition(n, k)}");
                    break;
                case "2":
                    Console.WriteLine($"Prawdopodobieństwo: {NativeMethods.BernoulliTrial((int)ReadDouble("Liczba prób n"), (int)ReadDouble("Sukcesy k"), ReadDouble("Prawd. p"))}");
                    break;
                case "3":
                    double[] data = ReadArray("Podaj dane (po spacji)");
                    Console.WriteLine($"Średnia: {NativeMethods.CalculateMean(data.Length, data)}");
                    Console.WriteLine($"Odchylenie: {NativeMethods.CalculateStandardDeviation(data.Length, data)}");
                    break;
            }
            Pause();
        }

        // =============================================================
        //                      DZIAŁ 7: STUDIA
        // =============================================================
        static void MenuUniversity()
        {
            Console.Clear();
            DrawHeader("POZIOM STUDENCKI (UNIVERSITY)", ConsoleColor.Red);

            Console.WriteLine("1. Liczby Zespolone (Działania podstawowe)");
            Console.WriteLine("2. Regresja Liniowa (Dopasowanie prostej)");
            Console.WriteLine("3. Całkowanie Numeryczne (Metoda Trapezów)");

            Console.Write("\nWybór: ");
            switch (Console.ReadLine())
            {
                case "1":
                    Console.WriteLine("Z1 = a+bi, Z2 = c+di");
                    double r1 = ReadDouble("a"), i1 = ReadDouble("b"), r2 = ReadDouble("c"), i2 = ReadDouble("d");
                    NativeMethods.Complex_Multiply(r1, i1, r2, i2, out double rm, out double im);
                    Console.WriteLine($"Mnożenie: {rm} + {im}i");
                    NativeMethods.Complex_Divide(r1, i1, r2, i2, out double rd, out double id);
                    Console.WriteLine($"Dzielenie: {rd} + {id}i");
                    Console.WriteLine($"Moduł Z1: {NativeMethods.Complex_Modulus(r1, i1)}");
                    break;
                case "2":
                    int n = (int)ReadDouble("Ile punktów?");
                    double[] x = new double[n], y = new double[n];
                    for (int i = 0; i < n; i++) { x[i] = ReadDouble($"x{i}"); y[i] = ReadDouble($"y{i}"); }
                    NativeMethods.LinearRegression(n, x, y, out double a, out double b, out double rsq);
                    Console.WriteLine($"Model: y = {a:F4}x + {b:F4}");
                    Console.WriteLine($"R^2: {rsq:F4}");
                    break;
                case "3":
                    double[] vals = ReadArray("Wartości funkcji y[]");
                    double h = ReadDouble("Krok h");
                    Console.WriteLine($"Całka: {NativeMethods.Integrate_Trapezoidal(vals.Length, vals, h)}");
                    break;
            }
            Pause();
        }

        // --- POMOCNICZE ---
        static void DrawHeader(string title, ConsoleColor color = ConsoleColor.Cyan)
        {
            Console.ForegroundColor = color;
            Console.WriteLine("=============================================");
            Console.WriteLine($"   {title}");
            Console.WriteLine("=============================================");
            Console.ResetColor();
        }

        static double ReadDouble(string prompt)
        {
            Console.Write($" {prompt}: ");
            string s = Console.ReadLine().Replace('.', ',');
            if (double.TryParse(s, out double res)) return res;
            return 0;
        }

        static double[] ReadArray(string prompt)
        {
            Console.WriteLine($" {prompt}:");
            try
            {
                var parts = Console.ReadLine().Replace('.', ',').Split(new[] { ' ', ',' }, StringSplitOptions.RemoveEmptyEntries);
                double[] arr = new double[parts.Length];
                for (int i = 0; i < parts.Length; i++) arr[i] = double.Parse(parts[i]);
                return arr;
            }
            catch { return new double[0]; }
        }

        static void PrintMatrix(double[] m)
        {
            for (int i = 0; i < 3; i++)
                Console.WriteLine($"| {m[i * 3]:F2} {m[i * 3 + 1]:F2} {m[i * 3 + 2]:F2} |");
        }

        static void Pause()
        {
            Console.ForegroundColor = ConsoleColor.DarkGray;
            Console.WriteLine("\n[Enter] aby wrócić...");
            Console.ResetColor();
            Console.ReadKey();
        }
    }
}