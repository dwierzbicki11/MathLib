using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

namespace MathApp
{
    internal class Program
    {
        [DllImport("MathCore.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern int Dodaj(int a, int b);
        static void Main(string[] args)
        {
            Console.WriteLine(Dodaj(5, 7));
        }
    }
}
