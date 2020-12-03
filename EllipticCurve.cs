using LongArithmetics;
using static LongArithmetics.LongNumber;
using System;
using System.Collections.Generic;
using System.Text;

namespace Cryptography {
    public class EllipticCurve {
        // y^2 = x^3 + Ax + B (mod P)
        public LongNumber P { get; set; }
        public LongNumber A { get; set; }
        public LongNumber B { get; set; }
        public LongNumber G { get; set; }   // Generator (base point)
        public LongNumber N { get; set; }   // Order of the curve

        public EllipticCurve(LongNumber p, LongNumber a, LongNumber b, LongNumber g, LongNumber n) {
            P = p;
            A = a;
            B = b;
            G = g;
            N = n;
        }

        private bool IsPointOnCurve((LongNumber x, LongNumber y) p) {
            return p.x >= 0 && p.y >= 0 && p.x < P && p.y < P &&
                   p.y * p.y == (Pow(p.x, 3) + A * p.x + B) % P;
        }


    }
}
