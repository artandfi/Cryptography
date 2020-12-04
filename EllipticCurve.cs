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


        public EllipticCurve(ref LongNumber p, ref LongNumber a, ref LongNumber b, ref LongNumber g, ref LongNumber n) {
            P = p;
            A = a;
            B = b;
            G = g;
            N = n;
        }

        public (LongNumber, LongNumber) AddPoints(LongNumber x1, LongNumber y1, LongNumber x2, LongNumber y2) {
            if (!IsPointOnCurve(ref x1, ref y1) || !IsPointOnCurve(ref x2, ref y2)) {
                return (null, null);
            }

            bool pointsEqual = PointsEqual(ref x1, ref y1, ref x2, ref y2);

            if (pointsEqual) {
                if (y1 == 0) {
                    return (-1, -1);
                }

                var m = MulMod(y2 - y1, MulInverse(x2 - x1, P), P);
                var x3 = (m * m - x1 - x2) % P;
                return (x3, (-y1 + m * (x1 - x3)) % P);
            }
            else {
                if (x1 == x2) {
                    return (-1, 1);
                }

                var m = MulMod(3 * x1 * x1 + A, MulInverse(2 * y1, P), P);
                var x3 = (m * m - 2 * x1) % P;
                return (x3, (-y1 + m * (x1 - x3)) % P);
            }
        }

        private bool IsPointOnCurve(ref LongNumber x, ref LongNumber y) {
            return x >= 0 && y >= 0 && x < P && y < P &&
                   y * y == (Pow(x, 3) + A * x + B) % P;
        }

        private bool PointsEqual(ref LongNumber x1, ref LongNumber y1, ref LongNumber x2, ref LongNumber y2) {
            return x1 == x2 && y1 == y2;
        }
    }
}
