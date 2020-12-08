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
        public LongNumber N { get; set; }                     // Order of the curve
        public (LongNumber x, LongNumber y) G { get; set; }   // Generator (base point)


        public EllipticCurve(LongNumber p, LongNumber a, LongNumber b, LongNumber n, LongNumber gx, LongNumber gy) {
            P = p;
            A = a;
            B = b;
            G = (gx, gy);
            N = n;
        }

        public (LongNumber, LongNumber) AddPoints(LongNumber x1, LongNumber y1, LongNumber x2, LongNumber y2) {
            if (!IsPointOnCurve(x1, y1) || !IsPointOnCurve(x2, y2)) {
                return (null, null);
            }

            if (PointsEqual(ref x1, ref y1, ref x2, ref y2)) {
                if (y1 == 0 || y1 == -1) {
                    return (-1, -1);
                }

                var m = MulMod(3 * x1 * x1 + A, MulInverse(2 * y1, P), P);
                var x3 = (m * m - 2 * x1) % P;
                return (x3, (-y1 + m * (x1 - x3)) % P);
            }
            else {
                if (x1 == x2) {
                    return (-1, 1);
                }

                if (x1 == -1) {
                    return (x2, y2);
                } 

                if (x2 == -1) {
                    return (x1, y1);
                }

                var m = MulMod(y2 - y1, MulInverse(x2 - x1, P), P);
                var x3 = (m * m - x1 - x2) % P;
                return (x3, (-y1 + m * (x1 - x3)) % P);
            }
        }

        public (LongNumber, LongNumber) AddPoints((LongNumber, LongNumber) p1, (LongNumber, LongNumber) p2) {
            return AddPoints(p1.Item1, p1.Item2, p2.Item1, p2.Item2);
        }

        private bool IsPointOnCurve(LongNumber x, LongNumber y) {
            return x == -1 && y == -1 ||
                   x >= 0 && y >= 0 && x < P && y < P &&
                   (y * y) % P == (Pow(x, 3) + A * x + B) % P;
        }

        private bool PointsEqual(ref LongNumber x1, ref LongNumber y1, ref LongNumber x2, ref LongNumber y2) {
            return x1 == x2 && y1 == y2;
        }

        public (LongNumber, LongNumber) PointSelfSum(LongNumber k, (LongNumber, LongNumber) p) {
            var res = p;
            for (LongNumber i = 1; i < k; i++) {
                res = AddPoints(res, p);
            }

            return res;
        }

        public override string ToString() {
            var aStr = A == 0 ? "" : A == 1 ? " + x" : " + " + A + "x";
            var bStr = B == 0 ? "" : " + " + B;
            return "y^2 = x^3" + aStr + bStr + " (mod " + P + ")";
        }
    }
}
