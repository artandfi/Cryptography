using System;
using System.Collections;
using LongArithmetics;
using static LongArithmetics.LongNumber;

namespace Cryptography {
    public static class Cryptography {


        public static (LongNumber, LongNumber) FactorizePollard(LongNumber n) {
            var rnd = new Random();
            LongNumber x = n < int.MaxValue ? rnd.Next(0, n) : rnd.Next(0, int.MaxValue);
            LongNumber y = x;
            LongNumber d = new LongNumber(1);

            while (d == 1) {
                x = F(x, n);
                y = F(F(y, n), n);
                d = Gcd(n, Abs(x - y));
            }

            return d == n ? (null, null) : (d, n / d);
        }

        public static LongNumber LogBabyStepGiantStep(LongNumber a, LongNumber b, LongNumber n) {
            var m = Sqrt(n) + 1;
            var g0 = PowMod(a, m, n);
            var g = g0;
            var t = new Hashtable();

            for (var i = new LongNumber(1); i <= m; i++) {
                t.Add(i, g);
                g = MulMod(g, g0, n);
            }

            for (var j = new LongNumber(1); j <= m; j++) {
                foreach (LongNumber key in t.Keys) {
                    if ((LongNumber)t[key] == MulMod(b, PowMod(a, j, n), n)) {
                        return m * key - j;
                    }
                }
            }

            return null;
        }

        public static LongNumber Phi(LongNumber n) {
            var res = n;
            
            for (var i = new LongNumber(2); i * i <= n; ++i) {
                if (n % i == 0) {
                    while (n % i == 0) {
                        n /= i;
                    }
                    res -= res / i;
                }
            }

            return n > 1 ? res - res / n : res;
        }

        public static int Mu(LongNumber n) {
            if (n == 1)
                return 1;

            var p = new LongNumber(0);
            for (var i = new LongNumber(2); i <= n; i++) {
                if (n % i == 0 && IsPrime(i)) {
                    if (n % (i * i) == 0)
                        return 0;

                    p++;
                }
            }

            return p % 2 == 0 ? 1 : -1;
        }


        private static bool IsPrime(LongNumber n) {
            if (n == 2)
                return true;

            if (n % 2 == 0)
                return false;

            var t = Sqrt(n);
            for (var k = new LongNumber(3); k <= t; k += 2) {
                if (n % k == 0)
                    return false;
            }

            return true;
        }
        private static LongNumber F(LongNumber x, LongNumber n) => (x * x + 1) % n;
    }
}
