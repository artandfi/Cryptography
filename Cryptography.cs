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

        public static LongNumber LogPollard(LongNumber a, LongNumber b, LongNumber p) {
            var partition = MakePartition(p);
            var ai = new LongNumber(0);
            var bi = ai;
            var xi = new LongNumber(1);
            var a2i = ai;
            var b2i = bi;
            var x2i = xi;

            while (true) {
                xi = F(xi, a, b, partition);
                ai = G(xi, ai, partition);
                bi = H(xi, bi, partition);

                var f = F(x2i, a, b, partition);
                x2i = F(f, a, b, partition);
                a2i = G(f, G(x2i, a2i, partition), partition);
                b2i = H(f, H(x2i, b2i, partition), partition);

                if (xi == x2i) {
                    var r = bi - b2i;
                    if (r == 0) {
                        return null;
                    }

                    return (MulInverse(r, p) * (a2i - ai)) % p;
                }
            }
        }

        public static LongNumber LogBabyStepGiantStep(LongNumber a, LongNumber b, LongNumber n) {
            a %= n;
            b %= n;
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
            if (n < 1) {
                return null;
            }
            
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

        public static int? Legendre(LongNumber a, LongNumber p) {
            if (p < 3 || !IsPrime(p))
                return null;

            if (a % p == 0)
                return 0;

            return PowMod(a, (p - 1) / 2, p) == 1 ? 1 : -1;
        }

        public static int? Jacobi(LongNumber a, LongNumber b) {
            if (b < 1 || b % 2 == 0)
                return null;

            if (Gcd(a, b) != 1)
                return 0;

            a %= b;
            var t = new LongNumber(1);
            while (a != 0) {
                while (a % 2 == 0) {
                    a /= 2;
                    var r = b % 8;
                    if (r == 3 || r == 5)
                        t = -t;
                }
                Swap(ref a, ref b);

                if (a % 4 == 3 && b % 4 == 3)
                    t = -t;
                a %= b;
            }
            return b == 1 ? t : new LongNumber(0);
        }

        public static (LongNumber, LongNumber) SqrtCipolla(LongNumber n, LongNumber p) {
            if (Legendre(n, p) != 1) {
                return (null, 0);
            }

            LongNumber a = 0;
            LongNumber w2;
            while (true) {
                w2 = (a * a + p - n) % p;
                if (Legendre(w2, p) != 1)
                    break;
                a++;
            }

            var finalW = w2;
            (LongNumber, LongNumber) MulExtended((LongNumber, LongNumber) aa, (LongNumber, LongNumber) bb) {
                return ((aa.Item1 * bb.Item1 + aa.Item2 * bb.Item2 * finalW) % p,
                        (aa.Item1 * bb.Item2 + bb.Item1 * aa.Item2) % p);
            }

            var r = (new LongNumber(1), new LongNumber(0));
            var s = (a, new LongNumber(1));
            var nn = (p + 1) / 2;
            while (nn > 0) {
                if (nn % 2 != 0) {
                    r = MulExtended(r, s);
                }
                s = MulExtended(s, s);
                nn /= 2;
            }

            if (r.Item2 != 0 || r.Item1 * r.Item1 % p != n) {
                return (0, null);
            }

            return (r.Item1, p - r.Item1);
        }

        public static string IsPrimeSolovayStrassen(LongNumber n, int k) {
            if (n < 2) {
                return "The number " + n + " is less than 2, hence isn't prime.";
            }

            if (n % 2 == 0) {
                return "The number " + n + " is even, hence isn't prime.";
            }

            for (int i = 0; i < k; i++) {
                var a = Rand(2, n);
                if (Gcd(a,n) > 1 || PowMod(a, (n - 1) / 2, n) != Legendre(a, n) % n) {
                    return "The number " + n + " is odd and isn't prime.";
                }
            }

            double prob = 1 - Math.Pow(2, -k);
            return "The number " + n + " is prime with propability " + prob + ".";
        }

        #region Inner methods
        public static bool IsPrime(LongNumber n) {
            if (n == 2)
                return true;

            if (n % 2 == 0)
                return false;

            var t = Sqrt(n);
            for (LongNumber k = 3; k <= t; k += 2) {
                if (n % k == 0)
                    return false;
            }

            return true;
        }
        private static LongNumber F(LongNumber x, LongNumber n) => (x * x + 1) % n;

        private static LongNumber[] MakePartition(LongNumber p) {
            var r = p % 3;
            var d = p / 3;

            if (r == 0) {
                return new LongNumber[] { d, d * 2, d * 3 };
            }
            if (r == 1) {
                return new LongNumber[] { d, d * 2 + 1, d * 3 + 1 };
            }

            return new LongNumber[] { d + 1, d * 2 + 1, d * 3 + 2 };
        }

        private static LongNumber F(LongNumber x, LongNumber a, LongNumber b, LongNumber[] partition) {
            if (x >= 0 && x < partition[0]) {
                return (b * x) % partition[2];
            }

            if (x >= partition[0] && x < partition[1]) {
                return (x * x) % partition[2];
            }

            if (x >= partition[1] && x < partition[2]) {
                return (a * x) % partition[2];
            }

            return null;
        }
        
        private static LongNumber G(LongNumber x, LongNumber n, LongNumber[] partition) {
            if (x >= 0 && x < partition[0]) {
                return n % partition[2];
            }

            if (x >= partition[0] && x < partition[1]) {
                return (n * 2) % partition[2];
            }

            if (x >= partition[1] && x < partition[2]) {
                return (n + 1) % partition[2];
            }

            return null;
        }

        private static LongNumber H(LongNumber x, LongNumber n, LongNumber[] partition) {
            if (x >= 0 && x < partition[0]) {
                return (n + 1) % partition[2];
            }

            if (x >= partition[0] && x < partition[1]) {
                return (n * 2) % partition[2];
            }

            if (x >= partition[1] && x < partition[2]) {
                return n % partition[2];
            }

            return null;
        }

        #endregion
    }
}