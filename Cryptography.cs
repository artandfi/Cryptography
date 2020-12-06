using System;
using System.Collections;
using System.Collections.Generic;
using LongArithmetics;
using static LongArithmetics.LongNumber;

namespace Cryptography {
    public static class Cryptography {
        private const int ALPH_BASE = 28;
        private const int ALPH_DIFF = 64;
        
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
            a %= n;
            b %= n;
            var m = Sqrt(n) + 1;
            var g0 = PowMod(a, m, n);
            var g = g0;
            var t = new Dictionary<string, LongNumber>();

            for (var i = new LongNumber(1); i <= m; i++) {
                t.Add(g.ToString(), i);
                g = MulMod(g, g0, n);
            }

            for (var j = new LongNumber(0); j < m; j++) {
                var y = MulMod(b, PowMod(a, j, n), n);
                if (t.ContainsKey(y.ToString())) {
                    return m * t[y.ToString()] - j;
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

        public static void ElGamal(EllipticCurve curve) {
            // Step 1. Bob chooses a random number k = 1, ..., N-1.
            var k = Rand(1, curve.N);
            var Y = curve.PointSelfSum(k, curve.G);

            // Step 2. Alice's message
            (LongNumber x, LongNumber y) M;
            Console.Write("Enter Alice's message: ");
            M.x = MessageToLongNumber(Console.ReadLine());
            M.y = SqrtCipolla(Pow(M.x, 3) + curve.A * M.x + curve.B, curve.P).Item1;
            
            // Step 3. Encryption
            var r = Rand(1, curve.N);
            var D = curve.PointSelfSum(r, Y);
            var G = curve.PointSelfSum(r, curve.G);
            var H = curve.AddPoints(M, D);

            // Step 4. Decryption
            var S = curve.PointSelfSum(k, G);
            var S1 = (S.Item1, (S.Item1 + S.Item2) % curve.P);
            var M1 = curve.AddPoints(S1, H);
            var res = LongNumberToMessage(M1.Item1);
            Console.WriteLine("Message decrypted: " + res);
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

        public static LongNumber MessageToLongNumber(string msg) {
            LongNumber res = 0;
            for (int k = msg.Length - 1; k >= 0; k--) {
                var cur = msg[msg.Length - k - 1];
                var c = cur == ' ' ? 0 : cur - ALPH_DIFF;
                res += c * Pow(ALPH_BASE, k);
            }

            return res;
        }

        public static string LongNumberToMessage(LongNumber n) {
            string res = "";
            int r;
            do {
                r = n % ALPH_BASE;
                res += r == 0 ? ' ' : (char)(r + ALPH_DIFF);
                n /= ALPH_BASE;
            } while (n >= ALPH_BASE);
            r = n;
            res += r == 0 ? ' ' : (char)(r + ALPH_DIFF);

            var arr = res.ToCharArray();
            Array.Reverse(arr);
            return new string(arr);
        }
        #endregion
    }
}