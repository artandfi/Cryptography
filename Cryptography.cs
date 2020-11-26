using System;
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

        private static LongNumber F(LongNumber x, LongNumber n) => (x * x + 1) % n;
    }
}
