import math
import random
import sympy

"""
Construct a DSA prime p with Pocklington primality test and give a proof that p is a DSA prime.
The prime p should be a 512-bit integer. The prime q such that q divides p−1 should be a
    160-bit integer.
Construct a residue g modulo p of order q.
    1. Start the construction with given prime q1 as a seed.
    Find a prime q2 = 2h1q1+ 1,where q2 < q1^2 for some h1,
    with Pocklington test.
    Find a prime q3= 2h2q2+ 1, where q3 < q22 for some h2,
    etc., up to a 160-bit prime q is reached.2.

    2. Similarly, construct a (>256)-bit primes with using given
    prime s1 as a seed.

    3. Find a 512-bit prime p= 2hqs+ 1 for some h with Pocklington test.
    Construct g.
    Keep track of successful Pocklington tests and present that in the report as a proof
"""


def BinExp(a: int, b: int, mod: int):
    # Modulo function to perform binary exponentiation
    temp, base_number = 1, a
    while b > 0:  # exponent larger than zero
        if b % 2 == 1:  # if not even
            temp = (temp * base_number) % mod  # pow(x*y, 1, mod)
        base_number, b = (base_number * base_number) % mod, b // 2
    return temp % mod


def pocklington(N: int) -> bool:
    """
    :param N: n > 1
    :return:
    """
    assert N > 1 and N % 2 == 1

    # a^(n-1) congruent 1 mod n
    # p | n-1 and p > sqrt(n) -1
    # gcd(a ^(n-1)/p - 1 == 1
    # then n prime-
    primes_ = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103,
               107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223,
               227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347,
               349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
               467, 479, 487, 491, 499, 503]

    for p in primes_:
        for a in range(1, N - 1):
            if (N - 1) % p == 0:
                if BinExp(a, (N - 1), N) == 1:
                    # if ((a ** (N - 1)) % N) == 1:
                    tmp = BinExp(a, ((N - 1) // p), N)
                    if math.gcd(tmp, N) == 1:  # a ** ((N - 1) // q) - 1
                        return True
    return False


def driverCodePocklingtonTask(q, s, bit_length=512):
    p = (s * q + 1)
    h = 1
    escape_local_minima = 1
    stuck_in_bit_size = (0, 0)
    while p.bit_length() < bit_length:
        tmp = int(2 * h * q * s + 1)
        if pocklington(tmp) and sympy.isprime(tmp):
            h = 1
            stuck_in_bit_size = tmp.bit_length(), stuck_in_bit_size[1]
            if stuck_in_bit_size[0] == tmp.bit_length():
                stuck_in_bit_size = tmp.bit_length(), stuck_in_bit_size[1] + 1
            if stuck_in_bit_size[1] > 1000:
                escape_local_minima += escape_local_minima
                h = escape_local_minima
            p = tmp
            continue
        h += 1
    return p


def recoverG(p, q):
    while True:
        r = random.randrange(1, p)
        g = BinExp(r, (p - 1) // q, p)
        if g != 1:
            return g


if __name__ == '__main__':
    # my values are:
    q1 = 1299709
    s1 = 4994460703

    _find160BitPrime = driverCodePocklingtonTask(q1, 1, 160)  # driverCodeV2(q1, 160)
    _find256BitPrime = driverCodePocklingtonTask(1, s1, 256)

    if _find160BitPrime is not None:
        print(f"Successfully found 160 - bit prime: {_find160BitPrime}, with len {_find160BitPrime.bit_length()}")

    if _find256BitPrime is not None:
        print(f"Successfully found 256 - bit prime: {_find256BitPrime}, with len {_find256BitPrime.bit_length()}")

    _find512BitPrime = driverCodePocklingtonTask(_find160BitPrime, _find256BitPrime)

    if _find512BitPrime is not None:
        print(f"Successfully found 512 - bit prime: {_find512BitPrime}, with len {_find512BitPrime.bit_length()}")
        g = recoverG(_find512BitPrime, _find160BitPrime)
        print(f"Also found g:  {g}")
