import math
import sys

import numpy as np
from MandatoryAssignment1.A1.Mandatory1 import BinaryExponentiationWithoutRecursion as binexp
from MandatoryAssignment3.Code import createBsmooth


def posA(p: int):
    a = 2
    while pow(a, p - 1, p) != 1:  # mulig binexp !
        a += 1
    return a


def prime_factors(n):
    i = 2
    factors = []
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1:
        factors.append(n)
    if factors:
        return sorted(factors)
    else:
        return [0]


def pocklington(N: int) -> bool:
    """
    :param N: n > 1
    :return:
    """
    assert N > 1 and N % 2 == 1

    # a^(n-1) congruent 1 mod n
    # p | n-1 and p > sqrt(n) -1
    # gcd(a ^(n-1)/p - 1 == 1
    # then n prime

    primes_ = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]

    for q in primes_:
        for a in range(0, N - 1):
            if (N - 1) % q == 0:
                if ((a ** (N - 1)) % N) == 1:
                    if math.gcd(a ** ((N - 1) // q) - 1, N) == 1:
                        return True
    return False


print(pocklington(147))
# print(prime_factors(107, 9999))
