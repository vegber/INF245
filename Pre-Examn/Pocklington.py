import math
import random
import sys
import time
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

    # primes_ = [int(x) for x in open("primes.txt",'r+').read().split(",")]
    primes_ = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]

    for q in primes_:
        for a in range(0, N - 1):
            if (N - 1) % q == 0:
                if binexp(a, (N - 1), N) == 1:
                    # if ((a ** (N - 1)) % N) == 1:
                    tmp = binexp(a, ((N - 1) // q) - 1, N)
                    if math.gcd(tmp, N) == 1:  # a ** ((N - 1) // q) - 1
                        return True
    return False


start = time.time()
for x in range(201, 1000, 2):
    pocklington(x)
print(f"time: {time.time() - start}")


def driverCode(q1=1299709, bit_length=160):
    qi = q1
    h = 0
    while qi.bit_length() != bit_length:
        # find qi+1 == 2*hi*qi + 1
        for hi in range(1, qi // 2):
            qi2 = int(2 * hi * qi) + 1  # 2*h*qi + 1
            isprime = pocklington(qi2)
            if isprime and qi2 < pow(qi, 2):
                print(qi2.bit_length())
                qi = qi2
                h = hi
                break
        # print(f"Found new candidate {qi}")
    print(h)
    return qi


def driverAlt(q1=1299709, bit_length=160):
    q = q1
    h = 1
    while q.bit_length() < bit_length:
        tmp = int(2 * h * q + 1)
        if pocklington(tmp) and tmp < pow(q, 2):
            q = tmp
            h = 1
            continue
        h += 1
    return q


def drivercodetask3(q=905765673186944881352668387379493810404306452479,
                    s=67325197105509819051461517808547342938036883286575830551267881971120479404031,
                    bit_length=512):
    p = (s * q + 1) * 99
    h = 1
    escape_local_mimima = 1
    stuck_in_bit_size = (0,0)
    while p.bit_length() < bit_length:
        tmp = int(2 * h * q * s + 1)
        if random.randint(0, 10) > 2:
            print(tmp.bit_length())
        if pocklington(tmp):
            h = 1
            stuck_in_bit_size = tmp.bit_length(), stuck_in_bit_size[1]
            if stuck_in_bit_size[0] == tmp.bit_length():
                stuck_in_bit_size = tmp.bit_length(), stuck_in_bit_size[1] + 1
            if stuck_in_bit_size[1] > 1000:
                escape_local_mimima += escape_local_mimima
                h = escape_local_mimima
            p = tmp
            continue
        h += 1
    return p


a = drivercodetask3()
print(f"number  was: {a} with length: {a.bit_length()}")
# print(f"Found 160 bit prime: {driverCode()}")
# print(f"Found 256 bit prime:  {driverCode(q1=4994460703, bit_length=256)}")
# print(f"Found 512 bit prime: {drivercodetask3()}")
