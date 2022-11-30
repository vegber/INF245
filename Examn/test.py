import math
def isGermainPrime( q):
    """Determines if q is a Sophie Germain prime"""
    p = 2 * q + 1  # p satisfies q | p-1
    a = candidateTo(p)  # s satisfies a**(p-1) % p == 1 mod p

    return math.gcd(a * a - 1, p) == 1


def candidateTo(p):
    """Return an a such that a**(p-1) % p == 1"""
    a = 2
    while pow(a, p - 1, p) != 1:
        a += 1
    return a


def powmod(a, exponent, modulus):
    """Return a**exponent % modulus"""
    result, power = 1, a
    while exponent > 0:
        if exponent % 2 != 0:
            result = (result * power) % modulus
        exponent, power = exponent / 2, power ** 2 % modulus
    return result


def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a


isGermainPrime(536322828007183377480253606803869765986172446936943)
