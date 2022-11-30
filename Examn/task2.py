import math
import random

"""
Factor RSA modulus N given RSa public and private exponents
    :parameter: e = 3
    :parameter: d = variable from exam
"""

# my values are:
N = 15541351133097319658694005067857282074739726828480643400789860273748542845849754037903
d = 10360900755398213105796003378571521383159812310467463479085539107171089161979092336827
e = 3


def BinaryExponentiation(a: int, b: int, mod: int):
    # Modulo function to perform binary exponentiation
    temp, base_number = 1, a
    while b > 0:  # exponent larger than zero
        if b % 2 == 1:  # if not even
            temp = (temp * base_number) % mod  # pow(x*y, 1, mod)
        base_number, b = (base_number * base_number) % mod, b // 2
    return temp % mod


def factorN(N, e, d):
    """
    :param N:
    :param e:
    :param d:
    :return: pq = N
    """
    k = (d * e) - 1
    t = k
    while True:
        g = random.randint(2, N - 1)
        if t % 2 == 0:
            t = t // 2
            x = BinaryExponentiation(g, t, N)
            y = math.gcd(x - 1, N)
            if x > 1 and y > 1:
                p = y
                q = N // y
                return p, q


def verify(p, q, N):
    return p * q == N


# Solution Found:
# (p) 5575186298788107701075327939401920743688443
# (q) 2787593149394053850537663969700960371844221

if __name__ == '__main__':
    p, q = factorN(N, e, d)
    print(f"Found possible p {p} and q {q}")
    print(f"Is it correct?\np*q=N ? {verify(p, q, N)}")
