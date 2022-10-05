import math
from typing import Any
import time
from fractions import Fraction

import gmpy2 as gmpy2


def get_Quotients(a: int, b: int) -> list:
    """
    Perform Extended Euclidean Algorithm, to find quotients

    :param a: int a
    :param b: int b
    :return: quotients
    """
    s, old_s = 0, 1
    r, old_r = b, a

    q_n = []
    while r != 0:
        q_i = old_r // r
        q_n.append(q_i)
        old_r, r = r, old_r - q_i * r
        old_s, s = s, old_s - q_i * s

    return q_n


def find_roots(a, b, c):
    """
    Solving linear equation
    :param phi:
    :param N:
    :return:
    """
    # print(f"a: {a}, b: {b}, c: {c}")
    dis = (b ** 2) - 4 * a * c
    sqrt_val = gmpy2.isqrt(dis)  # math.sqrt(abs(dis))

    pluss = -b + sqrt_val

    negative = -b - sqrt_val

    x_1 = int(pluss // 2)  # Fraction(-b + sqrt_val, 2)
    x_2 = int(negative // 2)  # Fraction(-b - sqrt_val, 2)
    return x_1, x_2


def valid_roots(phi, N):
    b = ((N + 1 - phi) ** 2)
    d = (4 * 1 * N)
    return b > d and math.sqrt(abs(b - d)).is_integer()


def ContinuedFractionAlgorithm(N: int, e: int):
    """
    Let N, e be an RSA pub key, where ed congruent with 1 mod phi(N)
    d is RSA secret exponent.
    We know d is rel. small.

    Task: Factor N = pq
    :return:
    """

    q_n = get_Quotients(e, N)
    P_n = [1] + [q_n[1]] + [0] * (len(q_n) - 2)
    Q_n = [0, 1] + [0] * (len(q_n) - 2)
    # P_n = q_n * P_(n-1) + P_(n-2)
    # Q_n = q_n * Q_(n-1) + P_(n_2)

    for n in range(2, len(q_n)):
        P_n[n] = q_n[n] * P_n[n - 1] + P_n[n - 2]
        Q_n[n] = q_n[n] * Q_n[n - 1] + Q_n[n - 2]

    # phi equation
    #   1. not integer, try another
    #   2. integer
    #       x^2 - (N+1 - phi)x + N
    P_n, Q_n = P_n[1:], Q_n[1:]
    P_n, Q_n = Q_n, P_n  # now: P = Numerator, Q = Denominator

    for P_i, Q_i in zip(P_n, Q_n):
        if P_i == 0 or Q_i == 0 or Q_i % 2 == 0 or e * Q_i % P_i != 1:
            continue
        phi = (e * Q_i - 1) // P_i
        if valid_roots(phi, N):
            x_1, x_2 = find_roots(1, (N - int(phi) + 1), N)
            print(f"factor: {x_1} , {x_2}")
            # return abs(x_1), abs(x_1)
    return 0, 0


def g_X(x, P: int):
    return (x ** 2 + x) % P


def FactorByP_Method(N: int) -> str | int | Any:
    x = 2
    y = 2
    d = 1
    while d == 1:
        x = g_X(x, N)
        y = g_X(g_X(y, N), N)
        d = math.gcd(abs(x - y), N)  # md1.GCD(abs(x - y), N)

    if d == N:
        return "failed"
    else:
        return d


def DixonsMethod():
    pass


if __name__ == '__main__':
    N = 1098667602555738997359355609545017255443451483195436981354791765341639135658156206242197992115989996829728203054347117299
    e = 415884005149779743103130951097941968793336745983309971301432658775763996247677181243042840232106535367251782466233724389

    # e = 42667  # 388511
    # N = 64741  # 331879
    a, b = ContinuedFractionAlgorithm(N, e)
    print(a, b)
    p_method_N = 10862216162096506735513546937

    # var = ContinuedFractionAlgorithm(int(N), int(e))
    # print(continued_fraction(int(N) / int(e)))
    # start = time.time()
    # p_method_factor = FactorByP_Method(p_method_N)
    # print(f"Time took: {time.time() - start} and result was: {p_method_factor}")
