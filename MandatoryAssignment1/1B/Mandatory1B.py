"""
Implement arithmetic operations with polynomials g(x) & f(x) under
modulus a prime number p.
"""

import re
from turtle import end_fill
from typing import Any


def Polynomial_addition(f_x: list, g_x: list, p: int) -> list:
    f_x, g_x = EqualizeLengthOnOrder(f_x, g_x)

    return [(x + y) % p for x, y in zip(f_x, g_x)]


def EqualizeLengthOnOrder(f_x, g_x):
    # this is naive approach! todo
    # if deg(f) != deg(g)
    if len(f_x) < len(g_x):
        f_x = [0] + f_x

    elif len(g_x) < len(f_x):
        g_x = [0] + g_x
    return f_x, g_x


def Polynomial_multiplication(f_x: list, g_x: list, p: int) -> list:
    init = [0] * (len(f_x) + len(g_x) - 1)
    for f in range(len(f_x)):
        for g in range(len(g_x)):
            init[f + g] += (f_x[f] * g_x[g]) % p
    return init


def normalize(poly):
    while poly and poly[-1] == 0:
        poly.pop()
    if not poly:
        poly.append(0)


def degree(poly):
    while poly and poly[-1] == 0:
        poly.pop()  # normalize
    return len(poly) - 1


def poly_div(N: list, D: list, p: int):
    dD = degree(D)
    dN = degree(N)
    if dD < 0: raise ZeroDivisionError
    if dN >= dD:
        q = [0] * (dN - 1)
        while dN >= dD:
            d = [0] * (dN - dD) + D
            mult = q[dN - dD] = modDivide(N[-1], d[-1], p)  # N[-1] / float(d[-1])
            d = [(coeff * mult) % p for coeff in d]
            N = [(coeffN - coeffd) % p for coeffN, coeffd in zip(N, d)]
            dN = degree(N)
        r = N
    else:
        q = [0]
        r = N
    return q, r


# Function to compute a/b under modulo m
def modDivide(a: int, b: int, m: int) -> Exception | int | float | Any:
    a = a % m
    inv = modInverse(b, m)
    if inv == -1:
        return 0
    else:
        return (inv * a) % m


def modInverse(b: int, m: int) -> int:
    import math
    g = math.gcd(b, m)
    if g != 1:
        # print("Inverse doesn't exist")
        return -1
    else:
        # If b and m are relatively prime,
        # then modulo inverse is b^(m-2) mode m
        return pow(b, m - 2, m)


def PolyDiv(f_x, g_x, p: int):
    import copy
    # make sure deg(f_x) > deg(g_x)
    assert len(f_x) >= len(g_x)
    g_x = [0] * (len(f_x) - len(g_x)) + g_x
    # Copy input arr. so we can change them
    f_x = copy.deepcopy(f_x)
    g_x = copy.deepcopy(g_x)

    i_th_elem = 0
    # find index of element with the highest order != 0
    highest_order_elemt_index = [i for i, e in enumerate(f_x) if e != 0][0]
    order_of_highest_elemt = len(f_x) - highest_order_elemt_index
    while True:
        if order_of_highest_elemt < len(g_x):
            # can't divide this element, rest is remainder
            return
        # find a s.t a*g[0] == inv(f_x(0)) = 0
        theta = modDivide(f_x[highest_order_elemt_index], g_x[1], p)

        f_x = [(f + (modInverse((theta * g), p))) % p for f, g in zip(f_x, g_x)]

        print(f_x)

        break


def PrintPolynomial(h_x: list):
    out = ""
    for i, e in enumerate(h_x):
        if i == len(h_x) - 1:
            out += str(e)
        else:
            out += f"{e}x^{(len(h_x) - 1) - i} + "
    print(out)


if __name__ == '__main__':
    PRIME = 3

    a = [1, 0, 3, 1]
    b = [1, 2, 1]

    print("POLYNOMIAL LONG DIVISION")
    N = [4, 2, 1]
    D = [0, 2, 2]
    print(PolyDiv(a, b, PRIME))
    # print(" %s / %s =" % (N, D))
    # print(" %s remainder %s" % poly_div(N, D, 5))
    # print(Polynomial_addition(list_1, list_2, PRIME))
    # PrintPolynomial(Polynomial_multiplication(list_1, list_2, PRIME))
