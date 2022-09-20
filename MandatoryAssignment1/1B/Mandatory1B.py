"""
Implement arithmetic operations with polynomials g(x) & f(x) under
modulus a prime number p.
"""

import re
from turtle import end_fill


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


def Polynomial_division(num, den):
    # Create normalized copies of the args
    num = num[:]
    normalize(num)
    den = den[:]
    normalize(den)

    if len(num) >= len(den):
        # Shift den towards right so it's the same degree as num
        shiftlen = len(num) - len(den)
        den = [0] * shiftlen + den
    else:
        return [0], num

    quot = []
    divisor = float(den[-1])
    for i in range(shiftlen + 1):
        # Get the next coefficient of the quotient.
        mult = num[-1] / divisor
        quot = [mult] + quot

        # Subtract mult * den from num, but don't bother if mult == 0
        # Note that when i==0, mult!=0; so quot is automatically normalized.
        if mult != 0:
            d = [mult * u for u in den]
            num = [u - v for u, v in zip(num, d)]

        num.pop()
        den.pop(0)

    normalize(num)
    return quot, num


def PolyDiv(f_x, g_x, p: int):
    import copy
    # make sure deg(f_x) > deg(g_x)
    assert len(f_x) >= len(g_x)

    # Copy input arr. so we can change them
    f_x = copy.deepcopy(f_x)
    g_x = copy.deepcopy(g_x)

    i_th_elem = 0
    while True:
        dividator = f_x[i_th_elem]

        divident = g_x[0]

        # mult. divident with x s.t dividator == 0 mod p
        print(divident)
        faktor = (divident % p) * pow(divident, -1, p)
        print(faktor)
        print(dividator * faktor % p)
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

    PolyDiv(a, b, PRIME)

    # print(Polynomial_addition(list_1, list_2, PRIME))
    # PrintPolynomial(Polynomial_multiplication(list_1, list_2, PRIME))
