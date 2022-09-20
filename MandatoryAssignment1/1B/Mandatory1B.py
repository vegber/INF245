"""
Implement arithmetic operations with polynomials g(x) & f(x) under
modulus a prime number p.
"""


def Polynomial_addition(f_x: list, g_x: list, p: int) -> list:
    f_x, g_x = EqualizeLengthOnOrder(f_x, g_x)

    return [(x + y) % p for x, y in zip(f_x, g_x)]


def EqualizeLengthOnOrder(f_x, g_x):
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


def test(num, den):
    print("%s / %s ->" % (num, den))
    q, r = Polynomial_division(num, den)
    print("quot: %s, rem: %s\n" % (q, r))
    return q, r


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

    list_1 = [1, 2, 0]
    list_2 = [2, 1, 1, 2]

    num = [1, 5, 10, 10, 5, 1]
    den = [1, 2, 1]
    a, b = test(num, den)
    PrintPolynomial(a)
    PrintPolynomial(b)

    num = [5, 16, 10, 22, 7, 11, 1, 3]
    den = [1, 2, 1, 3]

    quot = [5, 1, 3, 0, 1]
    rem = [0, 5]

    q, r = test(num, den)
    PrintPolynomial(q)
    PrintPolynomial(r)
    assert quot == q
    assert rem == r

    # print(Polynomial_addition(list_1, list_2, PRIME))
    # PrintPolynomial(Polynomial_multiplication(list_1, list_2, PRIME))
