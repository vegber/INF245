"""
Implement arithmetic operations with polynomials g(x) & f(x) under
modulus a prime number p.
"""


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
            init[f + g] = (init[f + g] + (f_x[f] * g_x[g]) % p) % p
    return init


def PolynomialDivision(f_x: list, g_x: list, p: int):
    import copy
    # make sure deg(f_x) > deg(g_x)
    assert len(f_x) >= len(g_x)

    g_x = g_x[get_highest_order_elem(g_x):]
    f_x = f_x[get_highest_order_elem(f_x):]
    tmp_G = g_x
    g_x = [0] * (len(f_x) - len(g_x)) + g_x

    # Copy input arr. so we can change them
    f_x = copy.deepcopy(f_x)
    g_x = copy.deepcopy(g_x)

    quotient = []
    _round = 0
    while True:
        # find index of element with the highest order != 0
        highest_order_elemt_index = get_highest_order_elem(f_x)
        order_of_highest_elemt = len(f_x) - highest_order_elemt_index

        if order_of_highest_elemt < len(g_x) - 1:
            # can't divide this element, rest is remainder
            return quotient, f_x
        # find a s.t a*g[0] == inv(f_x(0)) = 0
        factor = [x for x in range(p) if ((x * g_x[get_highest_order_elem(g_x)]) % p) == f_x[
            highest_order_elemt_index]][0]  # modDivide(f_x[highest_order_elemt_index], get_highest_order_elem(g_x), p)

        quotient.append((factor, abs(get_highest_order_elem(f_x) - get_highest_order_elem(g_x))))

        factor_order = quotient[_round][1]
        y = 0
        for g_i in range(len(tmp_G)):
            g_i_z = -(tmp_G[g_i] * factor) % p
            order_g_i = len(tmp_G) - g_i
            f_x[len(f_x) - (order_g_i + factor_order)] = (f_x[len(f_x) - (order_g_i + factor_order)] + g_i_z) % p
            y += 1
        _round += 1


def Polynomial_GCD(f_x: list, g_x: list, p: int):
    import copy
    f = copy.deepcopy(f_x)
    g = copy.deepcopy(g_x)

    q, r = PolynomialDivision(f, g_x, p)
    g = [0] * (len(r) - len(g)) + g
    while True:
        if getDegre(r) == 0 and r[-1] == 0:
            break
        q, r_prime = PolynomialDivision(g, r, p)
        g = r
        r = r_prime
    # make gcd Monic
    inv_largest = pow(g[get_highest_order_elem(g)], -1, p)
    return [(g * inv_largest) % p for g in g]


def Polynomial_BinaryExponentiation(f_x, a, g_x, p):
    res = [1]

    if a == 0:
        return res

    while a > 0:

        f_x = f_x[get_highest_order_elem(f_x):]
        if a % 2 == 1:
            res = Polynomial_multiplication(res, f_x, p)
            q, res = PolynomialDivision(res, g_x, p)
        a = a // 2
        h = Polynomial_multiplication(h, h, p)
        q, h = PolynomialDivision(h, g_x, p)

    return res


def get_highest_order_elem(l: list) -> int:
    if len([x for x in l if x == 0]) == len(l):
        return len(l)
    else:
        return [i for i, e in enumerate(l) if e != 0][0]


def getDegre(l: list):
    x = 0
    for i, e in enumerate(l):
        if e != 0:
            x = (len(l)) - i - 1
            break
    return x


def PrintPolynomial(h_x: list):
    out = ""
    for i, e in enumerate(h_x):
        if i == len(h_x) - 1:
            out += str(e)
        else:
            out += f"{e}x^{(len(h_x) - 1) - i} + "
    return out


def PrintQuotient(q_x: list):
    out = ""
    for i, e in enumerate(q_x):
        out += f"{e[0]}x^{(i)}"
    return out


if __name__ == '__main__':
    PRIME = 3

    # a = [1, 0, 3, 1]
    # b = [1, 2, 1]

    f = [2, 1, 3, 0]
    g = [1, 1]
    a = 2
    test_1 = [2, 2, 2, 1]

    test_2 = [2, 0, 5, 1]

    # print(Polynomial_GCD(test_1, test_2, p=PRIME))
    # print("POLYNOMIAL LONG DIVISION")
    # a, b = PolynomialDivision(test_2, test_1, 5)
    # print(a)
    # print(f"Dividing {PrintPolynomial(test_2)} with {PrintPolynomial(test_1)}")
    # out = f"quotient: {PrintQuotient(a)} with remainder: {PrintPolynomial(b)}"
    # print(out)
    # print(f"quotient: {a} remainder", end=" "), PrintPolynomial(b)
    # print(PolyDiv(a, b, PRIME))
    # print(" %s / %s =" % (N, D))
    # print(" %s remainder %s" % poly_div(N, D, 5))
    # print("POLYNOMAIAL ADDITION")
    # PrintPolynomial(Polynomial_addition(test_1, test_2, PRIME))
    # print("MULTIPLICATION: ")
    # print(PrintPolynomial(Polynomial_multiplication(test_1, test_1, 3)))
    # print(pow(2, -2, PRIME))

    print("Exponentiation under modulo: ")
    print(Polynomial_BinaryExponentiation(f, a, g, PRIME))
