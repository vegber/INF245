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


def Polynomial_division_remainder(f_x: list, g_x: list) -> list:
    return []


def DivsionUnderModulo(a: int, b: int, p: int) -> int:
    return (a % p) * (pow(b, -1, p)) % p


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

    # print(Polynomial_addition(list_1, list_2, PRIME))
    # PrintPolynomial(Polynomial_multiplication(list_1, list_2, PRIME))
    print(DivsionUnderModulo(8, 3, 5))
