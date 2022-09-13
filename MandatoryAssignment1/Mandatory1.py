import random
from decimal import Decimal, localcontext

import sympy


def ExtendedEuclideanAlgorithm(a, b):
    if a == 0:
        return b, 0, 1

    gcd, c_1, c_2 = ExtendedEuclideanAlgorithm(b % a, a)
    c1_prime = c_2 - (b // a) * c_1
    c_2_prime = c_1
    return gcd, c1_prime, c_2_prime


def GCD(a, b):
    if b == 0:
        return abs(a)
    else:
        return GCD(b, a % b)


def BinaryExponentiation(a: int, b: int):
    """
    This is my own implementation of Binary Exponentiation.
    However, Python is very bad at recursion, and the native
    implementation of pow(a, b) is more effective due to
    RAM/Memory concern.
    :param a:
    :param b:
    :return:
    """
    if b == 0:
        return 1
    a_prime = a * BinaryExponentiation(a, b - 1)

    return int(a_prime)


def VerifyECA(a, b, u, v):
    return GCD(a, b) == u * a + v * b


def ReduceToEchelonForm(matrix: list, modulus: int):
    r, t = 0, 0
    for j in range(len(matrix)):
        t = t + 1
        for i in range(len(matrix[j])):
            r = r + 1
            b = matrix[j][i]
            if b % modulus != 0:
                d = GCD(b, modulus)
                if d > 1:
                    # n = d * (n/d) => terminate
                    return
                # if d == 1, z * b congruent with  1 mod n
            # swap Mr+1 <=> Mi , Mr+1 <= r* Mr+1
            matrix[r], matrix[i] = matrix[i], matrix[r]
            matrix[r] = [(x * matrix[r]) % modulus for x in matrix[r]]
            print(matrix)


# modulo function to perform binary
# exponentiation
def modulo(base, exponent, mod):
    x = 1
    y = base
    while exponent > 0:
        if exponent % 2 == 1:
            x = (x * y) % mod
        y = (y * y) % mod
        exponent = exponent // 2
    return x % mod


def JacobiSymbol(alpha: int, n: int) -> int:
    # alpha should be int, and n should be odd positive int
    # assert (n > 0 and n % 2 == 1)
    alpha = pow(alpha, 1, n)
    t = 1
    while alpha != 0:
        while alpha % 2 == 0:
            alpha = alpha // 2
            r = n % 8
            if r == 3 or r == 5:
                t = -t
        alpha, n = n, alpha
        if alpha % 4 == 3 and n % 4 == 3:
            t = -t
        alpha = alpha % n

        if alpha > n // 2:
            alpha = alpha - n
    if n == 1:
        return t
    else:
        return 0


def Solovay_Strassen_Test(n, k=1) -> str:
    """
    :input: n, a value to test primality
    :out: composite if test fails, probably prime else
    :rtype: str
    """
    # check n odd prime > 1
    import math
    assert (n > 1)

    for i in range(k):
        a = random.randint(2, int(math.sqrt(n - 1)))  # random is ge && le
        x = JacobiSymbol(a, n)
        print(modulo(a, (n - 1) / 2, n))
        mod = a ** ((n - 1) / 2)
        if (x == 0) or mod != modulo(a, (n - 1) / 2, n):
            # if (x == 0) or (pow(a, (n - 1) // 2, n) != x % n):
            return "Composite"
    return "Probably Prime"


if __name__ == '__main__':
    matrix = [
        [1, -2, -2, -2, -1, 2],
        [0, 3, -2, -3, 1, 3],
        [3, 0, 0, 1, -1, -3],
        [3, -3, -2, 0, 1, 1],
        [0, -3, 3, -3 - 3, 2]]
    MODULO = 456995412589

    test_matrix = [
        [1, 2, -2, 0, 1],
        [-1, 0, 0, 1, 1],
        [0, 2, -1, 0, 1]]
    TEST_MODULO = 15
    ########
    # Task one - Extended euclidean algorithm
    # a = 620709603821307061
    # b = 390156375246520685
    # R, u, v = ExtendedEuclideanAlgorithm(a, b)
    # print(f"Extended Euc. Algorithm of\nA  == {a}\nB  == {b}\GCD is {R}")
    # print(f"Coefficient for bigger int (u) u == {u}\nCoefficient for smaller int (v) == {v}")
    # print(BinaryExponentiation(2, 4))
    # print(5 // 2)
    # ReduceToEchelonForm(matrix=test_matrix, modulus=TEST_MODULO)
    # print(ExtendedEuclideanAlgorithm(2, 15))

    # var_1 = -776439811
    # var_2 = 50556018318800449023
    # print(JacobiSymbol(var_1, var_2))
    print(sympy.isprime(2 ** 127 - 1))
    print(Solovay_Strassen_Test(51))
