import random
from typing import Type

from test import calculateJacobian as JacobiSymbol
import copy


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


def BinaryExponentiationWithoutRecursion(a: int, b: int, mod: int):
    """
    Modulo function to perform binary exponentiation - without recursion
    :param a:
    :param b:
    :param mod:
    :return:
    """
    temp, base_number = 1, a
    while b > 0:  # exponent larger than zero
        if b % 2 == 1:  # if not even
            temp = (temp * base_number) % mod  # pow(x*y, 1, mod)
        base_number, b = (base_number * base_number) % mod, b // 2
    return temp % mod


def VerifyECA(a, b, u, v):
    """
    Verification method to check that:
        For a ≥ b > 0, Extended euclidean alg.
        should find integers u & v
        such that ua + vb = gcd(a, b)
    :param a:
    :param b:
    :param u:
    :param v:
    :return:
    """
    return GCD(a, b) == u * a + v * b


def RowReduceEchelonForm(m: list, modulus: int) -> Type[list] | list:
    # A is concat of matrix A and vector
    # Our matrix M = m x (t + 1)
    MATRIX = copy.deepcopy(m)
    c, r = 0, 0
    while r != len(MATRIX):
        # find pivot el
        Mij = MATRIX[r][c]
        z = 0
        d = GCD(Mij, modulus)
        if d > 1:  # if GCD(Mij, N) > 1, terminate
            return list
        if d == 1:  # if GCD is one, find Z * b congruent with 1 mod N
            z = pow(Mij, -1, modulus)
            # apply z to all elements of row
            MATRIX[r] = [MATRIX[r][x] * z % modulus for x in range(len(m[0]))]

        if Mij % modulus == 0:
            # switch
            for x in range(r + 1, len(MATRIX)):
                if MATRIX[x][Mij] != 0 % modulus:
                    MATRIX[r], MATRIX[r + 1] = MATRIX[r + 1], MATRIX[r]
                # already switched, now reduce

        # make zero space under pivot element
        for e in range(r + 1, len(MATRIX)):
            # apply zero element mult. to all elements of row r + 1
            MATRIX[e] = [(MATRIX[e][x] - MATRIX[e][r] % modulus * MATRIX[r][x]) % modulus for x in range(len(m[0]))]
        Mij, c, r = Mij + 1, c + 1, r + 1
    return MATRIX


"""def JacobiSymbol(alpha: int, n: int) -> int:
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

        # if alpha > n // 2:
        #     alpha = alpha - n
    if n == 1:
        return t
    else:
        return 0"""


def Solovay_Strassen_Test(n, k=1) -> str:
    """
    :input: n, a value to test primality
    :out: composite if test fails, probably prime else
    :rtype: str
    """
    # check n odd prime > 1
    assert (n > 1)
    import math
    for i in range(k):
        # Solovay strassen test:
        # 1) gcd(a, n) != 1 => composite
        # 2) a ** (n-1) / 2 congruent with jacobi(a/n) mod n
        a = random.randint(2, int(math.sqrt(n)))  # random is ge && le
        if GCD(a, n) != 1:
            print(f"gcd not != 1")
            return "Composite"
        x = (JacobiSymbol(a, n))
        # a ** ((n - 1) / 2) can be re - written to our bin. exp method as:
        mod = BinaryExponentiationWithoutRecursion(a, (n - 1) / 2, n)
        if (x == 0) or mod != x:
            return "Composite"
    return "Probably Prime"


if __name__ == '__main__':
    matrix = [
        [1, -2, -2, -2, -1, 2],
        [0, 3, -2, -3, 1, 3],
        [3, 0, 0, 1, -1, 2],
        [3, -3, -2, 0, 1, 1],
        [0, -3, 3, -3, -3, 2]]
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
    # print(sympy.isprime(2 ** 127 - 1))
    # print(Solovay_Strassen_Test(170141183460469231731687303715884105727, 10))

    M = RowReduceEchelonForm(matrix, MODULO)

    s = [[str(e) for e in row] for row in M]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print('\n'.join(table))
