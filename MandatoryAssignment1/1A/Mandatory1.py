import copy
import math
import random
from typing import Type


def ExtendedEuclideanAlgorithm(a, b):
    if a == 0:
        return b, 0, 1
    gcd, u_1, v_1 = ExtendedEuclideanAlgorithm(b % a, a)
    return gcd, v_1 - (b // a) * u_1, u_1


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
    """
    if b == 0:
        return 1
    a_prime = a * BinaryExponentiation(a, b - 1)
    return int(a_prime)


def BinaryExponentiationWithoutRecursion(a: int, b: int, mod: int):
    # Modulo function to perform binary exponentiation - without recursion
    temp, base_number = 1, a
    while b > 0:  # exponent larger than zero
        if b % 2 == 1:  # if not even
            temp = (temp * base_number) % mod  # pow(x*y, 1, mod)
        base_number, b = (base_number * base_number) % mod, b // 2
    return temp % mod


def VerifyECA(a, b, u, v) -> bool:
    """
    Verification method to check that:
        For a ≥ b > 0, Extended euclidean alg.
        should find integers u & v
        such that ua + vb = gcd(a, b)
    """
    return GCD(a, b) == u * a + v * b


def RowReduceEchelonForm(m: list, modulus: int) -> Type[list] | list:
    # A is concat of matrix A and vector
    # Our matrix M = m x (t + 1)
    MATRIX = copy.deepcopy(m)
    r = 0
    while r != len(MATRIX):
        # find pivot el
        Mij = MATRIX[r][r]
        d = GCD(Mij, modulus)
        if d > 1:  # if GCD(Mij, N) > 1, terminate
            return list
        if d == 1:  # if GCD is one, find Z * b congruent with 1 mod N
            z = pow(Mij, -1, modulus)
            # apply z to all elements of row
            MATRIX[r] = [MATRIX[r][x] * z % modulus for x in range(len(m[0]))]
        # pivot element completely divides modulo
        if Mij % modulus == 0:
            # switch rows r + 1
            for x in range(r + 1, len(MATRIX)):
                if MATRIX[x][Mij] != 0 % modulus:
                    MATRIX[r], MATRIX[r + 1] = MATRIX[r + 1], MATRIX[r]
                # already switched, now reduce

        # make zero space for r + 1 under pivot element
        for e in range(r + 1, len(MATRIX)):
            # apply zero element mult. to all elements of row r + 1
            MATRIX[e] = [(MATRIX[e][x] - MATRIX[e][r] % modulus * MATRIX[r][x]) % modulus for x in range(len(m[0]))]
        # if r = m, terminate, else r <- r+1 and
        Mij, r = Mij + 1, r + 1
    return MATRIX


# This is a helping function that will reduce the matrix upwards from the bottom left pivot
# which in the end will solve the matrix. We simply take each pivot, then take each position above
# the pivot, and multiply this number with the pivot row. then we subtract the pivot row with the
# current row above. this will create zeros above the pivot. when we are done with all the pivots
# we will be left with a solved matrix if the matrix is solvable.
def SolveRowEchelonForm(m, p):
    kpivot = 0
    for x in range(len(m) - 1, 0, -1):
        for y in range(len(m) - kpivot - 1, 0, -1):
            num = m[y - 1][len(m[0]) - 1 - kpivot - 1]
            for z in range(len(m[0]) - 1, 0, -1):
                m[y - 1][z] = (m[y - 1][z] - m[x][z] * num) % p
        kpivot += 1

    return m


def JacobiSymbol(a, n, d=1):
    if a == 0:
        return 0  # (0/n) = 0
    if a == 1:
        return d  # (1/n) = 1

    if a < 0:
        # property of Jacobi
        # (a/n) = (-a/n)*(-1/n)
        a = -a
        if n % 4 == 3:
            # (-1/n) = -1 if n = 3 (mod 4)
            d = -d

    while a:
        if a < 0:
            # (a/n) = (-a/n)*(-1/n)
            a = -a
            if n % 4 == 3:
                # (-1/n) = -1 if n = 3 (mod 4)
                d = -d
        while a % 2 == 0:
            a = a // 2
            if n % 8 == 3 or n % 8 == 5:
                d = -d
        # swap
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            d = -d
        a = a % n
        if a > n // 2:
            a = a - n
    if n == 1:
        return d
    return 0


def Solovay_Strassen_Test(n, k=20) -> str:
    """
    :input: n, a value to test primality
    :out: composite if test fails, probably prime else
    :rtype: str
    """
    # check n odd prime > 1
    assert (n > 1)
    for i in range(k):
        # Solovay strassen test:
        # 1) gcd(a, n) != 1 => composite
        # 2) a ** (n-1) / 2 congruent with jacobi(a/n) mod n
        a = random.randint(2, n)  # random is ge && le
        if GCD(a, n) != 1:
            return "Composite"
        x = (n + JacobiSymbol(a, n)) % n
        # a ** ((n - 1) / 2) can be re - written to our bin. exp method as:
        mod = BinaryExponentiationWithoutRecursion(a, (n - 1) // 2, n)
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
    a = 620709603821307061
    b = 390156375246520685
    R, u, v = ExtendedEuclideanAlgorithm(a, b)
    print(f"Extended Euc. Algorithm of\nA  == {a}\nB  == {b}\GCD is {R}")
    print(f"Coefficient for bigger int (u) u == {u}\nCoefficient for smaller int (v) == {v}")
    print(BinaryExponentiation(2, 4))

    var_1 = -776439811
    var_2 = 50556018318800449023
    print(f"Jacobi of {var_1} and {var_2} is: {JacobiSymbol(var_1, var_2)}")
    print(f"Solovay Strassen test: 2**127 -1 is {Solovay_Strassen_Test((2 ** 127) - 1, 20)}")
    a = 620709603821307061
    b = 390156375246520685
    #
    k, u, v = (ExtendedEuclideanAlgorithm(620709603821307061, 390156375246520685))
    print(f"EEA: {k} u is {u} v is {v}")
    print(f"Does u * {a} + v * {b} = GCD(a, b)?: {VerifyECA(a, b, u, v)}")
    M = RowReduceEchelonForm(matrix, MODULO)
    M = SolveRowEchelonForm(M, MODULO)

    s = [[str(e) for e in row] for row in M]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print('\n'.join(table))
    print("\n")


    print(
        f"Binary exponentiation modulo n: {BinaryExponentiationWithoutRecursion(393492946341, 103587276991, 72447943125)}")
