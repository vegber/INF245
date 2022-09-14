import random
from test import calculateJacobian as JacobiSymbol


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


def ReduceToEchelonForm(matrix: list, modulus: int):
    r, t = 0, 0
    for j in range(len(matrix) - 1):  #
        t = t + 1
        for i in range(len(matrix[j]) - 1):
            r = r + 1
            b = matrix[j][i]
            z = 0
            if b % modulus != 0:
                d = GCD(b, modulus)
                if d > 1:
                    # n = d * (n/d) => terminate
                    return
                if d == 1:
                    z = pow(b, -1, modulus)
                # swap Mr+1 <=> Mi , Mr+1 <= r* Mr+1
                matrix[j][r] = z * matrix[j][r]
                matrix[j][r] = matrix[i]

                # matrix[r] = [(x * matrix[r]) for x in matrix[r]]

            # gå gjennom rad
            # Mu =
            continue
        continue


def volume2(matrix: list, modulus: int):
    # A is concat of matrix A and vector a
    # Our matrix M = m x (t + 1)

    M = matrix
    col_length, rad_length, pivot_counter = 0, 0, 0

    while rad_length < (len(M[0]) - 1) and col_length < len(M):
        # find pivot el
        Pivot_elem = M[pivot_counter][rad_length]
        z = 0
        if Pivot_elem % modulus != 0:
            # Reduce
            d = GCD(Pivot_elem, modulus)
            if d > 1: return "Does not work!"
            if d == 1: z = pow(Pivot_elem, -1, modulus)
            M[pivot_counter] = [M[pivot_counter][x] * z % modulus for x in range(len(matrix[0]))]

            # make zero space under pivot element

            for a in range(pivot_counter + 1, len(M)):
                row = M[a]
                special_element = M[a][pivot_counter] % modulus
                M[a] = [(row[x] - special_element * M[pivot_counter][x]) % modulus for x in range(len(matrix[0]))]
        else:
            # switch
            for a in range(pivot_counter + 1, len(M)):
                row = M[a]
                if row[Pivot_elem] != 0 % 5:
                    M[pivot_counter], M[pivot_counter + 1] = M[pivot_counter + 1], M[pivot_counter]

                # already switched, now reduce
                if d == 1: z = pow(Pivot_elem, -1, modulus)
                M[pivot_counter] = [M[pivot_counter][x] * z % modulus for x in range(len(matrix[0]))]

                # make zero space under pivot element

                for a in range(pivot_counter + 1, len(M)):
                    row = M[a]
                    special_element = M[a][pivot_counter] % modulus
                    M[a] = [(row[x] - special_element * M[pivot_counter][x]) % modulus for x in range(len(matrix[0]))]

        Pivot_elem += 1
        rad_length += 1
        col_length += 1
        pivot_counter += 1
    return M


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

    M = volume2(matrix, MODULO)

    s = [[str(e) for e in row] for row in M]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print('\n'.join(table))
