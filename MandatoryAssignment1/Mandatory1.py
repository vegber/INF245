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
    if b == 0:
        return 1

    a_prime = a * BinaryExponentiation(a, b - 1)

    return a_prime


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
    print(ExtendedEuclideanAlgorithm(2, 15))
