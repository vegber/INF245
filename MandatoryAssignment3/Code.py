"""
Digital Signature Algorithm and Discrete Logarithms
"""
import copy
import math
import random

import numpy as np

from MandatoryAssignment1.A1.Mandatory1 import BinaryExponentiationWithoutRecursion as binExp
from MandatoryAssignment3.MatrixHelper import RowReduceMatrixRowEchelonForm


def findJ(p):
    """
    j should be gcd(j, p-1) = 1
    :param p:
    :return:
    """
    j = 2
    while math.gcd(j, (p - 1)) != 1:
        j = random.randint(1, p - 2)
    return j


def forgeElGamalSigAlgorithm(p: int, g: int, y: int):
    """
    Forge ElGamal signature without knowing the private key
        Construct triple: m, a, b
            m: int
            a, b = Signature
    Forge signature for the message m, find a,b s.t
        0 < a < p,
        0 < b < p - 1
    and
        g^m ≡ y^a * a^b mod p

    :param p: prime p
    :param g: primitive root mod p
    :param y: y ≡ 5^x mod p
    :return:
    """
    # system private key: x mod p - 1
    # public key: y ≡ 5^x mod p

    i, j = 2, findJ(p)
    print(f"\nFound valid co - prime j mod p - 1: {j} ")
    g_i, y_j = binExp(g, i, p), binExp(y, j, p)
    a = (g_i * y_j) % p
    j_1 = pow(j, -1, p - 1)
    b = pow((-a * j_1), 1, p - 1)
    m = pow((-a * i * j_1), 1, p - 1)
    print(f"Forged signature with values: \n\t(m): {m}\n\t(a): {a}\n\t(b): {b}\n")
    return m, a, b


def verifyElGamalSig(p, g, y, m, a, b):
    print("Verify ElGamalSig")
    return binExp(g, m, p) == (binExp(y, a, p) * (binExp(a, b, p)) % p)


def runTaskOne():
    # forge ElGamal
    print(f"Forge Elgamal Signature")
    P = 593831971123607
    G = 13
    Y = 239505966643112
    print(f"With param: (p) {P},  \t (g) {G}, ", end=" ")
    print(f"(y) ≡ g^x mod p =  {Y}")
    m, a, b = forgeElGamalSigAlgorithm(P, G, Y)
    print(f"Is g^m mod p == y^a * a^b mod p ? {verifyElGamalSig(P, G, Y, m, a, b)}")


def recoverXFromSameKDSA(q, r_1, s_1, s_2, m_1, m_2):
    """
    When the same k is used to generate the signatures of two messages: m1 & m2, we can easily recover the
    secret x.
    The following should hold:
    Given the two signatures M1 & M2, we solve for k
        s1 = inv(k)*(M1 + x*r) mod q
        s2 = inv(k)*(M2 + x*r) mod q
        # subtract both signatures:
        (**All is under mod q)
            s1 - s2 = inv(k)(M1 + x*r) - inv(k)(M2 + x*r)
            s1 - s2 = inv(k)(M1 + x*r - M2 + x*r) # reduce exp.
            s1 - s2 = inv(k)(M1 - M2) # reduce exp.
        Thus,
            k = (M1 - M2) / (s1 - s2)
            recover x from the expression
                x = r^-1 * (s * k - M) mod q
    :return: Secret x
    """
    delta = pow((s_1 - s_2) % q, -1, q)
    m_ = (m_1 - m_2) % q
    k = (m_ * delta) % q
    rINV = pow(r_1, -1, q)
    x1 = (rINV * (s_1 * k - m_1)) % q
    x2 = (rINV * (s_2 * k - m_2)) % q
    assert x1 == x2, "Recover x failed"
    return x1


def verifyXDSA(g, x, p, y):
    return y == binExp(g, x, p) and y == pow(g, x, p)


def runTaskTwo():
    """"""
    p = 949772751547464211
    q = 4748626326421
    g = 314668439607541235

    # y = g^x mod p
    y = 254337213994578435

    m_1, m_2 = 2393923168611338985551149, 9330804276406639874387938

    sig_1, r_1, s_1 = 2393923168611338985551149, 2381790971040, 3757634198511
    sig_2, r_2, s_2 = 9330804276406639874387938, 2381790971040, 4492765251707

    # print(rhoMethodForLogarithms(y, g, q))

    x = recoverXFromSameKDSA(q, r_1, s_1, s_2, m_1, m_2)
    print(f"Found X: {x}")
    print(f"Now test if correct answer by: y ≡ g^x mod p")
    print(f"\ty is {y}")
    print(f"\tg^x mod p is: {binExp(g, x, p)}\ny≡g^x mod p == {verifyXDSA(g, x, p, y)}")

    print("_" * 20)
    print("task two B: Rho method")

    f()


def fab(x, a, b, alfa, beta, N, n):
    if x % 3 == 0:
        x = pow(x, 2, N)  # x * x % N
        a = a * 2 % n
        b = b * 2 % n

    elif x % 3 == 1:
        x = x * alfa % N
        a = (a + 1) % n
    else:
        x = x * beta % N
        b = (b + 1) % n

    return x, a, b


def solveForX(a, b, A, B, n):
    print(f"{(A - a)} * {pow(b - B, -1, n)} % {n}")
    return ((A - a) * pow((b - B), -1, n)) % n


def f(N=949772751547464211, n=4748626326421, alfa=314668439607541235, beta=254337213994578435, y=254337213994578435):
    x, a, b, = 1, 0, 0
    X, A, B = x, a, b
    print("%8s %15s %15s %15s  %18s %15s %15s\n" % ("i", "x", "a", "b", "X", "A", "B"))
    for i in range(1, n):
        x, a, b = fab(x, a, b, alfa, beta, N, n)
        X, A, B = fab(X, A, B, alfa, beta, N, n)
        X, A, B = fab(X, A, B, alfa, beta, N, n)
        if x == X and math.gcd(b - B, n) == 1:
            print("%8d %15d %15d %15d   %18d %15d %15d\n" % (i, x, a, b, X, A, B))
            x = solveForX(a, b, A, B, n)
            print(f"x is: {x}")
            break


def is_prime(n):
    """
    :param n:
    :return: True/False n is prime
    """
    for i in range(2, n):
        if (n % i) == 0:
            return False
    return True


def createBsmooth(B: int):
    """
    Find primes in B
    :param B:
    :return: list(primes_n)
    """
    return [p for p in range(2, B + 1) if is_prime(p)]


def prime_factors(n):
    i = 2
    factors = []
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1:
        factors.append(n)
    if factors:
        return sorted(factors)
    else:
        return [0]


def createMatrix(random_X_s):
    org_m = copy.deepcopy(random_X_s)
    return [x[1] for x in random_X_s], org_m


def factor(n, B):
    # Factor base for the given number
    base = createBsmooth(B)

    # Starting from the ceil of the root
    # of the given number N
    start = int(math.sqrt(n))

    # Storing the related squares
    pairs = []

    # For every number from the square root
    # Till N
    for i in range(start, n):

        # Finding the related squares
        for j in range(len(base)):
            lhs = i ** 2 % n
            rhs = base[j] ** 2 % n

            # If the two numbers are the
            # related squares, then append
            # them to the array
            if lhs == rhs:
                pairs.append([i, base[j]])

    new = []

    # For every pair in the array, compute the
    # GCD such that
    for i in range(len(pairs)):
        factor = math.gcd(pairs[i][0] - pairs[i][1], n)

        # If we find a factor other than 1, then
        # appending it to the final factor array
        if factor != 1:
            new.append(factor)

    x = np.array(new)

    # Returning the unique factors in the array
    return np.unique(x)


def is_b_smooth(b: list, S_b: list) -> bool:
    return max(b) <= max(S_b)


def RowReduceEchelonForm(m: list, modulus: int):
    # A is concat of matrix A and vector
    # Our matrix M = m x (t + 1)
    MATRIX = copy.deepcopy(m)
    r = 0
    while r != len(MATRIX):
        # find pivot el
        Mij = MATRIX[r][r]
        d = math.gcd(Mij, modulus)
        # if d > 1:  # if GCD(Mij, N) > 1, terminate
        #     return list
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
    return list(MATRIX)


def SolveRowEchelonForm(m: list, p: int):
    """
    Tries to reduce matrix from pivot elements. Start in reversed order,  take each row[pivot] multiply
    under modulo with this row. Subtract row with current row. Return solved matrix
    :param M: reduced echenom form matrix
    :param p: prime
    :return: solved matrix
    """
    PIV = 0
    M = copy.deepcopy(m)
    print(type(M))
    length = sum(1 for _ in M)
    for j in range(length - 1, 0, -1):
        for i in range(length - PIV - 1, 0, -1):
            # Take each pivot row
            num = M[i - 1][len(M[0]) - 1 - PIV - 1]
            for z in range(len(M[0]) - 1, 0, -1):
                # subtract pivot row with row - 1 to create null elements
                M[i - 1][z] = (M[i - 1][z] - M[j][z] * num) % p
                # reduce under modulo
        PIV += 1
    return M


def printMatrix(b):
    for lmao in b:
        print(list(lmao), end="\n")


def IndexCalculusAlgorithm(p=2602163, alfa=1535637, B=30, g=2):
    # Algorithm divided into two parts: finding linear congruences
    # and trying to row reduce row echelon form
    # output from "FindRowsBsmooth" returns rows random_X_s
    # that we insert into the already implemented RowReduceMatrixRowEchelonForm
    #  from Assignment 1
    # Since not all systems is possible to rowreduce, we wrap the return into
    # a try catch loop, this triggers a re - find new parameters if it fails.

    # if unique system is found, we solve system according to algorihtm
    # spesification
    S_b, m, n, random_X_s = FindRowsBsmooth(B, p)
    M = None
    while M is None:
        try:
            while True:
                # Row Reduce Matrix to Row Echelon Form
                # This was implemented in assignment 1
                M = RowReduceMatrixRowEchelonForm(random_X_s, p - 1)
                # Convert to numpy arr. for easy list comprehension later
                M_numpyarr = np.array(M, dtype=int)
                rowReduced = M_numpyarr.reshape(m, n + 1)  # Reshape 1D list to 2D

                # Check if system has unique solution -- else find new/more numbers b smooth
                if rowReduced[0][-1] == 0:
                    S_b, m, n, random_X_s = FindRowsBsmooth(B, p)
                else:
                    # print matrix before row reduced
                    print(f"Linear system to solve: ")
                    printMatrix(random_X_s)
                    print()
                    # Print matrix with unique solution
                    print(f"Linear system on row reduced by gaussian elm. ")
                    printMatrix(rowReduced)
                    while True:
                        # find random valid y
                        y = random.randint(0, p - 1)
                        b = (pow(g, y, p) * alfa) % p
                        fac_ = prime_factors(b)
                        # if b smooth, break => solve
                        if is_b_smooth(fac_, S_b):
                            print(
                                f"Found valid y: ({y}),\n\twhich is bsmooth with factors: \n\t\t({' + '.join([f'{y}^{x}' for x, y in zip(fac_, S_b)])}) mod {p}", end='\n')
                            print()
                            break
                        # Found valid y
                    # solve:
                    l = [0] * n
                    for i_s in range(len(S_b)):
                        l[i_s] = fac_.count(S_b[i_s])
                    x_vals = findXi(rowReduced)
                    x = 0
                    print(f"Now, we want to sum and add each solution from the solved matrix, \n\tand multiply with "
                          f"the corr. index from S_b :: minus the y value closed under mod p-1"
                          f"\n\t\t({' + '.join([f'{x} * {y}' for x, y in zip(x_vals, l)])}) - {y}   mod {p-1}")
                    for x_i, l_i in zip(x_vals, l):
                        x += x_i * l_i  # (x_i**l_i)
                    x = (x - y) % (p - 1)
                    print(f"After calc. x is found to be: {x}")
                    return x
                # exceptions: row reduce failed, try with new parameters
        except:
            S_b, m, n, random_X_s = FindRowsBsmooth(B, p)


def FindRowsBsmooth(B, p):
    S_b = createBsmooth(B)
    n = len(S_b)
    random_X_s = []
    # find m + c rows
    m = n
    while len(random_X_s) != m:
        # find random b
        x_s = random.randint(math.floor(math.sqrt(p)), p - 1)
        # compute congruence
        b = pow(2, x_s, p)  # (x_s ** 2) % p
        if b == 1:
            continue  # if prime and >= B ??
        # find factors
        factors_of_b = prime_factors(b)
        if max(factors_of_b) > B: continue
        # check if valid
        if is_b_smooth(factors_of_b, S_b):
            l = [0] * n
            for i_s in range(len(S_b)):
                l[i_s] = factors_of_b.count(S_b[i_s])
            temp = l
            temp.append(x_s)
            random_X_s.append(temp)
    return S_b, m, n, random_X_s


def findXi(m: list):
    xi, piv = [], 0
    for _ in m:
        xi.append(m[piv][-1])
        piv += 1
    return xi


def runTaskThree():
    # should be 2116767
    while True:

        x = IndexCalculusAlgorithm()
        if x is not None:
            x = int(x)
            if pow(2, x, 2602163) == 1535637:
                print(f"Now, check if correct solution: 2^x ≡ 1535637 ?=> 2^{x} mod 2602163 is: {pow(2, x, 2602163)}")
                print(f"We found correct X: {x}")
                break
            else:
                print(f"Wrong: got x ==  {x}", end="\n")
                print()


if __name__ == '__main__':
    # runTaskOne()
    # runTaskTwo()
    runTaskThree()
