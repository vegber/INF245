import copy
import math
import random
from typing import Any
import time
from fractions import Fraction
import MandatoryAssignment1.A1.Mandatory1 as md1
import gmpy2 as gmpy2
import numpy as np


def get_Quotients(a: int, b: int) -> list:
    """
    Perform Extended Euclidean Algorithm, to find quotients

    :param a: int a
    :param b: int b
    :return: quotients
    """
    s, old_s = 0, 1
    r, old_r = b, a

    q_n = []
    while r != 0:
        q_i = old_r // r
        q_n.append(q_i)
        old_r, r = r, old_r - q_i * r
        old_s, s = s, old_s - q_i * s

    return q_n


def find_roots(a, b, c):
    """
    Solving linear equation
    :param phi:
    :param N:
    :return:
    """
    # print(f"a: {a}, b: {b}, c: {c}")
    sqrt_val = gmpy2.isqrt((b ** 2) - (4 * a * c))  # math.sqrt(abs(dis))

    pluss = -b + sqrt_val

    negative = -b - sqrt_val

    x_1 = int(pluss // 2)  # Fraction(-b + sqrt_val, 2)
    x_2 = int(negative // 2)  # Fraction(-b - sqrt_val, 2)
    return x_1, x_2


def valid_roots(phi, N):
    b = ((-(N + 1 - phi)) ** 2)
    d = (4 * 1 * N)
    return gmpy2.is_square(b - d)


def ContinuedFractionAlgorithm(N: int, e: int):
    """
    Let N, e be an RSA pub key, where ed congruent with 1 mod phi(N)
    d is RSA secret exponent.
    We know d is rel. small.

    Task: Factor N = pq
    :return:
    """

    q_n = get_Quotients(e, N)
    P_n = [1] + [q_n[1]] + [0] * (len(q_n) - 2)
    Q_n = [0, 1] + [0] * (len(q_n) - 2)
    # P_n = q_n * P_(n-1) + P_(n-2)
    # Q_n = q_n * Q_(n-1) + P_(n_2)

    for n in range(2, len(q_n)):
        P_n[n] = q_n[n] * P_n[n - 1] + P_n[n - 2]
        Q_n[n] = q_n[n] * Q_n[n - 1] + Q_n[n - 2]

    # phi equation
    #   1. not integer, try another
    #   2. integer
    #       x^2 - (N+1 - phi)x + N
    P_n, Q_n = P_n[1:], Q_n[1:]
    P_n, Q_n = Q_n, P_n  # now: P = Numerator, Q = Denominator

    for P_i, Q_i in zip(P_n, Q_n):
        if P_i == 0 or Q_i == 0 or Q_i % 2 == 0 or e * Q_i % P_i != 1:
            continue
        phi = (e * Q_i - 1) // P_i
        if valid_roots(phi, N):
            x_1, x_2 = find_roots(1, -(N - int(phi) + 1), N)
            print(f"factor: {x_1} , {x_2}")
            # return abs(x_1), abs(x_1)
    return 0, 0


def g_X(x, P: int):
    return (x ** 2 + 1) % P


def FactorByP_Method(N: int) -> str | int | Any:
    x = 2
    y = 2
    d = 1
    while d == 1:
        x = g_X(x, N)
        y = g_X(g_X(y, N), N)
        d = math.gcd(abs(x - y), N)  # md1.GCD(abs(x - y), N)

    if d == N:
        return "failed"
    else:
        return d


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


def is_b_smooth(b: list, S_b: list) -> bool:
    return max(b) <= max(S_b)


def RowReduceV2(matrix, MODULO: int):
    for row in range(len(matrix)):
        for col in range(len(matrix[row])):
            matrix[row][col] = matrix[row][col] % MODULO
    return matrix
    A = np.array(matrix, dtype=int)

    i = 0  # row
    j = 0  # column
    while True:
        # find next nonzero column
        while all(A.T[j] == 0):
            j += 1
            # if reached the end, break
            if j == len(A[0]) - 1: break
        # if a_ij == 0 find first row i_>=i with a
        # nonzero entry in column j and swap rows i and i_
        if A[i][j] == 0:
            i_ = i
            while A[i_][j] == 0:
                i_ += 1
                # if reached the end, break
                if i_ == len(A) - 1: break
            A[[i, i_]] = A[[i_, i]]
        # divide ith row a_ij to make it a_ij == 1
        A[i] = A[i] / A[i][j]
        # eliminate all other entries in the jth column by subtracting
        # multiples of the ith row from the others
        for i_ in range(len(A)):
            if i_ != i:
                A[i_] = A[i_] - A[i] * A[i_][j] / A[i][j]
        # if reached the end, break
        if (i == len(A) - 1) or (j == len(A[0]) - 1): break
        # otherwise, we continue
        i += 1
        j += 1

    return A


def DixonsMethod(N: int, B: int):
    """
    Dixon's algorithm, Random squares:

    :param N: N is the number to factor
    :param B: B is the smoothness / complexity
    :return: Factors of N / exl. trivial solutions
    """

    ### STEP 1 ###
    # Random 1 ≤ x < N
    # compute b ≡ x^2 mod N, 1 ≤ b < N
    # trial division by primes |S_b|
    # |S_b| = {2, 3 .... ≤ B}

    # if b i s B - smooth
    # Yes: store factor of trial div. given by b = PI from j=1 to n, q_j ^(lj), q_j is in set S_b
    # else repeat new x

    S_b = createBsmooth(B)
    n = len(S_b)
    random_X_s = []
    while len(random_X_s) != B:
        # find random b
        x_s = random.randint(math.floor(math.sqrt(N)), N)
        # compute congruence
        b = (x_s ** 2) % N
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
            random_X_s.append((b, l))

    ### STEP 2 ###
    # collect m = n * c, rows. ( create matrix)
    M = []

    for x in random_X_s:
        print(x[1])
        M.append(x[1])
    org_m = copy.deepcopy(M)
    ### STEP 3 ###
    # Contruct system of lin. equations, i.e matrix
    print()
    M_r = RowReduceV2(M, 2)  # row reduced mod m
    for x in M_r:
        print(x)

    ### STEP 4 ###
    # Solve, system of lin. equations: i.e matrix - mod 2
    match_1 = 0
    match_2 = 0
    M_r = np.array(M_r, dtype=bool)
    print(M_r)
    for x_ in range(len(M_r)):
        for y in range(len(M_r)):
            if x_ == y: continue
            elif (M_r[x_] & M_r[y]).any():
                match_1 = x_
                match_2 = y
                break

    if match_1 == 0 and match_2 == 0:
        print("repeat")

    x_m_1 = [(z**x) for x, z in zip(org_m[match_1], S_b) if x != 0]
    x_m_2 = [(z**x) for x, z in zip(org_m[match_2], S_b) if x != 0]
    Y = np.prod(x_m_1) * np.prod(x_m_2)
    X = random_X_s[match_1][0] * random_X_s[match_2][0]

    ### STEP 5 ###
    # set X = x_1^(Z_1), ... x_m^(Z_m) mod N
    # set Y = q_1^(k_1), ... q_n^(K_m) mod N

    ### STEP 6 ###
    N_1 = math.gcd(X - Y, N)
    if N_1 == 1 or N_1 == N:
        DixonsMethod(N, B)
    else:
        print(f"Factor of {N} is {N_1}")
    # GCD(X - Y, N) == N_1
    # if N_1 == 1 or N
    # Try another solution of step 3
    #
    # N_1 != N: N_1 is factor of N!

    pass


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


if __name__ == '__main__':
    # N = 1098667602555738997359355609545017255443451483195436981354791765341639135658156206242197992115989996829728203054347117299
    # e = 415884005149779743103130951097941968793336745983309971301432658775763996247677181243042840232106535367251782466233724389

    # e = 42667  # 388511
    # N = 64741  # 331879
    # a, b = ContinuedFractionAlgorithm(N, e)
    # print(a, b)
    # p_method_N = 10862216162096506735513546937

    # var = ContinuedFractionAlgorithm(int(N), int(e))
    # print(continued_fraction(int(N) / int(e)))
    # start = time.time()
    # p_method_factor = FactorByP_Method(p_method_N)
    # print(f"Time took: {time.time() - start} and result was: {p_method_factor}")

    # N = 12
    DixonsMethod(899, 7)
    # print(((2**4 )* 3 * 7) % 629)
