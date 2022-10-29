"""
Digital Signature Algorithm and Discrete Logarithms
"""
import math
import random
from MandatoryAssignment1.A1.Mandatory1 import BinaryExponentiationWithoutRecursion as binExp
from MandatoryAssignment1.A1.Mandatory1 import ExtendedEuclideanAlgorithm as eea
from sympy.ntheory.residue_ntheory import _discrete_log_pollard_rho


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
    print(f"a {(A - a)} * {pow(b - B, -1, n)} % {n}")
    return ((A - a) * pow((b - B), -1, n)) % n


def checkAnswer(generator, x, p, y):
    return y == pow(generator, x, p)


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


if __name__ == '__main__':
    # runTaskOne()
    runTaskTwo()
