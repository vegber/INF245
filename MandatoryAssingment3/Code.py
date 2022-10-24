"""
Digital Signature Algorithm and Discrete Logarithms
"""
import math
import random
from MandatoryAssignment1.A1.Mandatory1 import BinaryExponentiationWithoutRecursion as binExp
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
    P = 274742860459763
    G = 5
    Y = 262274678376340
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
    return y == binExp(g, x, p)


def runTaskTwo():
    """"""
    p = 3986625417249813809
    q = 19928344283621
    g = 2890026512265626020

    # y = g^x mod p
    y = 1621561995432343084

    m_1, m_2 = 1115209791959069177061830, 2151657259407048953791701

    sig_1, r_1, s_1 = 1115209791959069177061830, 12880312906177, 14706957637905
    sig_2, r_2, s_2 = 2151657259407048953791701, 12880312906177, 16242187205965

    x = recoverXFromSameKDSA(q, r_1, s_1, s_2, m_1, m_2)

    print(f"Found X: {x}")
    print(f"Now test if correct answer by: y ≡ g^x mod p")
    print(f"\ty is {y}")
    print(f"\tg^x mod p is: {binExp(g, x, p)}\ny≡g^x mod p == {verifyXDSA(g, x, p, y)}")


def f(x, a, b, G, H, P, Q):
    """
    Pollard function, three subsets equal size
    :param x:
    :param a:
    :param b:
    :param G:
    :param H:
    :param P:
    :param Q:
    :return:
    """
    subset = x % 3
    if subset == 0:
        x = x * G % P
        a = (a + 1) % Q
    elif subset == 1:
        x = x * H % P
        b = (b + 1) % Q
    elif subset == 2:
        x = x ** 2 % P
        a = a * 2 % Q
        b = b * 2 % Q
    return x, a, b


def rhoMethodForLogarithms(G, H, P):
    """

    :param a: a generator of G
    :param b: an element of G
    :return: An integer x such that a^x = b, or failure
    """
    # a_i, b_i, x_i = 0, 0, 1

    # P: prime
    # H:
    # G: generator
    Q = (P - 1) / 2  # sub group

    x = G * H
    a = 1
    b = 1

    X = x
    A = a
    B = b

    # Do not use range() here. It makes the algorithm amazingly slow.
    for i in range(1, P):
        # Who needs pass-by reference when you have Python!!! ;)

        # Hedgehog
        x, a, b = f(x, a, b, G, H, P, Q)

        # Hare
        X, A, B = f(X, A, B, G, H, P, Q)
        X, A, B = f(X, A, B, G, H, P, Q)

        if x == X:
            break

    nom = a - A
    denom = B - b

    res = pow(denom, -1, Q) * nom % Q
    if pow(G, res, P) == H:
        return res

    return res + Q



def test(p, a, b):
    print(_discrete_log_pollard_rho(p, a, b))


if __name__ == '__main__':
    # runTaskOne()
    runTaskTwo()
