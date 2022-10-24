"""
Digital Signature Algorithm and Discrete Logarithms
"""
import math
import random
from MandatoryAssignment1.A1.Mandatory1 import BinaryExponentiationWithoutRecursion as binExp


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


if __name__ == '__main__':
    runTaskOne()
