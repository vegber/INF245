"""
Digital Signature Algorithm and Discrete Logarithms
"""
import math
import random


def find_k_coprime_p(p):
    k_ = random.randint(1, p - 1)
    if math.gcd(k_, (p - 1)) == 1:
        return k_
    else:
        find_k_coprime_p(p)


def forgeElGamalSigAlgorithm(p: int, g: int, y: int, m: int):
    """
    Forge ElGamal signature without knowing the private key
        Construct triple: m, a, b
            m: int
            a, b = Signature
    :param p: prime p
    :param g: primitive root mod p
    :param y: y ≡ 5^x mod p
    :return:
    """
    # system private key: x mod p - 1
    # public key: y ≡ 5^x mod p

    # take random secret residue k mod p - 1 coprime to p - 1.
    # Compute a ≡ g^k mod p, 0 < a < p
    k = find_k_coprime_p(p)
    a = pow(g, k, p)

if __name__ == '__main__':
    # forge ElGamal
    P = 274742860459763
    G = 5
    Y = 262274678376340
    # choose M
    M = 313
    forgeElGamalSigAlgorithm(P, G, Y, M)
