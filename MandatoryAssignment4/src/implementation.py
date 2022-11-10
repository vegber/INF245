###
# NTRU Assignment INF 245
import math
import sys
from math import log

import numpy as np
from sympy import Poly, symbols, GF, invert, isprime, gcdex, div
from sympy.abc import x


def ntruEncDec():
    """
    1.1)
        Compute the inversion of
            f mod (x^N - 1)
            p mod (x^N - 1)
            q ??

    1.2)
        Compute the NTRU public key

    1.3)
        Compute the NTRU cipher - text
    """

    # (N, p, q) = (11, 3, 32)
    N, p, q = 11, 3, 32
    # Public parameters:
    #   L_f   = L(4, 3)
    #   L_g   = L(3,3)
    #   L_phi = L(3,3)

    print()
    print("=" * 40)
    print(f"\t\t\t\tTASK ONE")
    print("=" * 40)
    print()

    f = Poly(np.array([-1, 1, -1, 0, 1, 0, 1, 0, 1, 0, -1][::-1], dtype=int), x)
    g = Poly(np.array([-1, 1, 0, 0, 1, 1, 0, 0, -1, 0, -1][::-1], dtype=int), x)
    m = Poly(np.array([-1, 0, 0, -1, 1, 0, 0, 0, 1, 1, -1][::-1], dtype=int), x)
    phi = Poly(np.array([1, 0, 1, -1, 1, 0, -1, -1, 0, 0, 0][::-1], dtype=int), x)
    # residue class, x^n - 1
    resClass = Poly(np.array([-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1][::-1], dtype=int), x)

    f_inverse_mod_p = Poly(polyInverse(f, resClass, p), x)
    f_inverse_mod_q = Poly(polyInverse(f, resClass, q), x)

    print(f"fp invers {f_inverse_mod_p}")
    print(f"fq invers {f_inverse_mod_q}")
    print()
    #####
    # Create public key
    # h = p*f_q * g
    pubPoly = createPublicKey(f_inverse_mod_q, p, g, q, resClass)
    print(f"Public key: {pubPoly.all_coeffs()}")

    # Ciphertext
    # e = phi * h + m (mod q)
    ciphertext = createCipherText(pubPoly, m, phi, q, resClass)
    print(f"Ciphertext: {ciphertext}")

    print()
    print("=" * 40)
    print(f"\t\t\t\tTASK TWO")
    print("=" * 40)
    print()

    ciphertext = Poly(np.array([9, 28, 18, 20, 3, 24, 25, 28, 10, 1, 26][::-1], dtype=int), x)
    plaintext = decryptMessage(f, ciphertext, f_inverse_mod_p, q, p, resClass)
    print(f"plaintext: {Poly(plaintext, x)}")


def decryptMessage(f, c, fInvP, q, p, resClass):
    a = f.__mul__(c).__mod__(resClass)
    a = Poly(closeUnderMod(np.array(a.all_coeffs()), q), x)
    plaintext = a.__mul__(fInvP).__mod__(resClass)
    plaintext = Poly((closeUnderMod(np.array(plaintext.all_coeffs(), dtype=int), p)), x)
    return plaintext


def createCipherText(pubPol, message, phi, q, residueClass):
    """
    :param pubPol:
    :param message:
    :param phi:
    :param q:
    :param residueClass:
    :return: e = r * h + m mod(q)
    """
    rh = phi.__mul__(pubPol).__mod__(residueClass)
    return Poly(closeUnderMod(np.array(rh.__add__(message).all_coeffs(), dtype=int), q, False), x)


def createPublicKey(fInverseQ, p, g, q, residueClass):
    p_f_q = fInverseQ.__mul__(p)
    h = p_f_q.__mul__(g)
    var = h.__mod__(residueClass)
    test = closeUnderMod(np.array(var.all_coeffs(), dtype=int), q)
    return Poly(test, x)


def closeUnderMod(arr, n: int, intervall=True):
    if intervall:
        lower_bound, upper_bound = -n // 2, n // 2
        coeffs_under_mod = []
        for org in arr:
            x = org % n
            if x <= upper_bound:
                coeffs_under_mod.append(x)
            else:
                x = x - n
                coeffs_under_mod.append(x)
    else:
        return [x % n for x in arr]
    return np.array(coeffs_under_mod)


def printPolynomial(p, mod):
    for el in range(len(p)):
        print(f"{p[el]}x^{el} + ", end=" ")
    print(f" mod {mod}", end=" ")


def polyInverse(poly_in, poly_I, poly_mod):
    """
    Find the inverse of the polynomial poly_in in the Galois filed GF(poly_mod)
    i.e. the inverse in
        Z/poly_mod[X]/poly_I
    Inputs and outputs are given as an array of coefficients where
        x^4 + 5x^2 + 3 == [1,0,5,0,3]
    Returns
    =======
    Either an empty array if the inverse cannot be found, or the inverse of the
    polynomial poly_in as an array of coefficients.
    """
    Ppoly_I = Poly(poly_I, x)
    Npoly_I = len(Ppoly_I.all_coeffs())
    if isprime(poly_mod):
        # For prime poly_mod we only need use the sympy invert routine, we then pull out
        # all the coefficients for the inverse and return (not all_coeffs() also includes
        # zeros in the array
        try:
            inv = invert(Poly(poly_in, x).as_expr(), Ppoly_I.as_expr(), domain=GF(poly_mod, symmetric=False))
        except "Not able to invert!":
            return np.array([])
    elif log(poly_mod, 2).is_integer():
        try:
            inv = invert(Poly(poly_in, x).as_expr(), Ppoly_I.as_expr(), domain=GF(2, symmetric=False))
            ex = int(log(poly_mod, 2))
            for a in range(1, ex):
                inv = ((2 * Poly(inv, x) - Poly(poly_in, x) * Poly(inv, x) ** 2) % Ppoly_I).trunc(poly_mod)
            inv = Poly(inv, domain=GF(poly_mod, symmetric=False))
        except "Failed to invert":
            return np.array([])
    else:
        # Otherwise we cannot find the inverse
        return np.array([])

    # If we have got this far we have calculated an inverse, double-check the inverse via poly mult
    tmpCheck = np.array(Poly((Poly(inv, x) * Poly(poly_in, x)) % Ppoly_I,
                             domain=GF(poly_mod, symmetric=False)).all_coeffs(), dtype=int)
    if len(tmpCheck) > 1 or tmpCheck[0] != 1:
        sys.exit("ERROR : Error in calculation of polynomial inverse")

    # Passed the error check so return polynomial coefficients as array
    return padArr(np.array(Poly(inv, x).all_coeffs(), dtype=int), Npoly_I - 1)


def padArr(A_in, A_out_size):
    """
    :param A_in: Take an input numpy integer array A_in and pad with leading zeros.
    :param A_out_size:
    :return: Return the numpy array of size A_out_size with leading zeros
    """
    return np.pad(A_in, (A_out_size - len(A_in), 0), constant_values=0)


def rationalsToZZ(u, modolus):
    # denominator must be co - prime with Q!
    assert math.gcd(u.all_coeffs()[0].as_numer_denom()[1], modolus) == 1
    liz = []
    for el in u.all_coeffs():
        numerator = el.as_numer_denom()[0]
        denominator = el.as_numer_denom()[1]
        liz.append(numerator * pow(denominator, -1, modolus) % modolus)
    return liz


def computePolynomials(cipher_texts, u, q, resClass):
    cij = []
    for ci in range(len(cipher_texts)):
        for cj in range(ci + 1, len(cipher_texts)):
            sub = cipher_texts[ci].__sub__(cipher_texts[cj])
            mul = sub.__mul__(u).__mod__(resClass)
            cij.append(Poly(closeUnderMod(np.array(mul.all_coeffs(), dtype=int), q), x))
    return cij


def computeC(cij, N):
    X = Poly([1] * N, x)
    candidates = []
    for c in range(len(cij)):
        for i in range(-2, 3, 1):
            tmp = div(cij[c], Poly([1, -1], x), domain='ZZ')[0]
            tmp = tmp.add(X.__mul__(i))
            candidates.append(tmp)
    for i in range(len(candidates)):
        for j in range(i + 1, len(candidates)):
            res = candidates[i].sub(candidates[j]).all_coeffs()
            if len([*filter(lambda x: abs(x) <= 2, res)]) == len(res):
                print(buildPhi(candidates[i], candidates[j]))

    return 0, 0


def multipleEncSameMessage():
    N, p, q = 11, 3, 32

    cipher_texts = [
        Poly(np.array([15, 24, 0, 30, 5, 24, 5, 28, 20, 29, 10], dtype=int), x),
        Poly(np.array([20, 24, 12, 27, 19, 3, 3, 23, 28, 2, 29], dtype=int), x),
        Poly(np.array([0, 20, 4, 0, 30, 30, 15, 7, 4, 19, 29], dtype=int), x),
        Poly(np.array([14, 8, 20, 3, 11, 29, 29, 9, 6, 28, 1], dtype=int), x)]

    f = Poly(np.array([-1, 1, -1, 0, 1, 0, 1, 0, 1, 0, -1][::-1], dtype=int), x)
    g = Poly(np.array([-1, 1, 0, 0, 1, 1, 0, 0, -1, 0, -1][::-1], dtype=int), x)
    resClass = Poly(np.array([-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1][::-1], dtype=int), x)
    fInvP = Poly(polyInverse(f, resClass, p), x)
    fInvQ = Poly(polyInverse(f, resClass, q), x)
    publicKey = createPublicKey(fInvQ, p, g, q, resClass)
    u, v, gcd = gcdex(publicKey, resClass)

    u = Poly(rationalsToZZ(u, q), x)

    cij = computePolynomials(cipher_texts, u, q, resClass)

    # find c
    phi_1, phi_2 = computeC(cij, N)
    pass


def buildPhi(phi_1, phi_2):
    phi_1 = phi_1.all_coeffs()
    print(phi_2)
    phi_2 = phi_2.all_coeffs()
    lookup = {
        2: [1, -1],
        1: [1, 0],  # or 0, -1
        0: [0, 0],
        -1: [-1, 0],  # 0, 1
        -2: [-1, 1]}
    for i in range(len(phi_2)):
        lookup_val = (phi_1[i] - phi_2[i])
        if lookup_val < -2 or lookup_val > 2: return
        if lookup_val == 0: continue
        tmp_A = lookup[lookup_val][0]
        tmp_B = lookup[lookup_val][1]

        phi_1[i] = tmp_A
        phi_2[i] = tmp_B
    print("-"*40)
    print(phi_1)
    print(phi_2)


if __name__ == '__main__':
    global x
    x = symbols('x')
    # ntruEncDec()
    multipleEncSameMessage()
