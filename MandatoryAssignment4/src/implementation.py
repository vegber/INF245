###
# NTRU Assignment INF 245

import sys
from math import log
import numpy as np
from sympy import Poly, symbols, GF, invert, isprime
from sympy.abc import x


def taskOne():
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
    print("=" * 120)
    print(f"TASK ONE")
    print("=" * 120)
    print()

    f = np.array([-1, 1, -1, 0, 1, 0, 1, 0, 1, 0, -1][::-1], dtype=int)
    g = np.array([-1, 1, 0, 0, 1, 1, 0, 0, -1, 0, -1][::-1], dtype=int)
    m = np.array([-1, 0, 0, -1, 1, 0, 0, 0, 1, 1, -1][::-1], dtype=int)
    phi = np.array([1, 0, 1, -1, 1, 0, -1, -1, 0, 0, 0][::-1], dtype=int)

    irr_poly = np.array([-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1][::-1], dtype=int)
    resClass = Poly(irr_poly, x)
    # f_inv = poly_inv(f, irr_poly, p)
    # print(f_inv)
    # f & g private key

    # Test
    # wikipedia
    forelesning_f = np.array([-1, 1, 1, 0, -1, 0, 1, 0, 0, 1, -1][::-1], dtype=int)
    forelesning_g = np.array([-1, 0, 1, 1, 0, 1, 0, 0, -1, 0, -1][::-1], dtype=int)
    forelesning_m = np.array([-1, 0, 0, 1, -1, 0, 0, 0, -1, 1, 1][::-1], dtype=int)
    forelesning_phi = np.array([-1, 0, 1, 1, 1, -1, 0, -1, 0, 0, 0][::-1], dtype=int)

    forelesning_f = Poly(f, x)
    forelesning_g = Poly(g, x)
    f_inverse_mod_p = Poly(poly_inv(f, irr_poly, p), x)
    f_inverse_mod_q = Poly(poly_inv(f, irr_poly, q), x)
    forelesning_m = Poly(m, x)
    forelesning_phi = Poly(phi, x)

    print(f"fp invers {f_inverse_mod_p}")
    print(f"fq invers {f_inverse_mod_q}")
    print()
    #####
    # Create public key
    # h = p*f_q * g
    pubPoly = createPublicKey(f_inverse_mod_q, p, forelesning_g, q, resClass)
    print(f"Public key: {pubPoly}")

    # Ciphertext
    # e = phi * h + m (mod q)
    ciphertext = createCipherText(pubPoly, forelesning_m, forelesning_phi, q, resClass)
    print(f"Ciphertext: {ciphertext}")

    print()
    print("=" * 120)
    print(f"TASK TWO")
    print("=" * 120)
    print()

    ciphertext = Poly(np.array([9, 28, 18, 20, 3, 24, 25, 28, 10, 1, 26][::-1], dtype=int), x)
    plaintext = decryptMessage(forelesning_f, ciphertext, f_inverse_mod_p, q, p, resClass)
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


def print_polynomial(p, mod):
    for x in range(len(p)):
        print(f"{p[x]}x^{x} + ", end=" ")
    print(f" mod {mod}", end=" ")


def poly_inv(poly_in, poly_I, poly_mod):
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


if __name__ == '__main__':
    global x
    x = symbols('x')
    taskOne()
