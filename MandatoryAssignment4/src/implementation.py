###
# NTRU Assignment INF 245

import numpy as np
from math import log, gcd

import sympy
from sympy import Poly, symbols, GF, invert, isprime
from sympy.abc import x
import sys


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
    x = symbols('x')

    f = np.array([-1, 1, -1, 0, 1, 0, 1, 0, 1, 0, -1][::-1], dtype=int)
    g = np.array([-1, 1, 0, 0, 1, 1, 0, 0, -1, 0, -1][::-1], dtype=int)
    m = np.array([-1, 0, 0, -1, 1, 0, 0, 0, 1, 1, -1][::-1], dtype=int)
    phi = np.array([1, 0, 1, -1, 1, 0, -1, -1, 0, 0, 0][::-1], dtype=int)

    irr_poly = np.array([-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1][::-1], dtype=int)
    test = Poly(irr_poly, x)
    # f_inv = poly_inv(f, irr_poly, p)
    # print(f_inv)
    # f & g private key

    # Test
    # wikipedia
    forelesning_f = np.array([-1, 1, 1, 0, -1, 0, 1, 0, 0, 1, -1][::-1], dtype=int)
    forelesning_g = np.array([-1, 0, 1, 1, 0, 1, 0, 0, -1, 0, -1][::-1], dtype=int)
    forelesning_g = Poly(forelesning_g, x)
    f_inverse_mod_p = np.flip(poly_inv(forelesning_f, irr_poly, p))
    print_polynomial(f_inverse_mod_p, p)
    print()
    f_inverse_mod_q = Poly(poly_inv(forelesning_f, irr_poly, q), x)

    print(f_inverse_mod_q)
    #####
    # Create public key
    # h = p*f_q * g

    p_f_q = f_inverse_mod_q.__mul__(p)

    h = p_f_q.__mul__(forelesning_g)
    var = h.__mod__(test)
    print(var)


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
    x = symbols('x')
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
    taskOne()
