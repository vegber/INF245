###
# NTRU Assignment INF 245
import math
import sys
from math import log

import numpy as np
from sympy import Poly, symbols, GF, invert, isprime, gcdex, div, pdiv, pretty
from sympy.abc import x
from MandatoryAssignment4.src.implementation import *

N, p, q = 19, 3, 32
# Lf = 5, 4  #positive/negative
# Lg = 4, 4
# Lø = 4, 4

public_key = Poly(np.array([25, 14, 18, 29, 1, 27, 21, 3, 17, 21, 25, 4, 30, 1, 24, 10, 16, 0, 2][::-1], dtype=int), x)

cipher_texts = [
    Poly(np.array([20, 29, 31, 28, 24, 12, 25, 29, 30, 13, 9, 25, 15, 6, 5, 24, 12, 16, 4][::-1], dtype=int), x),
    Poly(np.array([6, 24, 17, 10, 13, 16, 4, 6, 5, 29, 8, 18, 7, 15, 13, 31, 6, 31, 4][::-1], dtype=int), x),
    Poly(np.array([13, 21, 0, 2, 18, 4, 25, 2, 18, 6, 13, 25, 1, 3, 14, 19, 11, 31, 30][::-1], dtype=int), x),
    Poly(np.array([14, 19, 0, 0, 18, 2, 24, 2, 16, 5, 13, 26, 0, 3, 13, 21, 11, 1, 30][::-1], dtype=int), x),
    Poly(np.array([13, 20, 1, 0, 19, 2, 24, 3, 17, 6, 13, 27, 1, 3, 13, 21, 12, 31, 31][::-1], dtype=int), x),
    Poly(np.array([20, 16, 4, 2, 16, 4, 7, 5, 27, 0, 5, 31, 21, 2, 30, 2, 2, 4, 26][::-1], dtype=int), x),
    Poly(np.array([13, 21, 0, 1, 19, 3, 25, 3, 17, 5, 13, 25, 31, 1, 14, 20, 11, 0, 0][::-1], dtype=int), x),
    Poly(np.array([9, 1, 0, 4, 10, 8, 29, 23, 17, 8, 2, 29, 19, 21, 17, 7, 25, 28, 27][::-1], dtype=int), x)]

# find public inverse

res_clas = Poly(np.array([-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1][::-1], dtype=int), x)

public_key = public_key.__sub__(288) # Trekker fra, h(1) - siden pub blir co prime under eea...

u, v, gcd = gcdex(public_key, res_clas)
u = Poly(rationalsToZZ(u, q), x)

cij = computePolynomials(cipher_texts, u, q, res_clas)
_ = [print(x) for x in cij]
phi_s = computeC(cij, N, q, res_clas, d=4)
