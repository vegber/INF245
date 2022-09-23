import operator
from random import randint

import numpy as np


def pop_zeros(lst):
    lst = lst.copy()
    while lst[-1] == 0:
        lst.pop()
    return lst


def EEA(a, b):
    '''
    For integers a >= b
    Finds the integer d = gcd(a, b) along with the integers u, v s.t.
    a*u + b*v = d
    :param a: int
    :param b: int
    :return: u, v, gcd(a, b)
    '''
    assert (isinstance(a, int))
    assert (isinstance(b, int))
    assert (a >= b)
    A1 = np.array([1, 0, a])
    B = np.array([0, 1, b])
    A2 = B
    while A2[2] != 0:
        q = int(np.floor(A1[2] / A2[2]))
        B = A1 - q * A2
        A1 = A2
        A2 = B

    A1 = [int(k) for k in A1]

    return A1


def binExpMod(a, m, n):
    '''
    iteratively calculates a^m mod n
    :param a: int
    :param m: int
    :param n: int
    :return: b = a^m mod n
    '''
    res = 1  # Initialize result

    # Update a modulo n
    a = a % n

    if a == 0:
        return 0

    while m > 0:
        # If m is odd,
        # multiply result with a
        if m % 2 == 1:
            res = (res * a) % n

        # now m is even, half it
        # and update a by squaring it
        m = m // 2
        a = (a * a) % n

    return res


class PolyRing:
    '''
    A class to define polynomials over the finite ring of residues Z/n"
    '''

    def __init__(self, coeff_lst, modulo):
        '''
        Instatiate a polynomial by its coefficients, given in increasing
        :param coeff_lst: list of coeffiecients, increasing order i.e
        [1, 2, 4, 0, 1] = 1*x^0 + 2*x^1 + 4*x^2 + 0*x^3 + 1*x^4
        :param coeff_lst: list of coefficients, in order of increasing degree
        :param modulo: the modulo n of the relevant residue ring
        '''
        assert (isinstance(modulo, int))
        assert (modulo > 1)
        for c in coeff_lst:
            assert (isinstance(c, int))

        coeff_lst = [c % modulo for c in coeff_lst]

        self.modulo = modulo

        degree = 0
        for d, c in list(enumerate(coeff_lst))[::-1]:  # find the degree of the polynomial
            if c != 0:
                degree = d
                break

        self.degree = degree
        self.coeffs = coeff_lst

    def __str__(self):
        '''
        Display the polynomial without the zero coefficients
        :return: str of the polynomial
        '''
        res = ""
        if self.degree == 0:
            if self.__lead__()[0] == 0:
                return "0"
            return str(self.__lead__()[0])

        for d, c in enumerate(self.coeffs):
            if c != 0 and d == 0:
                res += str(c) + " + "
            if c != 0 and d == 1:
                res += str(c) + "*x" + " + "
            if c != 0 and d > 1:
                res += str(c) + "*x^" + str(d) + " + "

        res = res[:-3]
        # res = res[::-1]

        return res

    def __eq__(self, other):
        '''
        Determine if this polynomial is equal to the other polynomial
        :param other:
        :return: boolean true or false
        '''
        assert (isinstance(other, PolyRing))
        assert (self.modulo == other.modulo)

        if self.degree != other.degree:
            return False

        m_d = max(len(self.coeffs), len(other.coeffs))
        # make the coefficient lists the same length
        self.coeffs += [0] * (m_d - len(self.coeffs))
        other.coeffs += [0] * (m_d - len(other.coeffs))

        assert (len(self.coeffs) == len(other.coeffs))

        for i in range(m_d - 1):
            if self.coeffs[i] % self.modulo != other.coeffs[i] % other.modulo:
                return False
        return True

    def __copy__(self):
        '''
        make a shallow copy of this polynomial f
        :return: h, a shallow copy of f
        '''
        h_coeffs = self.coeffs.copy()
        h_modulo = 0 + self.modulo

        return PolyRing(h_coeffs, h_modulo)

    def __neg__(self):
        '''
        Returns the additive inverse of the polynomial
        :return:
        '''
        h_coeffs = self.coeffs.copy()
        h_modulo = 0
        h_modulo += self.modulo

        for i in range(len(h_coeffs)):
            h_coeffs[i] = -h_coeffs[i]
        return PolyRing(h_coeffs, h_modulo)

    def __add__(self, other):
        '''
        Compute the sum h of this polynomial f with another entry polynomial g modulo n.
        :param other polynomial g
        :return: polynomial h
        '''
        assert (isinstance(other, PolyRing))
        assert (self.modulo == other.modulo)
        m_d = max(len(self.coeffs), len(other.coeffs))

        # make the coefficient lists the same length
        self.coeffs += [0] * (m_d - len(self.coeffs))
        other.coeffs += [0] * (m_d - len(other.coeffs))

        # make shallow copies of the polynomial f and its modulo
        h_coeffs = self.coeffs.copy()
        h_modulo = 0
        h_modulo += self.modulo

        for i in range(m_d - 1):
            h_coeffs[i] += other.coeffs[i]
            h_coeffs[i] %= h_modulo

        return PolyRing(h_coeffs, h_modulo)

    def __lead__(self):
        '''
        Find the leading term of the polynomial, e.g. f(x) = 3*x^2 + 1 -> c=3, d=2
        :return: c, d
        '''
        d = self.degree
        c = self.coeffs[self.degree]
        return c, d

    def __scale__(self, a):
        '''
        Computes the product a * f of the polynomial f and integer a
        :param a: int
        :return: a * f
        '''
        assert (isinstance(a, int))
        hx = self.__copy__()
        hx.coeffs = [(a * i) % hx.modulo for i in hx.coeffs]
        return hx

    def __compress__(self):
        zero = PolyRing([0], self.modulo)
        if self == zero:
            return zero
        return PolyRing(pop_zeros(self.coeffs), 0 + self.modulo)

    def __polyprod__(self, other):
        '''
        Compute the product h of this polynomial f with the other polynomial g modulo n
        :param other:
        :return:
        '''
        assert (isinstance(other, PolyRing))
        assert (self.modulo == other.modulo)
        m_d = max(len(self.coeffs), len(other.coeffs))

        # make the coefficient lists the same length
        self.coeffs += [0] * (m_d - len(self.coeffs))
        other.coeffs += [0] * (m_d - len(other.coeffs))

        # instantiate the product h
        h_modulo = 0 + self.modulo
        hx = PolyRing([0] * 2 * m_d, h_modulo)

        for i in range(0, self.degree + 1):
            res = [0] * 2 * m_d
            for j in range(0, other.degree + 1):
                res[i + j] = (self.coeffs[i] * other.coeffs[j]) % h_modulo
            hx += PolyRing(res, h_modulo)

        return hx

    def __longdiv__(self, other):
        '''
        Compute the quotient q and remainder r of this polynomial f with divided by entry polynomial g!=0 modulo n,
        if deg(f) < deg(g) it returns q=0, r=f, NOTE: the modulo n must be prime
        :param other polynomial g
        :return: polynomials q, r s.t. f = q*g + r
        '''
        assert (isinstance(other, PolyRing))
        assert (self.modulo == other.modulo)
        assert (other != PolyRing([0], other.modulo))
        m_d = max(len(self.coeffs), len(other.coeffs))

        # make the coefficient lists the same length
        self.coeffs += [0] * (m_d - len(self.coeffs))
        other.coeffs += [0] * (m_d - len(other.coeffs))

        gc = other.__copy__()

        # instantiate polynomial r
        rx = self.__copy__()

        # instantiate polynomial q
        qx = PolyRing([0] * m_d, 0 + self.modulo)

        if self.degree < other.degree:  # trivial case, but unnecessary as the while loop still works
            return qx, rx

        while rx != PolyRing([0], rx.modulo) and rx.degree >= other.degree:
            t_coeffs = [0] * m_d  # divide away the leading terms
            cr, dr = rx.__lead__()
            cg, dg = other.__lead__()
            z = EEA(rx.modulo, cg)[1]  # find the multiplicative inverse of cg modulo n
            # print(z)
            ct = (cr * z) % rx.modulo  # calculate the quotient of the lead terms
            dt = dr - dg
            t_coeffs[dt] = ct
            tx = PolyRing(t_coeffs, rx.modulo)
            qx += tx
            rx = rx + (-(tx.__polyprod__(gc)))
            # print(self, "/", gc, "=", qx, "with remainder: ", rx)

        assert (qx.__polyprod__(gc) + rx == self)

        return qx, rx

    def __evaluate__(self, x0):
        '''
        Calculate the value of the polynomial f at x0 modulo n
        :param x0:
        :return: f(x0) modulo n
        '''
        assert (isinstance(x0, int))
        res = 0
        for d, c in list(enumerate(self.coeffs[0:self.degree + 1]))[::-1]:
            base = pow(x0, d, self.modulo)
            base *= c
            base %= self.modulo
            res += base
            res %= self.modulo
        return res

    def __gcd__(self, other):
        '''
        For polynomial f. Returns the gcd of the two polynomials.
        :param other: polynomial g
        :return: gcd(f, g)
        '''
        assert (isinstance(other, PolyRing))
        # assert(self.degree >= other.degree)

        poly_ordered = sorted([self, other], key=operator.attrgetter('degree'))
        b = poly_ordered[0].__copy__()  # g
        a = poly_ordered[1].__copy__()  # f>g

        q, r = a.__longdiv__(b)

        while r != PolyRing([0], self.modulo):
            a = b
            b = r
            q, r = a.__longdiv__(b)

        return b

    def __EEA__(self, other, monic=False):
        '''
        Calculate the polynomial d = gcd(f, g) and the polynomials u, v s.t. this u*f + v*g = d modulo n
        :param other: polynomial g
        :return: gcd(f, g), u, v polynomials
        '''
        assert (isinstance(other, PolyRing))
        assert (self.modulo == other.modulo)
        if self == PolyRing([0], self.modulo) or other == PolyRing([0], self.modulo):
            return 0, 0, 0

        poly_ordered = sorted([self, other], key=operator.attrgetter('degree'))
        r0 = poly_ordered[0].__copy__()  # g
        r1 = poly_ordered[1].__copy__()  # f>g

        s0 = PolyRing([1], self.modulo)
        s1 = PolyRing([0], self.modulo)
        t0 = PolyRing([0], self.modulo)
        t1 = PolyRing([1], self.modulo)

        while r1 != PolyRing([0], self.modulo):
            # print(r1)
            q = r0.__longdiv__(r1)[0]
            # print("quotient q=r0/r1 :", q)
            r0 = r1
            r1 = r0 + -(q.__polyprod__(r1))
            s0 = s1
            s1 = s0 + -(q.__polyprod__(s1))
            t0 = t1
            t1 = t0 + -(q.__polyprod__(t1))

        r0LeadInv = EEA(self.modulo, r0.__lead__()[0])[1]

        assert (self.__longdiv__(r0)[1] == PolyRing([0], self.modulo))
        assert (other.__longdiv__(r0)[1] == PolyRing([0], self.modulo))

        if monic:
            r0 = r0.__scale__(r0LeadInv)
            s0 = s0.__scale__(r0LeadInv)
            t0 = t0.__scale__(r0LeadInv)
            assert (r0.__lead__()[0] == 1)

        return r0, s0, t0

    def __invert__(self, other):
        '''
        Calculates the inverse of a polynomial with respect to modulo polynomial other
        and the integer ring
        :param: other - the modulo polynomial
        :return: f^-1 with resp. to modulo polynomial other in the integer ring
        '''
        assert (isinstance(other, PolyRing))
        assert (self.modulo == other.modulo)
        assert (other != PolyRing([0], self.modulo))

        gcd, inv, quo = self.__EEA__(other=other, monic=True)
        if gcd.degree > 0:
            return inv
        inv = inv.__scale__(gcd)
        assert (inv * self == PolyRing([1], self.modulo))
        return inv

    def __polyMod__(self, other):
        '''
        Calculates the residue of this polynomial f modulo a polynomial g
        :param other: polynomial g
        :return: residue r
        '''
        assert (isinstance(other, PolyRing))
        assert (self.modulo == other.modulo)
        assert (other != PolyRing([0], self.modulo))

        return self.__longdiv__(other)[1]

    def __binExpPolyMod__(self, other, a):
        '''
        Calculates for this polynomial f, some modulo polynomial g and a positive integer a
         (f(x))^a mod g(x)
        :param other: g(x) polynomial
        :param a: positive integer
        :return: polynomial h(x) = (f(x))^a mod g(x)
        '''
        assert (isinstance(other, PolyRing))
        assert (isinstance(a, int))
        assert (a >= 0)

        res = PolyRing([1], self.modulo)

        if a == 0:
            return res.__scale__(0)
        if a == 1:
            return res.__polyMod__(other)

        g = other.__copy__()
        h = self.__copy__().__polyMod__(g)

        while a > 0:
            # print(a)
            # If a is odd,
            # multiply result with h
            if a % 2 == 1:
                res = (res.__polyprod__(h)).__polyMod__(g)
                res = res.__compress__()  # the computations on res makes the coefficient_lst unnecessarily long

            # now a is even, half it and round down
            # and update h by squaring it
            a = a // 2
            h = (h.__polyprod__(h)).__polyMod__(g)
            h = h.__compress__()  # the computations on h makes the coefficient_lst unnecessarily long

        return res

    def __findRoot__(self):
        '''
        For a given finite univariate polynomial f(x) in Z/p, where p prime, calculate a root r s.t. f(r) = 0.
        :return: r, root of f or None if f has no root
        '''

        done = False
        res = 0

        fx = self.__copy__()

        hx = PolyRing([0, 1], fx.modulo)
        hx = hx.__binExpPolyMod__(fx, fx.modulo)
        hx = hx + -(PolyRing([0, 1], fx.modulo))
        # hx = x^p - x
        # print("hx =", hx)

        gx = fx.__gcd__(hx)  # gx = gcd(f, x^p - x)
        gx = gx.__compress__()
        gx = gx.__scale__(EEA(fx.modulo, gx.coeffs[0])[1])  # monic

        if gx.degree == 0:
            done = True
            res = None

        if gx.degree == 1:
            done = True
            res = (res - EEA(gx.modulo, gx.coeffs[1])[1] * gx.coeffs[0]) % gx.modulo

        while not done:
            # print("gx =", gx)
            assert (gx.degree >= 2)  # this should be the case we are in

            b = randint(2, fx.modulo)  # take a random b mod(p)

            bx = PolyRing([b, 1], fx.modulo)
            one = PolyRing([1, 0], fx.modulo)
            exp = int((fx.modulo - 1) / 2)

            bx = (bx + -one).__binExpPolyMod__(fx, exp).__compress__()
            bx = bx.__compress__()
            # print("bx =", bx)
            vx = gx.__gcd__(bx).__compress__()
            # print("vx =", vx)

            if vx == PolyRing([1, 0], fx.modulo) or vx == gx:
                continue

            elif vx.degree == 1:
                done = True
                res = (res - EEA(gx.modulo, vx.coeffs[1])[1] * vx.coeffs[0]) % fx.modulo

            elif 2 <= vx.degree < gx.degree:
                kx, rx = gx.__longdiv__(vx)
                assert (rx == PolyRing([0], fx.modulo))  # just to be sure
                if kx.degree < vx.degree:
                    gx = kx.__compress__()
                else:
                    gx = vx.__compress__()

        return res

    def __findAllRoots__(self):
        '''
        For a given finite univariate polynomial f(x) in Z/p, where p prime, calculate all roots r_i s.t. f(r_i) = 0.
        For a polynomial of degree d, there are at least d roots counting multiplicity.
        :return: roots, a list of roots for f
        '''
        assert (self != PolyRing([0], self.modulo))

        roots = []

        fx = self.__copy__()
        ntrys = 0
        while fx.degree > 0 and ntrys < fx.degree * 20:
            root = fx.__findRoot__()
            ntrys += 1
            if root:
                # print(f"f({root}) =", fx.__evaluate__(root))
                fx = fx.__scale__(EEA(fx.modulo, fx.__lead__()[0])[1])
                roots.append(root)
                root = (fx.modulo - root) % fx.modulo
                factor = PolyRing([root, 1], fx.modulo)
                fx, zero = fx.__longdiv__(factor)
                assert (zero == PolyRing([0], fx.modulo))

        return roots


f_x = [3, 6, 3, 2, 3]
g_x = [2, 0, 5, 1]

PRIME = 7

obj_1 = PolyRing(f_x, modulo=PRIME)
obj_2 = PolyRing(g_x, modulo=PRIME)
a, b =(obj_1.__longdiv__(obj_2))
print(b)
