"""
Task ONE - exam 28.11.22
CRT
Use a variation of the Chinese Remainder Theorem  to
    compute an integer x, s.t.
        x ≡ a_1 mod n_1
        x ≡ a_2 mod n_2
    n_1 & n_2 is not co - prime

By the Chinese Remainder Theorem solving a congruence modulo a 
composite number is reducible to solving the same congruence 
modulo its prime power factors
"""


def lcm(a, b, g):
    """
    least common multiplier.
    :param a:
    :param b:
    :param g:
    :return:
    """
    return (a * b) // g


def ext_gcd(a, b):
    if a == 0:
        return b, 0, 1
    gcd, u_1, v_1 = ext_gcd(b % a, a)
    return gcd, v_1 - (b // a) * u_1, u_1


# general chinese remainder theorem.
def GCRT(a, n, b, m):
    g, p, q = ext_gcd(n, m)
    # set K to the least common multiplier.
    K = lcm(n, m, g)
    # no solution.
    if not a % g == b % g:
        return -1, -1

    # use bezout's identity:
    # (n/g)*p + (m/g)*q = 1

    res1 = (a * (m // g) * q) % K
    res2 = (b * (n // g) * p) % K
    res = res1 + res2

    return K, res


if __name__ == '__main__':
    a1 = 970215631
    a2 = 1569366205
    n1 = 1973671601
    n2 = 1973672789

    K, x = GCRT(a1, n1, a2, n2)
    assert x != -1, "No solution"
    print(f"Found x = {x}")

    print(x % n1 == a1 % n1)
    print(x % n2 == a2 % n2)
