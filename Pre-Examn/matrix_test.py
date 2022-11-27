from sympy import Matrix, pprint, mod_inverse

A = Matrix([
    [1, -2, -2, -2, -1, 2],
    [0, 3, -2, -3, 1, 3],
    [3, 0, 0, 1, -1, 2],
    [3, -3, -2, 0, 1, 1],
    [0, -3, 3, -3, -3, 2]
])

B = Matrix(
    [
        [1, 0, 1],
        [0, 0, 0],
        [2, 2, 0]
    ]
)

def mod(x, modulos):
    n, d = x.as_numer_denom()
    return n * mod_inverse(d, modulos) % modulos


def Row_Reduce(M, modulos):
    # Find row-reduced echelon form of B modulo 5:
    m_ref = M.rref(iszerofunc=lambda x: x % modulos == 0)
    return m_ref[0].applyfunc(lambda x: mod(x, modulos))


n = 456995412589
pprint(Row_Reduce(B, 3))
