import numpy as np

from MandatoryAssignment1.A1.Mandatory1 import SolveRowEchelonForm

def reduceMod2(matrix, MODULO: int):
    A = np.array(matrix, dtype=int)

    i = 0  # row
    j = 0  # column
    while True:
        # find next nonzero column
        while all(A.T[j] == 0):
            j += 1
            # if reached the end, break
            if j == len(A[0]) - 1: break
        # if a_ij == 0 find first row i_>=i with a
        # nonzero entry in column j and swap rows i and i_
        if A[i][j] == 0:
            i_ = i
            while A[i_][j] == 0:
                i_ += 1
                # if reached the end, break
                if i_ == len(A) - 1: break
            A[[i, i_]] = A[[i_, i]]
        # divide ith row a_ij to make it a_ij == 1
        A[i] = A[i] / A[i][j]
        # eliminate all other entries in the jth column by subtracting
        # multiples of the ith row from the others
        for i_ in range(len(A)):
            if i_ != i:
                A[i_] = A[i_] - A[i] * A[i_][j] / A[i][j]
        # if reached the end, break
        if (i == len(A) - 1) or (j == len(A[0]) - 1): break
        # otherwise, we continue
        i += 1
        j += 1

    return A


xi = [
    [2, 2, 1, 0, 0, 100],
    [4, 0, 0, 0, 0, 18],
    [0, 1, 1, 0, 1, 12],
    [1, 0, 0, 1, 1, 62],
    [1, 2, 0, 0, 1, 143],
    [1, 1, 1, 1, 1, 206]
]
# out = rref_mod_n(xi, 228)
# a = np.array(out)
# b = a.reshape(6, 6)

# print(b)

out = reduceMod2(xi, 228)

for x in out:
    print(x)