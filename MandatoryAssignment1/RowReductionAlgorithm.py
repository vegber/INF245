import numpy as np


def rref(matrix):
    if not matrix: return

    numRows = len(matrix)
    numCols = len(matrix[0])

    i, j = 0, 0
    while True:
        if i >= numRows or j >= numCols:
            break

        if matrix[i][j] == 0:
            nonzeroRow = i
            while nonzeroRow < numRows and matrix[nonzeroRow][j] == 0:
                nonzeroRow += 1

            if nonzeroRow == numRows:
                j += 1
                continue

            temp = matrix[i]
            matrix[i] = matrix[nonzeroRow]
            matrix[nonzeroRow] = temp

        pivot = matrix[i][j]
        matrix[i] = [elem_row / pivot for elem_row in matrix[i]]

        for otherRow in range(0, numRows):
            if otherRow == i:
                continue
            if matrix[otherRow][j] != 0:
                matrix[otherRow] = [y - matrix[otherRow][j] * x
                                    for (x, y) in zip(matrix[i], matrix[otherRow])]
        i += 1
        j += 1



    return matrix


# write rows in row echelon form
def upper_triangular(M):
    # move all zeros to buttom of matrix
    M = np.concatenate((M[np.any(M != 0, axis=1)], M[np.all(M == 0, axis=1)]), axis=0)

    # iterate over matrix rows
    for i in range(0, M.shape[0]):

        # initialize row-swap iterator
        j = 1

        # select pivot value
        pivot = M[i][i]

        # find next non-zero leading coefficient
        while pivot == 0 and i + j < M.shape[0]:
            # perform row swap operation
            M[[i, i + j]] = M[[i + j, i]]

            # incrememnt row-swap iterator
            j += 1

            # get new pivot
            pivot = M[i][i]

        # if pivot is zero, remaining rows are all zeros
        if pivot == 0:
            # return upper triangular matrix
            return M

        # extract row
        row = M[i]

        # get 1 along the diagonal
        M[i] = row / pivot

        # iterate over remaining rows
        for j in range(i + 1, M.shape[0]):
            # subtract current row from remaining rows
            M[j] = M[j] - M[i] * M[j][i]

    # return upper triangular matrix
    return M


if __name__ == '__main__':
    test_matrix = (
        (1, 2, -2, 0, 1),
        (-1, 0, 0, 1, 1),
        (0, 2, -1, 0, 1))
    TEST_MODULO = 15

    out = (upper_triangular(test_matrix))

