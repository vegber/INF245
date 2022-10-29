import random


def f(y, residueAlfa, generator, alfa, beta, P):
    subset = y % 3

    if subset == 0:
        y = residueAlfa * y
        alfa = alfa + 1
        beta = beta
    elif subset == 1:
        y *= y
        alfa = alfa * 2
        beta = beta * 2
    else:
        y = generator * y
        alfa = alfa
        beta = beta + 1
    return y % P, alfa % P, beta % P


def SolveForX():
    print("got here")


def rhomethod(residue, generator, P):

    alfa, beta = 0, random.randint(1, P -1)
    y = (generator ** beta) % P

    for i in range(P):

        y_i, alfa_i, beta_i = f(y, residue, generator, alfa, beta, P)

        y_2i, alfa_2i, beta_2i = f(y, residue, generator, alfa, beta, P)
        y_2i, alfa_2i, beta_2i = f(y_2i, residue, generator, alfa, beta, P)

        if y_i == y_2i:
            SolveForX()
            break
        y, alfa, beta = y_i, alfa_i, beta_i

rhomethod(7, 6, 41)