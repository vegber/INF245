import math
import random
import sys


def fab(x, a, b, alfa, beta, N, n):
    if x % 3 == 0:
        x = pow(x, 2, N)  # x * x % N
        a = a * 2 % n
        b = b * 2 % n

    elif x % 3 == 1:
        x = x * alfa % N
        a = (a + 1) % n
    else:
        x = x * beta % N
        b = (b + 1) % n

    return x, a, b


def solveForX(a, b, A, B, n):
    print(f"a {(A-a)} * {pow(b-B, -1, n)} % {n}")
    return ((A - a) * pow((b - B), -1, n)) % n


def checkAnswer(generator, x, p, y):
    return y == pow(generator, x, p)


def f(N=949772751547464211, n=4748626326421, alfa=314668439607541235, beta=254337213994578435, y=254337213994578435):
    x, a, b, = 1, 0, 0
    X, A, B = x, a, b
    print("%8s %15s %15s %15s  %20s %20s %15s\n" % ("i", "x", "a", "b", "X", "A", "B"))
    for i in range(1, n):
        x, a, b = fab(x, a, b, alfa, beta, N, n)
        X, A, B = fab(X, A, B, alfa, beta, N, n)
        X, A, B = fab(X, A, B, alfa, beta, N, n)
        if x == X and math.gcd(b - B, n) == 1:
            print("%8d %15d %15d %15d   %20d %20d %15d\n" % (i, x, a, b, X, A, B))
            x = solveForX(a, b, A, B, n)
            print(f"x is: {x}")
            break


# f(41, 40, 6, 7, 39)
f()