import copy
import math


def ContinuedFractionAlgorithm(N: int, e: int):
    """
    Let N, e be an RSA pub key, where ed congruent with 1 mod phi(N)
    d is RSA secret exponent.
    We know d is rel. small.

    Task: Factor N = pq
    :return:
    """
    N = copy.deepcopy(N)
    cf = []
    q = math.floor(N)
    cf.append(q)
    N = N - 1
    i = 0
    while N != 0 and i < int(e):
        q = math.floor(1 / N)
        cf.append(q)
        N = 1 / N - q
        i += 1
    return cf


if __name__ == '__main__':
    N = "1098667602555738997359355609545017255443451483195436981354791765341639135658156206242197992115989996829728203054347117299"
    e = "415884005149779743103130951097941968793336745983309971301432658775763996247677181243042840232106535367251782466233724389"

    var = ContinuedFractionAlgorithm(int(N), int(e))
    print(var)
