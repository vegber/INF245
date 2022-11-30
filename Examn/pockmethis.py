
def pockAlgorihtm(q, bit_length):
    for j in range(5):
        p1 = q

        k = 1
        primeFound = False
        while not primeFound:
            p2 = (2 * k * p1) + 1
            kp1 = k * p1

            if p2.bit_length() == bit_length