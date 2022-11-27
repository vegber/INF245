import math


def writeprimes():
    with open("primes.txt", "w+") as file:
        n = 10000  # first 10000 primes
        tmp_n = 1
        p = 3
        primes = [2]

        while tmp_n < n:

            is_prime = True
            for i in range(3, int(math.sqrt(p) + 1), 2):
                # range with step 2

                if p % i == 0:
                    is_prime = False

            if is_prime:
                primes += [p]
                tmp_n += 1

            p += 2

        for prim in primes:
            file.write(f"{str(prim)},")

writeprimes()