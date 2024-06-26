# Schoof's algorithm and other helper algorithms in Python
# by David Horvát, 2024

import functools
import math
import typing
import random
from math import gcd

import sympy


def mod_sqrt(number: int, p: int) -> int:
    """Uses the Tonelli-Shanks-Algorithm to find the modular square root of a given number"""
    s = (p - 1 & -(p - 1)).bit_length() - 1
    q = (p - 1) // pow(2, s)

    # check for shortcut
    if s == 1:
        return pow(number, (p + 1) // 4, p)

    # full tonelli shanks
    else:
        z = 0
        while pow(z, (p - 1) // 2, p) != p - 1:
            z = random.randint(0, p)

        c = pow(z, q, p)
        y = pow(number, (q + 1) // 2)
        t = pow(number, q)
        m = s

        while t % p != 1:
            i = 1
            while pow(t, pow(2, i), p) != 1:
                i += 1
                if i > m:
                    return 0

            # raises TypeError if no square-root is possible
            d = pow(c, pow(2, m - i - 1), p)

            c = pow(d, 2, p)
            y = (y * d) % p
            t = (t * pow(d, 2)) % p
            m = i

        return y


def get_primes_list(q: int):
    """returns a list of primes with prod(primes) > 4 * sqrt(q)"""

    def get_next_prime(base: int) -> int:
        """returns the highest prime using simple division checking including self"""
        for x in range(2, math.isqrt(base) + 1):
            if base % x == 0:
                return get_next_prime(base + 1)
        return base

    primes = [2]
    product = 2
    limit = math.ceil(4 * math.sqrt(q))

    while product < limit:
        next_prime = get_next_prime(primes[-1] + 1)

        primes.append(next_prime)
        product *= next_prime

    return primes


def is_quadratic_residue(k: int, p: int) -> bool:
    """using Legendre Symbol"""
    return pow(k, (p - 1) // 2, p) == 1


def f_m(m: int, a: int, b: int, x: int, y: int) -> int:
    """function call interface for division polynomials"""
    if m % 2 == 0:
        return div_pol(m, a, b, x, y) // y
    else:
        return div_pol(m, a, b, x, y)


@functools.cache
def div_pol(order: int, a: int, b: int, x: int, y: int) -> int:
    """division polynomial using only x as parameter"""

    if order in [0, 1]:
        return order

    elif order == 2:
        return 2 * y

    elif order == 3:
        return (
                3 * pow(x, 4) +
                6 * a * pow(x, 2) +
                12 * b * x -
                pow(a, 2)
        )

    elif order == 4:
        return (4 * y * (
                pow(x, 6) +
                5 * a * pow(x, 4) +
                20 * b * pow(x, 3) -
                5 * pow(a, 2) * pow(x, 2) -
                4 * a * b * x -
                8 * pow(b, 2) -
                pow(a, 3)
        ))

    m = order // 2
    if order % 2 == 0:
        return div_pol(m, a, b, x, y) // (2 * y) * (
                div_pol(m + 2, a, b, x, y) * pow(div_pol(m - 1, a, b, x, y), 2) -
                div_pol(m - 2, a, b, x, y) * pow(div_pol(m + 1, a, b, x, y), 2)
        )

    else:
        return (
                div_pol(m + 2, a, b, x, y) * pow(div_pol(m, a, b, x, y), 3) -
                div_pol(m - 1, a, b, x, y) * pow(div_pol(m + 1, a, b, x, y), 3)
        )


if __name__ == '__main__':
    # define your curve parameters with E: y² = x³ + ax + b and q = p^b
    q = 19
    a = 11
    b = 5

    primes = get_primes_list(q)
    t = dict()

    x = sympy.symbols('x')
    t[2] = 1 if sympy.gcd(pow(x, q) - x, pow(x, 3) + a * x + b) == 1 else 0
