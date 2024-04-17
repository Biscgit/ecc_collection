# Lenstra elliptic-curve factorization in Python
# by David Horvát, 2024

from __future__ import annotations
import math
import random

# set number to be factorized
number = <set value>

# adjust if too slow
max_iterations = 10
max_factor = 1_000


class InvalidCurve(Exception):
    pass


class WeierStrassEC:
    """Creates a Weierstrass elliptic curve from a, b and p"""

    def __init__(self, a: int, b: int, p: int):
        if (4 * pow(a, 3) + 27 * pow(b, 2)) % p == 0:
            raise InvalidCurve("Invalid curve parameters `a` and `b`")

        self.a = a
        self.b = b
        self.p = p

    def __eq__(self, other: WeierStrassEC) -> bool:
        return self.a == other.a and self.b == other.b and self.p == other.p


class Point:
    """Defines a Point on an elliptic curve"""

    def __init__(self, x: int, y: int | math.inf, ecc: WeierStrassEC):
        self.curve = ecc
        self.x = x % ecc.p
        self.y = int(y % ecc.p) if math.isfinite(y) else y

    def is_infinite(self) -> bool:
        return math.isinf(self.y)

    def __repr__(self) -> str:
        curve = self.curve

        if self.is_infinite():
            return f"on `EC[a={curve.a}, b={curve.b}, p={curve.p}]` with `Point[x={self.x}, y=∞]`"

        return f"on `EC[a={curve.a}, b={curve.b}, p={curve.p}]` with `Point[x={self.x}, y={self.y}]`"

    def __eq__(self, other: Point) -> bool:
        return self.x == other.x and self.y == other.y

    def __add__(self, other: Point) -> Point:
        """Adds two points on the same elliptic curve
        Add the same point for point-doubling"""

        curve = self.curve

        # points must be on the same curve
        if curve != other.curve:
            raise Exception("Adding points from different curves")

        # match for infinite point
        if self.is_infinite():
            return Point(self.x, math.inf, curve)
        if other.is_infinite():
            return Point(other.x, math.inf, curve)

        # calculate slope s
        try:
            s = self.get_slope(other)

        except (ValueError, ZeroDivisionError):
            return Point(other.x, math.inf, curve)

        # calculate new coordinates x and y
        x = (pow(s, 2) - self.x - other.x) % curve.p
        y = (s * (self.x - x) - self.y) % curve.p

        return Point(x, y, curve)

    def get_slope(self, other: Point) -> int:
        """Calculates the slope of adding two points/doubling.
        Raises ValueError or ZeroDivisionError for an infinite slope"""

        p = self.curve.p

        # point doubling
        if self == other:
            denominator = (2 * self.y)
            numerator = (3 * pow(self.x, 2) + self.curve.a) % (p * denominator)

        # point addition
        else:
            denominator = (other.x - self.x)
            numerator = (other.y - self.y) % (p * denominator)

        return numerator * pow(denominator, -1, p) % p

    def __mul__(self, scalar: int) -> Point:
        """Multiply a point by a provided scalar
        Uses the double-and-add algorithm"""

        res = self

        for position in range(scalar.bit_length() - 1, 0, -1):
            bit = scalar >> (position - 1) & 0b1
            res += res

            # add
            if bit:
                res += self

        return res

    __rmul__ = __mul__

    def lenstra(self) -> int | None:
        """Runs the lenstra algorithm to find p in n = p * q
        Returns an integer if found and None if nothing has been found."""

        point = self
        next_point = self

        p = self.curve.p

        def lenstra_mul(scalar: int) -> bool | None:
            """Multiplies point by a scalar.
            On finding a point in infinity return true"""

            nonlocal point, next_point

            for position in range(scalar.bit_length() - 1, 0, -1):
                bit = scalar >> (position - 1) & 0b1

                # double
                if (next_point := point + point).is_infinite():
                    return True

                point = next_point

                # add
                if bit:
                    if (next_point := point + self).is_infinite():
                        return True

                    point = next_point

        # check points with [k!]G
        for factorial in range(2, max_factor):
            if lenstra_mul(factorial):
                factor = math.gcd((point.x - next_point.x) % p, p)

                # ensure the result is not 1 or p
                return factor if p > factor > 1 else None


def run_lenstra(n: int, stdout=False) -> int | None:
    """Runs the full lenstra algorithm. This can be used by other modules.
    Adjust the parameters at the top of this file"""

    # filter even numbers
    if (n - 1) & 0b1:
        return 2

    for i in range(max_iterations):
        # limit random number size
        x = random.randint(0, math.isqrt(n))
        y = random.randint(0, math.isqrt(n))
        a = random.randint(0, math.isqrt(n))
        b = (pow(y, 2) - pow(x, 3) - a * x) % n

        try:
            # guaranteed point on the curve
            elliptic_curve = WeierStrassEC(a, b, n)
            start_point = Point(x, y, elliptic_curve)

        except InvalidCurve:
            continue

        if factor := start_point.lenstra():
            if stdout:
                print(f"Found factors `p={factor}` and `q={n // factor}` for `{n=}`\n"
                      f"using Weierstrass (`{a=}`, `{b=}`)\nand Point `({x=}, {y=})` "
                      f"on {i + 1}-{dict({1: 'st', 2: 'nd', 3: 'rd'}).get(i + 1, 'th')} iteration.\n")

            return factor


if __name__ == '__main__':
    print(f"---\nRunning Lenstra elliptic-curve factorization for `{number}` with max `{max_iterations}` iterations\n")
    start_time = time.perf_counter()

    if run_lenstra(number, stdout=True) is None:
        print("No factors found!")

    print(f"Factorization took {(time.perf_counter() - start_time) * 1e3:.2f}ms")
