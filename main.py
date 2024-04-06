# Tonelli Shanks Algorithm
import random
import typing


# generates an elliptic curve function which returns y²
def get_weier_strassen(a_val: int, b_val: int, mod_val: int) -> typing.Callable[[int], int]:
    def pow_y_2(x_val: int) -> int:
        return (pow(x_val, 3) + a_val * x_val + b_val) % mod_val

    return pow_y_2


def execute(rhs: int, p: int) -> typing.Optional[tuple[int, int]]:
    # split p_val- 1
    s = (p - 1 & -(p - 1)).bit_length() - 1
    q = (p - 1) // pow(2, s)

    # check for shortcut
    if s == 1:
        print(f"Shortcut possible because {s=}")

        y_pos = pow(rhs, (p + 1) // 4, p)
        y_neg = pow(-rhs, (p + 1) // 4, p)

        return y_pos, y_neg

    else:
        # find a random z with Legendre-Symbol = -1
        z = 0
        while pow(z, (p - 1) // 2, p) != p - 1:
            z = random.randint(0, p)

        # determine common variables
        c = pow(z, q, p)
        y = pow(rhs, (q + 1) // 2)
        t = pow(rhs, q)
        m = s

        while t % p != 1:
            # determine smallest i
            i = 1
            while pow(t, pow(2, i, p), p) != 1:
                i += 1

                # check for point at 0
                if i > m:
                    return 0, 0

            # if not all values are int, no result is possible
            try:
                d = pow(c, pow(2, m - i - 1, p), p)

            except TypeError:
                return None

            # set common variables
            c = pow(d, 2, p)
            y = (y * d) % p
            t = (t * pow(d, 2, p)) % p
            m = i

        return y, (-y) % p


if __name__ == '__main__':
    # define Weierstraß and searched x-value
    a =
    b =
    mod =
    x =

    curve = get_weier_strassen(a, b, mod)
    result = execute(curve(x), mod)

    if result is None:
        print("No matching coordinates found")

    else:
        y_1, y_2 = sorted(result)
        if y_1 == y_2:
            print(f"For x={x} corresponding point is P\u2081({x % mod}|{y_1})")
        else:
            print(f"For x={x} corresponding points are P\u2081({x % mod}|{y_1}) and P\u2082({x % mod}|{y_2})")
