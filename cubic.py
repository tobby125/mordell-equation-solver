def sol(k, n):
    d = (81 * k * 2**(n + 1)
         + ((81 * k * 2**(n + 1) - 162 * k)**2 + 4*(18*k - 3)**3)**(1/2)
         - 162*k)**(1/3)

    a = 2**(1/3) * (18*k - 3)
    return -a/(3*d) + d/(3*2**(1/3)) + 1

def solve(a, b, c, d):
    A = -b ** 3 / (27 * a ** 3) + b * c / (6 * a ** 2) - d / (2 * a)
    B = (c / (3 * a) - b ** 2 / (9 * a ** 2)) ** 3
    C = (A ** 2 + B) ** (1 / 2)

    s = (A + C) ** (1 / 3) + (A - C) ** (1 / 3) - b / (3 * a)

    assert s.imag == 0
    assert s.real >= 0
    return s.real

def quadratic(a, b, c):
    D = b**2 - 4*a*c
    if D < 0:
        return set()
    return {(-b + D**(1/2)) / (2*a), (-b - D**(1/2)) / (2*a)}
