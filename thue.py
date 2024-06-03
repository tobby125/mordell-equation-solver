from math import sqrt, floor, ceil
import cubic
from sage.interfaces.gp import gp

def disc(a, b, c, d):
    return -27*(a**2*d**2-6*a*b*c*d-3*b**2*c**2+4*a*c**3+4*b**3*d)

def H(a, b, c, d):
    return lambda x, y: 9*(b**2-a*c) * x**2 + 9*(b*c-a*d) * x*y + 9*(c**2-b*d) * y**2

def G(a, b, c, d):
    a1 = -a**2*d+3*a*b*c-2*b**3
    b1 = -b**2*c-a*b*d+2*a*c**2
    c1 = b*c**2-2*b**2*d+a*c*d
    d1 = -3*b*c*d+2*c**3+a*d**2

    return lambda x, y: 27 * (a1*x**3 + 3*b1*x**2*y + 3*c1*x*y**2 + d1*y**3)

def thue(a, b, c, d):
    return lambda x, y: a*x**3 + 3*b*x**2*y + 3*c*x*y**2 + d*y**3

def thue_solve(a, b, c, d):
    f = f'{a}*x^3 + {3*b}*x^2 + {3*c}*x + {d}'
    for x, y in gp(f'thue(thueinit({f}, 1), 1)'):
        yield int(x), int(y)

def thue_pos_irr(K):
    a_lower = 1
    a_upper = 2 * K ** (1 / 4) / (3 * sqrt(3))

    b_lower = lambda a: 0
    b_upper = lambda a: a / 2 + 1 / 3 * sqrt(sqrt(K) - 27 * a ** 2 / 4)

    P2 = lambda a, b: cubic.solve(-4, (3 * a + 6 * b) ** 2, 0, 27 * a ** 2 * K)

    c_lower = lambda a, b: (9 * b ** 2 - P2(a, b)) / (9 * a)
    c_upper = lambda a, b: b - a

    for a in range(ceil(a_lower), floor(a_upper + 1)):
        for b in range(ceil(b_lower(a)), floor(b_upper(a) + 1)):
            for c in range(ceil(c_lower(a, b)), floor(c_upper(a, b) + 1)):
                for d in cubic.quadratic(-27 * a ** 2, -27 * (4 * b ** 3 - 6 * a * b * c),
                                         -27 * (4 * a * c ** 3 - 3 * b ** 2 * c ** 2) - K):
                    d = round(d)
                    if disc(a, b, c, d) == K:
                        yield a, b, c, d

def thue_pos_red(K):
    C_lower = -(K/108)**(1/3)
    C_upper = sqrt(K/27)

    for C in range(ceil(C_lower), floor(C_upper + 1)):
        if C == 0:
            continue
        s = (K + 108*C**3) / (81*C**2)
        if s < 0:
            continue
        for B in sqrt(s), -sqrt(s):
            B = round(B)
            a, b, c, d = 1, B, C, 0
            if disc(a, b, c, d) == K:
                yield a, b, c, d

def thue_neg_irr(K):
    a_lower = 1
    a_upper = (-16*K/27)**(1/4)

    b_lower = lambda a: 0
    b_upper = lambda a: a / 2 + 1 / 3 * sqrt(sqrt(-K/3) - 3 * a ** 2 / 4)

    c_lower = lambda a, b: (1-3*b) / 3
    c_upper = lambda a, b: (-K/(4*a))**(1/3) / 3 + (b**2/a if a >= 2*b else b - a/4)

    for a in range(ceil(a_lower), floor(a_upper + 1)):
        for b in range(ceil(b_lower(a)), floor(b_upper(a) + 1)):
            for c in range(ceil(c_lower(a, b)), floor(c_upper(a, b) + 1)):
                for d in cubic.quadratic(-27 * a ** 2, -27 * (4 * b ** 3 - 6 * a * b * c),
                                         -27 * (4 * a * c ** 3 - 3 * b ** 2 * c ** 2) - K):
                    d = round(d)
                    if disc(a, b, c, d) == K:
                        yield a, b, c, d

def thue_neg_red(K):
    C_lower = 1
    C_upper = sqrt(-K / 27)

    for C in range(ceil(C_lower), floor(C_upper + 1)):
        if C == 0:
            continue
        s = (K + 108 * C ** 3) / (81 * C ** 2)
        if s < 0:
            continue
        for B in sqrt(s), -sqrt(s):
            B = round(B)
            a, b, c, d = 1, B, C, 0
            if disc(a, b, c, d) == K:
                yield a, b, c, d

def mordell_solve(k):
    K = -108 * k
    if K > 0:
        thue_iters = thue_pos_irr, thue_pos_red
    else:
        thue_iters = thue_neg_irr, thue_neg_red
    for thue_iter in thue_iters:
        for a, b, c, d in thue_iter(K):
            for x0, y0 in thue_solve(a, b, c, d):
                X = H(a, b, c, d)(x0, y0) // 9
                Y = G(a, b, c, d)(x0, y0) // 54
                assert Y**2 == X**3 + k
                yield X, Y
                yield X, -Y

if __name__ == '__main__':
    for k in range(-1000, 1000):
        if k == 0:
            continue
        mordell_solve(k)
        solutions = list(set(mordell_solve(k)))
        solutions.sort()
        print(k, solutions)
