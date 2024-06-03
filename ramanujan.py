import os

'''d = {i: [] for i in range(-10**7, 10**7 + 1)}

d[-12044272] = []
d[-15190000] = []

for file in os.listdir('mordell'):
    with open(f'mordell/{file}') as f:
        i = None
        for line in f:
            x = line.split()
            if not x:
                continue
            try:
                i = int(x[0])
                x = x[2:]
            except ValueError:
                pass
            for j in x:
                d[i].append(tuple(int(k) for k in j[1:-1].split(',')))'''

import thue
# y^3 + k = x^2
def mordell(k):
    r = set(thue.mordell_solve(k))
    return r

# A*2^n + B = x^2
def ramanujan(A, B):
    # get rid of powers of 2 in A and powers of 4 in B
    n_off = 0
    while A % 2 == 0:
        A //= 2
        n_off += 1
    x_mult = 1
    while B % 4 == 0:
        B //= 4
        x_mult *= 2
        n_off -= 2

    sols = []
    # y^3 + k = x^2

    # Case 1: 2^n = y^3
    # (A*y)^3 + A^2*B = (A*x)^2

    for r, s in mordell(A**2 * B):
        en = (r // A) ** 3
        if r % A == 0 and s % A == 0 and s >= 0 and en & (en-1) == 0 and en != 0:
            n = en.bit_length() - 1 - n_off
            x = s // A * x_mult
            sols.append((n, x))

    # Case 2: 2^n = 2*y^3
    # (2A*y)^3 + 4*A^2*B = (2A*x)^2

    for r, s in mordell(4 * A**2 * B):
        en = 2 * (r // (2 * A)) ** 3
        if r % (2*A) == 0 and s % (2*A) == 0 and s >= 0 and en & (en-1) == 0 and en != 0:
            n = en.bit_length() - 1 - n_off
            x = s // (2*A) * x_mult
            sols.append((n, x))

    # Case 3: 2^n = 4*y^3
    # (4A*y)^3 + 16*A^2*B = (4A*x)^2

    for r, s in mordell(16 * A**2 * B):
        en = 4 * (r // (4 * A)) ** 3
        if r % (4*A) == 0 and s % (4*A) == 0 and s >= 0 and en & (en-1) == 0 and en != 0:
            n = en.bit_length() - 1 - n_off
            x = s // (4*A) * x_mult
            sols.append((n, x))

    return sorted(sols)

def test_ramanujan(A, B):
    sols = []

    for n in range(-30, 30):
        S = A*2**n + B
        if S < 0:
            continue
        x = round(S**(1/2))
        if S == x**2:
            sols.append((n, x))

    return sorted(sols)

if __name__ == '__main__':
    A = 16
    B = 28
    print(ramanujan(A, B))
    print(test_ramanujan(A, B))