from ramanujan import ramanujan

def solution(m, j):
    A = 2**(m+3)*j
    B = 2**(2*m) - 3*2**(m+2)*j + 4*j**2

    for n, x in ramanujan(A, B):
        p = 9*2**m*j - 6*j**2 + 3*j*x
        q = 18*2**(2*m)
        if p % q == 0:
            k = p // q
            if (3*k) % j == 0:
                s = 3*k // j * 2**m
                print((k, n, s))

for m in range(1, 10):
    for j in range(1, 2**m, 2):
        print(f'm = {m}, j = {j}')
        solution(m, j)