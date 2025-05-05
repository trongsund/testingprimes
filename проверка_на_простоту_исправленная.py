import math
import random
import numpy as np

def powmod(a, p, N):
    pow2 = 1
    result = 1
    while (pow2 <= p):
        pow2 *= 2
    pow2 //= 2
    while (pow2 > 0):
        if ((p // pow2) % 2 == 1):
            result *= a
        if (pow2 > 1):
            result **= 2
        pow2 //= 2
        result %= N
    return result

def matrix_powmod(a, p, N):
    pow2 = 1
    result = np.eye(a.shape[0], dtype=int)
    while (pow2 <= p):
        pow2 *= 2
    pow2 //= 2
    while (pow2 > 0):
        if ((p // pow2) % 2 == 1):
            result = np.dot(result, a)
        if (pow2 > 1):
            result = np.dot(result, result)
        pow2 //= 2
        result %= N
    return result


def jacobi(P, Q):
    P %= Q
    R = 0
    result = 1
    if (Q % 2 == 0):
        return 0
    while (P > 1):
        while (P % 2 == 0):
            if ((Q ** 2 - 1) % 16 != 0):
                result *= -1
            P //= 2
        R = Q % P
        if ((Q - 1) % 2 == 1 and (P - 1) % 2 == 1):
            result *= -1
        Q = P
        P = R
    return result

def is_prime(N):
    prime = True
    for i in range(2, min([N, 10000000])):
        prime = prime and (N % i != 0)
        if not prime:
            print("Small divisors")
            return prime
    if (N < 10 ** 14):
        return prime
    divisors_f1 = []
    f1 = 1
    r1 = N - 1
    for i in range(2, 10000000):
        while (r1 % i == 0):
            r1 //= i
            f1 *= i
            if (i != 2):
                divisors_f1.append(i)
    random.seed()
    a_example = random.randrange(2, N - 1)
    prime = prime and (powmod(a_example, N - 1, N) == 1) 
    if not prime:
        print("Not Fermat-pseudoprime")
        return prime
    print(divisors_f1)
    suits_a = False
    exists_a = False
    is_primitive_root = True
    for p in divisors_f1:
        for i in range(7):
            a1_example = random.randrange(2, N - 1)
            suits_a = (powmod(a1_example, N-1, N) == 1) and (powmod(a1_example, (N - 1)//p, N) != 1)
            exists_a = exists_a or suits_a
        is_primitive_root = is_primitive_root and exists_a
        exists_a = False
    prime = prime and is_primitive_root
    for i in range(12):
        a1_example = random.randrange(2, N - 1)
        suits_a = (powmod(a1_example, N-1, N) == 1) and (powmod(a1_example, (N - 1)//r1, N) != 1)
        exists_a = exists_a or suits_a
    prime = prime and exists_a
    if not prime:
        return prime
    print("F1-pseudoprime")
    if (N < 10 ** 21):
    	return prime
    divisors_f2 = []
    f2 = 1
    r2 = N + 1
    for i in range(2, 10000000):
        while (r2 % i == 0):
            r2 //= i
            f2 *= i
            if (i != 2):
                divisors_f2.append(i)
    print(divisors_f2)
    suits_a = False
    exists_a = False
    is_primitive_root = True
    iters = 0
    for q in divisors_f2:
        while (iters < 288 and not exists_a):
            a2_example = random.randrange(2, N - 1)
            b2_example = random.randrange(2, N - 1)
            if (jacobi(a2_example ** 2 + 4 * b2_example, N) == -1):
                U_example = np.array([[0, 1], [b2_example, a2_example]])
                vec = np.array([0, 1])
                vecres = np.dot(matrix_powmod(U_example, N+1, N), vec)
                vecresq = np.dot(matrix_powmod(U_example, (N+1)//q, N), vec)
                suits_a = (vecres[0] == 0) and (vecresq[0] % N != 0)
                exists_a = exists_a or suits_a
                if (vecres[0] == 0):
                    print([q, iter])
            iters += 1
        if (not exists_a):
            print(q)
        is_primitive_root = is_primitive_root and exists_a
        exists_a = False
        iters = 0
    prime = prime and is_primitive_root
    exists_a = False
    while (iters < 288 and not exists_a):
        a2_example = random.randrange(2, N - 1)
        b2_example = random.randrange(2, N - 1)
        if (jacobi(a2_example ** 2 + 4 * b2_example, N) == -1):
            U_example = np.array([[0, 1], [b2_example, a2_example]])
            vec = np.array([0, 1])
            vecres = np.dot(matrix_powmod(U_example, N+1, N), vec)
            vecresq = np.dot(matrix_powmod(U_example, (N+1)//r2, N), vec)
            suits_a = (vecres[0] % N == 0 and vecresq[0] % N != 0)
            exists_a = exists_a or suits_a
        iters += 1
    iters = 0
    if (not exists_a):
        print(r2)
    prime = prime and exists_a
    if not prime:
        return prime
    print("F2-pseudoprime")
    divisors_f4 = []
    f4 = 1
    r4 = N ** 2 + 1
    for i in range(2, 10000000):
        while (r4 % i == 0):
            r4 //= i
            f4 *= i
            if (i != 2):
                divisors_f4.append(i)
    print(divisors_f4)
    suits_a = False
    exists_a = False
    is_primitive_root = True
    for q in divisors_f4:
        while (iters < 288 and not exists_a):
            c_c = random.randrange(2, N - 1)
            d_d = random.randrange(2, N - 1)
            h_h = random.randrange(2, N - 1)
            k_k = random.randrange(2, N - 1)
            if (jacobi(d_d, N) == -1 and jacobi(c_c ** 2 - 16 * d_d, N) == -1):
                p1_p1 = 4 * (h_h ** 2 * 2 + 2 * h_h * k_k * c_c + 2 * k_k * k_k * d_d)
                p2_p2 = (p1_p1 ** 2 - 16 * d_d) // 4
                q_q = (p1_p1 ** 2) // 16 + d_d - (h_h * h_h * c_c + k_k * k_k * c_c * d_d +
                                                        8 * h_h * k_k * d_d)
                U_example = np.array([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1],
                                      [-q_q ** 2, q_q * p1_p1, -p2_p2 - 2 * q_q, p1_p1]])
                vec_example = np.array([0, 1, p1_p1, p1_p1 ** 2 - p2_p2 - 3 * q_q])
                vecres = np.dot(matrix_powmod(U_example, N ** 2 + 1, N), vec_example)
                vecresq = np.dot(matrix_powmod(U_example, (N ** 2 + 1)//q, N), vec_example)
                suits_a = (vecres[0] % N == 0 and vecresq[0] % N != 0)
                exists_a = exists_a or suits_a
            iters += 1
        if (not exists_a):
            print(q)
        is_primitive_root = is_primitive_root and exists_a
        exists_a = False
        iters = 0
    exists_a = False
    prime = prime and is_primitive_root
    while (iters < 288 and not exists_a):
        c_c = random.randrange(2, N - 1)
        d_d = random.randrange(2, N - 1)
        h_h = random.randrange(2, N - 1)
        k_k = random.randrange(2, N - 1)
        if (jacobi(d_d, N) == -1 and jacobi(c_c ** 2 - 16 * d_d, N) == -1):
            p1_p1 = 4 * (h_h ** 2 * 2 + 2 * h_h * k_k * c_c + 2 * k_k * k_k * d_d)
            p2_p2 = (p1_p1 ** 2 - 16 * d_d) // 4
            q_q = (p1_p1 ** 2 + 16 * d_d - 16 * (h_h * h_h * c_c + k_k * k_k * c_c * d_d +
                                                8 * h_h * k_k * d_d)) // 16
            U_example = np.array([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1],
                                [-q_q ** 2, q_q * p1_p1, -p2_p2 - 2 * q_q, p1_p1]])
            vec_example = np.array([0, 1, p1_p1, p1_p1 ** 2 - p2_p2 - 3 * q_q])
            vecres = np.dot(matrix_powmod(U_example, N ** 2 + 1, N), vec_example)
            vecresq = np.dot(matrix_powmod(U_example, (N ** 2 + 1)//r4, N), vec_example)
            suits_a = (vecres[0] % N == 0 and vecresq[0] % N != 0)
            exists_a = exists_a or suits_a
        iters += 1
    prime = prime and exists_a
    if not prime:
        return prime
    print("F4-pseudoprime")
    return prime
