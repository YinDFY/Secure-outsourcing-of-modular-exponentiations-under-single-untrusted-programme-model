import random
import time
import gmpy2
from gmpy2 import mpz
from Crypto.Util.number import getPrime
gmpy2.get_context().precision = 2048




def calculate_average(numbers):
    if not numbers:
        return 0  # 避免除以零错误

    total = sum(numbers)
    average = total / len(numbers)
    return average

def RandN(g1,p ,M ,flag = False,table_beta = [],table_alpha = [],n = 64,h = 64):
    """
    :param p: a large prime number
    :param M: a large number
    :param n: integer
    :param h: integer  h>=1
    :return: pair(x,X), where g^x = X mod p
    """
    times = gmpy2.mpz(gmpy2.div(p - 1, M))
    g = gmpy2.powmod(g1, times, p)
    #transform data into big number library form
    p = gmpy2.mpz(p)
    M = gmpy2.mpz(M)
    k = gmpy2.mpz(random.randint(1,n))
    h = gmpy2.mpz(h)
    n = gmpy2.mpz(n)
    #times = gmpy2.mpz(gmpy2.div(p-1,M))
    if flag == False:
        for i in range(n):
            alpha = gmpy2.mpz(random.randint(0,M-1))
            beta = gmpy2.powmod(g, alpha, p)
            table_alpha.append(alpha)
            table_beta.append(beta)
        return table_alpha,table_beta
    else:
        #Pair Generation
        x = gmpy2.mpz(0)
        X = gmpy2.mpz(1)
        i = 0
        while i < k:
           s = random.randint(0,n-1)
           xj = gmpy2.mpz(random.randint(1,h-1))
           x = gmpy2.add(x,gmpy2.mul(table_alpha[s],xj))
           x = gmpy2.f_mod(x,M)

           X = gmpy2.mul(X,gmpy2.powmod(table_beta[s],xj,p))
           X = gmpy2.f_mod(X, p)
           i = i + 1
           if gmpy2.f_mod(x,M) == 0:
                x = gmpy2.mpz(0)
                X = gmpy2.mpz(1)
                i = 0

        x = gmpy2.f_mod(x,M)
        X = gmpy2.f_mod(X,p)
        g = gmpy2.mpz(g)
        x = gmpy2.mpz(x)
        X = gmpy2.mpz(X)
        return g,x,X


# generate u
def generate_u(p,number):
    roots = set()  # Use a set to avoid duplicates
    while len(roots) < number:
        u = random.randint(2, p - 1)
        if gmpy2.gcd(u, p) == 1:
            roots.add(u)

    return list(roots)

# generate a
def generate_a(q, number):
    roots = set()  # Use a set to avoid duplicates
    while len(roots) < number:
        u = random.randint(2, q - 1)
        if gmpy2.gcd(u, q) == 1:
            roots.add(u)

    return list(roots)



def generate_large_prime(bit):
    prime_candidate = getPrime(bit)
    return gmpy2.mpz(prime_candidate)

def generate_p_q(bit):
    q = generate_large_prime(bit)  # 生成 512 位的素数 q
    i = gmpy2.mpz(2)
    while True:
        p_candidate = gmpy2.mul(q, i) + 1  # p = 2q + 1
        i = i + 1
        if gmpy2.is_prime(p_candidate):
            return p_candidate, q


if __name__ == "__main__":

    start_time = time.time()
    p ,q = generate_p_q()
    print(p)
    print(q)

    N,t = generate_p_q()
    print(N)
    print(t)
    L  = gmpy2.mpz(gmpy2.mul(p,N))
    #k = random.randint(2,L)
    #kN = gmpy2.mul(k,N)
    number = 10000

    u = generate_u(L, 1)
    a = generate_a(gmpy2.mul(p-1,N-1),1)
    print(u)
    print(a)
    #a = generate_a()
    '''
    for i in range(number):
        y = gmpy2.add(u[i],kN)
    '''
    end_time = time.time()
    run_time = end_time - start_time
    print(f"程序运行时间：{run_time * 1000}秒")
