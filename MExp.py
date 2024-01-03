from outsourcing_utils import *
import gmpy2
from gmpy2 import *
import random
import time
gmpy2.get_context().precision = 2048

def fast_power(base, exponent):
    result = mpz(1)
    while exponent > 0:
        # 如果指数是奇数，则将结果乘以底数
        if f_mod(exponent,2) == 1:
            result = mul(result,base)
        # 将底数平方
        base = mul(base,base)
        # 将指数右移一位（相当于除以2）
        exponent //= 2
    return result

def Exp(u_list, a_list, p,q,g1,table_alpha,table_beta,r):
    s3_time = time.time()
    blinding_pairs = []

    for i in range(6):
        g, x, X = RandN(g1, p, q,True,table_beta,table_alpha)
        blinding_pairs.append([x, X])

    v1_inverse = invert(blinding_pairs[4][1], p)
    v2_inverse = invert(blinding_pairs[5][1], p)
    a_sum = mpz(0)

    w = []
    w_ = []
    x = []
    x_ = []
    i = 0
    hw_x = blinding_pairs[0][1]
    hw_y = blinding_pairs[1][1]

    for u in u_list:
        wx = f_mod(mul(u, v1_inverse), p)
        w.append(wx)
        hw_x = f_mod(mul(hw_x,wx),p)
        wy = f_mod(mul(u,v2_inverse),p)
        w_.append(wy)
        hw_y = f_mod(mul(hw_y,wy),p)
        a_sum = f_mod(add(a_sum, a_list[i]), q)
        i = i + 1


    k5_a_sum = f_mod(mul(blinding_pairs[4][0],a_sum),q)
    k1t1y1 = f_mod(sub(k5_a_sum,blinding_pairs[2][0]),q)
    k1_inverse = invert(blinding_pairs[0][0],q)
    t1y1 = f_mod(mul(k1t1y1,k1_inverse),q)

    k6r = f_mod(mul(r,blinding_pairs[5][0]),q)
    k6r_a_sum = f_mod(mul(k6r,a_sum),q)
    k2t2y2 = f_mod(sub(k6r_a_sum,blinding_pairs[3][0]),q)
    k2_inverse = invert(blinding_pairs[1][0],q)
    t2y2 = f_mod(mul(k2t2y2,k2_inverse),q)



    for a in a_list:
        x.append(f_mod(sub(a,t1y1),q))
        x_.append(f_mod(sub(f_mod(mul(r,a),q),t2y2),q))

    e3_time = time.time()
    keyGen_time = (e3_time - s3_time)*1000
    #print("KeyGen:",(e_time - s_time)*1000)

    #Compute
    s2_time = time.time()
    TK = blinding_pairs[2][1]
    VK = blinding_pairs[3][1]
    hw_t1y1 = powmod(hw_x,t1y1,p)
    hw_t2y2 = powmod(hw_y,t2y2,p)
    wx = []
    wy = []
    for k in range(len(w)):
        wx.append(powmod(w[k],x[k],p))
        wy.append(powmod(w_[k],x_[k],p))

    e2_time = time.time()
    compute_time = (e2_time - s2_time) * 1000
    #print("Compute time:",(e_time - s_time)*1000)
    #verify
    flag = False
    s1_time = time.time()
    TK = f_mod(mul(TK,hw_t1y1),p)
    VK = f_mod(mul(VK,hw_t2y2),p)
    for tmp in wx:
        TK = f_mod(mul(TK,tmp),p)
    for tmp in wy:
        VK = f_mod(mul(VK,tmp),p)
    if(powmod(TK,r,p) == VK):
        flag = True
    else:
        flag = False
    e1_time = time.time()
    verify_time = (e1_time - s1_time)*1000
    #print("Verification time:",(e1_time - s1_time)*1000)
    return flag,keyGen_time,compute_time,verify_time

if __name__ == "__main__":

    batch = [10000,20000,30000,40000,50000]
    for number in batch:
        setup = []
        key = []
        com = []
        ver = []
        for _ in range(50):
            start_time = time.time()
            p, q = generate_p_q(1024)
            u1 = generate_u(p, number)
            g1 = gmpy2.mpz(random.randint(2, p - 1))
            table_alpha, table_beta = RandN(g1, p, q, False)
            u = []
            for tmp in u1:
                #u1 = mpz(powmod(tmp,mpz(div(p-1,q)),p))
                k = mpz(random.randint(2,q))
                kN = mul(k,q)
                u1 = f_mod(add(tmp,kN),p)
                u.append(u1)
            end_time = time.time()
            setup.append((end_time - start_time)*1000)
            #print("Setup time",(end_time - start_time)*1000)

            #print("p",p)
            r = mpz(random.randint(2,11))
            a1 = generate_a(q, number)


            flag ,keygen,compute,verify= Exp(u,a1,p,q,g1,table_alpha, table_beta,r)
            key.append(keygen)
            com.append(compute)
            ver.append(verify)
            #print(flag)
        print("Group:",number)
        print("setup_ave:",calculate_average(setup))
        print("keygen_ave:",calculate_average(key))
        print("compute_ave:",calculate_average(com))
        print("ver_ave",calculate_average(ver))
