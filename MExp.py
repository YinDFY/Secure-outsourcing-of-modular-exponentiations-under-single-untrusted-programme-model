from outsourcing_utils import RandN
import gmpy2
from gmpy2 import *
import random
import time
gmpy2.get_context().precision = 2048


def Exp(u_list, a_list, p,q):
    blinding_pairs = []


    g1 = gmpy2.mpz(random.randint(2, p - 1))

    for i in range(6):
        g, x, X = RandN(g1, p, q)
        blinding_pairs.append([x, X])

    for kl in range(len(a_list)):
        u1 = u_list[kl]
        a1 = a_list[kl]
        w1 = gmpy2.invert(blinding_pairs[4][1], p)
        w1 = gmpy2.f_mod(gmpy2.mul(w1, u1), p)

        k1t1y1 = gmpy2.f_mod(gmpy2.mul(blinding_pairs[4][0], a1), q)
        k1t1y1 = gmpy2.f_mod(gmpy2.sub(k1t1y1, blinding_pairs[2][0]), q)
        t1y1 = gmpy2.invert(blinding_pairs[0][0], q)
        t1y1 = gmpy2.f_mod(gmpy2.mul(t1y1, k1t1y1), q)

        x1 = gmpy2.f_mod(gmpy2.sub(a1, t1y1), q)

        # Verification Key
        r = random.randint(2,11)

        v2_invert = gmpy2.invert(blinding_pairs[5][1], p)
        w2 = gmpy2.f_mod(gmpy2.mul(v2_invert, u1), p)

        k6r = gmpy2.f_mod(gmpy2.mul(blinding_pairs[5][0], r), q)
        k6ra = gmpy2.f_mod(gmpy2.mul(k6r, a1), q)
        k2t2y2 = gmpy2.f_mod(gmpy2.sub(k6ra, blinding_pairs[3][0]), q)
        k2_invert = gmpy2.invert(blinding_pairs[1][0], q)
        t2y2 = gmpy2.f_mod(gmpy2.mul(k2_invert, k2t2y2), q)

        ar = gmpy2.f_mod(gmpy2.mul(a1, r), q)
        x2 = gmpy2.f_mod(gmpy2.sub(ar, t2y2), q)

        h1w1 = gmpy2.f_mod(gmpy2.mul(blinding_pairs[0][1], w1), p)
        h2w2 = gmpy2.f_mod(gmpy2.mul(blinding_pairs[1][1], w2), p)

        # Cloud compute
        w1x1 = gmpy2.powmod(w1, x1, p)
        h1w1t1y1 = gmpy2.powmod(h1w1, t1y1, p)
        w2x2 = gmpy2.powmod(w2, x2, p)
        h2w2t2y2 = gmpy2.powmod(h2w2, t2y2, p)

        # Verification
        real_answer = gmpy2.f_mod(gmpy2.mul(w1x1, blinding_pairs[2][1]), p)
        real_answer = gmpy2.f_mod(gmpy2.mul(real_answer, h1w1t1y1), p)
        verification = gmpy2.f_mod(gmpy2.mul(w2x2, blinding_pairs[3][1]), p)
        verification = gmpy2.f_mod(gmpy2.mul(verification, h2w2t2y2), p)
        real = gmpy2.powmod(u1, a1, p)
        if real_answer == real:
            print("Real answer is right!")
        #print(real)
        ans = gmpy2.powmod(real_answer, r, p)
        #print(ans)
        #ans = real_answer
        if ans == verification:
            print("Right Answer!")
        else:
            print("Wrong Answer!")


if __name__ == "__main__":
    p = gmpy2.mpz(
        281878283546411134331955053915719078231298790298153670604738367255324835663025510192086892441535016481116338563638344784808982357501005264068784771242695427024340056890189080913241880334753875579302434068869289045538220269912981798151095757328017107872133600753431315640237281631536195165981838724850783448162913)
    q = gmpy2.mpz(
        179769313486231590772930519078902473361797697894230657273430081157732675805500963132708477322407536021120113879871393357658789768814416622492847430639474124377767893424865485276302219601246094119453082952085005768838150682342462881473913110540827237163350510684586298239947245938479716304835356329624224137859)
    u1 = gmpy2.mpz(
        93335961247381333622948086567469035919897422875673769285617485691513550631500202812931576004542073287485376630916307088361801383038384057561565819161703979573113800759009783073515401144515476530669485458019737020553072834988501503944918367822567461547462476061086861281813129544650265044766153907288416469786580)
    a1 = gmpy2.mpz(
        123754834094517264786251613465163955329118838291911092349429048166615726908463192259531999898195839643729638854959475797185234487594815622208075506633961477446217636996304078545523387538944174707028468897255828437774652167526383541712411125502362345275266360538349268189936548024732341143025040327088953139021)

    u = []
    u1 = mpz(pow(u1,mpz(div(p,q-1)),p))
    u.append(u1)
    a = []
    a.append(a1)

    Exp(u,a,p,q)
"""
r = 9
real:65198125551157851040047117901326940796668407819971663731705946559112112925840067372641764076162950836695193576852840348314062795821189854530698548759483882638886484278911895233062530758356385148002681132354806004168819764593678029393245165335116210157147502289690981173259310971850928721965837170448231042482018

verify:55543189065326500947025307327243992407738722204089820036392224037408874948993521145941883108762546040973541817168204711052714133107050541270420514374401344573165597870031832033436793146631392280849478196126594578834355872718368111345805382582120902544723223480539944919452292732705655741147672230574030483277920
"""
