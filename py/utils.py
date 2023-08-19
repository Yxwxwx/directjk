import numpy as np
import ctypes

# slots of atm
CHARGE_OF       = 0
PTR_COORD       = 1
NUC_MOD_OF      = 2
PTR_ZETA        = 3
PTR_FRAC_CHARGE = 3
RESERVE_ATMLOT1 = 4
RESERVE_ATMLOT2 = 5
ATM_SLOTS       = 6


# slots of bas
ATOM_OF         = 0
ANG_OF          = 1
NPRIM_OF        = 2
NCTR_OF         = 3
KAPPA_OF        = 4
PTR_EXP         = 5
PTR_COEFF       = 6
RESERVE_BASLOT  = 7
BAS_SLOTS       = 8

def get_nao(bas):
    nshl = len(bas)
    for i in range (nshl):
        nao += (bas[i][ANG_OF] * 2 + 1) * bas[i][NCTR_OF]
        return nao
    
def shell_slice(idx, bas):
    n = 0
    for i in range(idx):
        n += nbas_per_shell(bas, i)
    return n

def shell_pairs(nshls):
    pairs = np.zeros((2, nshls*(nshls+1)//2), dtype=int)
    k = 0
    for i in range(nshls):
        for j in range(i + 1):
            pairs[0, k] = i
            pairs[1, k] = j
            k += 1
    return pairs

def nbas_per_shell(bas, idx):
    n = bas[idx][ANG_OF] * 2 + 1
    return n

def merge(a, b, condition):
    return a if condition else b

def unique_factor(M, N, P, Q):
    factors = np.zeros(6, dtype=np.float64)
    flag1 = merge(0, 1, M == N)
    flag2 = merge(0, 1, P == Q)
    flag3 = merge(0, 1, (M == P) and (N == Q))
    flag4 = merge(1, 0, (flag1 == 1) and (flag2 == 1))
    flag5 = merge(1, 0, (flag1 == 1) and (flag3 == 1))
    flag6 = merge(1, 0, (flag2 == 1) and (flag3 == 1))
    flag7 = merge(1, 0, (flag4 == 1) and (flag3 == 1))
    factors[0] = 1.0 + flag1 + flag2 + flag4  # for J_MN
    factors[1] = flag3 + flag5 + flag6 + flag7  # for J_PQ
    factors[2] = 1.0 + flag3  # for K_MP
    factors[3] = flag1 + flag5  # for K_NP
    factors[4] = flag2 + flag6  # for K_MQ
    factors[5] = flag4 + flag7  # for K_NQ
    return factors
