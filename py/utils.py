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
        nao = 0
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

