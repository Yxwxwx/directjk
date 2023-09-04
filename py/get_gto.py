import numpy as np
import ctypes
from pyscf import gto
from utils import shell_pairs, get_nao


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

# Create a molecular object for H2O molecule
mol = gto.M(atom='''
            O  0.0  0.0  0.0
            O  0.0  0.0  1.5
            H  1.0  0.0  0.0
            H  0.0  0.7  1.0
            ''', basis='6-31g**')

# necessary parameters for Libcint
atm = mol._atm.astype(np.intc)
bas = mol._bas.astype(np.intc)
env = mol._env.astype(np.double)

'''
print("atm: ", atm)
print("bas: ", bas)
print("env: ", env)
'''

# Extract the necessary parameters from the molecular object
nao = mol.nao_nr().astype(np.intc)
#nao_yx = get_nao(bas)
#assert nao_yx == nao


nshl = len(bas)
natm = len(atm)

'''
print("nao: ", nao)
print("nshl: ", nshl)
print("natm: ", natm)
'''

# Necessary parameters for 2e integer
pairs = shell_pairs(nshl).astype(np.intc)

# ==>Get informt <==
# get atom coordinate
A_t = mol.atom_coords()
# get atom charge
Z_A = mol.atom_charges()
# get number of electron
ne = sum(mol.nelec)
# get the number of alpha orb
nalpha = mol.nelec[0]
# get the number of beta orb
nbeta = mol.nelec[1]
# get the number of doubly occupied orbitals
ndocc = min(nalpha, nbeta)
# get the number of single occupied orbitals
nsocc = abs(nalpha - nbeta)
# print
print("Number of AO: ", nao)
print("Number of alpha electrons: ", nalpha)
print("Number of beta electrons: ", nbeta)
print("Number of electrons: ", ne)
print("Number of doubly occupied orbitals: ", ndocc)
print("Number of single occupied orbitals: ", nsocc)

# get r_AB martix
r_AB = np.empty((mol.natm, mol.natm))
for A in range(mol.natm):
    for B in range(mol.natm):
        if A != B:
            r_AB[A, B] = np.linalg.norm(A_t[A] - A_t[B])
        else:
            r_AB[A, B] = np.infty

E_nuc = 0.5 * np.einsum("A, B, AB ->", Z_A, Z_A, 1 / r_AB)
