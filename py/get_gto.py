import numpy as np
import ctypes
from pyscf import gto, scf
from utils import shell_slice, shell_pairs, nbas_per_shell, merge, unique_factor

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
mol = gto.M(atom='O 0 0 0; H 0 -0.757 0.587; H 0 0.757 0.587', basis='sto-3g')

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
nshl = len(bas)
natm = len(atm)

'''
print("nao: ", nao)
print("nshl: ", nshl)
print("natm: ", natm)
'''

# Necessary parameters for 2e integer
pairs = shell_pairs(nshl).astype(np.intc)
