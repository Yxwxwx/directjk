import numpy as np
import ctypes
# 导入get_gto.py中的变量
from get_gto import atm, bas, env, nshl, nao, natm, pairs

'''
print("atm from get_gto:", atm)
print("bas from get_gto:", bas)
print("env from get_gto:", env)
print("nshl from get_gto:", nshl)
print("nao from get_gto:", nao)
print("natm from get_gto:", natm)

'''
print("pairs: ", pairs)




# Load the C library
lib = ctypes.CDLL('../c/libcalc_int2e.so')

# Define the function signature and arguments
lib.calculate_j.restype = None
lib.calculate_j.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.double, ndim=2),
    np.ctypeslib.ndpointer(dtype=np.intc, ndim=2),
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.intc, ndim=2),
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.double, ndim=1),
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.intc, ndim=2)
]

# Allocate memory for the overlap matrix
j_matrix = np.zeros((nao, nao), dtype=np.double)

# Call the C function
lib.calculate_j(
    j_matrix,
    atm,
    ctypes.c_int(natm),
    bas,
    ctypes.c_int(nshl),
    env,
    ctypes.c_int(nao),
    pairs
)

print("J matrix: ", j_matrix)