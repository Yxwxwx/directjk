import numpy as np
import ctypes
# 导入get_gto.py中的变量
from get_gto import atm, bas, env, nshl, nao, natm

'''
print("atm from get_gto:", atm)
print("bas from get_gto:", bas)
print("env from get_gto:", env)
print("nshl from get_gto:", nshl)
print("nao from get_gto:", nao)
print("natm from get_gto:", natm)

'''



# Load the C library
lib = ctypes.CDLL('../c/libcalc_int1e.so')

# Define the function signature and arguments
lib.calculate_overlap.restype = None
lib.calculate_overlap.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.double, ndim=2),
    np.ctypeslib.ndpointer(dtype=np.intc, ndim=2),
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.intc, ndim=2),
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.double, ndim=1),
    ctypes.c_int
]

lib.calculate_nuce.restype = None
lib.calculate_nuce.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.double, ndim=2),
    np.ctypeslib.ndpointer(dtype=np.intc, ndim=2),
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.intc, ndim=2),
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.double, ndim=1),
    ctypes.c_int
]

lib.calculate_kinetic.restype = None
lib.calculate_kinetic.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.double, ndim=2),
    np.ctypeslib.ndpointer(dtype=np.intc, ndim=2),
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.intc, ndim=2),
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.double, ndim=1),
    ctypes.c_int
]

# Allocate memory for the overlap matrix
overlap_matrix = np.zeros((nao, nao), dtype=np.double)
nuce_matrix =  np.zeros((nao, nao), dtype=np.double)
kinetic_matrix =  np.zeros((nao, nao), dtype=np.double)

# Call the C function
lib.calculate_overlap(
    overlap_matrix,
    atm,
    ctypes.c_int(natm),
    bas,
    ctypes.c_int(nshl),
    env,
    ctypes.c_int(nao)
)
lib.calculate_nuce(
    nuce_matrix,
    atm,
    ctypes.c_int(natm),
    bas,
    ctypes.c_int(nshl),
    env,
    ctypes.c_int(nao)
)
lib.calculate_kinetic(
    kinetic_matrix,
    atm,
    ctypes.c_int(natm),
    bas,
    ctypes.c_int(nshl),
    env,
    ctypes.c_int(nao)
)

kinetic_matrix = kinetic_matrix + kinetic_matrix.T - np.diag(kinetic_matrix.diagonal())
nuce_matrix = nuce_matrix + nuce_matrix.T - np.diag(nuce_matrix.diagonal())
overlap_matrix = overlap_matrix + overlap_matrix.T - np.diag(overlap_matrix.diagonal())
# Print the calculated overlap matrix

'''
print("Overlap matrix:")
print(overlap_matrix)
print("Nuce matrix:")
print(nuce_matrix)
print("Kinetic matrix:")
print(kinetic_matrix)
'''
# Calculate Hamiltion matrix
hamiltion_matrix = np.zeros((nao, nao), dtype=np.double)
hamiltion_matrix = kinetic_matrix + nuce_matrix 


print("Hamiltion matrix: ", hamiltion_matrix) 
