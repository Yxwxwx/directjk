import numpy as np
import scipy 
import ctypes
from get_gto import ndocc, natm, E_nuc, nao
import get_int1e
import get_int2e

def get_fock_matrix(dm):
    return get_int1e.get_hamilton_matrix() + get_int2e.get_g_matrix(dm)

def get_density_matrix(fock, nocc):
    dm = np.zeros((nao, nao))
    S = get_int1e.get_overlap_matrix()
    A = scipy.linalg.fractional_matrix_power(S, -0.5)
    #A = np.linalg.inv(np.linalg.cholesky(S).T)
    F_p = A.T @ fock @ A
    eigs,  coeffsm = np.linalg.eigh(F_p)
    
    c_occ = A @ coeffsm
    c_occ = c_occ[:, :nocc]
    dm =  np.einsum('pi,qi->pq', c_occ, c_occ, optimize=True)
    #dm = c_occ[:, nocc] @ c_occ[:, nocc].T
    return dm
'''
def get_density_matrix(fock, nocc):
    dm = np.zeros((nao, nao))
    S = get_int1e.get_overlap_matrix()
    e, c = scipy.linalg.eigh(fock, S)
    idx = np.argmax(abs(c.real), axis=0)
    c[:,c[idx,np.arange(len(e))].real<0] *= -1
    c_occ = c[:, :nocc]
    dm =  np.einsum('pi,qi->pq', c_occ, c_occ, optimize=True)
    return dm
'''
# ==> Build DIIS Extrapolation Function <==
def diis_xtrap(F_list, DIIS_RESID):
    # Build B matrix
    B_dim = len(F_list) + 1
    B = np.empty((B_dim, B_dim))
    B[-1, :] = -1
    B[:, -1] = -1
    B[-1, -1] = 0
    for i in range(len(F_list)):
        for j in range(len(F_list)):
            B[i, j] = np.einsum('ij,ij->', DIIS_RESID[i], DIIS_RESID[j], optimize=True)

    # Build RHS of Pulay equation
    rhs = np.zeros((B_dim))
    rhs[-1] = -1

    # Solve Pulay equation for c_i's with NumPy
    coeff = np.linalg.solve(B, rhs)

    # Build DIIS Fock matrix
    F_DIIS = np.zeros_like(F_list[0])
    for x in range(coeff.shape[0] - 1):
        F_DIIS += coeff[x] * F_list[x]

    return F_DIIS

def set_looper(max_iter, E_conv, D_conv):
    scf_e = 0.0
    e_old = 0.0
    H = get_int1e.get_hamilton_matrix()
    S = get_int1e.get_overlap_matrix()
    #dm = get_density_matrix(H, ndocc)
    dm = np.eye(nao)
    for scf_iter in range(1, max_iter + 1):
        F = get_fock_matrix(dm)
        #print("F: \n", F, "\n")
        diis_r = F.dot(dm).dot(S) - S.dot(dm).dot(F)
        scf_e = np.einsum('pq,pq->', (F + H), dm, optimize=True) + E_nuc
        dE = scf_e - e_old
        dRMS = 0.5 * np.mean(diis_r ** 2) ** 0.5
        print('SCF Iteration %3d: Energy = %4.16f dE = % 1.5E dRMS = %1.5E' % (scf_iter, scf_e, dE, dRMS))

        if(abs(dE) < E_conv) and (dRMS < D_conv):
            print("SCF convergence! Congrats")
            break
        
        e_old = scf_e
        
        dm = get_density_matrix(F, ndocc)
    
def set_looper_diis(max_iter, E_conv, D_conv):
    scf_e = 0.0
    e_old = 0.0
    H = get_int1e.get_hamilton_matrix()
    S = get_int1e.get_overlap_matrix()
    F_list = []
    R_list = []
    dm = get_density_matrix(H, ndocc)
    #dm = np.eye(nao)
    for scf_iter in range(1, max_iter + 1):
        F = get_fock_matrix(dm)
        #print("F: \n", F, "\n")
        diis_r = F.dot(dm).dot(S) - S.dot(dm).dot(F)
        # Trial & Residual Vector Lists -- one each for alpha & beta
        F_list.append(F)
        R_list.append(diis_r)
        scf_e = np.einsum('pq,pq->', (F + H), dm, optimize=True) + E_nuc
        dE = scf_e - e_old
        dRMS = 0.5 * np.mean(diis_r ** 2) ** 0.5
        print('SCF Iteration %3d: Energy = %4.16f dE = % 1.5E dRMS = %1.5E' % (scf_iter, scf_e, dE, dRMS))

        if(abs(dE) < E_conv) and (dRMS < D_conv):
            print("SCF convergence! Congrats")
            break
        
        e_old = scf_e
        if scf_iter >= 1:
            if scf_iter == 2:
                print("DIIS start!")
            F = diis_xtrap(F_list, R_list)
        
        dm = get_density_matrix(F, ndocc)