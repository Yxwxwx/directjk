#include <stdio.h>
#include <stdlib.h>
#include "cint.h"

double* unique_factor(int M, int N, int P, int Q);

void  calculate_j(double* j_matrix, int *atm, int natm, int *bas, int nshl, double *env, int nao, int* pairs, double * dm)
{
    int shls[4];
    int* shl_dim = malloc(sizeof(int)*nshl);
    int* shl_slices = malloc(sizeof(int)*nshl);
    int num_pairs = nshl * (nshl + 1) / 2;
    
    int dimi, dimj, dimk, diml;
    int ii0, ij0, ik0, il0;


    for (int i = 0; i < nshl; i++)
    {
        int n = 0;
        n = CINTcgto_spheric(i, bas);
        shl_dim[i] = n;

        int value = 0;
        for (int j = 0; j < i; j++)
            {
                
                value += CINTcgto_spheric(j, bas);
               
            }

        shl_slices[i] = value; 
       
    }
    
    
    CINTOpt *opt = NULL;
    cint2e_sph_optimizer(&opt, atm, natm, bas, nshl, env);
    

    for(int ipr = 0; ipr < num_pairs; ipr++)
    {
        int i = pairs[0 * num_pairs + ipr];
        int j = pairs[1 * num_pairs + ipr];

        shls[0] = i;
        shls[1] = j;

        dimi = shl_dim[i];
        ii0 = shl_slices[i];

        dimj = shl_dim[j];
        ij0 = shl_slices[j];

        for (int jpr = 0; jpr < ipr + 1; jpr++)
        {
            int k = pairs[0 * num_pairs + jpr];
            int l = pairs[1 * num_pairs + jpr];

            shls[2] = k;
            shls[3] = l;

            dimk = shl_dim[k];
            ik0 = shl_slices[k];

            diml = shl_dim[l];
            il0 = shl_slices[l];

            double *buf;
            buf = malloc(sizeof(double) * dimi * dimj * dimk * diml);
            cint2e_sph(buf, shls, atm, natm, bas, nshl, env, opt);          


            double* factors = unique_factor(i, j, k, l);
            
            for (int ii = 0; ii < dimi ; ii++)
            {
                for (int ij = 0; ij < dimj; ij++)
                {

                    double j_ij = 0;

                    for (int ik = 0; ik < dimk; ik++)
                    {
                        for (int il = 0; il < diml; il++)
                        {
                            double fac1_ij = factors[1] * dm[(ii + ii0) * nao + (ij + ij0)];
                            
                            //double eri = buf[ii * dimj * dimk * diml + ij * dimk * diml + ik * diml + il];

                            double eri = buf[il * dimi * dimj * dimk + ik * dimi * dimj + ij * dimi + ii];
                            
                            j_ij += eri * dm[(ik + ik0) * nao + (il + il0)];

                            j_matrix[(ik0 + ik) * nao + (il0 + il)] += eri * fac1_ij;
                       
                       }
                        
                    }
                    
                    j_matrix[(ii0 + ii) * nao + (ij0 + ij)] += factors[0] * j_ij;
                    
                }
                
            }
            
            free(buf);
            free(factors);

        }
        
        
    }

    CINTdel_optimizer(&opt);
/*

    for (int  i = 0; i < nshl; i++)
    {
        printf("shl_dim: %d\n ", shl_dim[i]);
        printf("shl_slices: %d\n", shl_slices[i]);
    }

*/
    

    free(shl_dim);
    free(shl_slices);
    
    /*
    for (int i = 0; i < nao; i++) {
		printf("j_matrix[%d]: ", i);
		for (int j = 0; j < nao; j++) {
			printf("%f ", j_matrix[i * nao + j]);
		}
		printf("\n");
	}
    */


}

void calculate_k(double* k_matrix, int *atm, int natm, int *bas, int nshl, double *env, int nao, int* pairs, double * dm)
{
    int shls[4];
    int* shl_dim = malloc(sizeof(int)*nshl);
    int* shl_slices = malloc(sizeof(int)*nshl);
    int num_pairs = nshl * (nshl + 1) / 2;
    
    int dimi, dimj, dimk, diml;
    int ii0, ij0, ik0, il0;


    for (int i = 0; i < nshl; i++)
    {
        int n = 0;
        n = CINTcgto_spheric(i, bas);
        shl_dim[i] = n;

        int value = 0;
        for (int j = 0; j < i; j++)
            {
                
                value += CINTcgto_spheric(j, bas);
               
            }

        shl_slices[i] = value; 
       
    }
    
    
    CINTOpt *opt = NULL;
    cint2e_sph_optimizer(&opt, atm, natm, bas, nshl, env);

    for(int ipr = 0; ipr < num_pairs; ipr++)
    {
        int i = pairs[0 * num_pairs + ipr];
        int j = pairs[1 * num_pairs + ipr];

        shls[0] = i;
        shls[1] = j;

        dimi = shl_dim[i];
        ii0 = shl_slices[i];

        dimj = shl_dim[j];
        ij0 = shl_slices[j];

        for (int jpr = 0; jpr < ipr + 1; jpr++)
        {
            int k = pairs[0 * num_pairs + jpr];
            int l = pairs[1 * num_pairs + jpr];

            shls[2] = k;
            shls[3] = l;

            dimk = shl_dim[k];
            ik0 = shl_slices[k];

            diml = shl_dim[l];
            il0 = shl_slices[l];

            double *buf;
            buf = malloc(sizeof(double) * dimi * dimj * dimk * diml);
            cint2e_sph(buf, shls, atm, natm, bas, nshl, env, opt);          


            double* factors = unique_factor(i, j, k, l);

            for (int ii = 0; ii < dimi; ii++)
            {
                for (int ij = 0; ij < dimj; ij++)
                {

                    for (int ik = 0; ik < dimk; ik++)
                    {
                        double fac4_jk = factors[4] * dm[(ij0 + ij) * nao + (ik0 + ik)];
                        double fac5_ik = factors[5] * dm[(ii0 + ii) * nao + (ik0 + ik)];

                        double k_jk = 0;
                        double k_ik = 0;

                        for (int il = 0; il < diml; il++)
                        {
                            double eri = buf[il * dimi * dimj * dimk + ik * dimi * dimj + ij * dimi + ii];
                            k_ik += eri * dm[(ij0 + ij) * nao + (il0 + il)];
                            k_jk += eri * dm[(ii0 + ii) * nao + (il0 + il)];
                            
                            k_matrix[(ii0 + ii) * nao + (il0 + il)] += eri * fac4_jk;
                            k_matrix[(ij0 + ij) * nao + (il0 + il)] += eri * fac5_ik;
                                                    
                            
                        }
                        
                        k_matrix[(ii0 + ii) * nao + (ik0 + ik)] += factors[2] * k_ik;
                        k_matrix[(ij0 + ij) * nao + (ik0 + ik)] += factors[3] * k_jk;

                    }
                    
                }
                
            }

            free(buf);
            free(factors);
            
        }
    }

    CINTdel_optimizer(&opt);
    free(shl_dim);
    free(shl_slices);

}

void calculate_g(double* g_matrix, int *atm, int natm, int *bas, int nshl, double *env, int nao, int* pairs, double * dm)
{
    int shls[4];
    int* shl_dim = malloc(sizeof(int)*nshl);
    int* shl_slices = malloc(sizeof(int)*nshl);
    int num_pairs = nshl * (nshl + 1) / 2;
    
    int dimi, dimj, dimk, diml;
    int ii0, ij0, ik0, il0;
    double * j_matrix;
    double* k_matrix;
    j_matrix = malloc(sizeof(double) * nao * nao);
    k_matrix = malloc(sizeof(double) * nao * nao);

    for (int i = 0; i < nshl; i++)
    {
        int n = 0;
        n = CINTcgto_spheric(i, bas);
        shl_dim[i] = n;

        int value = 0;
        for (int j = 0; j < i; j++)
            {
                
                value += CINTcgto_spheric(j, bas);
               
            }

        shl_slices[i] = value; 
       
    }
    
    
    CINTOpt *opt = NULL;
    cint2e_sph_optimizer(&opt, atm, natm, bas, nshl, env);
    

    for(int ipr = 0; ipr < num_pairs; ipr++)
    {
        int i = pairs[0 * num_pairs + ipr];
        int j = pairs[1 * num_pairs + ipr];

        shls[0] = i;
        shls[1] = j;

        dimi = shl_dim[i];
        ii0 = shl_slices[i];

        dimj = shl_dim[j];
        ij0 = shl_slices[j];

        for (int jpr = 0; jpr < ipr + 1; jpr++)
        {
            int k = pairs[0 * num_pairs + jpr];
            int l = pairs[1 * num_pairs + jpr];

            shls[2] = k;
            shls[3] = l;

            dimk = shl_dim[k];
            ik0 = shl_slices[k];

            diml = shl_dim[l];
            il0 = shl_slices[l];

            double *buf;
            buf = malloc(sizeof(double) * dimi * dimj * dimk * diml);
            cint2e_sph(buf, shls, atm, natm, bas, nshl, env, opt);          

            double* factors = unique_factor(i, j, k, l);

            for (int ii = 0; ii < dimi; ii++)
            {
                for (int ij = 0; ij < dimj; i++)
                {
                    double j_ij = 0;

                    for (int ik = 0; ik < dimk; ik++)
                    {
                        double fac4_jk = factors[4] * dm[(ij0 + ij) * nao + (ik0 + ik)];
                        double fac5_ik = factors[5] * dm[(ii0 + ii) * nao + (ik0 + ik)];

                        double k_jk = 0;
                        double k_ik = 0;
                        
                        for (int il = 0; il < diml; il++)
                        {
                            double fac1_ij = factors[1] * dm[(ii + ii0) * nao + (ij + ij0)];
                            double eri = buf[il * dimi * dimj * dimk + ik * dimi * dimj + ij * dimi + ii];
                            j_ij += eri * dm[(ik + ik0) * nao + (il + il0)];
                            j_matrix[(ik0 + ik) * nao + (il0 + il)] += eri * fac1_ij;

                            k_ik += eri * dm[(ij0 + ij) * nao + (il0 + il)];
                            k_jk += eri * dm[(ii0 + ii) * nao + (il0 + il)];
                            
                            k_matrix[(ii0 + ii) * nao + (il0 + il)] += eri * fac4_jk;
                            k_matrix[(ij0 + ij) * nao + (il0 + il)] += eri * fac5_ik;

                        }
                        
                        k_matrix[(ii0 + ii) * nao + (ik0 + ik)] += factors[2] * k_ik;
                        k_matrix[(ij0 + ij) * nao + (ik0 + ik)] += factors[3] * k_jk;

                    }
                    
                    j_matrix[(ii0 + ii) * nao + (ij0 + ij)] += factors[0] * j_ij;
                }
                
            }
            
            free(buf);
            free(factors);
            

        }
    }
    
    CINTdel_optimizer(&opt);
    free(shl_dim);
    free(shl_slices);
    
    for (int i = 0; i < nao; i++) {
        for (int j = 0; j < i + 1; j++)
        {
            g_matrix[i * nao + j] += j_matrix[i * nao + j] * 2.0;
            g_matrix[i * nao + j] -= k_matrix[i * nao + j];
        }
    }
    
    free(j_matrix);
    free(k_matrix);
    

}







double* unique_factor(int M, int N, int P, int Q) 
{
	double* factors = (double*)malloc(6 * sizeof(double));
	
	int flag1 = (M == N) ? 0 : 1;
	int flag2 = (P == Q) ? 0 : 1;
	int flag3 = ((M == P) && (N == Q)) ? 0 : 1;
	int flag4 = ((flag1 == 1) && (flag2 == 1)) ? 1 : 0;
	int flag5 = ((flag1 == 1) && (flag3 == 1)) ? 1 : 0;
	int flag6 = ((flag2 == 1) && (flag3 == 1)) ? 1 : 0;
	int flag7 = ((flag4 == 1) && (flag3 == 1)) ? 1 : 0;
	
	factors[0] = 1.0 + flag1 + flag2 + flag4;  // for J_MN
	factors[1] = flag3 + flag5 + flag6 + flag7;  // for J_PQ
	factors[2] = 1.0 + flag3;  // for K_MP
	factors[3] = flag1 + flag5;  // for K_NP
	factors[4] = flag2 + flag6;  // for K_MQ
	factors[5] = flag4 + flag7;  // for K_NQ
	
	return factors;
}

