#include <stdio.h>
#include <stdlib.h>
#include "cint.h"

double* unique_factor(int M, int N, int P, int Q);

void  calculate_j(double* j_matrix, int *atm, int natm, int *bas, int nshl, double *env, int nao, int* pairs)
{
    int shls[2];
    //int shl_dim[nshl], shl_slices[nshl];
    int* shl_dim = malloc(sizeof(int)*nshl);
    int* shl_slices = malloc(sizeof(int)*nshl);
    int num_pairs = nshl * (nshl + 1) / 2;
    
    int dimi, dimj, dimk, diml;
    int ii0, ij0, ik0, il0;
    double *buf;

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
        //printf("shl_dim[%d] = %d, shl_slices[%d] = %d\n", i, shl_dim[i], i, shl_slices[i]);
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
            /**/
            printf("i: %d\n", i);
            printf("j: %d\n", j);
            printf("k: %d\n", k);
            printf("l: %d\n", l);
            printf("\n");

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
                            double eri = buf[ii * dimj * dimk * diml + ij * dimk * diml + ik * diml + il];
                            //printf("buf[ii * dimj * dimk * diml + ij * dimk * diml + ik * diml + il]: %f\n", eri);

                            j_ij += eri;

                            j_matrix[(ik0 + ik) * nao + (il0 + il)] += eri * factors[1];
                            printf("ik0 + ik: %d\n", (ik0 + ik));
                            printf("il0 + il: %d\n", (il0 + il));

                            if (il0 + il > ik0 + ik)
                            {
                                printf("FUCK!\n");
                                printf("eir: %f\n", eri);
                            }
                            
                            printf("\n");
                            //printf("j_matrix[(ik0 + ik) * nao + (il0 + il)]: %f\n", eri);
                        }
                        
                    }
                    
                    j_matrix[(ii0 + ii) * nao + (ij0 + ij)] += factors[0] * j_ij;
                    //printf("j_matrix[(ii0 + ii) * nao + (ij0 + ij)]: %f\n", j_matrix[(ii0 + ii) * nao + (ij0 + ij)]);
                    printf("ii0 + ii: %d\n", (ii0 + ii));
                    printf("ij0 + ij: %d\n", (ij0 + ij));
                    
                    if (ij0 + ij > ii0 + ii)
                            {
                                printf("FUCK!\n");
                                //printf("eir: %d\n", eri);
                            }

                    printf("\n");
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
    
    for (int i = 0; i < nao; i++) {
		printf("j_matrix[%d]: ", i);
		for (int j = 0; j < nao; j++) {
			printf("%f ", j_matrix[i * nao + j]);
		}
		printf("\n");
	}

    //printf("num_pairs: %d\n", num_pairs);

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