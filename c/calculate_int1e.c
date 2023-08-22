#include <stdio.h>
#include <stdlib.h>
#include "cint.h"

void calculate_overlap(double *overlap_matrix, int *atm, int natm, int *bas, int nshl, double *env, int nao){
    int shls[2];
    int dimi = 0, dimj = 0 ;

    //overlap_matrix = malloc(sizeof(double) * nao * nao);
    for(int ipr = 0; ipr < nshl; ipr++)
    {
        shls[0] = ipr;
        dimi = CINTcgto_spheric(ipr, bas);
        int ii0 = 0;
        for (int i = 0; i < ipr; i++)
        {

            ii0 += CINTcgto_spheric(i, bas);
            
        }

        
        for (int jpr = 0; jpr < ipr + 1; jpr++)
        {
            shls[1] = jpr;
            dimj = CINTcgto_spheric(jpr, bas);
            int ij0 = 0;
            for (int j = 0; j < jpr; j++)
            {
                
                ij0 += CINTcgto_spheric(j, bas);
                
            }

                        
            double *buf;
            buf = malloc(sizeof(double) * dimi * dimj);

            cint1e_ovlp_sph(buf, shls, atm, natm, bas, nshl, env);

            for (int ii = 0; ii < dimi; ii++)
            {
                for (int jj = 0; jj < dimj; jj++)
                {
                    overlap_matrix[(ii0 + ii) * nao +  (ij0 + jj)] = buf[ii * dimj + jj];
                    /*
                    printf("ipr: %d\n", ipr);  // Output the value of ii0
                    printf("jpr: %d\n", jpr);  // Output the value of ii0
                    printf("ii0: %d\n", ii0);  // Output the value of ii0
                    printf("ij0: %d\n", ij0);  // Output the value of ii0
                    printf("dimi: %d\n", dimi);  // Output the value of ii0
                    printf("dimj: %d\n", dimj);  // Output the value of ii0
                    printf("ii: %d\n", ii);  // Output the value of ii0
                    printf("jj: %d\n", jj);  // Output the value of ii0
                    printf("buf[ii * dimj + jj]: %f\n", buf[ii * dimj + jj]);
                    printf("\n");
                    */
                }
                
            }
            
            free(buf);

        }
        
    }
    
}


void calculate_nuce(double *nuce_matrix, int *atm, int natm, int *bas, int nshl, double *env, int nao){
    int shls[2];
    int dimi = 0, dimj = 0 ;

    //overlap_matrix = malloc(sizeof(double) * nao * nao);
    for(int ipr = 0; ipr < nshl; ipr++)
    {
        shls[0] = ipr;
        dimi = CINTcgto_spheric(ipr, bas);
        int ii0 = 0;
        for (int i = 0; i < ipr; i++)
        {

            ii0 += CINTcgto_spheric(i, bas);
            
        }

        
        for (int jpr = 0; jpr < ipr + 1; jpr++)
        {
            shls[1] = jpr;
            dimj = CINTcgto_spheric(jpr, bas);
            int ij0 = 0;
            for (int j = 0; j < jpr; j++)
            {
                
                ij0 += CINTcgto_spheric(j, bas);
                
            }

                        
            double *buf;
            buf = malloc(sizeof(double) * dimi * dimj);

            cint1e_nuc_sph(buf, shls, atm, natm, bas, nshl, env);

            for (int ii = 0; ii < dimi; ii++)
            {
                for (int jj = 0; jj < dimj; jj++)
                {
                    nuce_matrix[(ii0 + ii) * nao +  (ij0 + jj)] = buf[ii * dimj + jj];
                    /*
                    printf("ipr: %d\n", ipr);  // Output the value of ii0
                    printf("jpr: %d\n", jpr);  // Output the value of ii0
                    printf("ii0: %d\n", ii0);  // Output the value of ii0
                    printf("ij0: %d\n", ij0);  // Output the value of ii0
                    printf("dimi: %d\n", dimi);  // Output the value of ii0
                    printf("dimj: %d\n", dimj);  // Output the value of ii0
                    printf("ii: %d\n", ii);  // Output the value of ii0
                    printf("jj: %d\n", jj);  // Output the value of ii0
                    printf("buf[ii * dimj + jj]: %f\n", buf[ii * dimj + jj]);
                    printf("\n");
                    */
                }
                
            }
            
            free(buf);

        }
        
    }
    
}


void calculate_kinetic(double *kinetic_matrix, int *atm, int natm, int *bas, int nshl, double *env, int nao){
    int shls[2];
    int dimi = 0, dimj = 0 ;

    //overlap_matrix = malloc(sizeof(double) * nao * nao);
    for(int ipr = 0; ipr < nshl; ipr++)
    {
        shls[0] = ipr;
        dimi = CINTcgto_spheric(ipr, bas);
        int ii0 = 0;
        for (int i = 0; i < ipr; i++)
        {

            ii0 += CINTcgto_spheric(i, bas);
            
        }

        
        for (int jpr = 0; jpr < ipr + 1; jpr++)
        {
            shls[1] = jpr;
            dimj = CINTcgto_spheric(jpr, bas);
            int ij0 = 0;
            for (int j = 0; j < jpr; j++)
            {
                
                ij0 += CINTcgto_spheric(j, bas);
                
            }

                        
            double *buf;
            buf = malloc(sizeof(double) * dimi * dimj);

            cint1e_kin_sph(buf, shls, atm, natm, bas, nshl, env);

            for (int ii = 0; ii < dimi; ii++)
            {
                for (int jj = 0; jj < dimj; jj++)
                {
                    kinetic_matrix[(ii0 + ii) * nao +  (ij0 + jj)] = buf[ii * dimj + jj];
                    /*
                    printf("ipr: %d\n", ipr);  // Output the value of ii0
                    printf("jpr: %d\n", jpr);  // Output the value of ii0
                    printf("ii0: %d\n", ii0);  // Output the value of ii0
                    printf("ij0: %d\n", ij0);  // Output the value of ii0
                    printf("dimi: %d\n", dimi);  // Output the value of ii0
                    printf("dimj: %d\n", dimj);  // Output the value of ii0
                    printf("ii: %d\n", ii);  // Output the value of ii0
                    printf("jj: %d\n", jj);  // Output the value of ii0
                    printf("buf[ii * dimj + jj]: %f\n", buf[ii * dimj + jj]);
                    printf("\n");
                    */
                }
                
            }

            free(buf);
            
        }
        
    }
    
}