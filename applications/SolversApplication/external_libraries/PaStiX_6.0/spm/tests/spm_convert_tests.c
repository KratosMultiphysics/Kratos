/**
 *
 * @file spm_convert_tests.c
 *
 * @copyright 2011-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * Test and validate the spmConvert routine.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 **/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <spm_tests.h>

#define PRINT_RES(_ret_)                        \
    if(_ret_ == -1) {                           \
        printf("UNDEFINED\n");                  \
    }                                           \
    else if(_ret_ > 0) {                        \
        printf("FAILED(%d)\n", _ret_);          \
        err++;                                  \
    }                                           \
    else {                                      \
        printf("SUCCESS\n");                    \
    }

char* fltnames[] = { "Pattern", "", "Float", "Double", "Complex32", "Complex64" };
char* mtxnames[] = { "General", "Symmetric", "Hermitian" };

int spmComp( const spmatrix_t *spm1,
             const spmatrix_t *spm2 )
{
    spm_int_t *colptr1, *colptr2;
    spm_int_t *rowptr1, *rowptr2;
    int       *valptr1, *valptr2;
    spm_int_t  i;

    if ( spm1->fmttype != SpmCSC ) {
        fprintf(stderr, "Function made to compare only two SPM matrices in CSC format\n");
        return -1;
    }

    if ((spm1->mtxtype != spm2->mtxtype) ||
        (spm1->flttype != spm2->flttype) ||
        (spm1->fmttype != spm2->fmttype) ||
        (spm1->gN      != spm2->gN     ) ||
        (spm1->n       != spm2->n      ) ||
        (spm1->gnnz    != spm2->gnnz   ) ||
        (spm1->nnz     != spm2->nnz    ) ||
        (spm1->dof     != spm2->dof    ) ||
        (spm1->gNexp   != spm2->gNexp  ) ||
        (spm1->nexp    != spm2->nexp   ) ||
        (spm1->gnnzexp != spm2->gnnzexp) ||
        (spm1->nnzexp  != spm2->nnzexp ) ||
        (spm1->layout  != spm2->layout ))
    {
        return 1;
    }

    colptr1 = spm1->colptr;
    colptr2 = spm2->colptr;
    for (i=0; i<=spm1->n; i++, colptr1++, colptr2++) {
        if (*colptr1 != *colptr2 ) {
            return 2;
        }
    }

    rowptr1 = spm1->rowptr;
    rowptr2 = spm2->rowptr;
    for (i=0; i<spm1->nnz; i++, rowptr1++, rowptr2++) {
        if (*rowptr1 != *rowptr2 ) {
            return 3;
        }
    }

    /* Check values */
    if (spm1->values != NULL) {
        spm_int_t size = spm1->nnzexp * (spm_size_of( spm1->flttype ) / sizeof(int));
        valptr1 = (int*)(spm1->values);
        valptr2 = (int*)(spm2->values);
        for (i=0; i<size; i++, valptr1++, valptr2++) {
            if (*valptr1 != *valptr2) {
                z_spmPrintElt( stderr, i, i, *valptr1 );
                z_spmPrintElt( stderr, i, i, *valptr2 );
                return 4;
            }
        }
    }

    return 0;
}

int main (int argc, char **argv)
{
    char *filename;
    spmatrix_t  spm, *spm2;
    spm_driver_t driver;
    spm_mtxtype_t mtxtype;
    int baseval;
    int ret = SPM_SUCCESS;
    int err = 0;
    FILE *f;
    int rc;

    spmGetOptions( argc, argv,
                   &driver, &filename );

    rc = spmReadDriver( driver, filename, &spm );
    free(filename);

    if ( rc != SPM_SUCCESS ) {
        fprintf(stderr, "ERROR: Could not read the file, stop the test !!!\n");
        return EXIT_FAILURE;
    }

    printf(" -- SPM Conversion Test --\n");
    spmConvert(SpmCSC, &spm);

    printf(" Datatype: %s\n", fltnames[spm.flttype] );
    for( baseval=0; baseval<2; baseval++ )
    {
        printf(" Baseval : %d\n", baseval );
        spmBase( &spm, baseval );

        /**
         * Backup the spm
         */
        spm2 = spmCopy( &spm );

        for( mtxtype=SpmGeneral; mtxtype<=SpmHermitian; mtxtype++ )
        {
            if ( (mtxtype == SpmHermitian) &&
                 ((spm.flttype != SpmComplex64) && (spm.flttype != SpmComplex32)) )
            {
                continue;
            }
            spm.mtxtype  = mtxtype;
            spm2->mtxtype = mtxtype;

            printf("   Matrix type : %s\n", mtxnames[mtxtype - SpmGeneral] );

            /**
             * Test cycle CSC -> CSR -> IJV -> CSC
             */
            rc = asprintf( &filename, "convert_b%d_%s_CSC_cycle1.dat",
                           baseval, mtxnames[mtxtype - SpmGeneral] );
            if ( (f = fopen( filename, "w" )) == NULL ) {
                perror("spm_convert_test:cycle1:csc");
                return EXIT_FAILURE;
            }
            spmPrint( &spm, f );
            fclose(f); free(filename);

            printf("   -- Test Conversion CSC -> CSR: ");
            ret = spmConvert( SpmCSR, &spm );
            ret = (ret != SPM_SUCCESS) || (spm.fmttype != SpmCSR );
            PRINT_RES(ret);

            rc = asprintf( &filename, "convert_b%d_%s_CSR_cycle1.dat",
                           baseval, mtxnames[mtxtype - SpmGeneral] );
            if ( (f = fopen( filename, "w" )) == NULL ) {
                perror("spm_convert_test:cycle1:csr");
                return EXIT_FAILURE;
            }
            spmPrint( &spm, f );
            fclose(f); free(filename);

            printf("   -- Test Conversion CSR -> IJV: ");
            ret = spmConvert( SpmIJV, &spm );
            ret = (ret != SPM_SUCCESS) || (spm.fmttype != SpmIJV );
            PRINT_RES(ret);

            rc = asprintf( &filename, "convert_b%d_%s_IJV_cycle1.dat",
                           baseval, mtxnames[mtxtype - SpmGeneral] );
            if ( (f = fopen( filename, "w" )) == NULL ) {
                perror("spm_convert_test:cycle1:ijv");
                return EXIT_FAILURE;
            }
            spmPrint( &spm, f );
            fclose(f); free(filename);

            printf("   -- Test Conversion IJV -> CSC: ");
            ret = spmConvert( SpmCSC, &spm );
            ret = (ret != SPM_SUCCESS) || (spm.fmttype != SpmCSC );
            PRINT_RES(ret);

            /**
             * Check that we came back to the initial state.
             * Do not check if Symmetric or Hermitian due to transposition made
             * in the function.
             */
            if (mtxtype == SpmGeneral) {
                printf("   -- Check the spm after cycle : ");
                ret = spmComp( spm2, &spm );
                PRINT_RES(ret);
            }

            rc = asprintf( &filename, "convert_b%d_%s_CSC_cycle2.dat",
                           baseval, mtxnames[mtxtype - SpmGeneral] );
            if ( (f = fopen( filename, "w" )) == NULL ) {
                perror("spm_convert_test:cycle2:csc");
                return EXIT_FAILURE;
            }
            spmPrint( &spm, f );
            fclose(f); free(filename);

            /**
             * Test second cycle CSC -> IJV -> CSR -> CSC
             */
            printf("   -- Test Conversion CSC -> IJV: ");
            ret = spmConvert( SpmIJV, &spm );
            ret = (ret != SPM_SUCCESS) || (spm.fmttype != SpmIJV );
            PRINT_RES(ret);

            rc = asprintf( &filename, "convert_b%d_%s_IJV_cycle2.dat",
                           baseval, mtxnames[mtxtype - SpmGeneral] );
            if ( (f = fopen( filename, "w" )) == NULL ) {
                perror("spm_convert_test:cycle2:ijv");
                return EXIT_FAILURE;
            }
            spmPrint( &spm, f );
            fclose(f); free(filename);

            printf("   -- Test Conversion IJV -> CSR: ");
            ret = spmConvert( SpmCSR, &spm );
            ret = (ret != SPM_SUCCESS) || (spm.fmttype != SpmCSR );
            PRINT_RES(ret);

            rc = asprintf( &filename, "convert_b%d_%s_CSR_cycle2.dat",
                           baseval, mtxnames[mtxtype - SpmGeneral] );
            if ( (f = fopen( filename, "w" )) == NULL ) {
                perror("spm_convert_test:cycle2:csr");
                return EXIT_FAILURE;
            }
            spmPrint( &spm, f );
            fclose(f); free(filename);

            printf("   -- Test Conversion CSR -> CSC: ");
            ret = spmConvert( SpmCSC, &spm );
            ret = (ret != SPM_SUCCESS) || (spm.fmttype != SpmCSC );
            PRINT_RES(ret);

            rc = asprintf( &filename, "convert_b%d_%s_CSC_end.dat",
                           baseval, mtxnames[mtxtype - SpmGeneral] );
            if ( (f = fopen( filename, "w" )) == NULL ) {
                perror("spm_convert_test:end");
                return EXIT_FAILURE;
            }
            spmPrint( &spm, f );
            fclose(f); free(filename);

            /* Check that we came back to the initial state */
            printf("   -- Check the spm after cycle : ");
            ret = spmComp( spm2, &spm );
            PRINT_RES(ret);
        }
        printf("\n");
        spmExit( spm2 );
        free( spm2 );
    }
    spmExit( &spm  );

    if( err == 0 ) {
        printf(" -- All tests PASSED --\n");
        return EXIT_SUCCESS;
    }
    else
    {
        printf(" -- %d tests FAILED --\n", err);
        return EXIT_FAILURE;
    }

    (void)rc;
}
