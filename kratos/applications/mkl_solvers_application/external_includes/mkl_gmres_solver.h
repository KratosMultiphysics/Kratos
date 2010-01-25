/*          
* =======================================================================*
* kkkk   kkkk  kkkkkkkkkk   kkkkk    kkkkkkkkkk kkkkkkkkkk kkkkkkkkkK    *
* kkkk  kkkk   kkkk   kkkk  kkkkkk   kkkkkkkkkk kkkkkkkkkk kkkkkkkkkK    *
* kkkkkkkkk    kkkk   kkkk  kkkkkkk     kkkk    kkk    kkk  kkkk         *
* kkkkkkkkk    kkkkkkkkkkk  kkkk kkk	kkkk    kkk    kkk    kkkk       *
* kkkk  kkkk   kkkk  kkkk   kkkk kkkk   kkkk    kkk    kkk      kkkk     *
* kkkk   kkkk  kkkk   kkkk  kkkk  kkkk  kkkk    kkkkkkkkkk  kkkkkkkkkk   *
* kkkk    kkkk kkkk    kkkk kkkk   kkkk kkkk    kkkkkkkkkk  kkkkkkkkkk 	 *
*                                                                        *
* krATos: a fREe opEN sOURce CoDE for mULti-pHysIC aDaptIVe SoLVErS,     *
* aN extEnsIBLe OBjeCt oRiEnTEd SOlutION fOR fInITe ELemEnt fORmULatIONs *
* Copyleft by 2003 ciMNe                                                 *
* Copyleft by 2003 originary authors Copyleft by 2003 your name          *
* This library is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License as         *
* published by the Free Software Foundation; either version 2.1 of       *
* the License, or any later version.                                     *
*                                                                        *
* This library is distributed in the hope that it will be useful, but    *
* WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                   *
* See the GNU Lesser General Public License for more details.            *
*                                                                        *
* You should have received a copy of the GNU Lesser General Public       *
* License along with this library; if not, write to International Centre *
* for Numerical Methods in Engineering (CIMNE),                          *
* Edifici C1 - Campus Nord UPC, Gran Capità s/n, 08034 Barcelona.        *
*                                                                        *
* You can also contact us to the following email address:                *
* kratos@cimne.upc.es                                                    *
* or fax number: +34 93 401 65 17                                        *
*                                                                        *
* Created at Institute for Structural Mechanics                          *
* Ruhr-University Bochum, Germany                                        *
* Last modified by:    $Author: janosch $  				 *
* Date:                $Date: 2009-01-14 09:40:28 $			 *
* Revision:            $Revision: 1.1 $ 				 *
*========================================================================*
* International Center of Numerical Methods in Engineering - CIMNE	 *
* Barcelona - Spain 							 *
*========================================================================*
*/

#if !defined(KRATOS_MKL_GMRES_SOLVER_H_INCLUDED )
#define  KRATOS_MKL_GMRES_SOLVER_H_INCLUDED

// #define BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT

// External includes 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <boost/timer.hpp>
#include <omp.h>

#include "boost/smart_ptr.hpp"
#include "includes/ublas_interface.h"

#include <boost/numeric/bindings/traits/sparse_traits.hpp>

#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>

// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"
                 
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"


namespace ublas = boost::numeric::ublas;

namespace Kratos
{
    template< class TSparseSpaceType, class TDenseSpaceType, 
              class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
        class MKLGMRESSolver : public DirectSolver< TSparseSpaceType, 
              TDenseSpaceType, TReordererType>
    {
        public:
            /**
             * Counted pointer of SuperLUSolver
             */
            typedef boost::shared_ptr<MKLGMRESSolver> Pointer;
            
            typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType; 
            
            typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
            
            typedef typename TSparseSpaceType::VectorType VectorType;
            
            typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
            
            /**
             * Default constructor
             */
            MKLGMRESSolver(){}
            
            /**
             * Destructor
             */
            virtual ~MKLGMRESSolver(){}
            
            /** 
             * Normal solve method.
             * Solves the linear system Ax=b and puts the result on SystemVector& rX.
             * rX is also th initial guess for iterative methods.
             * @param rA. System matrix
             * @param rX. Solution vector.
             * @param rB. Right hand side vector.
             */
            bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
            {
                std::cout << "======= ENTERING GMRES SOLVER =======" << std::endl;
                double start_solver = omp_get_wtime();
                typedef boost::numeric::bindings::traits::sparse_matrix_traits<SparseMatrixType> matraits;
                typedef boost::numeric::bindings::traits::vector_traits<VectorType> mbtraits;
                typedef typename matraits::value_type val_t;
                
                MKL_INT n = matraits::size1 (rA);
                assert (n == matraits::size2 (rA));
                assert (n == mbtraits::size1 (rB));
                assert (n == mbtraits::size1 (rX));
                
                /** 
                 * nonzeros in rA
                 */
                double* a = matraits::value_storage(rA);
                
                /** 
                 * manual index vector generation
                 */
                int *index1_vector = new (std::nothrow) int[rA.index1_data().size()];
                int *index2_vector = new (std::nothrow) int[rA.index2_data().size()];
                std::cout << "Size of the problem: " << n << std::endl;
//                 std::cout << "Size of index1_vector: " << rA.index1_data().size() << std::endl;
                std::cout << "number of nonzeros: " << rA.index2_data().size() << std::endl;
                for(unsigned int i = 0; i < rA.index1_data().size(); i++ )
                {
                    index1_vector[i] = (int)(rA.index1_data()[i])+1;
                } 
                for(unsigned int i = 0; i < rA.index2_data().size(); i++ )
                {
                    index2_vector[i] = (int)(rA.index2_data()[i])+1;
                }
                
                /* RHS and solution vectors. */
                double *b = mbtraits::storage(rB);
//                 KRATOS_WATCH(rB);

                //commented by Riccardo as it is unusd
//                double *x = mbtraits::storage(rX);
                
                
                

    //---------------------------------------------------------------------------
    // Define arrays for the upper triangle of the coefficient matrix
    // Compressed sparse row storage is used for sparse representation
    //---------------------------------------------------------------------------
//                 MKL_INT ia[6]={1,3,6,9,12,14};
//                 MKL_INT ja[13]={    1,        3,
//                     1,   2,        4,
//                     2,   3,        5,
//                     3,   4,   5,
//                     4,   5  };
//                     double A[13]={ 1.0,     -1.0,
//                         -1.0, 1.0,     -1.0,
//                         1.0,-2.0,      1.0,
//                         -1.0, 2.0,-1.0,
//                         -1.0,-3.0 };
    //---------------------------------------------------------------------------
    // Allocate storage for the ?par parameters and the solution/rhs/residual vectors
    //---------------------------------------------------------------------------
                        MKL_INT ipar[128];
                        int size_tmp = fmin(n,150);
                        double tmp[((2*size_tmp+1)*n+size_tmp*size_tmp+9)/2 + 1];
                        double dpar[128];
                        double expected_solution[n];
                        for( int i=0; i<n; i++ )
                            expected_solution[i] = b[i];
//                         std::cout << "expected solution generated" << std::endl;
                        double rhs[n];
                        double computed_solution[n];
                        double residual[n];
                        //---------------------------------------------------------------------------
                        // Some additional variables to use with the RCI (P)FGMRES solver
                        //---------------------------------------------------------------------------*/
                        MKL_INT itercount;
                        MKL_INT RCI_request, i, ivar;
                        double dvar;
                        char cvar;

//                         std::cout << "--------------------------------------------------------" << std::endl;
//                         std::cout << "The FULLY ADVANCED example of usage of RCI FGMRES solver" << std::endl;
//                         std::cout << "   to solve a non-symmetric indefinite non-degenerate" << std::endl;
//                         std::cout << "          algebraic system of linear equations" << std::endl;
//                         std::cout << "--------------------------------------------------------" << std::endl;;
                        //---------------------------------------------------------------------------
                        // Initialize variables and the right hand side through matrix-vector product
                        //---------------------------------------------------------------------------*/
                        ivar=n;
                        cvar='N';
//                         std::cout << "before:" << std::endl;
//                         std::cout << "B: ";
//                         for( int i=0; i<n; i++ )
//                             std::cout << b[i] << "\t";
//                         std::cout << std::endl;
//                         mkl_dcsrgemv(&cvar, &ivar, a, index1_vector, index2_vector, expected_solution, b);
//                         std::cout << "after:" << std::endl;
//                         std::cout << "B: ";
//                         for( int i=0; i<n; i++ )
//                             std::cout << b[i] << "\t";
//                         std::cout << std::endl;
//                         std::cout << "===============================" << std::endl;
                        //---------------------------------------------------------------------------
                        // Save the right-hand side in vector rhs for future use
                        //---------------------------------------------------------------------------*/
                        i=1;
                        dcopy(&ivar, b, &i, rhs, &i);
                        //---------------------------------------------------------------------------
                        // Initialize the initial guess
                        //---------------------------------------------------------------------------*/
                        for(i=0;i<n;i++)
                        {
                            computed_solution[i]=0.0;
                        }
                        //---------------------------------------------------------------------------
                        // Initialize the solver
                        //---------------------------------------------------------------------------*/
//                         std::cout << "initializing the solver...";
                        dfgmres_init(&ivar, computed_solution, b, &RCI_request, ipar, dpar, tmp);
//                         std::cout << " done. Exit with " << RCI_request << std::endl;
                        if (RCI_request!=0) Error( RCI_request );
                        //---------------------------------------------------------------------------
//                         / Set the desired parameters:
//                         * do the restart after 2 iterations
//                         * LOGICAL parameters:
//                         * do not do the stopping test for the maximal number of iterations
//                         * do the Preconditioned iterations of FGMRES method
//                         * DOUBLE PRECISION parameters
//                         * set the relative tolerance to 1.0D-3 instead of default value 1.0D-6 */
//                        /*---------------------------------------------------------------------------*/
                        ipar[14]=2;
                        ipar[7]=1;
                        ipar[8]=1;
                        ipar[9] = 0;
                        ipar[10]=0;
                        dpar[0]=1.0E-3;
//                        /*---------------------------------------------------------------------------
//                         * Check the correctness and consistency of the newly set parameters
//                         *---------------------------------------------------------------------------*/
//                         std::cout << "checking consistency...";
                        dfgmres_check(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp);
//                         std::cout << " done. Exit with " << RCI_request << std::endl;
                        if (RCI_request!=0) Error( RCI_request );
                        
                        /*---------------------------------------------------------------------------
                         * Print the info about the RCI FGMRES method 
                         *---------------------------------------------------------------------------*/
//                         std::cout << "Some info about the current run of RCI FGMRES method:" << std::endl;
//                         if (ipar[7]) 
//                         {
//                            std::cout << "As ipar[7]=" << ipar[7] <<", the automatic test for the maximal number of iterations will be performed" << std::endl;
//                         }
//                         else 
//                         {
//                             std::cout << "As ipar[7]=" << ipar[7] <<", the automatic test for the maximal number of iterations will be skipped" << std::endl;
//                         }
//                         std::cout << "+++" << std::endl;
//                         if (ipar[8]) 
//                         {
//                             std::cout << "As ipar[8]=" << ipar[8] <<", the automatic residual test will be performed" << std::endl;
//                         }
//                         else 
//                         {
//                             std::cout << "As ipar[8]=" << ipar[8] <<", the automatic residual test will be skipped" << std::endl;
//                         }
//                         std::cout << "+++" << std::endl;
//                         if (ipar[9]) 
//                         {
//                             std::cout << "As ipar[9]=" << ipar[9] <<", the user-defined stopping test will be requested via RCI_request=2" << std::endl;
//                         }
//                         else 
//                         {
//                             std::cout << "As ipar[9]=" << ipar[9] <<", the user-defined stopping test will not be requested, thus, RCI_request will not take the value 2" << std::endl;
//                         }
//                         std::cout << "+++" << std::endl;
//                         if (ipar[10]) 
//                         {
//                             std::cout << "As ipar[10]=" << ipar[10] <<", the Preconditioned FGMRES iterations will be performed, thus, the preconditioner action will be requested via RCI_request=3" << std::endl;
//                         }
//                         else 
//                         {
//                             std::cout << "As ipar[10]=" << ipar[10] <<", the Preconditioned FGMRES iterations will not be performed, thus, RCI_request will not take the value" << std::endl;
//                         }
//                         std::cout << "+++" << std::endl;
//                         if (ipar[11]) 
//                         {
//                             std::cout << "As ipar[11]=" << ipar[11] <<", the automatic test for the norm of the next generated vector is not equal to zero up to rounding and computational errors will be performed, thus, RCI_request will not take the value 4" << std::endl;
//                         }
//                         else 
//                         {
//                             std::cout << "As ipar[11]=" << ipar[11] <<", the automatic test for the norm of the next generated vector is not equal to zero up to rounding and computational errors will be skipped, thus, the user-defined test will be requested via RCI_request=4" << std::endl;
//                         }
//                         std::cout << "+++" << std::endl;
                        //---------------------------------------------------------------------------
//                         * Compute the solution by RCI (P)FGMRES solver with preconditioning
//                         * Reverse Communication starts here
                        //---------------------------------------------------------------------------
//                         std::cout << "computing solution...";
//                         std::cout << "before: " << std::endl;
//                         std::cout << "b: ";
//                         for( int i=0; i<n; i++ )
//                             std::cout << b[i] << "\t";
//                         std::cout << std::endl;
//                         std::cout << "computed_solution: ";
//                         for( int i=0; i<n; i++ )
//                             std::cout << computed_solution[i] << "\t";
//                         std::cout << std::endl;
                        
                        ONE:  dfgmres(&ivar, computed_solution, b, &RCI_request, ipar, dpar, tmp);
//                         std::cout << " done. Exit with " << RCI_request << std::endl;
                        
//                         std::cout << "after:" << std::endl;
//                         std::cout << "b: ";
//                         for( int i=0; i<n; i++ )
//                             std::cout << b[i] << "\t";
//                         std::cout << std::endl;
//                         std::cout << "computed_solution: ";
//                         for( int i=0; i<n; i++ )
//                             std::cout << computed_solution[i] << "\t";
//                         std::cout << std::endl;

                        
                        
                        
                        
                        
                        /*---------------------------------------------------------------------------
                         * If RCI_request=0, then the solution was found with the required precision
                         *---------------------------------------------------------------------------*/
                        if (RCI_request==0)
                        {
//                             std::cout << "COMPLETE!!!" << std::endl;
                            goto COMPLETE;
                        }
                        /*---------------------------------------------------------------------------
                         * If RCI_request=1, then compute the vector A*tmp[ipar[21]-1] 
                         * and put the result in vector tmp[ipar[22]-1]
                         *---------------------------------------------------------------------------
                         * NOTE that ipar[21] and ipar[22] contain FORTRAN style addresses,
                         * therefore, in C code it is required to subtract 1 from them to get C style
                         * addresses
                         *---------------------------------------------------------------------------*/
                        if (RCI_request==1)
                        {
//                             std::cout << "computing A*tmp...";
                            mkl_dcsrgemv(&cvar, &ivar, a, index1_vector, index2_vector, &tmp[ipar[21]-1], &tmp[ipar[22]-1]);
//                             std::cout << " recomputing... ";
                            goto ONE;
                        }
                        
//                        /*---------------------------------------------------------------------------
//                        /* If RCI_request=2, then do the user-defined stopping test
//                        /* The residual stopping test for the computed solution is performed here
//                        /*---------------------------------------------------------------------------
//                        /* NOTE: from this point vector b[N] is no longer containing the right-hand
//                        /* side of the problem! It contains the current FGMRES approximation to the
//                        /* solution. If you need to keep the right-hand side, save it in some other
//                        /* vector before the call to dfgmres routine. Here we saved it in vector
//                        /* rhs[N]. The vector b is used instead of rhs to preserve the
//                        /* original right-hand side of the problem and guarantee the proper
//                        /* restart of FGMRES method. Vector b will be altered when computing the
//                        /* residual stopping criterion!
//                        /*---------------------------------------------------------------------------*/
                        if (RCI_request==2)
                        {
//                            /* Request to the dfgmres_get routine to put the solution into b[N] via ipar[12]
//                            /*---------------------------------------------------------------------------
//                            /* WARNING: beware that the call to dfgmres_get routine with ipar[12]=0 at this stage may
//                            /* destroy the convergence of the FGMRES method, therefore, only advanced users should
//                            /* exploit this option with care */
                            ipar[12]=1;
                            /* Get the current FGMRES solution in the vector b[N] */
                            dfgmres_get(&ivar, computed_solution, b, &RCI_request, ipar, dpar, tmp, &itercount);
                            /* Compute the current true residual via MKL (Sparse) BLAS routines */
//                             std::cout << "computing residual...";
                            mkl_dcsrgemv(&cvar, &ivar, a, index1_vector, index2_vector, b, residual);
                            dvar=-1.0E0;
                            i=1;
                            daxpy(&ivar, &dvar, rhs, &i, residual, &i);
                            dvar=dnrm2(&ivar,residual,&i);
                            std::cout << " done. dvar = " << dvar << std::endl;
                            if (dvar<1.0E-3) goto COMPLETE;
                            else
                            {
//                                 std::cout << "residual too large: recomputing...";
                                goto ONE;
                            }
                        }
                        
//                        /*---------------------------------------------------------------------------
//                        /* If RCI_request=3, then apply the preconditioner on the vector
//                        /* tmp[ipar[21]-1] and put the result in vector tmp[ipar[22]-1]
//                        /*---------------------------------------------------------------------------
//                        /* NOTE that ipar[21] and ipar[22] contain FORTRAN style addresses,
//                        /* therefore, in C code it is required to subtract 1 from them to get C style
//                        /* addresses
//                        /*---------------------------------------------------------------------------*/
//                         if (RCI_request==3)
//                         {
//                             std::cout << "applying preconditioner...";
//                             if (ipar[3]==3)
//                             {
//                                 tmp[ipar[22]-1+0]=-2.0;
//                                 tmp[ipar[22]-1+1]=0.08519601586107672;
//                                 tmp[ipar[22]-1+2]=-1.1590871369607090;
//                                 tmp[ipar[22]-1+3]=-0.65791939687456980;
//                                 tmp[ipar[22]-1+4]=0.75660051476696133;
//                             }
//                             else 
//                             {
//                                 if(ipar[3]==4) 
//                                 {
//                                     tmp[ipar[22]-1+0]=0.0;
//                                     tmp[ipar[22]-1+1]=0.0;
//                                     tmp[ipar[22]-1+2]=0.0;
//                                     tmp[ipar[22]-1+3]=1.0;
//                                     tmp[ipar[22]-1+4]=-1.0;
//                                 }
//                                 else
//                                 {
//                                     for(i=0;i<n;i++)
//                                     {
//                                         tmp[ipar[22]-1+i]=(double)(i)*tmp[ipar[21]-1+i];
//                                     }
//                                 }
//                             }
//                             std::cout << " done. recomputing...";
//                             goto ONE;
//                         }
//                        /*---------------------------------------------------------------------------
//                        /* If RCI_request=4, then check if the norm of the next generated vector is
//                        /* not zero up to rounding and computational errors. The norm is contained
//                        /* in dpar[6] parameter
//                        /*---------------------------------------------------------------------------*/
                        if (RCI_request==4)
                        {
//                             std::cout << "checking norm of next vector: " << dpar[6];
                            if (dpar[6]<1.0E-12)
                            {
//                                 std::cout << " All OK --> Exiting." << std::endl;
                                goto COMPLETE;
                            }
                            else
                            {
//                                 std::cout << "recomputing...";
                                goto ONE;
                            }
                        }
//                        /*---------------------------------------------------------------------------
//                        /* If RCI_request=anything else, then dfgmres subroutine failed
//                        /* to compute the solution vector: computed_solution[N]
//                        /*---------------------------------------------------------------------------*/
                        else
                        {
                            Error( RCI_request );
                        }
//                        /*---------------------------------------------------------------------------
//                        /* Reverse Communication ends here
//                        /* Get the current iteration number and the FGMRES solution (DO NOT FORGET to
//                        /* call dfgmres_get routine as computed_solution is still containing
//                        /* the initial guess!). Request to dfgmres_get to put the solution
//                        /* into vector computed_solution[N] via ipar[12]
//                        /*---------------------------------------------------------------------------*/
                        COMPLETE:   ipar[12]=0;
                        std::cout << "completing...";
                        dfgmres_get(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp, &itercount);
                        std::cout << "done." << std::endl;
//                        /*---------------------------------------------------------------------------
//                        /* Print solution vector: computed_solution[N] and the number of iterations: itercount
//                        /*---------------------------------------------------------------------------*/
//                         std::cout << "The system has been SUCCESSFULLY solved" << std::endl;
//                         std::cout << "The following solution has been obtained:" << std::endl;
//                         std::cout << "[ ";
//                         for (i=0;i<n;i++) 
//                         {
//                             std::cout << computed_solution[i] << " ";
//                         }
//                         std::cout << std::endl;
//                         printf("\n The expected solution is: \n");
//                         for (i=0;i<n;i++) 
//                         {
//                             printf("expected_solution[%d]=",i);
//                             printf("%e\n",expected_solution[i]);
//                         }
                        
                        goto SUCCEDED;
                        
                        SUCCEDED: std::cout << "GMRES solver succeeded, number of iterations: " << itercount << std::endl;
                        double stop_solver = omp_get_wtime();
                        std::cout << "Solver time: " << stop_solver-start_solver << std::endl;
//                        /*---------------------------------------------------------------------------
//                        /* write back solution to X vector
//                        /*---------------------------------------------------------------------------*/
                        for( int i=0; i<n; i++ )
                            rX[i] = computed_solution[i];
//                         KRATOS_WATCH(rX);
//                         delete [] ipar;
                        delete [] index1_vector;
                        delete [] index2_vector;
//                         delete [] computed_solution;
//                         delete [] rhs;
                        std::cout << "======= EXITING GMRES SOLVER ========" << std::endl;
                        return true;
                        
            }
            
            void Error( int error_code )
            {
                KRATOS_ERROR( std::logic_error, "The solver has returned the ERROR code ", error_code );
            }
            
            /** 
             * Multi solve method for solving a set of linear systems with same coefficient matrix.
             * Solves the linear system Ax=b and puts the result on SystemVector& rX. 
             * rX is also th initial guess for iterative methods.
             * @param rA. System matrix
             * @param rX. Solution vector.
             * @param rB. Right hand side vector.
             */
            bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
            {
                KRATOS_ERROR(std::logic_error,"ERROR: This solver can be used for single RHS only", "");
                return false;
            }
            
            /**
             * Print information about this object.
             */
            void  PrintInfo(std::ostream& rOStream) const
            {
                rOStream << "GMRES solver finished.";
            }
            
            /**
             * Print object's data.
             */
            void  PrintData(std::ostream& rOStream) const 
            {
            }
        
        private:
            
            /**
             * Assignment operator.
             */
            MKLGMRESSolver& operator=(const MKLGMRESSolver& Other);
            
            /**
             * Copy constructor.
             */
//             ParallelSuperLUSolver(const ParallelSuperLUSolver& Other);
    
    }; // Class ParallelSuperLUSolver

    
    /**
     * input stream function
     */
    template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
    inline std::istream& operator >> (std::istream& rIStream, MKLGMRESSolver< TSparseSpaceType, 
                                      TDenseSpaceType, TReordererType>& rThis)
    {
        return rIStream;
    }
    
    /**
     * output stream function
     */
    template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
    inline std::ostream& operator << (std::ostream& rOStream, 
                                      const MKLGMRESSolver<TSparseSpaceType, 
                                      TDenseSpaceType, TReordererType>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
 
  
}  // namespace Kratos.

#endif // KRATOS_MKL_GMRES_SOLVER_H_INCLUDED  defined 


