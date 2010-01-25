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
* Revision:            $Revision: 1.3 $ 				 *
*========================================================================*
* International Center of Numerical Methods in Engineering - CIMNE	 *
* Barcelona - Spain 							 *
*========================================================================*
*/

#if !defined(KRATOS_MKL_PARDISO_SOLVER_H_INCLUDED )
#define  KRATOS_MKL_PARDISO_SOLVER_H_INCLUDED

// #define BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT

// External includes 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* PARDISO prototype. */
#if defined(_WIN32) || defined(_WIN64)
#define pardiso_ PARDISO
#else
#define PARDISO pardiso_
#endif
#if defined(MKL_ILP64)
#define MKL_INT long long
#else
#define MKL_INT int
#endif
extern "C"
{
    extern MKL_INT PARDISO
            (void *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
             double *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
             MKL_INT *, double *, double *, MKL_INT *);
}


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
#include "linear_solvers/direct_solver.h"

namespace ublas = boost::numeric::ublas;

namespace Kratos
{
    template< class TSparseSpaceType, class TDenseSpaceType, 
              class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
        class MKLPardisoSolver : public DirectSolver< TSparseSpaceType, 
              TDenseSpaceType, TReordererType>
    {
        public:
            /**
             * Counted pointer of SuperLUSolver
             */
            typedef boost::shared_ptr<MKLPardisoSolver> Pointer;
            
            typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType; 
            
            typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
            
            typedef typename TSparseSpaceType::VectorType VectorType;
            
            typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
            
            /**
             * @param niter number of iterative refinements allowed
             */
            MKLPardisoSolver(unsigned int niter)
            {
                mRefinements = niter;
            }
            
            MKLPardisoSolver()
            {
                mRefinements = 0;
            }
            
            /**
             * Destructor
             */
            virtual ~MKLPardisoSolver(){}
            
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
                std::cout << "Size of index1_vector: " << rA.index1_data().size() << std::endl;
                std::cout << "Size of index2_vector: " << rA.index2_data().size() << std::endl;
                for(unsigned int i = 0; i < rA.index1_data().size(); i++ )
                {
                    index1_vector[i] = (int)(rA.index1_data()[i])+1;
                } 
                for(unsigned int i = 0; i < rA.index2_data().size(); i++ )
                {
                    index2_vector[i] = (int)(rA.index2_data()[i])+1;
                }
                
                MKL_INT mtype = 11; /* Real nonsymmetric matrix */
                /* RHS and solution vectors. */
                double *b = mbtraits::storage(rB);
                double *x = mbtraits::storage(rX);
                
                MKL_INT nrhs = 1; /* Number of right hand sides. */
                /* Internal solver memory pointer pt, */
                /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
                /* or void *pt[64] should be OK on both architectures */
                void *pt[64];
                /* Pardiso control parameters. */
                MKL_INT iparm[64];
                MKL_INT maxfct, mnum, phase, error, msglvl;
                /* Auxiliary variables. */
                MKL_INT i;
                double ddum; /* Double dummy */
                MKL_INT idum; /* Integer dummy. */
                /* -------------------------------------------------------------------- */
                /* .. Setup Pardiso control parameters. */
                /* -------------------------------------------------------------------- */
                for (i = 0; i < 64; i++) {
                    iparm[i] = 0;
                }
                iparm[0] = 1; /* No solver default */
                iparm[1] = 2; /* Fill-in reordering from METIS */
                /* Numbers of processors, value of OMP_NUM_THREADS */
                iparm[2] = omp_get_max_threads();
                std::cout << "number of threads: " << iparm[2] << std::endl;
                if( mRefinements > 0 )
                    iparm[3] = 1; /* iterative-direct algorithm */
                else
                    iparm[3] = 0; /* no iterative-direct algorithm */
                iparm[4] = 0; /* No user fill-in reducing permutation */
                iparm[5] = 0; /* Write solution into x */
                iparm[6] = 0; /* Not in use */
                iparm[7] = mRefinements; /* Max numbers of iterative refinement steps */
                iparm[8] = 0; /* Not in use */
                iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
                iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
                iparm[11] = 0; /* Not in use */
                iparm[12] = 0; /* Not in use */
                iparm[13] = 0; /* Output: Number of perturbed pivots */
                iparm[14] = 0; /* Not in use */
                iparm[15] = 0; /* Not in use */
                iparm[16] = 0; /* Not in use */
                iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
                iparm[18] = -1; /* Output: Mflops for LU factorization */
                iparm[19] = 0; /* Output: Numbers of CG Iterations */
                maxfct = 1; /* Maximum number of numerical factorizations. */
                mnum = 1; /* Which factorization to use. */
                msglvl = 0; /* Print statistical information in file */
                error = 0; /* Initialize error flag */
                /* -------------------------------------------------------------------- */
                /* .. Initialize the internal solver memory pointer. This is only */
                /* necessary for the FIRST call of the PARDISO solver. */
                /* -------------------------------------------------------------------- */
                for (i = 0; i < 64; i++) {
                    pt[i] = 0;
                }
                /* -------------------------------------------------------------------- */
                /* .. Reordering and Symbolic Factorization. This step also allocates */
                /* all memory that is necessary for the factorization. */
                /* -------------------------------------------------------------------- */
                
                phase = 11;
                PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                         &n, a, index1_vector, index2_vector, &idum, &nrhs,
                         iparm, &msglvl, &ddum, &ddum, &error);
                
                if (error != 0) {
                    std::cout << "ERROR during symbolic factorization: " << error << std::endl;
                    exit(1);
                }
                std::cout << "Reordering completed ... " << std::endl;
                //printf("\nNumber of nonzeros in factors = %d", iparm[17]);
                //printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
                /* -------------------------------------------------------------------- */
                /* .. Numerical factorization. */
                /* -------------------------------------------------------------------- */
                phase = 22;
                PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                         &n, a, index1_vector, index2_vector, &idum, &nrhs,
                         iparm, &msglvl, &ddum, &ddum, &error);
                if (error != 0) {
                    std::cout << "ERROR during numerical factorization: " << error << std::endl;
                    exit(2);
                }
                std::cout << "Factorization completed ... " << std::endl;
                /* -------------------------------------------------------------------- */
                /* .. Back substitution and iterative refinement. */
                /* -------------------------------------------------------------------- */
                phase = 33;
                iparm[7] = 2; /* Max numbers of iterative refinement steps. */
                /* Set right hand side to one. */
                //for (i = 0; i < n; i++) {
                //    b[i] = 1;
                //}
                PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                         &n, a, index1_vector, index2_vector, &idum, &nrhs,
                         iparm, &msglvl, b, x, &error);
                if (error != 0) {
                    std::cout << "ERROR during solution: " << error << std::endl;
                    exit(3);
                }
                /* -------------------------------------------------------------------- */
                /* .. Termination and release of memory. */
                /* -------------------------------------------------------------------- */
                phase = -1; /* Release internal memory. */
                PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                         &n, &ddum, index1_vector, index2_vector, &idum, &nrhs,
                         iparm, &msglvl, &ddum, &ddum, &error);
                delete [] index1_vector;
                delete [] index2_vector;
                
                std::cout << "#### SOLVER TIME: " << omp_get_wtime()-start_solver << " ####" << std::endl;
                return true;
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
                rOStream << "PARDISO solver finished.";
            }
            
            /**
             * Print object's data.
             */
            void  PrintData(std::ostream& rOStream) const 
            {
            }
        
        private:
            
            int mRefinements;
            
            /**
             * Assignment operator.
             */
            MKLPardisoSolver& operator=(const MKLPardisoSolver& Other);
            
            /**
             * Copy constructor.
             */
//             ParallelSuperLUSolver(const ParallelSuperLUSolver& Other);
    
    }; // Class ParallelSuperLUSolver

    
    /**
     * input stream function
     */
    template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
    inline std::istream& operator >> (std::istream& rIStream, MKLPardisoSolver< TSparseSpaceType, 
                                      TDenseSpaceType, TReordererType>& rThis)
    {
        return rIStream;
    }
    
    /**
     * output stream function
     */
    template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
    inline std::ostream& operator << (std::ostream& rOStream, 
                                      const MKLPardisoSolver<TSparseSpaceType, 
                                      TDenseSpaceType, TReordererType>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
 
  
}  // namespace Kratos.

#endif // KRATOS_MKL_PARDISO_SOLVER_H_INCLUDED  defined 


