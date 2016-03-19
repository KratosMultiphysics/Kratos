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
#include <mkl_types.h>
/* PARDISO prototype. */
extern "C"
{
    extern MKL_INT PARDISO
    (void *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
     double *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
     MKL_INT *, double *, double *, MKL_INT *);
}


#include <boost/timer.hpp>

#include "utilities/openmp_utils.h"

#include "boost/smart_ptr.hpp"
#include "includes/ublas_interface.h"

#include <boost/numeric/bindings/traits/sparse_traits.hpp>
#include <boost/numeric/bindings/traits/matrix_traits.hpp>
#include <boost/numeric/bindings/traits/vector_traits.hpp>
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
    KRATOS_CLASS_POINTER_DEFINITION( MKLPardisoSolver );

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
    virtual ~MKLPardisoSolver() {}

    virtual bool AdditionalPhysicalDataIsNeeded()
    {
        return false;
    }

    virtual void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    )
    {}

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
        double start_solver = OpenMPUtils::GetCurrentTime();
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
        MKL_INT *index1_vector = new (std::nothrow) MKL_INT[rA.index1_data().size()];
        MKL_INT *index2_vector = new (std::nothrow) MKL_INT[rA.index2_data().size()];
        std::cout << "Size of the problem: " << n << std::endl;
        std::cout << "Size of index1_vector: " << rA.index1_data().size() << std::endl;
        std::cout << "Size of index2_vector: " << rA.index2_data().size() << std::endl;
//                 std::cout << "pardiso_solver: line 156" << std::endl;
        for(unsigned int i = 0; i < rA.index1_data().size(); i++ )
        {
            index1_vector[i] = (MKL_INT)(rA.index1_data()[i])+1;
        }
//                 std::cout << "pardiso_solver: line 161" << std::endl;
        for(unsigned int i = 0; i < rA.index2_data().size(); i++ )
        {
            index2_vector[i] = (MKL_INT)(rA.index2_data()[i])+1;
        }
        /**
         *  Matrix type flag:
         * 1    real and structurally symmetric
         * 2    real and symmetic positive definite
         * -2   real and symmetric indefinite
         * 3    complex and structurally symmetric
         * 4    complex and Hermitian positive definite
         * -4   complex and Hermitian indefinite
         * 6    complex and symmetic
         * 11   real and nonsymmetric
         * 13   complex and nonsymmetric
         */
        MKL_INT mtype = 11;
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
        for (i = 0; i < 64; i++)
        {
            iparm[i] = 0;
        }
        iparm[0] = 1; /* No solver default */
        iparm[1] = 2; /* Fill-in reordering from METIS */
        /* Numbers of processors, value of OMP_NUM_THREADS */
//        iparm[2] = OpenMPUtils::GetNumThreads(); //omp_get_max_threads();
        iparm[2] = OpenMPUtils::GetNumThreads(); //omp_get_num_procs();
        std::cout << "Number of threads/procs (for MKL): " << iparm[2] << std::endl;
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
        iparm[12] = 1; /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
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
        for (i = 0; i < 64; i++)
        {
            pt[i] = 0;
        }
        /* -------------------------------------------------------------------- */
        /* .. Reordering and Symbolic Factorization. This step also allocates */
        /* all memory that is necessary for the factorization. */
        /* -------------------------------------------------------------------- */
//                 std::cout << "pardiso_solver: line 241" << std::endl;
        phase = 11;
        PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, a, index1_vector, index2_vector, &idum, &nrhs,
                 iparm, &msglvl, &ddum, &ddum, &error);

        if (error != 0)
        {
            std::cout << "ERROR during symbolic factorization: " << error << std::endl;
            ErrorCheck(error);
            exit(1);
        }
//                 std::cout << "pardiso_solver: line 251" << std::endl;
        std::cout << "Reordering completed ... " << std::endl;
        //printf("\nNumber of nonzeros in factors = %d", iparm[17]);
        //printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
        /* -------------------------------------------------------------------- */
        KRATOS_WATCH(iparm[63]);

        /* .. Numerical factorization. */
        /* -------------------------------------------------------------------- */
        phase = 22;
        PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, a, index1_vector, index2_vector, &idum, &nrhs,
                 iparm, &msglvl, &ddum, &ddum, &error);
        if (error != 0)
        {
            std::cout << "ERROR during numerical factorization: " << error << std::endl;
            ErrorCheck(error);
            exit(2);
        }
        std::cout << "Factorization completed ... " << std::endl;
//                 std::cout << "pardiso_solver: line 267" << std::endl;
        /* -------------------------------------------------------------------- */
        /* .. Back substitution and iterative refinement. */
        /* -------------------------------------------------------------------- */
        phase = 33;
        iparm[7] = 2; /* Max numbers of iterative refinement steps. */
        /* Set right hand side to one. */
        //for (i = 0; i < n; i++) {
        //    b[i] = 1;
        //}
//                 std::cout << "pardiso_solver: line 277" << std::endl;
        PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, a, index1_vector, index2_vector, &idum, &nrhs,
                 iparm, &msglvl, b, x, &error);
        if (error != 0)
        {
            std::cout << "ERROR during solution: " << error << std::endl;
            ErrorCheck(error);
            exit(3);
        }
//                 std::cout << "pardiso_solver: line 285" << std::endl;
        /* -------------------------------------------------------------------- */
        /* .. Termination and release of memory. */
        /* -------------------------------------------------------------------- */
//                 std::cout << "pardiso_solver: line 289" << std::endl;
        phase = -1; /* Release internal memory. */
        PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, &ddum, index1_vector, index2_vector, &idum, &nrhs,
                 iparm, &msglvl, &ddum, &ddum, &error);
//                 std::cout << "pardiso_solver: line 294" << std::endl;
        delete [] index1_vector;
//                 std::cout << "pardiso_solver: line 296" << std::endl;
        delete [] index2_vector;
//                 std::cout << "pardiso_solver: line 298" << std::endl;

        std::cout << "#### SOLVER TIME: " << OpenMPUtils::GetCurrentTime()-start_solver << " ####" << std::endl;
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
        double start_solver = OpenMPUtils::GetCurrentTime();
        typedef boost::numeric::bindings::traits::sparse_matrix_traits<SparseMatrixType> matraits;
        typedef boost::numeric::bindings::traits::matrix_traits<DenseMatrixType> mbtraits;
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
        MKL_INT *index1_vector = new (std::nothrow) MKL_INT[rA.index1_data().size()];
        MKL_INT *index2_vector = new (std::nothrow) MKL_INT[rA.index2_data().size()];
        std::cout << "Size of the problem: " << n << std::endl;
        std::cout << "Size of index1_vector: " << rA.index1_data().size() << std::endl;
        std::cout << "Size of index2_vector: " << rA.index2_data().size() << std::endl;
//                 std::cout << "pardiso_solver: line 156" << std::endl;
        for(unsigned int i = 0; i < rA.index1_data().size(); i++ )
        {
            index1_vector[i] = (MKL_INT)(rA.index1_data()[i])+1;
        }
//                 std::cout << "pardiso_solver: line 161" << std::endl;
        for(unsigned int i = 0; i < rA.index2_data().size(); i++ )
        {
            index2_vector[i] = (MKL_INT)(rA.index2_data()[i])+1;
        }
        /**
         *  Matrix type flag:
         * 1    real and structurally symmetric
         * 2    real and symmetic positive definite
         * -2   real and symmetric indefinite
         * 3    complex and structurally symmetric
         * 4    complex and Hermitian positive definite
         * -4   complex and Hermitian indefinite
         * 6    complex and symmetic
         * 11   real and nonsymmetric
         * 13   complex and nonsymmetric
         */
        MKL_INT mtype = 11;
        MKL_INT nrhs = mbtraits::size2(rB); /* Number of right hand sides. */

        /* RHS and solution vectors. */
        DenseMatrixType Bt = trans(rB);
        DenseMatrixType Xt = ZeroMatrix(nrhs, n);
        double *b = mbtraits::storage(Bt);
        double *x = mbtraits::storage(Xt);

        // inefficient copy
//        double b[nrhs * n];
//        double x[nrhs * n];
//        for(int i = 0; i < nrhs; ++i)
//        {
//            std::copy(column(rB, i).begin(), column(rB, i).end(), b + i * n);
//        }

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
        for (i = 0; i < 64; i++)
        {
            iparm[i] = 0;
        }
        iparm[0] = 1; /* No solver default */
        iparm[1] = 2; /* Fill-in reordering from METIS */
        /* Numbers of processors, value of OMP_NUM_THREADS */
        iparm[2] = OpenMPUtils::GetNumThreads(); //omp_get_max_threads();
//        iparm[2] = OpenMPUtils::GetNumProcs(); //omp_get_num_procs();
        std::cout << "Number of threads/procs (for MKL): " << iparm[2] << std::endl;
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
        iparm[12] = 1; /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
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
        for (i = 0; i < 64; i++)
        {
            pt[i] = 0;
        }
        /* -------------------------------------------------------------------- */
        /* .. Reordering and Symbolic Factorization. This step also allocates */
        /* all memory that is necessary for the factorization. */
        /* -------------------------------------------------------------------- */
//                 std::cout << "pardiso_solver: line 241" << std::endl;
        phase = 11;
        PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, a, index1_vector, index2_vector, &idum, &nrhs,
                 iparm, &msglvl, &ddum, &ddum, &error);

        if (error != 0)
        {
            std::cout << "ERROR during symbolic factorization: " << error << std::endl;
            ErrorCheck(error);
            exit(1);
        }
//                 std::cout << "pardiso_solver: line 251" << std::endl;
        std::cout << "Reordering completed ... " << std::endl;
        //printf("\nNumber of nonzeros in factors = %d", iparm[17]);
        //printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
        /* -------------------------------------------------------------------- */
        KRATOS_WATCH(iparm[63]);

        /* .. Numerical factorization. */
        /* -------------------------------------------------------------------- */
        phase = 22;
        PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, a, index1_vector, index2_vector, &idum, &nrhs,
                 iparm, &msglvl, &ddum, &ddum, &error);
        if (error != 0)
        {
            std::cout << "ERROR during numerical factorization: " << error << std::endl;
            ErrorCheck(error);
            exit(2);
        }
        std::cout << "Factorization completed ... " << std::endl;
//                 std::cout << "pardiso_solver: line 267" << std::endl;
        /* -------------------------------------------------------------------- */
        /* .. Back substitution and iterative refinement. */
        /* -------------------------------------------------------------------- */
        phase = 33;
        iparm[7] = 2; /* Max numbers of iterative refinement steps. */
        /* Set right hand side to one. */
        //for (i = 0; i < n; i++) {
        //    b[i] = 1;
        //}
//                 std::cout << "pardiso_solver: line 277" << std::endl;
        PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, a, index1_vector, index2_vector, &idum, &nrhs,
                 iparm, &msglvl, b, x, &error);
        if (error != 0)
        {
            std::cout << "ERROR during solution: " << error << std::endl;
            ErrorCheck(error);
            exit(3);
        }
//                 std::cout << "pardiso_solver: line 285" << std::endl;
        /* -------------------------------------------------------------------- */
        /* .. Termination and release of memory. */
        /* -------------------------------------------------------------------- */
//                 std::cout << "pardiso_solver: line 289" << std::endl;
        phase = -1; /* Release internal memory. */
        PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, &ddum, index1_vector, index2_vector, &idum, &nrhs,
                 iparm, &msglvl, &ddum, &ddum, &error);
//                 std::cout << "pardiso_solver: line 294" << std::endl;
        delete [] index1_vector;
//                 std::cout << "pardiso_solver: line 296" << std::endl;
        delete [] index2_vector;
//                 std::cout << "pardiso_solver: line 298" << std::endl;

        // inefficient copy
//        for(int i = 0; i < nrhs; ++i)
//        {
//            std::copy(x + i * n, x + (i + 1) * n, column(rX, i).begin());
//        }
        
        noalias(rX) = trans(Xt);
        
        std::cout << "#### SOLVER TIME: " << OpenMPUtils::GetCurrentTime()-start_solver << " ####" << std::endl;
        return true;
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

    void ErrorCheck(MKL_INT error)
    {
        switch(error)
        {
            case -1:
                std::cout << "Input inconsistent" << std::endl;
                break;
            case -2:
                std::cout << "Not enough memory" << std::endl;
                break;
            case -3:
                std::cout << "Reordering problem" << std::endl;
                break;
            case -4:
                std::cout << "Zero pivot, numerical factorization or iterative refinement problem" << std::endl;
                break;
            case -5:
                std::cout << "Unclassified (internal) error" << std::endl;
                break;
            case -6:
                std::cout << "Reordering failed (matrix types 11, 13 only)" << std::endl;
                break;
            case -7:
                std::cout << "Diagonal matrix problem" << std::endl;
                break;
            case -8:
                std::cout << "32-bit integer overflow problem" << std::endl;
                break;
            case -9:
                std::cout << "Not enough memory for OOC" << std::endl;
                break;
            case -10:
                std::cout << "Problems with opening OOC temporary files" << std::endl;
                break;
            case -11:
                std::cout << "Read/write problems with the OOC data file" << std::endl;
                break;
        }
    }

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


