//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Denis Demidov
//                   Riccardo Rossi
//

#if !defined(KRATOS_AMGCL_MPI_SOLVER_H_INCLUDED )
#define  KRATOS_AMGCL_MPI_SOLVER_H_INCLUDED

// #ifndef AMGCL_PARAM_MISSING
// #define AMGCL_PARAM_MISSING(name) std::cout << "unset AMGCL parameter with name " << name <<std::endl;
// #endif
// KRATOS_ERROR << , #name)
// Unknown parameter action
#ifndef AMGCL_PARAM_UNKNOWN
#  define AMGCL_PARAM_UNKNOWN(name)                                            \
      std::cerr << "AMGCL WARNING: unknown parameter " << name << std::endl
#endif

// System includes
#include <iostream>
#include <fstream>
#include <utility>

// External includes
/* BOOST */
#include <boost/range/iterator_range.hpp>
#include <boost/property_tree/json_parser.hpp>

/* AMGCL */
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/ublas.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/adapter/block_matrix.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/amgcl_solver.h"


namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{


///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class AmgclMPISolver
 * @ingroup KratosCore
 * @brief This is a multigrid solver based on the AMGCL library
 * @details Created by Denis Deminov: https://github.com/ddemidov/amgcl
 * AMGCL is a header-only C++ library for solving large sparse linear systems with algebraic multigrid (AMG) method. AMG is one of the most effective iterative methods for solution of equation systems arising, for example, from discretizing PDEs on unstructured grids. The method can be used as a black-box solver for various computational problems, since it does not require any information about the underlying geometry. AMG is often used not as a standalone solver but as a preconditioner within an iterative solver (e.g. Conjugate Gradients, BiCGStab, or GMRES).
 * AMGCL builds the AMG hierarchy on a CPU and then transfers it to one of the provided backends. This allows for transparent acceleration of the solution phase with help of OpenCL, CUDA, or OpenMP technologies. Users may provide their own backends which enables tight integration between AMGCL and the user code.
 * @author Denis Demidov
 * @author Riccardo Rossi
 */
template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class AmgclMPISolver : public AMGCLSolver< TSparseSpaceType,TDenseSpaceType, TReordererType >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AmgclMPISolver
    KRATOS_CLASS_POINTER_DEFINITION( AmgclMPISolver );

    /// The base class definition
    typedef AMGCLSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    /// The sparse matric type
    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    /// Vector type definition
    typedef typename TSparseSpaceType::VectorType VectorType;

    /// Dense matrix type
    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    /// DofArray type
    typedef ModelPart::DofsArrayType DofsArrayType;

    /// The index type definition
    typedef std::size_t IndexType;

    /// The size type definition
    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief This is the default constructor
     * @param ThisParameters The configuration parameters
     */
    AmgclMPISolver(Parameters ThisParameters = Parameters(R"({})"))
        :
        AMGCLSolver< TSparseSpaceType,TDenseSpaceType, TReordererType>(ThisParameters)
    {

    }


    /**
     * Destructor
     */
    ~AmgclMPISolver() override {};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Normal solve method.
     * @details Solves the linear system Ax=b and puts the result on SystemVector& rX. rVectorx is also th initial guess for iterative methods.
     * @param rA System matrix
     * @param rX Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB Right hand side vector.
     */
    bool Solve(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB
        ) override
    {
/*        // Initial checks
        KRATOS_ERROR_IF(TSparseSpaceType::Size1(rA) != rA.size2() ) << "matrix A is not square! sizes are " << rA.size1() << " and " << rA.size2() << std::endl;
        KRATOS_ERROR_IF(rX.size() != rA.size1()) << "size of x does not match the size of A. x size is " << rX.size() << " matrix size is " << rA.size1() << std::endl;
        KRATOS_ERROR_IF(rB.size() != rA.size1()) << "size of b does not match the size of A. b size is " << rB.size() << " matrix size is " << rA.size1() << std::endl;

        // Set block size
        if((this->mAMGCLParameters).get<std::string>("precond.coarsening.type") != std::string("ruge_stuben")) {
            this->mAMGCLParameters.put("precond.coarsening.aggr.eps_strong",0.0);
            this->mAMGCLParameters.put("precond.coarsening.aggr.block_size",this->mBlockSize);
        }
        this->mAMGCLParameters.put("solver.tol", this->mTolerance);
        this->mAMGCLParameters.put("solver.maxiter", this->mMaxIterationsNumber);
        this->mAMGCLParameters.put("precond.coarse_enough",this->mCoarseEnough/this->mBlockSize);

        Matrix B;
        if(this->mProvideCoordinates) {
            B = ZeroMatrix(  rA.size1(), this->mBlockSize*4  );
            for(IndexType i=0; i<rA.size1(); i+=this->mBlockSize) {
                for( IndexType j=0; j<static_cast<IndexType>(this->mBlockSize); j++) {
                    B(i+j,  j) = 1.0;

                    IndexType inode = i/this->mBlockSize;

                    B(i+j, this->mBlockSize +j*3 + 0) = this->mCoordinates[inode][0];
                    B(i+j, this->mBlockSize +j*3 + 1) = this->mCoordinates[inode][1];
                    B(i+j, this->mBlockSize +j*3 + 2) = this->mCoordinates[inode][2];
                }
            }
            this->mAMGCLParameters.put("precond.coarsening.nullspace.cols", B.size2());
            this->mAMGCLParameters.put("precond.coarsening.nullspace.rows", B.size1());
            this->mAMGCLParameters.put("precond.coarsening.nullspace.B",    &(B.data()[0]));
        }

        if(this->mVerbosity > 1)
            write_json(std::cout, this->mAMGCLParameters);

        if(this->mVerbosity == 4) {
            KRATOS_ERROR << " Verbosity = 4 for AMGCL not implemented in mpi" << std::endl;
        }

        size_t iters;
        double resid;
        {
            if(this->mFallbackToGMRES) this->mAMGCLParameters.put("solver.type", "bicgstab"); //first we need to try with bicgstab

            if(this->mUseBlockMatricesIfPossible) {
                KRATOS_ERROR_IF(rA.size1()%this->mBlockSize != 0) << "The block size employed " << this->mBlockSize << " is not an exact multiple of the matrix size " << rA.size1() << std::endl;
                if(this->mBlockSize == 1) ScalarSolve(rA,rX,rB, iters, resid);
                else if(this->mBlockSize == 2) BlockSolve<2>(rA,rX,rB, iters, resid);
                else if(this->mBlockSize == 3) BlockSolve<3>(rA,rX,rB, iters, resid);
                else if(this->mBlockSize == 4) BlockSolve<4>(rA,rX,rB, iters, resid);
                else
                    ScalarSolve(rA,rX,rB, iters, resid);
            } else {
                ScalarSolve(rA,rX,rB, iters, resid);
            }
        } //please do not remove this parenthesis!

        if(this->mFallbackToGMRES && resid > this->mTolerance ) {
            this->mAMGCLParameters.put("solver.type", "gmres");
            ScalarSolve(rA,rX,rB, iters, resid);
        }

        KRATOS_WARNING_IF("AMGCL Linear Solver", this->mTolerance < resid)<<"Non converged linear solution. ["<< resid << " > "<< this->mTolerance << "]" << std::endl;

        KRATOS_INFO_IF("AMGCL Linear Solver", this->mVerbosity > 1)
                    << "Iterations: " << iters << std::endl
                    << "Error: " << resid << std::endl << std::endl;

        // Setting values
        this->SetResidualNorm(resid);
        this->SetIterationsNumber(iters);

        // We check the convergence
        if(resid > this->mTolerance)
            return false;
*/
        return true;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /**
     * Print information about this object.
     */
    void  PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AMGCL MPI solver finished.";
    }

    /**
     * Print object's data.
     */
    void  PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    /**
     * @brief This method computes a scalar solve
     * @param rA System matrix
     * @param rX Solution vector. It's also the initial guess for iterative linear solvers.
     * @param rB Right hand side vector.
     * @param rIterationNumber The current number of iterations
     * @param rResidual The current residual of the problem
     */
    void ScalarSolve(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        IndexType& rIterationNumber,
        double& rResidual
        )
    {
        typedef amgcl::backend::builtin<double> Backend;

        amgcl::make_solver<
            amgcl::amg<
                Backend,
                amgcl::runtime::coarsening::wrapper,
                amgcl::runtime::relaxation::wrapper
                >,
            amgcl::runtime::solver::wrapper<Backend>
            > solve(amgcl::adapter::zero_copy(rA.size1(), rA.index1_data().begin(), rA.index2_data().begin(), rA.value_data().begin()), this->mAMGCLParameters);

        // Compute preconditioner
//         KRATOS_INFO_IF("AMGCL Linear Solver", mVerbosity > 0) << solve.precond() << std::endl;
//         else solve.precond();

        std::tie(rIterationNumber, rResidual) = solve(rB, rX);
    }

    /**
     * @brief This method solves by block a Ax=b system
     * @param rA System matrix
     * @param rX Solution vector. It's also the initial guess for iterative linear solvers.
     * @param rB Right hand side vector.
     * @param rIterationNumber The current number of iterations
     * @param rResidual The current residual of the problem
     */
    template< SizeType TBlockSize>
    void BlockSolve(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        IndexType& rIterationNumber,
        double& rResidual
        )
    {
        this->mAMGCLParameters.put("precond.coarsening.aggr.block_size",1);

        typedef amgcl::static_matrix<double, TBlockSize, TBlockSize> value_type;
        typedef amgcl::static_matrix<double, TBlockSize, 1> rhs_type;
        typedef amgcl::backend::builtin<value_type> Backend;

        SizeType n = rA.size1();

        amgcl::make_solver<
            amgcl::amg<
                Backend,
                amgcl::runtime::coarsening::wrapper,
                amgcl::runtime::relaxation::wrapper
                >,
            amgcl::runtime::solver::wrapper<Backend>
            > solve( amgcl::adapter::block_matrix<value_type>(std::tie(n,rA.index1_data(),rA.index2_data(),rA.value_data() )), this->mAMGCLParameters);

//         //compute preconditioner
//         if(mverbosity > 0) std::cout << solve.precond() << std::endl;
//         else solve.precond();

        rhs_type* x_begin = reinterpret_cast<rhs_type*>(&rX[0]);
        boost::iterator_range<rhs_type*> x_range = boost::make_iterator_range(x_begin, x_begin + n / TBlockSize);

        const rhs_type* b_begin = reinterpret_cast<const rhs_type*>(&rB[0]);
        boost::iterator_range<const rhs_type*> b_range = boost::make_iterator_range(b_begin, b_begin + n / TBlockSize);

        std::tie(rIterationNumber, rResidual) = solve(b_range, x_range);
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{
    /**
     * Assignment operator.
     */
    AmgclMPISolver& operator=(const AmgclMPISolver& Other);

    /**
     * Copy constructor.
     */
    AmgclMPISolver(const AmgclMPISolver& Other);

}; // Class AmgclMPISolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, AmgclMPISolver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AmgclMPISolver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

//#undef MPI_COMM_WORLD

}  // namespace Kratos.


#endif // KRATOS_AMGCL_MPI_SOLVER_H_INCLUDED  defined
