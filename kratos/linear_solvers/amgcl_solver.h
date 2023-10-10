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

#if !defined(KRATOS_AMGCL_SOLVER )
#define  KRATOS_AMGCL_SOLVER

#ifndef AMGCL_PARAM_UNKNOWN
#include "input_output/logger.h"
#  define AMGCL_PARAM_UNKNOWN(name)                                            \
    Kratos::Logger("AMGCL") << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::WARNING << "Unknown parameter " << name << std::endl
#endif

// System includes
#include <iostream>
#include <fstream>
#include <utility>

// External includes
/* BOOST */
#include <boost/range/iterator_range.hpp>
#include <boost/property_tree/json_parser.hpp>

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/iterative_solver.h"
#include "includes/ublas_interface.h"
#include "spaces/ublas_space.h"

#include <amgcl/coarsening/rigid_body_modes.hpp>

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
enum AMGCLSmoother
{
    SPAI0,SPAI1,ILU0,DAMPED_JACOBI,GAUSS_SEIDEL,CHEBYSHEV
};

enum AMGCLIterativeSolverType
{
    LGMRES,FGMRES,GMRES,BICGSTAB,CG,BICGSTAB_WITH_GMRES_FALLBACK,BICGSTAB2
};

enum AMGCLCoarseningType
{
    RUGE_STUBEN,AGGREGATION,SA,SA_EMIN
};

///@}
///@name  Functions
///@{


void Parameters2PTree(Parameters input,
                      boost::property_tree::ptree& r_output);


boost::property_tree::ptree Parameters2PTree(Parameters parameters);

/**
 * @brief This function solves with Ublas Matrix type
 * @param block_size Block size
 * @param rA System matrix
 * @param rX Solution vector. It's also the initial guess for iterative linear solvers.
 * @param rB Right hand side vector.
 * @param rIterationNumber The current number of iterations
 * @param rResidual The current residual of the problem
 */
void KRATOS_API(KRATOS_CORE) AMGCLSolve(
    int block_size,
    TUblasSparseSpace<double>::MatrixType& rA,
    TUblasSparseSpace<double>::VectorType& rX,
    TUblasSparseSpace<double>::VectorType& rB,
    TUblasSparseSpace<double>::IndexType& rIterationNumber,
    double& rResidual,
    Parameters amgclParams,
    int verbosity_level,
    bool use_gpgpu
    );

///@}
///@name Kratos Classes
///@{

/**
 * @class AMGCLSolver
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
class AMGCLSolver : public LinearSolver< TSparseSpaceType,
    TDenseSpaceType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AMGCLSolver
    KRATOS_CLASS_POINTER_DEFINITION( AMGCLSolver );

    /// The base class definition
    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    /// The sparse matric type
    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    /// Vector type definition
    typedef typename TSparseSpaceType::VectorType VectorType;

    /// Dense matrix type
    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    /// DofArray type
    typedef ModelPart::DofsArrayType DofsArrayType;

    /// The index type definition to be consistent
    typedef typename TSparseSpaceType::IndexType IndexType;

    /// The size type definition
    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief This is the default constructor
     * @param ThisParameters The configuration parameters
     */
    AMGCLSolver(Parameters ThisParameters = Parameters(R"({})"))
    {
        ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
        mAMGCLParameters = ThisParameters["solver_settings"];
    }

    /**
     * Destructor
     */
    ~AMGCLSolver() override {};

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
        // Initial checks
        KRATOS_ERROR_IF(TSparseSpaceType::Size1(rA) != TSparseSpaceType::Size2(rA) ) << "matrix A is not square! sizes are "
            << TSparseSpaceType::Size1(rA) << " and " << TSparseSpaceType::Size2(rA) << std::endl;
        KRATOS_ERROR_IF(TSparseSpaceType::Size(rX)  != TSparseSpaceType::Size1(rA)) << "size of x does not match the size of A. x size is " << TSparseSpaceType::Size(rX)
            << " matrix size is " << TSparseSpaceType::Size1(rA) << std::endl;
        KRATOS_ERROR_IF(TSparseSpaceType::Size(rB) != TSparseSpaceType::Size1(rA)) << "size of b does not match the size of A. b size is " << TSparseSpaceType::Size(rB)
            << " matrix size is " << TSparseSpaceType::Size1(rA) << std::endl;

        const unsigned static_block_size = this->GetBlockSize();
        size_t iters;
        double resid;
        AMGCLSolve(static_block_size, rA,rX,rB, iters, resid, mAMGCLParameters, 1, false);

        KRATOS_INFO("AMGCL Linear Solver")
            << "Iterations: " << iters << std::endl
            << "Residual:   " << resid << std::endl;

        return true;
    }

    /**
     * @brief Multi solve method for solving a set of linear systems with same coefficient matrix.
     * @details Solves the linear system Ax=b and puts the result on SystemVector& rX. rVectorx is also th initial guess for iterative methods.
     * @param rA System matrix
     * @param rX Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {
        return false;
    }

    /**
     * @brief Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * @details Another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function tells if the solver requires such data
     */
    bool AdditionalPhysicalDataIsNeeded() override
    {
        return true;
    }

    /**
     * @brief Some solvers may require a minimum degree of knowledge of the structure of the matrix.
     * @details To make an example when solving a mixed u-p problem, it is important to identify the row associated to v and p. Another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers which require knowledge on the spatial position of the nodes associated to a given dof. This function is the place to eventually provide such data
     * @param rA System matrix
     * @param rX Solution vector. It's also the initial guess for iterative linear solvers.
     * @param rB Right hand side vector.
     */
    void ProvideAdditionalData (
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        DofsArrayType& rDofSet,
        ModelPart& rModelPart
        ) override
    {
        int old_ndof = -1;
        int ndof=0;

        if (!rModelPart.IsDistributed())
        {
            unsigned int old_node_id = rDofSet.size() ? rDofSet.begin()->Id() : 0;
            for (auto it = rDofSet.begin(); it!=rDofSet.end(); it++) {
                if(it->EquationId() < TSparseSpaceType::Size1(rA) ) {
                    IndexType id = it->Id();
                    if(id != old_node_id) {
                        old_node_id = id;
                        if(old_ndof == -1) old_ndof = ndof;
                        else if(old_ndof != ndof) { //if it is different than the block size is 1
                            old_ndof = -1;
                            break;
                        }

                        ndof=1;
                    } else {
                        ndof++;
                    }
                }
            }

            //if(old_ndof == -1)
            //    mBlockSize = 1;
            //else
            //    mBlockSize = ndof;

        }
        else //distribute
        {
            const std::size_t system_size = TSparseSpaceType::Size1(rA);
            int current_rank = rModelPart.GetCommunicator().GetDataCommunicator().Rank();
            unsigned int old_node_id = rDofSet.size() ? rDofSet.begin()->Id() : 0;
            for (auto it = rDofSet.begin(); it!=rDofSet.end(); it++) {
                if(it->EquationId() < system_size  && it->GetSolutionStepValue(PARTITION_INDEX) == current_rank) {
                    IndexType id = it->Id();
                    if(id != old_node_id) {
                        old_node_id = id;
                        if(old_ndof == -1) old_ndof = ndof;
                        else if(old_ndof != ndof) { //if it is different than the block size is 1
                            old_ndof = -1;
                            break;
                        }

                        ndof=1;
                    } else {
                        ndof++;
                    }
                }
            }

            //if(old_ndof != -1)
            //    mBlockSize = ndof;
            //int max_block_size = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(mBlockSize);
            //if( old_ndof == -1) {
            //    mBlockSize = max_block_size;
            //}
            //KRATOS_ERROR_IF(mBlockSize != max_block_size) << "Block size is not consistent. Local: " << mBlockSize  << " Max: " << max_block_size << std::endl;
        }
    }

    static Parameters GetDefaultParameters()
    {
        Parameters output(R"({
            "solver_type" : "amgcl",
            "solver_settings" : {
                "precond" : {
                    "class" : "amg",
                    "relax" : {
                        "type" : "ilu0"
                    },
                    "coarsening" : {
                        "type" : "aggregation",
                        "aggr" : {
                            "eps_strong" : 0.0,
                            "block_size" : 1
                        }
                    },
                    "coarse_enough" : 1000,
                    "npre" : 1,
                    "npost" : 1
                },
                "solver" : {
                    "type" : "gmres",
                    "maxiter" : 100,
                    "M" : 100,
                    "tol" : 1e-6
                }
            }
        })");
        return output;
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
        rOStream << "AMGCL solver:";
    }

    /**
     * Print object's data.
     */
    void  PrintData(std::ostream& rOStream) const override
    {
        rOStream << "Settings: " << mAMGCLParameters;
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

    unsigned GetBlockSize() const noexcept
    {
        unsigned static_block_size = 1;
        if (mAMGCLParameters.Has("precond")) {
            Parameters precond = mAMGCLParameters["precond"];
            if (precond.Has("coarsening")) {
                Parameters coarsening = precond["coarsening"];
                if (coarsening.Has("aggr")) {
                    Parameters aggr = coarsening["aggr"];
                    if (aggr.Has("block_size")) {
                        static_block_size = aggr["block_size"].Get<int>();
                    }
                }
            }
        }
        return static_block_size;
    }

    // Helper function for checking if a selected option is available
    // and printing the available options
    void CheckIfSelectedOptionIsAvailable(
        const Parameters ThisParameters,
        const std::string& rOptionName,
        const std::set<std::string>& rAvailableOptions)
    {
        if (rAvailableOptions.find(ThisParameters[rOptionName].GetString()) == rAvailableOptions.end()) {
            std::stringstream msg;
            msg << "Currently prescribed " << rOptionName << " : " << ThisParameters[rOptionName].GetString() << std::endl;
            msg << "Admissible values are :";
            for (const auto& r_name : rAvailableOptions) {
                msg << std::endl << "    " << r_name;
            }
            KRATOS_ERROR << "AMGCL Linear Solver : " << rOptionName << " is invalid!" << std::endl << msg.str() << std::endl;
        }
    }

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

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    Parameters mAMGCLParameters;  /// The configuration parameters of the AMGCl

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
    AMGCLSolver& operator=(const AMGCLSolver& Other);

    /**
     * Copy constructor.
     */
    AMGCLSolver(const AMGCLSolver& Other);

}; // Class AMGCLSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, AMGCLSolver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AMGCLSolver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos.


#endif // KRATOS_AMGCL_SOLVER  defined
