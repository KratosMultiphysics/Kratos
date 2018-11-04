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
#include <amgcl/preconditioner/runtime.hpp>

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/iterative_solver.h"
#include "includes/ublas_interface.h"

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

/**
 * @brief This function computes a scalar solve for Ublas Matrix type
 * @param rA System matrix
 * @param rX Solution vector. It's also the initial guess for iterative linear solvers.
 * @param rB Right hand side vector.
 * @param rIterationNumber The current number of iterations
 * @param rResidual The current residual of the problem
 */
template <class TSparseSpaceType>
typename std::enable_if<!TSparseSpaceType::IsDistributed(), void>::type
AMGCLScalarSolve(
    typename TSparseSpaceType::MatrixType& rA,
    typename TSparseSpaceType::VectorType& rX,
    typename TSparseSpaceType::VectorType& rB,
    typename TSparseSpaceType::IndexType& rIterationNumber,
    double& rResidual,
    const boost::property_tree::ptree &amgclParams,
    int verbosity_level
    )
{
    typedef amgcl::backend::builtin<double> Backend;

    amgcl::make_solver<
        amgcl::runtime::preconditioner<Backend>,
        amgcl::runtime::solver::wrapper<Backend>
        > solve(amgcl::adapter::zero_copy(TSparseSpaceType::Size1(rA), rA.index1_data().begin(), rA.index2_data().begin(), rA.value_data().begin()), amgclParams);

    std::tie(rIterationNumber, rResidual) = solve(rB, rX);

    if(verbosity_level > 1 )
        std::cout << "AMGCL Memory Occupation : " << amgcl::human_readable_memory(amgcl::backend::bytes(solve)) << std::endl;
}

/**
 * @brief This method solves by block a Ax=b system
 * @param rA System matrix
 * @param rX Solution vector. It's also the initial guess for iterative linear solvers.
 * @param rB Right hand side vector.
 * @param rIterationNumber The current number of iterations
 * @param rResidual The current residual of the problem
 */
template <int TBlockSize, class TSparseSpaceType>
typename std::enable_if<!TSparseSpaceType::IsDistributed(), void>::type
AMGCLBlockSolve(
    typename TSparseSpaceType::MatrixType & rA,
    typename TSparseSpaceType::VectorType& rX,
    typename TSparseSpaceType::VectorType& rB,
    typename TSparseSpaceType::IndexType& rIterationNumber,
    double& rResidual,
    boost::property_tree::ptree amgclParams,
    int verbosity_level
    )
{
    amgclParams.put("precond.coarsening.aggr.block_size",1);

    typedef amgcl::static_matrix<double, TBlockSize, TBlockSize> value_type;
    typedef amgcl::static_matrix<double, TBlockSize, 1> rhs_type;
    typedef amgcl::backend::builtin<value_type> Backend;

    std::size_t n = TSparseSpaceType::Size1(rA);

    amgcl::make_solver<
        amgcl::runtime::preconditioner<Backend>,
        amgcl::runtime::solver::wrapper<Backend>
        > solve( amgcl::adapter::block_matrix<value_type>(std::tie(n,rA.index1_data(),rA.index2_data(),rA.value_data() )), amgclParams);

    rhs_type* x_begin = reinterpret_cast<rhs_type*>(&rX[0]);
    boost::iterator_range<rhs_type*> x_range = boost::make_iterator_range(x_begin, x_begin + n / TBlockSize);

    const rhs_type* b_begin = reinterpret_cast<const rhs_type*>(&rB[0]);
    boost::iterator_range<const rhs_type*> b_range = boost::make_iterator_range(b_begin, b_begin + n / TBlockSize);

    std::tie(rIterationNumber, rResidual) = solve(b_range, x_range);

    if(verbosity_level > 1 )
        std::cout << "AMGCL Memory Occupation : " << amgcl::human_readable_memory(amgcl::backend::bytes(solve)) << std::endl;

}

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
    AMGCLSolver(Parameters ThisParameters = Parameters(R"({})"))
    {
        Parameters default_parameters( R"(
        {
            "preconditioner_type"            : "amg",
            "solver_type"                    : "AMGCL",
            "smoother_type"                  : "ilu0",
            "krylov_type"                    : "gmres",
            "coarsening_type"                : "aggregation",
            "max_iteration"                  : 100,
            "provide_coordinates"            : false,
            "gmres_krylov_space_dimension"   : 100,
            "verbosity"                      : 1,
            "tolerance"                      : 1e-6,
            "scaling"                        : false,
            "block_size"                     : 1,
            "use_block_matrices_if_possible" : true,
            "coarse_enough"                  : 1000,
            "max_levels"                     : -1,
            "pre_sweeps"                     : 1,
            "post_sweeps"                    : 1
        }  )" );

        // Now validate agains defaults -- this also ensures no type mismatch
        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        std::set<std::string> available_preconditioner = {"amg","relaxation","dummy"};

        //selecting preconditioner type - default is AMG
        mAMGCLParameters.put("precond.class", ThisParameters["preconditioner_type"].GetString());
        if(ThisParameters["preconditioner_type"].GetString() == "relaxation") //this implies not using. Use a relaxation sweep as preconditioning. Relaxation type is taken from smoother_type 
        {
            mAMGCLParameters.put("precond.type", ThisParameters["smoother_type"].GetString());
        }
        

        // Validate if values are admissible
        std::set<std::string> available_smoothers = {"spai0","spai1","ilu0","ilut","iluk","damped_jacobi","gauss_seidel","chebyshev"};
        std::set<std::string> available_solvers = {"gmres","bicgstab","cg","bicgstabl","lgmres","fgmres", "bicgstab_with_gmres_fallback","idrs"};
        std::set<std::string> available_coarsening = {"ruge_stuben","aggregation","smoothed_aggregation","smoothed_aggr_emin"};
        


        std::stringstream msg;

        if(available_smoothers.find(ThisParameters["smoother_type"].GetString()) == available_smoothers.end()) {
            msg << "Currently prescribed smoother_type : " << ThisParameters["smoother_type"].GetString() << std::endl;
            msg << "Admissible values are : spai0,spai1,ilu0,ilut,iluk,damped_jacobi,gauss_seidel,chebyshev" << std::endl;
            KRATOS_ERROR << " smoother_type is invalid: " << msg.str() << std::endl;
        }
        if(available_solvers.find(ThisParameters["krylov_type"].GetString()) == available_solvers.end()) {
            msg << "Currently prescribed krylov_type : " << ThisParameters["krylov_type"].GetString() << std::endl;
            msg << "Admissible values are : gmres,bicgstab,cg,bicgstabl,lgmres,fgmres, bicgstab_with_gmres_fallback,idrs" << std::endl;
            KRATOS_ERROR << " krylov_type is invalid: available possibilities are : " << msg.str() << std::endl;
        }
        if(available_coarsening.find(ThisParameters["coarsening_type"].GetString()) == available_coarsening.end()) {
            msg << "Currently prescribed krylov_type : " << ThisParameters["coarsening_type"].GetString() << std::endl;
            msg << "Admissible values are : ruge_stuben,aggregation,smoothed_aggregation,smoothed_aggr_emin" << std::endl;
            KRATOS_ERROR << " coarsening_type is invalid: available possibilities are : " << msg.str() << std::endl;
        }
        if(available_preconditioner.find(ThisParameters["preconditioner_type"].GetString()) == available_preconditioner.end()) {
            msg << "Currently prescribed preconditioner_type : " << ThisParameters["preconditioner_type"].GetString() << std::endl;
            msg << "Admissible values are : amg, relaxation, dummy" << std::endl;
            KRATOS_ERROR << " preconditioner_type is invalid: available possibilities are : " << msg.str() << std::endl;
        }

        mProvideCoordinates = ThisParameters["provide_coordinates"].GetBool();
        mCoarseEnough = ThisParameters["coarse_enough"].GetInt();

        mBlockSize = ThisParameters["block_size"].GetInt(); //set the mndof to an inital number
        mTolerance = ThisParameters["tolerance"].GetDouble();
        mMaxIterationsNumber = ThisParameters["max_iteration"].GetInt();
        mVerbosity=ThisParameters["verbosity"].GetInt();
        mGMRESSize = ThisParameters["gmres_krylov_space_dimension"].GetInt();

        const std::string& solver_type = ThisParameters["krylov_type"].GetString();
        if(solver_type == "gmres" || solver_type == "lgmres" || solver_type == "fgmres") {
            //KRATOS_ERROR << "------------------------  aaaaaaa";
            mAMGCLParameters.put("solver.M",  mGMRESSize);
            mAMGCLParameters.put("solver.type", solver_type);
        } else if(solver_type == "bicgstab_with_gmres_fallback") {
            mAMGCLParameters.put("solver.M",  mGMRESSize);
            mFallbackToGMRES = true;
            mAMGCLParameters.put("solver.type", "bicgstab");
        } else {
            mFallbackToGMRES = false;
            mAMGCLParameters.put("solver.type", solver_type);
        }

        //settings only needed if full AMG is used
        // if(ThisParameters["preconditioner_type"].GetString() == "amg")
        {
            mAMGCLParameters.put("precond.relax.type", ThisParameters["smoother_type"].GetString());
            mAMGCLParameters.put("precond.coarsening.type",  ThisParameters["coarsening_type"].GetString());

            

            int max_levels = ThisParameters["max_levels"].GetInt();
            if(max_levels >= 0)
                mAMGCLParameters.put("precond.max_levels",  max_levels); 

            mAMGCLParameters.put("precond.npre",  ThisParameters["pre_sweeps"].GetInt());
            mAMGCLParameters.put("precond.npost",  ThisParameters["post_sweeps"].GetInt());
        }

        mUseBlockMatricesIfPossible = ThisParameters["use_block_matrices_if_possible"].GetBool();

        if(mProvideCoordinates && mUseBlockMatricesIfPossible) {
            KRATOS_WARNING("AMGCL Linear Solver") << "Sorry coordinates can not be provided when using block matrices, hence setting muse_block_matrices_if_possible to false" << std::endl;
            mUseBlockMatricesIfPossible = false;
            ThisParameters["use_block_matrices_if_possible"].SetBool(false);
        }
        


    }

    /**
     * @brief Default constructor - uses ILU+GMRES
     * @param Smoother The smoother type considered
     * @param Solver The solver type considered
     * @param Tolerance tolerance that will be achieved by the iterative solver
     * @param MaxIterationsNumber this number represents both the number of iterations AND the size of the krylov space
     * @param Verbosity, a number from 0 (no output) to 2 (maximal output)
     * @param GMRESSize The size of the GMRES
     */
    AMGCLSolver(
        AMGCLSmoother Smoother,
        AMGCLIterativeSolverType Solver,
        double Tolerance,
        int MaxIterationsNumber,
        int Verbosity,
        int GMRESSize = 50
        ) : mTolerance(Tolerance),
            mMaxIterationsNumber(MaxIterationsNumber),
            mVerbosity(Verbosity),
            mBlockSize(1),
            mGMRESSize(GMRESSize),
            mCoarseEnough(1000),
            mFallbackToGMRES(false),
            mProvideCoordinates(false),
            mUseBlockMatricesIfPossible(false)
    {
        KRATOS_INFO_IF("AMGCL Linear Solver", mVerbosity > 0) << "Setting up AMGCL for iterative solve " << std::endl;

        // Choose smoother in the list "gauss_seidel, multicolor_gauss_seidel, ilu0, parallel_ilu0, ilut, damped_jacobi, spai0, chebyshev"
        SetSmootherType(Smoother);
        // Setting iterative solver
        SetIterativeSolverType(Solver);

        mAMGCLParameters.put("precond.coarsening.type", "aggregation");
    }

    /**
     * Default constructor - uses ILU+GMRES
     * @param Smoother The smoother type considered
     * @param Solver The solver type considered
     * @param Coarsening The coarsening type considered
     * @param Tolerance tolerance that will be achieved by the iterative solver
     * @param MaxIterationsNumber this number represents both the number of iterations AND the size of the krylov space
     * @param Verbosity, a number from 0 (no output) to 2 (maximal output)
     * @param GMRESSize The size of the GMRES
     */
    AMGCLSolver(
        AMGCLSmoother Smoother,
        AMGCLIterativeSolverType Solver,
        AMGCLCoarseningType Coarsening,
        double Tolerance,
        int MaxIterationsNumber,
        int Verbosity,
        int GMRESSize = 50,
        bool ProvideCoordinates = false
        ) : mTolerance(Tolerance),
            mMaxIterationsNumber(MaxIterationsNumber),
            mVerbosity(Verbosity),
            mBlockSize(1),
            mGMRESSize(GMRESSize),
            mCoarseEnough(1000),
            mFallbackToGMRES(false),
            mProvideCoordinates(ProvideCoordinates),
            mUseBlockMatricesIfPossible(false)
    {
        KRATOS_INFO_IF("AMGCL Linear Solver", mVerbosity > 0) << "Setting up AMGCL for iterative solve " << std::endl;

        // Choose smoother in the list "gauss_seidel, multicolor_gauss_seidel, ilu0, parallel_ilu0, ilut, damped_jacobi, spai0, chebyshev"
        SetSmootherType(Smoother);
        // Setting iterative solver
        SetIterativeSolverType(Solver);
        // Setting coarsening type
        SetCoarseningType(Coarsening);
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

        // Set block size
        if(mAMGCLParameters.get<std::string>("precond.coarsening.type") != std::string("ruge_stuben")) {
            mAMGCLParameters.put("precond.coarsening.aggr.eps_strong",0.0);
            mAMGCLParameters.put("precond.coarsening.aggr.block_size",mBlockSize);
        }
        mAMGCLParameters.put("solver.tol", mTolerance);
        mAMGCLParameters.put("solver.maxiter", mMaxIterationsNumber);

        mAMGCLParameters.put("precond.coarse_enough",mCoarseEnough/mBlockSize);

        Matrix B;
        if(mProvideCoordinates) {
            B = ZeroMatrix(  TSparseSpaceType::Size1(rA), mBlockSize*4  );
            for(IndexType i=0; i<TSparseSpaceType::Size1(rA); i+=mBlockSize) {
                for( IndexType j=0; j<static_cast<IndexType>(mBlockSize); j++) {
                    B(i+j,  j) = 1.0;

                    IndexType inode = i/mBlockSize;

                    B(i+j, mBlockSize +j*3 + 0) = mCoordinates[inode][0];
                    B(i+j, mBlockSize +j*3 + 1) = mCoordinates[inode][1];
                    B(i+j, mBlockSize +j*3 + 2) = mCoordinates[inode][2];
                }
            }
            mAMGCLParameters.put("precond.coarsening.nullspace.cols", B.size2());
            mAMGCLParameters.put("precond.coarsening.nullspace.rows", B.size1());
            mAMGCLParameters.put("precond.coarsening.nullspace.B",    &(B.data()[0]));
        }

        if(mVerbosity > 1)
            write_json(std::cout, mAMGCLParameters);

        if(mVerbosity == 4) {
            //output to matrix market
            std::stringstream matrix_market_name;
            matrix_market_name << "A" <<  ".mm";
            TSparseSpaceType::WriteMatrixMarketMatrix((char*) (matrix_market_name.str()).c_str(), rA, false);

            std::stringstream matrix_market_vectname;
            matrix_market_vectname << "b" << ".mm.rhs";
            TSparseSpaceType::WriteMatrixMarketVector((char*) (matrix_market_vectname.str()).c_str(), rB);

            if(mProvideCoordinates) {
                //output of coordinates
                std::ofstream coordsfile;
                coordsfile.open ("coordinates.txt");
                for(unsigned int i=0; i<mCoordinates.size(); i++) {
                    coordsfile << mCoordinates[i][0] << " " << mCoordinates[i][1] << " " << mCoordinates[i][2] << std::endl;
                }
                coordsfile.close();
            }

            KRATOS_ERROR << " Verbosity = 4 prints the matrix and exits" << std::endl;
        }

        size_t iters;
        double resid;
        {
            if(mFallbackToGMRES) mAMGCLParameters.put("solver.type", "bicgstab"); //first we need to try with bicgstab

            if(mUseBlockMatricesIfPossible) {
                KRATOS_ERROR_IF(TSparseSpaceType::Size1(rA)%mBlockSize != 0) << "The block size employed " << mBlockSize << " is not an exact multiple of the matrix size "
                    << TSparseSpaceType::Size1(rA) << std::endl;
                if(mBlockSize == 1) AMGCLScalarSolve<TSparseSpaceType>(rA,rX,rB, iters, resid, mAMGCLParameters, mVerbosity);
                else if(mBlockSize == 2) AMGCLBlockSolve<2, TSparseSpaceType>(rA,rX,rB, iters, resid, mAMGCLParameters, mVerbosity);
                else if(mBlockSize == 3) AMGCLBlockSolve<3, TSparseSpaceType>(rA,rX,rB, iters, resid, mAMGCLParameters, mVerbosity);
                else if(mBlockSize == 4) AMGCLBlockSolve<4, TSparseSpaceType>(rA,rX,rB, iters, resid, mAMGCLParameters, mVerbosity);
                else
                    AMGCLScalarSolve<TSparseSpaceType>(rA,rX,rB, iters, resid, mAMGCLParameters, mVerbosity);
            } else {
                AMGCLScalarSolve<TSparseSpaceType>(rA,rX,rB, iters, resid, mAMGCLParameters, mVerbosity);
            }
        } //please do not remove this parenthesis!

        if(mFallbackToGMRES && resid > mTolerance ) {
            mAMGCLParameters.put("solver.type", "gmres");
            AMGCLScalarSolve<TSparseSpaceType>(rA,rX,rB, iters, resid, mAMGCLParameters, mVerbosity);
        }

        KRATOS_WARNING_IF("AMGCL Linear Solver", mTolerance < resid)<<"Non converged linear solution. ["<< resid << " > "<< mTolerance << "]" << std::endl;

        KRATOS_INFO_IF("AMGCL Linear Solver", mVerbosity > 1)
                    << "Iterations: " << iters << std::endl
                    << "Error: " << resid << std::endl << std::endl;

        // Setting values
        SetResidualNorm(resid);
        SetIterationsNumber(iters);

        // We check the convergence
        if(resid > mTolerance)
            return false;

        return true;
    }

    /**
     * @brief This method sets the current iteration number
     * @param IterationsNumber The current iteration number
     */
    virtual void SetIterationsNumber(const IndexType IterationsNumber)
    {
        mIterationsNumber = IterationsNumber;
    }

    /**
     * @brief This method sets the current residual norm
     * @param ResidualNorm The current residual norm
     */
    virtual void SetResidualNorm(const double ResidualNorm)
    {
        mResidualNorm = ResidualNorm;
    }

    /**
     * @brief This method returns the current iteration number
     * @return mIterationsNumber The current iteration number
     */
    virtual IndexType GetIterationsNumber()
    {
        return mIterationsNumber;
    }

    /**
     * @brief This method returns the current residual norm
     * @return mResidualNorm The current residual norm
     */
    virtual double GetResidualNorm()
    {
        return mResidualNorm;
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
        unsigned int old_node_id = rDofSet.begin()->Id();
        int ndof=0;
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

        if(old_ndof == -1)
            mBlockSize = 1;
        else
            mBlockSize = ndof;

        KRATOS_INFO_IF("AMGCL Linear Solver", mVerbosity > 1) << "mndof: " << mBlockSize << std::endl;

        if(mProvideCoordinates) {
            mCoordinates.resize(TSparseSpaceType::Size1(rA)/mBlockSize);
            unsigned int i=0;
            for (auto it_dof = rDofSet.begin(); it_dof!=rDofSet.end(); it_dof+=mBlockSize) {
                if(it_dof->EquationId() < TSparseSpaceType::Size1(rA) ) {
                    auto it_node = rModelPart.Nodes().find(it_dof->Id());
                    mCoordinates[ i ] = it_node->Coordinates();
                    i++;
                }
            }
        }
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
        rOStream << "AMGCL solver finished.";
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

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    double mTolerance;                /// The tolerance considered
    IndexType mMaxIterationsNumber;   /// The maximum number of iterations considered
    int mVerbosity;                   /// The versoisty level
    int mBlockSize;                   /// The size of the dof block
    SizeType mGMRESSize;              /// The size of the GMRES
    SizeType mCoarseEnough;           /// The level of coarsening allowed
    bool mFallbackToGMRES;            /// Of consider GMRES as fallback (TODO: Local flag?)
    bool mProvideCoordinates;         /// If the coordinates are provided (TODO: Local flag?)
    bool mUseBlockMatricesIfPossible; /// If use the bloack matrices if possible  (TODO: Local flag?)

    std::vector<array_1d<double,3> > mCoordinates; /// The vector containing the local coordinates

    amgcl::runtime::coarsening::type mCoarsening;  /// The coarsening type considered
    amgcl::runtime::relaxation::type mRelaxation;  /// The relaxation type considered
    amgcl::runtime::solver::type mIterativeSolver; /// The iterative solver considered
    boost::property_tree::ptree mAMGCLParameters;  /// The configuration parameters of the AMGCl

    double mResidualNorm = 0.0;      /// The current residual norm
    IndexType mIterationsNumber = 0; /// The current iteration number

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method sets the smother type to be considered
     * @param SmootherType The smother type to be considered
     */
    void SetSmootherType(const AMGCLSmoother SmootherType)
    {
        switch(SmootherType)
        {
            case SPAI0:
            {
                mAMGCLParameters.put("precond.relax.type","spai0");
                mRelaxation = amgcl::runtime::relaxation::spai0;
                break;
            }
            case SPAI1:
            {
                mAMGCLParameters.put("precond.relax.type","spai1");
                mRelaxation = amgcl::runtime::relaxation::spai1;
                break;
            }
            case ILU0:
            {
                mAMGCLParameters.put("precond.relax.type","ilu0");
                mRelaxation = amgcl::runtime::relaxation::ilu0;
                break;
            }
            case DAMPED_JACOBI:
            {
                mAMGCLParameters.put("precond.relax.type","damped_jacobi");
                mRelaxation = amgcl::runtime::relaxation::damped_jacobi;
                break;
            }
            case GAUSS_SEIDEL:
            {
                mAMGCLParameters.put("precond.relax.type","gauss_seidel");
                mRelaxation = amgcl::runtime::relaxation::gauss_seidel;
                break;
            }
            case CHEBYSHEV:
            {
                mAMGCLParameters.put("precond.relax.type","chebyshev");
                mRelaxation = amgcl::runtime::relaxation::chebyshev;
                break;
            }
        };
    }

    /**
     * @brief This method sets the iterative solver to be considered
     * @param SolverType The iterative solver to be considered
     */
    void SetIterativeSolverType(const AMGCLIterativeSolverType SolverType)
    {
        switch(SolverType)
        {
            case GMRES:
            {
                mAMGCLParameters.put("solver.M",  mGMRESSize);
                mAMGCLParameters.put("solver.type", "gmres");
                break;
            }
            case FGMRES:
            {
                mAMGCLParameters.put("solver.M",  mGMRESSize);
                mAMGCLParameters.put("solver.type", "fgmres");
                break;
            }
            case LGMRES:
            {
                mAMGCLParameters.put("solver.M",  mGMRESSize);
                mAMGCLParameters.put("solver.type", "lgmres");
                break;
            }
            case BICGSTAB:
            {
                mAMGCLParameters.put("solver.type", "bicgstab");
                break;
            }
            case CG:
            {
                mAMGCLParameters.put("solver.type", "cg");
                break;
            }
            case BICGSTAB2:
            {
                mAMGCLParameters.put("solver.type", "bicgstabl");
                break;
            }
            case BICGSTAB_WITH_GMRES_FALLBACK:
            {
                mAMGCLParameters.put("solver.M",  mGMRESSize);
                mAMGCLParameters.put("solver.type", "bicgstab");
                mFallbackToGMRES=true;
                break;
            }
        };
    }

    /**
     * @brief This method sets the coarsening type to be considered
     * @param CoarseningType The coarsening type to be considered
     */
    void SetCoarseningType(const AMGCLCoarseningType CoarseningType)
    {
        switch(CoarseningType)
        {
            case RUGE_STUBEN:
            {
                mAMGCLParameters.put("precond.coarsening.type", "ruge_stuben");
                break;
            }
            case AGGREGATION:
            {
                mAMGCLParameters.put("precond.coarsening.type", "aggregation");
                break;
            }
            case SA:
            {
                mAMGCLParameters.put("precond.coarsening.type", "smoothed_aggregation");
                break;
            }
            case SA_EMIN:
            {
                mAMGCLParameters.put("precond.coarsening.type", "smoothed_aggr_emin");
                break;
            }
        };
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

//#undef MPI_COMM_WORLD

}  // namespace Kratos.


#endif // KRATOS_AMGCL_SOLVER  defined
