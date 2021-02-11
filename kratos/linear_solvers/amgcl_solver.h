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
    boost::property_tree::ptree amgclParams,
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
            "post_sweeps"                    : 1,
            "use_gpgpu"                      : false
        }  )" );

        // Now validate agains defaults -- this also ensures no type mismatch
        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // specify available options
        std::set<std::string> available_smoothers = {"spai0","spai1","ilu0","ilut","iluk","damped_jacobi","gauss_seidel","chebyshev"};
        std::set<std::string> available_solvers = {"gmres","bicgstab","cg","bicgstabl","lgmres","fgmres", "bicgstab_with_gmres_fallback","idrs"};
        std::set<std::string> available_coarsening = {"ruge_stuben","aggregation","smoothed_aggregation","smoothed_aggr_emin"};
        std::set<std::string> available_preconditioner = {"amg","relaxation","dummy"};

        // Validate if values are admissible
        CheckIfSelectedOptionIsAvailable(ThisParameters, "smoother_type",       available_smoothers);
        CheckIfSelectedOptionIsAvailable(ThisParameters, "krylov_type",         available_solvers);
        CheckIfSelectedOptionIsAvailable(ThisParameters, "coarsening_type",     available_coarsening);
        CheckIfSelectedOptionIsAvailable(ThisParameters, "preconditioner_type", available_preconditioner);

        //selecting preconditioner type - default is AMG
        mAMGCLParameters.put("precond.class", ThisParameters["preconditioner_type"].GetString());
        if(ThisParameters["preconditioner_type"].GetString() != "amg"){
            mUseAMGPreconditioning = false;
        }

        if(ThisParameters["preconditioner_type"].GetString() == "relaxation") //this implies not using. Use a relaxation sweep as preconditioning. Relaxation type is taken from smoother_type
        {
            mAMGCLParameters.put("precond.type", ThisParameters["smoother_type"].GetString());
        }

        mProvideCoordinates = ThisParameters["provide_coordinates"].GetBool();
        mCoarseEnough = ThisParameters["coarse_enough"].GetInt();

        mBlockSize = ThisParameters["block_size"].GetInt(); //set the mndof to an inital number
        mTolerance = ThisParameters["tolerance"].GetDouble();
        mMaxIterationsNumber = ThisParameters["max_iteration"].GetInt();
        mVerbosity=ThisParameters["verbosity"].GetInt();
        mGMRESSize = ThisParameters["gmres_krylov_space_dimension"].GetInt();

        const std::string& solver_type = ThisParameters["krylov_type"].GetString();
        mAMGCLParameters.put("solver.type", solver_type);
        mFallbackToGMRES = false;

        if(solver_type == "bicgstab_with_gmres_fallback") {
            mFallbackToGMRES = true;
            mAMGCLParameters.put("solver.type", "bicgstab");
        }

        //settings only needed if full AMG is used
        if(mUseAMGPreconditioning)
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

        mUseGPGPU = ThisParameters["use_gpgpu"].GetBool();
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

        mAMGCLParameters.put("solver.tol", mTolerance);
        mAMGCLParameters.put("solver.maxiter", mMaxIterationsNumber);

        if(mUseAMGPreconditioning)
            mAMGCLParameters.put("precond.coarse_enough",mCoarseEnough/mBlockSize);

        // Use rigid body modes or set block size
        int static_block_size = mUseBlockMatricesIfPossible ? mBlockSize : 1;
        if(mUseAMGPreconditioning && mProvideCoordinates && (mBlockSize == 2 || mBlockSize == 3)) {
            std::vector<double> B;
            int nmodes = amgcl::coarsening::rigid_body_modes(mBlockSize,
                    boost::make_iterator_range(
                        &mCoordinates[0][0],
                        &mCoordinates[0][0] + TSparseSpaceType::Size1(rA)),
                    B);

            static_block_size = 1;
            mAMGCLParameters.put("precond.coarsening.aggr.eps_strong", 0.0);
            mAMGCLParameters.put("precond.coarsening.aggr.block_size", 1);
            mAMGCLParameters.put("precond.coarsening.nullspace.cols",  nmodes);
            mAMGCLParameters.put("precond.coarsening.nullspace.rows",  TSparseSpaceType::Size1(rA));
            mAMGCLParameters.put("precond.coarsening.nullspace.B",     &B[0]);
        } else if(mUseAMGPreconditioning && mAMGCLParameters.get<std::string>("precond.coarsening.type") != std::string("ruge_stuben")) {
            mAMGCLParameters.put("precond.coarsening.aggr.eps_strong", 0.0);
            mAMGCLParameters.put("precond.coarsening.aggr.block_size", mBlockSize);
        }

        if (mVerbosity > 2) {
            write_json(std::cout, mAMGCLParameters);
        }

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
                    coordsfile << mCoordinates[i][0] << " " << mCoordinates[i][1] << " " << mCoordinates[i][2] << "\n";
                }
                coordsfile.close();
            }

            KRATOS_ERROR << " Verbosity = 4 prints the matrix and exits" << std::endl;
        }

        size_t iters;
        double resid;
        {
            if(mFallbackToGMRES) mAMGCLParameters.put("solver.type", "bicgstab"); //first we need to try with bicgstab

            if(mAMGCLParameters.get<std::string>("solver.type") == "gmres" ||
                mAMGCLParameters.get<std::string>("solver.type") == "lgmres" ||
                mAMGCLParameters.get<std::string>("solver.type") == "fgmres" )
                mAMGCLParameters.put("solver.M",  mGMRESSize);
            else
                mAMGCLParameters.erase("solver.M");

            if(mUseBlockMatricesIfPossible) {
                KRATOS_ERROR_IF(TSparseSpaceType::Size1(rA)%mBlockSize != 0) << "The block size employed " << mBlockSize << " is not an exact multiple of the matrix size "
                    << TSparseSpaceType::Size1(rA) << std::endl;
            }
            AMGCLSolve(static_block_size, rA,rX,rB, iters, resid, mAMGCLParameters, mVerbosity, mUseGPGPU);
        } //please do not remove this parenthesis!

        if(mFallbackToGMRES && resid > mTolerance ) {
            mAMGCLParameters.put("solver.type", "gmres");
            mAMGCLParameters.put("solver.M",  mGMRESSize);
            AMGCLSolve(1, rA,rX,rB, iters, resid, mAMGCLParameters, mVerbosity, mUseGPGPU);
        }

        KRATOS_WARNING_IF("AMGCL Linear Solver", mTolerance < resid)<<"Non converged linear solution. ["<< resid << " > "<< mTolerance << "]" << std::endl;

        KRATOS_INFO_IF("AMGCL Linear Solver", mVerbosity > 1)
            << "Iterations: " << iters << std::endl
            << "Error: "      << resid << std::endl;

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
    IndexType GetIterationsNumber() override
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

            if(old_ndof == -1)
                mBlockSize = 1;
            else
                mBlockSize = ndof;

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

            if(old_ndof != -1)
                mBlockSize = ndof;

            int max_block_size = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(mBlockSize);

            if( old_ndof == -1) {
                mBlockSize = max_block_size;
            }

            KRATOS_ERROR_IF(mBlockSize != max_block_size) << "Block size is not consistent. Local: " << mBlockSize  << " Max: " << max_block_size << std::endl;
        }

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
        rOStream << "AMGCL solver:";
    }

    /**
     * Print object's data.
     */
    void  PrintData(std::ostream& rOStream) const override
    {
        rOStream << "Settings: ";
        write_json(rOStream, mAMGCLParameters);
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

    double mTolerance;                /// The tolerance considered
    IndexType mMaxIterationsNumber;   /// The maximum number of iterations considered
    int mVerbosity;                   /// The versoisty level
    int mBlockSize;                   /// The size of the dof block
    SizeType mGMRESSize;              /// The size of the GMRES
    SizeType mCoarseEnough;           /// The level of coarsening allowed
    bool mFallbackToGMRES;            /// Of consider GMRES as fallback (TODO: Local flag?)
    bool mProvideCoordinates;         /// If the coordinates are provided (TODO: Local flag?)
    bool mUseBlockMatricesIfPossible; /// If use the bloack matrices if possible  (TODO: Local flag?)
    bool mUseGPGPU;                   /// Use GPGPU if available

    std::vector<array_1d<double,3> > mCoordinates; /// The vector containing the local coordinates

    boost::property_tree::ptree mAMGCLParameters;  /// The configuration parameters of the AMGCl

    double mResidualNorm = 0.0;      /// The current residual norm
    IndexType mIterationsNumber = 0; /// The current iteration number
    bool mUseAMGPreconditioning = true; ///by default this includes AMG preconditioning

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
                break;
            }
            case SPAI1:
            {
                mAMGCLParameters.put("precond.relax.type","spai1");
                break;
            }
            case ILU0:
            {
                mAMGCLParameters.put("precond.relax.type","ilu0");
                break;
            }
            case DAMPED_JACOBI:
            {
                mAMGCLParameters.put("precond.relax.type","damped_jacobi");
                break;
            }
            case GAUSS_SEIDEL:
            {
                mAMGCLParameters.put("precond.relax.type","gauss_seidel");
                break;
            }
            case CHEBYSHEV:
            {
                mAMGCLParameters.put("precond.relax.type","chebyshev");
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
