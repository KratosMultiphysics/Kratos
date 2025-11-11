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

#pragma once

#ifndef AMGCL_PARAM_UNKNOWN
#include "input_output/logger.h"
#  define AMGCL_PARAM_UNKNOWN(name)                                            \
    Kratos::Logger("AMGCL") << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::WARNING << "Unknown parameter " << name << std::endl
#endif

// External includes
#include <boost/range/iterator_range.hpp>
#include <boost/property_tree/json_parser.hpp>

// Project includes
#include "linear_solvers/linear_solver.h"
#include "includes/kratos_parameters.h"

// System includes
#include <iostream>
#include <optional>

// The implementation of AMGCLSolver is split between
// - an implementation header ("linear_solvers/amgcl_solver_impl.hpp")
// - implementation sources (e.g.: "linear_solvers/amgcl_solver_impl.cpp")
//
// The reason is twofold:
// - includes from the AMGCL library are extremely heavy, so they are
//   avoided in the class declaration ("linear_solvers/amgcl_solver.h").
//   Instead, the implementation header includes them and defines logic
//   common to any matrix/vector representations. Each source file that
//   defines an instantiation of AMGCLSolver includes the implementation
//   header.
// - Shared memory and distributed memory matrix/vector representations
//   are handled in separate source files to avoid adding a Trilinos
//   dependency to core.


namespace Kratos {


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
template< class TSparseSpaceType, class TDenseSpaceType>
class KRATOS_API(KRATOS_CORE) AMGCLSolver
    : public LinearSolver<TSparseSpaceType,TDenseSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AMGCLSolver
    KRATOS_CLASS_POINTER_DEFINITION( AMGCLSolver );

    /// The base class definition
    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType> BaseType;

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

    AMGCLSolver();

    AMGCLSolver(Parameters Settings);

    /**
     * @brief Default constructor - uses ILU+GMRES
     * @param rSmootherName The smoother type considered
     * @param rSolverName The solver type considered
     * @param Tolerance tolerance that will be achieved by the iterative solver
     * @param MaxIterationsNumber this number represents both the number of iterations AND the size of the krylov space
     * @param Verbosity, a number from 0 (no output) to 2 (maximal output)
     * @param GMRESSize The size of the GMRES
     */
    AMGCLSolver(
        const std::string& rSmootherName,
        const std::string& rSolverName,
        double Tolerance,
        int MaxIterationsNumber,
        int Verbosity,
        int GMRESSize = 50
        );

    /**
     * Default constructor - uses ILU+GMRES
     * @param rSmootherName The smoother type considered
     * @param rSolverName The solver type considered
     * @param rCoarseningName The coarsening type considered
     * @param Tolerance tolerance that will be achieved by the iterative solver
     * @param MaxIterationsNumber this number represents both the number of iterations AND the size of the krylov space
     * @param Verbosity, a number from 0 (no output) to 2 (maximal output)
     * @param GMRESSize The size of the GMRES
     */
    AMGCLSolver(
        const std::string& rSmootherName,
        const std::string& rSolverName,
        const std::string& rCoarseningName,
        double Tolerance,
        int MaxIterationsNumber,
        int Verbosity,
        int GMRESSize = 50,
        bool ProvideCoordinates = false
        );

    AMGCLSolver(AMGCLSolver&&) noexcept;

    AMGCLSolver(const AMGCLSolver&) = delete;

    ~AMGCLSolver() override;

    ///@}
    ///@name Operations
    ///@{

    /// @copydoc LinearSolver::InitializeSolutionStep
    void InitializeSolutionStep(SparseMatrixType& rLhs,
                                VectorType& rSolution,
                                VectorType& rRhs) override;

    /// @copydoc LinearSolver::PerformSolutionStep
    bool PerformSolutionStep(SparseMatrixType& rLhs,
                             VectorType& rSolution,
                             VectorType& rRhs) override;

    /// @copydoc LinearSolver::AdditionalPhysicalDataIsNeeded
    bool AdditionalPhysicalDataIsNeeded() override
    {
        return true;
    }

    /// @copydoc LinearSolver::ProvideAdditionalData
    void ProvideAdditionalData (
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        DofsArrayType& rDofSet,
        ModelPart& rModelPart
        ) override;

    /// @copydoc LinearSolver::Clear
    void Clear() override;

    static Parameters GetDefaultParameters()
    {
        return Parameters(R"(
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
            "block_size"                     : "auto",
            "use_block_matrices_if_possible" : true,
            "coarse_enough"                  : 1000,
            "max_levels"                     : -1,
            "pre_sweeps"                     : 1,
            "post_sweeps"                    : 1,
            "use_gpgpu"                      : false
        })");
    }

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

protected:
    ///@name Protected Operations
    ///@{

    // Helper function for checking if a selected option is available
    // and printing the available options
    void CheckIfSelectedOptionIsAvailable(
        const Parameters Settings,
        const std::string& rOptionName,
        const std::set<std::string>& rAvailableOptions)
    {
        if (rAvailableOptions.find(Settings[rOptionName].GetString()) == rAvailableOptions.end()) {
            std::stringstream msg;
            msg << "Currently prescribed " << rOptionName << " : " << Settings[rOptionName].GetString() << std::endl;
            msg << "Admissible values are :";
            for (const auto& r_name : rAvailableOptions) {
                msg << std::endl << "    " << r_name;
            }
            KRATOS_ERROR << "AMGCL Linear Solver : " << rOptionName << " is invalid!" << std::endl << msg.str() << std::endl;
        }
    }

    ///@}
    ///@name Member Variables
    ///@{

    double mTolerance;                              //< The tolerance considered
    IndexType mMaxIterationsNumber;                 //< The maximum number of iterations considered
    int mVerbosity;                                 //< The versoisty level
    std::optional<int> mBlockSize;                  //< The size of the dof block
    SizeType mGMRESSize;                            //< The size of the GMRES
    SizeType mCoarseEnough;                         //< The level of coarsening allowed
    bool mFallbackToGMRES;                          //< Of consider GMRES as fallback (TODO: Local flag?)
    bool mProvideCoordinates;                       //< If the coordinates are provided (TODO: Local flag?)
    bool mUseBlockMatricesIfPossible;               //< If use the bloack matrices if possible  (TODO: Local flag?)
    bool mUseGPGPU;                                 //< Use GPGPU if available
    std::vector<array_1d<double,3> > mCoordinates;  //< The vector containing the local coordinates
    boost::property_tree::ptree mAMGCLParameters;   //< The configuration parameters of the AMGCl
    bool mUseAMGPreconditioning = true;             //< by default this includes AMG preconditioning

    struct Impl;
    std::unique_ptr<Impl> mpImpl;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method sets the smother type to be considered
     * @param rSmootherName The smother type to be considered
     * @warning This method invalidates the existing hierarchy.
     */
    void SetSmootherType(const std::string& rSmootherName)
    {
        // Clear the existing hierarchy (but not other data).
        auto coordinates = std::move(mCoordinates);
        this->Clear();
        mCoordinates = std::move(coordinates);

        const std::vector<std::string> options {
            "spai0",
            "spai1",
            "ilu0",
            "ilut",
            "iluk",
            "damped_jacobi",
            "gauss_seidel",
            "chebyshev"
        };

        if (std::find(options.begin(), options.end(), rSmootherName) == options.end()) {
            std::stringstream message;
            message << "Invalid AMGCL smoother name: \"" << rSmootherName << "\". Options are:\n";
            for (const auto& r_option : options) message << "\t\"" << r_option << "\"\n";
            KRATOS_ERROR << message.str();
        }

        mAMGCLParameters.put("precond.relax.type", rSmootherName);
    }

    /**
     * @brief This method sets the iterative solver to be considered
     * @param rSolverName The iterative solver to be considered
     * @warning This method invalidates the existing hierarchy.
     */
    void SetIterativeSolverType(const std::string& rSolverName)
    {
        // Clear the existing hierarchy (but not other data).
        auto coordinates = std::move(mCoordinates);
        this->Clear();
        mCoordinates = std::move(coordinates);

        const std::vector<std::string> options {
            "gmres",
            "bicgstab",
            "cg",
            "bicgstabl",
            "lgmres",
            "fgmres",
            "idrs"
        };

        if (std::find(options.begin(), options.end(), rSolverName) == options.end()) {
            std::stringstream message;
            message << "Invalid AMGCL solver name: \"" << rSolverName << "\". Options are:\n";
            for (const auto& r_option : options) message << "\t\"" << r_option << "\"\n";
            KRATOS_ERROR << message.str();
        }

        mAMGCLParameters.put("solver.type", rSolverName);
        if (rSolverName == "gmres" || rSolverName == "fgmres" || rSolverName == "lgmres")
            mAMGCLParameters.put("solver.M",  mGMRESSize);
        else {
            mAMGCLParameters.erase("solver.M");
        }
    }

    /**
     * @brief This method sets the coarsening type to be considered
     * @param rCoarseningName The coarsening type to be considered
     * @warning This method invalidates the existing hierarchy.
     */
    void SetCoarseningType(const std::string& rCoarseningName)
    {
        // Clear the existing hierarchy (but not other data).
        auto coordinates = std::move(mCoordinates);
        this->Clear();
        mCoordinates = std::move(coordinates);

        const std::vector<std::string> options {
            "ruge_stuben",
            "aggregation",
            "smoothed_aggregation",
            "smoothed_aggr_emin"
        };

        if (std::find(options.begin(), options.end(), rCoarseningName) == options.end()) {
            std::stringstream message;
            message << "Invalid AMGCL coarsening name: \"" << rCoarseningName << "\". Options are:\n";
            for (const auto& r_option : options) message << "\t\"" << r_option << "\"\n";
            KRATOS_ERROR << message.str();
        }

        mAMGCLParameters.put("precond.coarsening.type", rCoarseningName);
    }

    void ApplySettings(Parameters Settings);

    ///@}
}; // class AMGCLSolver



template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
std::istream& operator >> (std::istream& rIStream, AMGCLSolver< TSparseSpaceType,
                           TDenseSpaceType>& rThis)
{
    return rIStream;
}


template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
std::ostream& operator << (std::ostream& rOStream,
                           const AMGCLSolver<TSparseSpaceType,
                           TDenseSpaceType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos.
