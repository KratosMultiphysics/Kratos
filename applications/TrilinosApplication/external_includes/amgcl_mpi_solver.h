//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
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

// System includes
#include <iostream>
#include <fstream>
#include <utility>

// Project includes
#include "includes/kratos_parameters.h"
#include "linear_solvers/amgcl_solver.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class AmgclMPISolver
 * @ingroup KratosTrilinosApplication
 * @brief This is a multigrid solver based on the AMGCL library
 * @details Created by Denis Deminov: https://github.com/ddemidov/amgcl
 * AMGCL is a header-only C++ library for solving large sparse linear systems with algebraic multigrid (AMG) method. AMG is one of the most effective iterative methods for solution of equation systems arising, for example, from discretizing PDEs on unstructured grids. The method can be used as a black-box solver for various computational problems, since it does not require any information about the underlying geometry. AMG is often used not as a standalone solver but as a preconditioner within an iterative solver (e.g. Conjugate Gradients, BiCGStab, or GMRES).
 * AMGCL builds the AMG hierarchy on a CPU and then transfers it to one of the provided backends. This allows for transparent acceleration of the solution phase with help of OpenCL, CUDA, or OpenMP technologies. Users may provide their own backends which enables tight integration between AMGCL and the user code.
 * @author Denis Demidov
 * @author Riccardo Rossi
 */
template< class TSparseSpaceType, class TDenseSpaceType>
class KRATOS_API(TRILINOS_APPLICATION) AmgclMPISolver
    : public AMGCLSolver< TSparseSpaceType,TDenseSpaceType >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AmgclMPISolver
    KRATOS_CLASS_POINTER_DEFINITION( AmgclMPISolver );

    /// The sparse matrix type
    using SparseMatrixType = typename TSparseSpaceType::MatrixType;

    /// Vector type definition
    using VectorType = typename TSparseSpaceType::VectorType;

    /// Dense matrix type
    using DenseMatrixType = typename TDenseSpaceType::MatrixType;

    /// DofArray type
    using DofsArrayType = ModelPart::DofsArrayType;

    /// The index type definition
    using IndexType = std::size_t;

    /// The size type definition
    using SizeType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief This is the default constructor
     * @param ThisParameters The configuration parameters
     */
    AmgclMPISolver(Parameters ThisParameters = Parameters(R"({})"))
        :
        AMGCLSolver< TSparseSpaceType,TDenseSpaceType>(ThisParameters) { }

    /// Copy constructor.
    AmgclMPISolver(const AmgclMPISolver& Other) = delete;

    /// Destructor.
    ~AmgclMPISolver() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AmgclMPISolver& operator=(const AmgclMPISolver& Other) = delete;

    ///@}
    ///@name Input and output
    ///@{

    /// Print information about this object.
    void  PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AMGCL-MPI-Solver";
    }

    ///@}

}; // Class AmgclMPISolver

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AmgclMPISolver<TSparseSpaceType,
                                  TDenseSpaceType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos.
