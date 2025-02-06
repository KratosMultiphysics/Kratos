//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/csr_matrix.h"
#include "containers/system_vector.h"
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class AssemblyUtilities
 * @ingroup KratosCore
 * @brief Utility class for handling the build
 * @details This static class collects the methods required to
 * handle the building of the sparse sytem matrices
 * @author Ruben Zorrilla
 */
template<class TDataType = double, class TIndexType=std::size_t, class TMatrixType=CsrMatrix<TDataType, TIndexType>, class TVectorType=SystemVector<TDataType, TIndexType>>
class AssemblyUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssemblyUtilities
    KRATOS_CLASS_POINTER_DEFINITION(AssemblyUtilities);

    /// Index type definition
    using IndexType = TIndexType;

    /// Size type definition
    using SizeType = std::size_t;

    /// Matrix type definition
    using SystemMatrixType = TMatrixType;

    /// Vector type definition
    using SystemVectorType = TVectorType;

    /// Dense space definition
    using DenseSpaceType = UblasSpace<double, Matrix, Vector>;

    /// Local system matrix type definition
    using LocalSystemMatrixType = DenseSpaceType::MatrixType;

    /// Local system vector type definition
    using LocalSystemVectorType = DenseSpaceType::VectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AssemblyUtilities() = delete;

    ///@}
    ///@name Operations
    ///@{

    template<class TAssemblyImplementation, class TThreadLocalStorage>
    static void AssembleLocalSystemElements(
        const ModelPart& rModelPart,
        const TAssemblyImplementation& rAssemblyImplementation,
        SystemMatrixType& rLHS,
        SystemVectorType& rRHS,
        TThreadLocalStorage& rTLS)
    {
        Assemble(rModelPart.Elements(), rModelPart.GetProcessInfo(), rAssemblyImplementation, rLHS, rRHS, rTLS);
    }

    template<class TAssemblyImplementation, class TThreadLocalStorage>
    static void AssembleLocalSystemConditions(
        const ModelPart& rModelPart,
        const TAssemblyImplementation& rAssemblyImplementation,
        SystemMatrixType& rLHS,
        SystemVectorType& rRHS,
        TThreadLocalStorage& rTLS)
    {
        Assemble(rModelPart.Conditions(), rModelPart.GetProcessInfo(), rAssemblyImplementation, rLHS, rRHS, rTLS);
    }

    ///@}
    ///@name Input and output
    ///@{


    ///@}
private:
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    template<class TEntityContainer, class TAssemblyImplementation, class TThreadLocalStorage>
    static void Assemble(
        const TEntityContainer& rEntities,
        const ProcessInfo& rProcessInfo,
        const TAssemblyImplementation& rAssemblyImplementation,
        SystemMatrixType& rLHS,
        SystemVectorType& rRHS,
        TThreadLocalStorage& rTLS)
    {
        // Getting entities container data
        auto ent_begin = rEntities.begin();
        const SizeType n_ent = rEntities.size();

        // Assemble entities
        #pragma omp parallel firstprivate(n_ent, rTLS, rProcessInfo)
        {
            // Assemble elements
            # pragma omp for schedule(guided, 512) nowait
            for (int k = 0; k < n_ent; ++k) {
                auto it_ent = ent_begin + k;
                if (it_ent->IsActive()) {
                    // Calculate local LHS and RHS contributions
                    rAssemblyImplementation(it_ent, rProcessInfo, rTLS);

                    // Get the positions in the global system
                    it_ent->EquationIdVector(rTLS.LocEqIds, rProcessInfo);

                    // Assemble the local contributions to the global system
                    rRHS.Assemble(rTLS.LocRhs, rTLS.LocEqIds);
                    rLHS.Assemble(rTLS.LocLhs, rTLS.LocEqIds);
                }
            }
        }
    }

    ///@}
}; // Class AssemblyUtilities

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
