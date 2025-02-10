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
#include "containers/sparse_graph.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h" //TODO: We use this for the scaling enums

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
 * @details This class collects the methods required to
 * handle the building of the sparse sytem matrices
 * @author Ruben Zorrilla
 */
class AssemblyUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Build type definition
    enum class BuildType
    {
        Block,
        Elimination
    };

    /// Pointer definition of AssemblyUtilities
    KRATOS_CLASS_POINTER_DEFINITION(AssemblyUtilities);

    /// DoF array type definition
    using DofsArrayType = ModelPart::DofsArrayType;

    /// Size type definition
    using SizeType = typename ModelPart::SizeType;

    /// Index type definition
    using IndexType = typename ModelPart::IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AssemblyUtilities() = delete;

    ///@}
    ///@name Operations
    ///@{

    static void SetUpSparseGraph(
        const ModelPart& rModelPart,
        SparseGraph<typename ModelPart::IndexType>& rSparseGraph);
    // {
    //     if (!rSparseGraph.IsEmpty()) {
    //         KRATOS_WARNING("AssemblyUtilities") << "Provided sparse graph is not empty and will be cleared." << std::endl;
    //         rSparseGraph.Clear();
    //     }

    //     Element::EquationIdVectorType eq_ids;
    //     for (auto& r_elem : rModelPart.Elements()) {
    //         r_elem.EquationIdVector(eq_ids, rModelPart.GetProcessInfo());
    //         rSparseGraph.AddEntries(eq_ids);
    //     }
    //     for (auto& r_cond : rModelPart.Conditions()) {
    //         r_cond.EquationIdVector(eq_ids, rModelPart.GetProcessInfo());
    //         rSparseGraph.AddEntries(eq_ids);
    //     }
    // }

    template<BuildType TBuildType, class TMatrixType, class TVectorType>
    static void ApplyDirichletConditions(
        const ModelPart& rModelPart,
        const DofsArrayType& rDofSet,
        TMatrixType& rLHS,
        TVectorType& rRHS)
    {
        if constexpr (TBuildType == BuildType::Block) {

            // Set the scaling factor vector
            const SizeType system_size = rLHS.size1();
            Vector scaling_factors (system_size);

            // Note that DOFs are assumed to be numbered consecutively in the block building
            const auto dof_begin = rDofSet.begin();
            IndexPartition<std::size_t>(rDofSet.size()).for_each([&](IndexType Index){
                auto it_dof = dof_begin + Index;
                scaling_factors[Index] = it_dof->IsFixed() ? 0.0 : 1.0;
            });

            // // Detect if there is a line of all zeros and set the diagonal to a certain number if this happens (1 if not scale, some norms values otherwise)
            // mScaleFactor = TSparseSpace::CheckAndCorrectZeroDiagonalValues(rModelPart.GetProcessInfo(), rA, rb, mScalingDiagonal);

            // double* Avalues = rA.value_data().begin();
            // std::size_t* Arow_indices = rA.index1_data().begin();
            // std::size_t* Acol_indices = rA.index2_data().begin();

            // IndexPartition<std::size_t>(system_size).for_each([&](std::size_t Index){
            //     const std::size_t col_begin = Arow_indices[Index];
            //     const std::size_t col_end = Arow_indices[Index+1];
            //     const double k_factor = scaling_factors[Index];
            //     if (k_factor == 0.0) {
            //         // Zero out the whole row, except the diagonal
            //         for (std::size_t j = col_begin; j < col_end; ++j)
            //             if (Acol_indices[j] != Index )
            //                 Avalues[j] = 0.0;

            //         // Zero out the RHS
            //         rb[Index] = 0.0;
            //     } else {
            //         // Zero out the column which is associated with the zero'ed row
            //         for (std::size_t j = col_begin; j < col_end; ++j)
            //             if(scaling_factors[ Acol_indices[j] ] == 0 )
            //                 Avalues[j] = 0.0;
            //     }
            // });



        } else if (TBuildType == BuildType::Elimination) {

        } else {
            static_assert(TBuildType == BuildType::Block || TBuildType == BuildType::Elimination, "Unssupported build type.");
        }
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


    ///@}
}; // Class AssemblyUtilities

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
