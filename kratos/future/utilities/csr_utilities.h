//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project include
#include "containers/nd_data.h"
#include "includes/model_part.h"

namespace Kratos
{

///@addtogroup KratosCore
///@{

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

class CsrUtilities final
{
public:
    ///@name Type Definitions
    ///@{

    using EquationIdVectorType = std::vector<std::size_t>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CsrUtilities() = delete;

    /// Destructor.
    ~CsrUtilities() = delete;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // template<class TContainerType, class TCsrMatrixType>
    // static void GetEquationIdCsrIndices(
    //     const TContainerType& rContainer,
    //     const ProcessInfo& rProcessInfo,
    //     const TCsrMatrixType& rCsrMatrix,
    //     NDData<int>& rNDData)
    // {
    //     // Get equation ids size from the first entity (assuming all entities have the same size)
    //     EquationIdVectorType equation_ids;
    //     rContainer.begin()->EquationIdVector(equation_ids, rProcessInfo);
    //     const std::size_t local_size = equation_ids.size();

    //     // Assign the input NDData to have shape: number of entities * local_size * local_size
    //     const std::size_t n_entities = std::distance(rContainer.begin(), rContainer.end());
    //     DenseVector<unsigned int> nd_data_shape(3);
    //     nd_data_shape[0] = n_entities;
    //     nd_data_shape[1] = local_size;
    //     nd_data_shape[2] = local_size;
    //     rNDData = NDData<int>(nd_data_shape);

    //     // Loop over the container
    //     auto data_view = rNDData.ViewData();
    //     IndexPartition<std::size_t>(n_entities).for_each([&](std::size_t i) {
    //         // Get current entity
    //         auto it = rContainer.begin() + i;
    //         const std::size_t it_pos = i * (local_size * local_size);

    //         // Get current entity equation ids
    //         EquationIdVectorType equation_ids;
    //         it->EquationIdVector(equation_ids, rProcessInfo);

    //         // Loop over the DOFs
    //         for (unsigned int i_local = 0; i_local < local_size; ++i_local) {
    //             const unsigned int i_global = equation_ids[i_local]; // Row global equation id
    //             for(unsigned int j_local = 0; j_local < local_size; ++j_local) {
    //                 const unsigned int j_global = equation_ids[j_local]; // Column global equation id
    //                 const unsigned int csr_index = rCsrMatrix.FindValueIndex(i_global,j_global); // Index in the CSR matrix values vector
    //                 data_view[it_pos + i_local * local_size + j_local] = csr_index;
    //             }
    //         }
    //     });
    // }

    // template<class TCsrMatrixType>
    // static void Assemble(
    //     const NDData<double>& rLeftHandSideContributions,
    //     const NDData<int>& rEqIdsCsrIndices,
    //     TCsrMatrixType& rLeftHandSide)
    // {
    //     // Get and check the provided data shapes
    //     const auto& shape = rEqIdsCsrIndices.Shape();
    //     const std::size_t n_entities = shape[0];
    //     const std::size_t local_size_1 = shape[1];
    //     const std::size_t local_size_2 = shape[2];
    //     KRATOS_ERROR_IF(shape.size() != rLeftHandSideContributions.Shape().size()) << "The shapes of the equation ids and left hand side contributions do not match." << std::endl;
    //     for (std::size_t i = 0; i < shape.size(); ++i) {
    //         KRATOS_ERROR_IF(shape[i] != rLeftHandSideContributions.Shape()[i]) << "The shapes of the equation ids and left hand side contributions do not match in component " << i << "." << std::endl;
    //     }

    //     // Get the left hand side data
    //     auto& r_lhs_data = rLeftHandSide.value_data();
    //     const auto& r_idx_data = rEqIdsCsrIndices.ViewData();
    //     const auto& r_lhs_contribution_data = rLeftHandSideContributions.ViewData();

    //     // Loop over the entities
    //     IndexPartition<std::size_t>(n_entities).for_each([&](std::size_t i) {
    //         const std::size_t entity_pos = i * (local_size_1 * local_size_2);
    //         for (std::size_t i_local = 0; i_local < local_size_1; ++i_local) {
    //             for (std::size_t j_local = 0; j_local < local_size_2; ++j_local) {
    //                 const std::size_t aux_idx = entity_pos + i_local * local_size_1 + j_local; // Position in the contributions and equation ids data
    //                 const int csr_index = r_idx_data[aux_idx]; // Index in the CSR matrix values vector
    //                 const double lhs_contribution = r_lhs_contribution_data[aux_idx]; // Scalar contribution to the left hand side
    //                 r_lhs_data[csr_index] += lhs_contribution;
    //             }
    //         }
    //     });
    // }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Friends
    ///@{


    ///@}
}; // Class CsrUtilities

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
