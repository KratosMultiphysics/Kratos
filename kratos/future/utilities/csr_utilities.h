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

    template<class TContainerType, class TCsrMatrixType>
    static void GetEquationIdCsrIndices(
        const TContainerType& rContainer,
        const ProcessInfo& rProcessInfo,
        const TCsrMatrixType& rCsrMatrix,
        NDData<int>& rNDData)
    {
        KRATOS_WATCH("B")
        // Get equation ids size from the first entity (assuming all entities have the same size)
        EquationIdVectorType equation_ids;
        rContainer.begin()->EquationIdVector(equation_ids, rProcessInfo);
        const std::size_t local_size = equation_ids.size();

        // Assign the input NDData to have shape: number of entities * local_size * local_size
        KRATOS_WATCH("C")
        DenseVector<unsigned int> nd_data_shape(3);
        nd_data_shape[0] = rContainer.size();
        nd_data_shape[1] = local_size;
        nd_data_shape[2] = local_size;
        KRATOS_WATCH(nd_data_shape)
        rNDData = NDData<int>(nd_data_shape);
        KRATOS_WATCH("D")

        //TODO: Parallelism?
        // Loop over the container
        auto data_view = rNDData.ViewData();
        for (auto it = rContainer.begin(); it != rContainer.end(); ++it) {
            // Get current entity equation ids
            EquationIdVectorType equation_ids;
            it->EquationIdVector(equation_ids, rProcessInfo);

            // Get local size
            const unsigned int local_size = equation_ids.size();

            // Loop over the local size
            for (unsigned int i_local = 0; i_local < local_size; ++i_local) {
                const unsigned int i_global = equation_ids[i_local];
                for(unsigned int j_local = 0; j_local < local_size; ++j_local) {
                    const unsigned int j_global = equation_ids[j_local];
                    const unsigned int k = rCsrMatrix.FindValueIndex(i_global,j_global);
                    data_view[it->Id() * (local_size * local_size) + i_local * local_size + j_local] = k;
                }
            }
        }
    }

    template<class TCsrMatrixType>
    static void Assemble(
        const NDData<double>& rLeftHandSideContributions,
        const NDData<int>& rEqIdsCsrIndices,
        TCsrMatrixType& rLeftHandSide)
    {
        const auto& shape = rEqIdsCsrIndices.Shape();
        const std::size_t n_entities = shape[0];
        const std::size_t local_size_1 = shape[1];
        const std::size_t local_size_2 = shape[2];

        //TODO: Parallelism?
        auto& r_lhs_data = rLeftHandSide.value_data();
        const auto& r_idx_data = rEqIdsCsrIndices.ViewData();
        const auto& r_lhs_contribution_data = rLeftHandSideContributions.ViewData();
        for (std::size_t i = 0; i < n_entities; ++i) {
            for (std::size_t j = 0; j < local_size_1; ++j) {
                for (std::size_t k = 0; k < local_size_2; ++k) {
                    std::size_t aux_idx = i * (local_size_1 * local_size_2) + j * local_size_1 + k;
                    const int csr_index = r_idx_data[aux_idx];
                    const double lhs_contribution = r_lhs_contribution_data[aux_idx];
                    r_lhs_data[csr_index] += lhs_contribution;
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
