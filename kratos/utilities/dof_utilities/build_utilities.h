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
#include "utilities/openmp_utils.h"

namespace Kratos
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class BuildUtilities
 * @ingroup KratosCore
 * @brief Utility class for handling the build
 * @details This static class collects the methods required to
 * handle the building of the sparse sytem matrices
 * @author Ruben Zorrilla
 */
template<class TDataType = double, class TIndexType=std::size_t, class TMatrixType=CsrMatrix<TDataType, TIndexType>, class TVectorType=SystemVector<TDataType, TIndexType>>
class BuildUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BuildUtilities
    KRATOS_CLASS_POINTER_DEFINITION(BuildUtilities);

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
    BuildUtilities() = delete;

    ///@}
    ///@name Operations
    ///@{

    static void Build(
        const ModelPart& rModelPart,
        SystemMatrixType& rLHS,
        SystemVectorType& rRHS)
    {

        // Getting model part data
        const SizeType n_elem = rModelPart.NumberOfElements();
        const SizeType n_cond = rModelPart.NumberOfConditions();
        auto el_begin = rModelPart.ElementsBegin();
        auto cond_begin = rModelPart.ConditionsBegin();
        const auto& r_process_info = rModelPart.GetProcessInfo();

        // Arrays for the local contributions to the system
        LocalSystemMatrixType loc_lhs = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType loc_rhs = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        Element::EquationIdVectorType eq_ids;

        KRATOS_WATCH("Before building elements and conditions")

        // Assemble elements and conditions
        #pragma omp parallel firstprivate(n_elem, n_cond, loc_lhs, loc_rhs, eq_ids, r_process_info)
        {
            // Assemble elements
            # pragma omp for schedule(guided, 512) nowait
            for (int k = 0; k < n_elem; ++k) {
                auto it_elem = el_begin + k;
                if (it_elem->IsActive()) {
                    // Calculate local LHS and RHS contributions
                    it_elem->CalculateLocalSystem(loc_lhs, loc_rhs, r_process_info);

                    // Get the positions in the global system
                    it_elem->EquationIdVector(eq_ids, r_process_info);

                    // Assemble the local contributions to the global system
                    rRHS.Assemble(loc_rhs, eq_ids);
                    rLHS.Assemble(loc_lhs, eq_ids);
                }
            }

            // Assemble conditions
            #pragma omp for schedule(guided, 512)
            for (int k = 0; k < n_cond; ++k) {
                auto it_cond = cond_begin + k;
                if (it_cond->IsActive()) {
                    // Calculate local LHS and RHS contributions
                    it_cond->CalculateLocalSystem(loc_lhs, loc_rhs, r_process_info);

                    // Get the positions in the global system
                    it_cond->EquationIdVector(eq_ids, r_process_info);

                    // Assemble the local contributions to the global system
                    rRHS.Assemble(loc_rhs, eq_ids);
                    rLHS.Assemble(loc_lhs, eq_ids);
                }
            }
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
}; // Class BuildUtilities

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
