//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// #if !defined(KRATOS_RESIDUALBASED_INCREMENTALUPDATE_STATIC_SCHEME_H )
// #define  KRATOS_RESIDUALBASED_INCREMENTALUPDATE_STATIC_SCHEME_H
#pragma once
/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "includes/variables.h"
#include "utilities/entities_utilities.h"

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
///@}
///@name  Functions
///@{
///@}
///@name Kratos Classes
///@{

/**
 * @class ResidualBasedIncrementalUpdateStaticScheme
 * @ingroup KratosCore
 * @brief This class provides the implementation of a static scheme
 * @details The only operation done in this  scheme is the update of the database, no predict is done
 * @tparam TSparseSpace The sparse space considered
 * @tparam TDenseSpace The dense space considered
 * @see Scheme
 * @author Riccardo Rossi
 */
template<class TSparseSpace, class TDenseSpace>
class MPMResidualBasedIncrementalUpdateStaticScheme
    : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>
{

public:
    KRATOS_CLASS_POINTER_DEFINITION( MPMResidualBasedIncrementalUpdateStaticScheme );

    /// Base class definition
    typedef Scheme<TSparseSpace,TDenseSpace>                                       BaseType;

    // The current class definition
    typedef MPMResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace> ClassType;

    /// DoF array type definition
    typedef typename BaseType::DofsArrayType                                  DofsArrayType;

    /// Data type definition
    typedef typename BaseType::TDataType                                          TDataType;
    /// Matrix type definition
    typedef typename BaseType::TSystemMatrixType                          TSystemMatrixType;
    /// Vector type definition
    typedef typename BaseType::TSystemVectorType                          TSystemVectorType;
    /// Local system matrix type definition
    typedef typename BaseType::LocalSystemVectorType                  LocalSystemVectorType;
    /// Local system vector type definition
    typedef typename BaseType::LocalSystemMatrixType                  LocalSystemMatrixType;

    /// Elements containers definition
    typedef ModelPart::ElementsContainerType                              ElementsArrayType;
    /// Conditions containers definition
    typedef ModelPart::ConditionsContainerType                          ConditionsArrayType;

    /// The definition of the vector containing the equation ids
    typedef Element::EquationIdVectorType                              EquationIdVectorType;


    /**
     * @brief Performing the prediction of the solution.
     * @param rModelPart The model part of the problem to solve
     * @param rA LHS matrix
     * @param rDx Incremental update of primary variables
     * @param rb RHS Vector
     */
    void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY
        
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        const auto &r_elements_array = rModelPart.Elements();
        const std::size_t n_elems = r_elements_array.size();
        IndexPartition<std::size_t>(n_elems).for_each([&](std::size_t i_elem) {
            auto it_elem = r_elements_array.begin() + i_elem;

            it_elem->AddExplicitContribution(r_current_process_info);
        });

        KRATOS_CATCH("")
    }

}; // Class ResidualBasedIncrementalUpdateStaticScheme
}  // namespace Kratos

// #endif /* KRATOS_RESIDUALBASED_INCREMENTALUPDATE_STATIC_SCHEME_H  defined */
