//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Andi Katili, Nicolò Crescenzio 
//

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
 * @class MPMResidualBasedIncrementalUpdateStaticScheme
 * @brief This class provides the implementation of a static scheme for MPM
 * @details The only additional operation done in this scheme is predict is done, since P2G mapping is moved to predict
 * @tparam TSparseSpace The sparse space considered
 * @tparam TDenseSpace The dense space considered
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

    /// Matrix type definition
    typedef typename BaseType::TSystemMatrixType                          TSystemMatrixType;
    /// Vector type definition
    typedef typename BaseType::TSystemVectorType                          TSystemVectorType;


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