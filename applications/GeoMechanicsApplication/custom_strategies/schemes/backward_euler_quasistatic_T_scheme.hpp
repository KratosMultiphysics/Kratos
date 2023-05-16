// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//

#pragma once

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/schemes/newmark_quasistatic_T_scheme.hpp"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class BackwardEulerQuasistaticTScheme : public NewmarkQuasistaticTScheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( BackwardEulerQuasistaticTScheme );

    using BaseType = Scheme<TSparseSpace, TDenseSpace>;
    using DofsArrayType = typename BaseType::DofsArrayType;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;
    using NewmarkQuasistaticTScheme<TSparseSpace,TDenseSpace>::mDeltaTime;

    ///Constructor
	// ============================================================================================
    // ============================================================================================
    BackwardEulerQuasistaticTScheme() :
        NewmarkQuasistaticTScheme<TSparseSpace,TDenseSpace>(1.0)
    {}

    ///Destructor
    // ============================================================================================
    // ============================================================================================
    ~BackwardEulerQuasistaticTScheme() override {}


protected:

    /// Member Variables
    // ============================================================================================
    // ============================================================================================
    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        //Update DtTemperature

        block_for_each(rModelPart.Nodes(), [&](Node& rNode){
            const double DeltaTemperature =  rNode.FastGetSolutionStepValue(TEMPERATURE)
                                           - rNode.FastGetSolutionStepValue(TEMPERATURE, 1);
            rNode.FastGetSolutionStepValue(DT_TEMPERATURE) = DeltaTemperature / mDeltaTime;
        });

        KRATOS_CATCH( "" )
    }

}; // Class BackwardEulerQuasistaticTScheme
}  // namespace Kratos
