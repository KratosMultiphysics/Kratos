// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#if !defined(KRATOS_BACKWARD_EULER_QUASISTATIC_PW_SCHEME )
#define  KRATOS_BACKWARD_EULER_QUASISTATIC_PW_SCHEME

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/schemes/newmark_quasistatic_Pw_scheme.hpp"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class BackwardEulerQuasistaticPwScheme : public NewmarkQuasistaticPwScheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( BackwardEulerQuasistaticPwScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>          BaseType;
    typedef typename BaseType::DofsArrayType          DofsArrayType;
    typedef typename BaseType::TSystemMatrixType      TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType      TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType  LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType  LocalSystemMatrixType;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mDeltaTime;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    BackwardEulerQuasistaticPwScheme() :
        NewmarkQuasistaticPwScheme<TSparseSpace,TDenseSpace>(1.0)
    {}

    //------------------------------------------------------------------------------------

    ///Destructor
    ~BackwardEulerQuasistaticPwScheme() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        //Update DtPressure

        block_for_each(rModelPart.Nodes(), [&](Node& rNode){
            const double DeltaPressure =  rNode.FastGetSolutionStepValue(WATER_PRESSURE)
                                        - rNode.FastGetSolutionStepValue(WATER_PRESSURE, 1);
            rNode.FastGetSolutionStepValue(DT_WATER_PRESSURE) = DeltaPressure / mDeltaTime;
        });

        KRATOS_CATCH( "" )
    }

}; // Class BackwardEulerQuasistaticPwScheme
}  // namespace Kratos

#endif // KRATOS_BACKWARD_EULER_QUASISTATIC_PW_SCHEME defined
