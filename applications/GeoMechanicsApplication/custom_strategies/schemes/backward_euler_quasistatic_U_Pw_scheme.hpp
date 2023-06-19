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

#if !defined(KRATOS_BACKWARD_EULER_QUASISTATIC_U_PW_SCHEME )
#define  KRATOS_BACKWARD_EULER_QUASISTATIC_U_PW_SCHEME

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "solving_strategies/schemes/scheme.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class BackwardEulerQuasistaticUPwScheme : public NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( BackwardEulerQuasistaticUPwScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>          BaseType;
    typedef typename BaseType::DofsArrayType          DofsArrayType;
    typedef typename BaseType::TSystemMatrixType      TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType      TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType  LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType  LocalSystemMatrixType;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mDeltaTime;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mBeta;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mGamma;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mTheta;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    BackwardEulerQuasistaticUPwScheme() :
        NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>(1.0, 1.0, 1.0)
    {
    }

    //------------------------------------------------------------------------------------

    ///Destructor
    ~BackwardEulerQuasistaticUPwScheme() override {}

protected:

    /// Member Variables
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    inline void SetTimeFactors(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        mDeltaTime = rModelPart.GetProcessInfo()[DELTA_TIME];
        rModelPart.GetProcessInfo()[VELOCITY_COEFFICIENT] = 1.0/mDeltaTime;
        rModelPart.GetProcessInfo()[DT_PRESSURE_COEFFICIENT] = 1.0/mDeltaTime;

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        //Update Acceleration, Velocity and DtPressure

        block_for_each(rModelPart.Nodes(), [&](Node& rNode) {

            noalias(rNode.FastGetSolutionStepValue(VELOCITY))     = (  rNode.FastGetSolutionStepValue(DISPLACEMENT)
                                                                     - rNode.FastGetSolutionStepValue(DISPLACEMENT, 1)) / mDeltaTime;

            noalias(rNode.FastGetSolutionStepValue(ACCELERATION)) = (  rNode.FastGetSolutionStepValue(VELOCITY)
                                                                     - rNode.FastGetSolutionStepValue(VELOCITY,1) ) / mDeltaTime;

            rNode.FastGetSolutionStepValue(DT_WATER_PRESSURE) = (  rNode.FastGetSolutionStepValue(WATER_PRESSURE)
                                                                 - rNode.FastGetSolutionStepValue(WATER_PRESSURE, 1)) / mDeltaTime;
        });

        KRATOS_CATCH( "" )
    }

}; // Class BackwardEulerQuasistaticUPwScheme
}  // namespace Kratos

#endif // KRATOS_BACKWARD_EULER_QUASISTATIC_U_PW_SCHEME defined
