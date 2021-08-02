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
        NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>()
    {
        mBeta = 1.0;
        mGamma = 1.0;
        mTheta = 1.0;
    }

    //------------------------------------------------------------------------------------

    ///Destructor
    ~BackwardEulerQuasistaticUPwScheme() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    inline void UpdateVariablesDerivatives(ModelPart& r_model_part) override
    {
        KRATOS_TRY

        //Update Acceleration, Velocity and DtPressure

        array_1d<double,3> DeltaDisplacement;
        double DeltaPressure;

        const int NNodes = static_cast<int>(r_model_part.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = r_model_part.NodesBegin();

        #pragma omp parallel for private(DeltaDisplacement,DeltaPressure)
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNode = node_begin + i;

            array_1d<double,3>& CurrentAcceleration = itNode->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double,3>& CurrentVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
            noalias(DeltaDisplacement) = itNode->FastGetSolutionStepValue(DISPLACEMENT) - itNode->FastGetSolutionStepValue(DISPLACEMENT, 1);
            const array_1d<double,3>& PreviousVelocity = itNode->FastGetSolutionStepValue(VELOCITY, 1);

            noalias(CurrentVelocity)     = DeltaDisplacement / mDeltaTime;
            noalias(CurrentAcceleration) = ( CurrentVelocity - PreviousVelocity ) / mDeltaTime;

            double& CurrentDtPressure = itNode->FastGetSolutionStepValue(DT_WATER_PRESSURE);
            DeltaPressure = itNode->FastGetSolutionStepValue(WATER_PRESSURE) - itNode->FastGetSolutionStepValue(WATER_PRESSURE, 1);
            CurrentDtPressure = DeltaPressure / mDeltaTime;
        }

        KRATOS_CATCH( "" )
    }

}; // Class BackwardEulerQuasistaticUPwScheme
}  // namespace Kratos

#endif // KRATOS_BACKWARD_EULER_QUASISTATIC_U_PW_SCHEME defined
