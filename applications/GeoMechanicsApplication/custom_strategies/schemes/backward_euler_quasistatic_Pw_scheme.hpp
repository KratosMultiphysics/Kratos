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
#include "solving_strategies/schemes/scheme.h"

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

            double& CurrentDtPressure = itNode->FastGetSolutionStepValue(DT_WATER_PRESSURE);
            DeltaPressure = itNode->FastGetSolutionStepValue(WATER_PRESSURE) - itNode->FastGetSolutionStepValue(WATER_PRESSURE, 1);
            CurrentDtPressure = DeltaPressure / mDeltaTime;
        }

        KRATOS_CATCH( "" )
    }

}; // Class BackwardEulerQuasistaticPwScheme
}  // namespace Kratos

#endif // KRATOS_BACKWARD_EULER_QUASISTATIC_PW_SCHEME defined
