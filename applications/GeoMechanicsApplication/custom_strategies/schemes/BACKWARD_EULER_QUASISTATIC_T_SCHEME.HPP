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

#if !defined(KRATOS_BACKWARD_EULER_QUASISTATIC_T_SCHEME )
#define  KRATOS_BACKWARD_EULER_QUASISTATIC_T_SCHEME

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

    typedef Scheme<TSparseSpace,TDenseSpace>          BaseType;
    typedef typename BaseType::DofsArrayType          DofsArrayType;
    typedef typename BaseType::TSystemMatrixType      TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType      TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType  LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType  LocalSystemMatrixType;
    //using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mDeltaTime;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    BackwardEulerQuasistaticTScheme() :
        NewmarkQuasistaticTScheme<TSparseSpace,TDenseSpace>(1.0)
    {}

    //------------------------------------------------------------------------------------

    ///Destructor
    ~BackwardEulerQuasistaticTScheme() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        //Update DtTemperature

        block_for_each(rModelPart.Nodes(), [&](Node<3>& rNode){
            const double DeltaTemperature =  rNode.FastGetSolutionStepValue(TEMPERATURE)
                                           - rNode.FastGetSolutionStepValue(TEMPERATURE, 1);
            rNode.FastGetSolutionStepValue(DT_TEMPERATURE) = DeltaTemperature / mDeltaTime;
        });

        KRATOS_CATCH( "" )
    }

}; // Class BackwardEulerQuasistaticTScheme
}  // namespace Kratos

#endif // KRATOS_BACKWARD_EULER_QUASISTATIC_T_SCHEME defined
