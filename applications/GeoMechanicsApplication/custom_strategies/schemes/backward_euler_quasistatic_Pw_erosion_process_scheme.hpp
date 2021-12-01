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

#if !defined(KRATOS_BACKWARD_EULER_QUASISTATIC_PW_EROSION_PROCESS_SCHEME )
#define  KRATOS_BACKWARD_EULER_QUASISTATIC_PW_EROSION_PROCESS_SCHEME

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

class BackwardEulerQuasistaticPwErosionProcessScheme : public BackwardEulerQuasistaticPwScheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( BackwardEulerQuasistaticPwErosionProcessScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>          BaseType;
    typedef typename BaseType::DofsArrayType          DofsArrayType;
    typedef typename BaseType::TSystemMatrixType      TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType      TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType  LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType  LocalSystemMatrixType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    BackwardEulerQuasistaticPwErosionProcessScheme() :
        BackwardEulerQuasistaticPwScheme<TSparseSpace,TDenseSpace>()
    {}

    //------------------------------------------------------------------------------------

    ///Destructor
    ~BackwardEulerQuasistaticPwErosionProcessScheme() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY
        // first the function from the mother class
        NewmarkQuasistaticUPwScheme::Initialize(rModelPart);

        // VG: TODO:
        // store the pointer/reference to the first pile element

        KRATOS_CATCH("")
    }

    void FinalizeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        // first the function from the mother class
        NewmarkQuasistaticUPwScheme::FinalizeNonLinIteration(rModelPart, A, Dx, b);

        // VG: TODO:
        // here there must be a loop over pile elements and check equilibrium in each element

        KRATOS_CATCH("")
    }

    /// Member Variables
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


}; // Class BackwardEulerQuasistaticPwErosionProcessScheme
}  // namespace Kratos

#endif // KRATOS_BACKWARD_EULER_QUASISTATIC_PW_EROSION_PROCESS_SCHEME defined
