//
//   Project Name:        KratosDamApplication   $
//   Last Modified by:    $Author:Lorenzo Gracia $
//   Date:                $Date:    October 2016 $
//   Revision:            $Revision:         1.0 $
//

#if !defined(KRATOS_INCREMENTAL_UPDATE_STATIC_SMOOTHING_SCHEME )
#define  KRATOS_INCREMENTAL_UPDATE_STATIC_SMOOTHING_SCHEME

// Application includes
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "custom_utilities/dam_initial_stress_utilities.hpp"
#include "custom_utilities/dam_nodal_stress_smoothing_utilities.hpp"
#include "custom_utilities/dam_nonlocal_damage_utilities.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class IncrementalUpdateStaticSmoothingScheme : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( IncrementalUpdateStaticSmoothingScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>     BaseType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    IncrementalUpdateStaticSmoothingScheme()
        : ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>() {}

    //------------------------------------------------------------------------------------

    ///Destructor
    virtual ~IncrementalUpdateStaticSmoothingScheme() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize(ModelPart& r_model_part) override
    {
        KRATOS_TRY

        // Initialize INITIAL_STRESS_TENSOR
        DamInitialStressUtilities::EnsureInitialStressTensorStorage(r_model_part);

        BaseType::mSchemeIsInitialized = true;

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        // 1. Normal element/condition nonlinear-iteration callbacks.
        BaseType::InitializeNonLinIteration(rModelPart, A, Dx, b);

        // 2. When the scheme owns the LOCAL equivalent-strain production
        //    (nonlocal damage active), recompute the current LOCAL quantities
        //    before control returns to the Poromechanics nonlocal strategy,
        //    which immediately performs the nonlocal averaging.
        if (DamNonlocalDamageUtilities::IsProcessBasedLocalOwnership(rModelPart)) {
            DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(rModelPart);
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void FinalizeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        // 1. Normal element/condition nonlinear-iteration callbacks.
        BaseType::FinalizeNonLinIteration(rModelPart, A, Dx, b);

        // 2. The strategy has already updated the displacement state (and
        //    optionally moved the mesh): recompute the current LOCAL quantities
        //    at the post-update state. The Poromechanics strategy performs the
        //    second nonlocal averaging immediately after this hook.
        if (DamNonlocalDamageUtilities::IsProcessBasedLocalOwnership(rModelPart)) {
            DamNonlocalDamageUtilities::CalculateLocalEquivalentStrain(rModelPart);
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        // Single-owner nodal Cauchy-stress smoothing lifecycle
        // (reset -> finalization -> accumulation -> normalization); see
        // DamNodalStressSmoothingUtilities.
        DamNodalStressSmoothingUtilities::FinalizeSolutionStep(*this, rModelPart, A, Dx, b);

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

}; // Class IncrementalUpdateStaticSmoothingScheme
}  // namespace Kratos

#endif // KRATOS_INCREMENTAL_UPDATE_STATIC_SMOOTHING_SCHEME defined
