//
//   Project Name:        KratosDamApplication   $
//   Last Modified by:    $Author:Lorenzo Gracia $
//   Date:                $Date:    October 2016 $
//   Revision:            $Revision:         1.0 $
//

#if !defined(KRATOS_BOSSAK_DISPLACEMENT_SMOOTHING_SCHEME )
#define  KRATOS_BOSSAK_DISPLACEMENT_SMOOTHING_SCHEME

// Application includes
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
#include "custom_utilities/dam_initial_stress_utilities.hpp"
#include "custom_utilities/dam_nodal_stress_smoothing_utilities.hpp"
#include "custom_utilities/dam_nonlocal_damage_utilities.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class BossakDisplacementSmoothingScheme : public ResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( BossakDisplacementSmoothingScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>     BaseType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    using Scheme<TSparseSpace,TDenseSpace>::mSchemeIsInitialized;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    BossakDisplacementSmoothingScheme(double rAlpham = 0.0, double rayleigh_m = 0.0, double rayleigh_k = 0.0)
        : ResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace>(rAlpham)

    {

        mRayleighAlpha = rayleigh_m;
        mRayleighBeta = rayleigh_k;

    }

    //------------------------------------------------------------------------------------

    ///Destructor
    virtual ~BossakDisplacementSmoothingScheme() {}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize(ModelPart& r_model_part) override
    {
        KRATOS_TRY

        r_model_part.GetProcessInfo()[RAYLEIGH_ALPHA] = mRayleighAlpha;
        r_model_part.GetProcessInfo()[RAYLEIGH_BETA] = mRayleighBeta;

        // Initialize INITIAL_STRESS_TENSOR
        DamInitialStressUtilities::EnsureInitialStressTensorStorage(r_model_part);

        mSchemeIsInitialized = true;

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

        // 2. Scheme-owned LOCAL equivalent-strain production (nonlocal damage).
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

        // 2. Post-update LOCAL equivalent-strain production (the dynamic scheme
        //    has already updated the state).
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

protected:

    //Member variables
    double mRayleighAlpha;
    double mRayleighBeta;


}; // Class BossakDisplacementSmoothingScheme
}  // namespace Kratos

#endif // KRATOS_BOSSAK_DISPLACEMENT_SMOOTHING_SCHEME defined
