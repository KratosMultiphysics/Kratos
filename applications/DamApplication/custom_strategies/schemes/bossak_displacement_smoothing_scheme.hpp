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
#include "custom_strategies/schemes/dam_nodal_stress_smoothing_utilities.hpp"
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
        block_for_each(r_model_part.Nodes(), [](Node& rNode){
            auto& r_initial_stress = rNode.FastGetSolutionStepValue(INITIAL_STRESS_TENSOR);
            if (r_initial_stress.size1() != 3 || r_initial_stress.size2() != 3) {
                r_initial_stress.resize(3,3,false);
            }
            r_initial_stress.clear();
        });

        mSchemeIsInitialized = true;

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

        // Single-owner nodal Cauchy-stress smoothing lifecycle:
        // reset -> element/condition finalization -> (process-based
        // accumulation when enabled) -> normalization. See
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
