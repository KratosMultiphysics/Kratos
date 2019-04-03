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

        unsigned int Dim = rModelPart.GetProcessInfo()[DOMAIN_SIZE];

        // Clear nodal variables
        #pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

            for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
            {
                itNode->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
                Matrix& rNodalStress = itNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                if(rNodalStress.size1() != Dim)
                    rNodalStress.resize(Dim,Dim,false);
                noalias(rNodalStress) = ZeroMatrix(Dim,Dim);
                itNode->FastGetSolutionStepValue(NODAL_JOINT_AREA) = 0.0;
                itNode->FastGetSolutionStepValue(NODAL_JOINT_WIDTH) = 0.0;
            }
        }

        BaseType::FinalizeSolutionStep(rModelPart,A,Dx,b);

        // Compute smoothed nodal variables
        #pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

            for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
            {
                const double& NodalArea = itNode->FastGetSolutionStepValue(NODAL_AREA);
                if (NodalArea>1.0e-20)
                {
                    const double InvNodalArea = 1.0/(NodalArea);
                    Matrix& rNodalStress = itNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                    for(unsigned int i = 0; i<Dim; i++)
                    {
                        for(unsigned int j = 0; j<Dim; j++)
                        {
                            rNodalStress(i,j) *= InvNodalArea;
                        }
                    }
                }

                const double& NodalJointArea = itNode->FastGetSolutionStepValue(NODAL_JOINT_AREA);
                if (NodalJointArea>1.0e-20)
                {
                    double& NodalJointWidth = itNode->FastGetSolutionStepValue(NODAL_JOINT_WIDTH);
                    NodalJointWidth = NodalJointWidth/NodalJointArea;
                }
            }
        }

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
