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

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    BossakDisplacementSmoothingScheme(double rAlpham = 0.0)
        : ResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace>(rAlpham) {}
    
    //------------------------------------------------------------------------------------
    
    ///Destructor
    virtual ~BossakDisplacementSmoothingScheme() {}
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY
        
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
                rNodalStress.clear();
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
                const double InvNodalArea = 1.0/(itNode->FastGetSolutionStepValue(NODAL_AREA));
                Matrix& rNodalStress = itNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                const unsigned int Dim = rNodalStress.size1();
                for(unsigned int i = 0; i<Dim; i++)
                {
                    for(unsigned int j = 0; j<Dim; j++)
                    {
                        rNodalStress(i,j) *= InvNodalArea;
                    }
                }
            }
        }
                
        KRATOS_CATCH("")
    }

}; // Class BossakDisplacementSmoothingScheme
}  // namespace Kratos

#endif // KRATOS_BOSSAK_DISPLACEMENT_SMOOTHING_SCHEME defined
