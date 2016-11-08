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
    typedef typename BaseType::DofsArrayType     DofsArrayType;
    typedef ModelPart::NodesContainerType        NodesArrayType;
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    IncrementalUpdateStaticSmoothingScheme()
        : ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>() {}
    
    //------------------------------------------------------------------------------------
    
    ///Destructor
    virtual ~IncrementalUpdateStaticSmoothingScheme() {}
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY;

        Scheme<TSparseSpace,TDenseSpace>::InitializeSolutionStep(rModelPart, A, Dx, b);

        double DeltaTime = rModelPart.GetProcessInfo()[DELTA_TIME];
        
        // Pure Newmark scheme 
        // Note: this is a quasistatic scheme. The velocity and acceleration are stored without affecting the displacement
        double beta = 0.25;
        double gamma = 0.5;

        // Initializing Newmark constants
        mNewmark0 = ( 1.0 / (beta * DeltaTime * DeltaTime) );
        mNewmark1 = ( gamma / (beta * DeltaTime) );
        mNewmark2 = ( 1.0 / (beta * DeltaTime) );
        mNewmark3 = ( 0.5 / (beta) - 1.0 );
        mNewmark4 = ( (gamma / beta) - 1.0  );
        mNewmark5 = ( DeltaTime * 0.5 * ( ( gamma / beta ) - 2.0 ) );

        KRATOS_CATCH( "" );
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Update(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b )
    {
        KRATOS_TRY;

        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

        // Update of displacement (by DOF)
        OpenMPUtils::PartitionVector DofPartition;
        OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, DofPartition);

        const int ndof = static_cast<int>(rDofSet.size());
        typename DofsArrayType::iterator DofBegin = rDofSet.begin();

        #pragma omp parallel for firstprivate(DofBegin)
        for(int i = 0;  i < ndof; i++)
        {
            typename DofsArrayType::iterator itDof = DofBegin + i;

            if (itDof->IsFree() )
            {
                itDof->GetSolutionStepValue() += Dx[itDof->EquationId()];
            }
        }

        // Updating time derivatives (nodally for efficiency)
        OpenMPUtils::PartitionVector NodePartition;
        OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

        const int nnodes = static_cast<int>(rModelPart.Nodes().size());
        NodesArrayType::iterator NodeBegin = rModelPart.Nodes().begin();

        #pragma omp parallel for firstprivate(NodeBegin)
        for(int i = 0;  i < nnodes; i++)
        {
            array_1d<double, 3 > DeltaDisplacement;

            NodesArrayType::iterator itNode = NodeBegin + i;

            noalias(DeltaDisplacement) = (itNode)->FastGetSolutionStepValue(DISPLACEMENT) - (itNode)->FastGetSolutionStepValue(DISPLACEMENT, 1);

            array_1d<double, 3 > & CurrentVelocity            = (itNode)->FastGetSolutionStepValue(VELOCITY, 0);
            const array_1d<double, 3 > & PreviousVelocity     = (itNode)->FastGetSolutionStepValue(VELOCITY, 1);

            array_1d<double, 3 > & CurrentAcceleration        = (itNode)->FastGetSolutionStepValue(ACCELERATION, 0);
            const array_1d<double, 3 > & PreviousAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 1);

            this->UpdateVelocity     (CurrentVelocity,     DeltaDisplacement, PreviousVelocity, PreviousAcceleration);

            this->UpdateAcceleration (CurrentAcceleration, DeltaDisplacement, PreviousVelocity, PreviousAcceleration);
        }

        KRATOS_CATCH( "" );
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY;

        //const double DeltaTime = rModelPart.GetProcessInfo()[DELTA_TIME];

        // Updating time derivatives (nodally for efficiency)
        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector NodePartition;
        OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

        const int nnodes = static_cast<int>( rModelPart.Nodes().size() );
        NodesArrayType::iterator NodeBegin = rModelPart.Nodes().begin();

        #pragma omp parallel for firstprivate(NodeBegin)
        for(int i = 0;  i< nnodes; i++)
        {
            array_1d<double, 3 > DeltaDisplacement;

            NodesArrayType::iterator itNode = NodeBegin + i;

            const array_1d<double, 3 > & PreviousAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 1);
            const array_1d<double, 3 > & PreviousVelocity     = (itNode)->FastGetSolutionStepValue(VELOCITY,     1);
            const array_1d<double, 3 > & PreviousDisplacement = (itNode)->FastGetSolutionStepValue(DISPLACEMENT, 1);
            array_1d<double, 3 > & CurrentAcceleration        = (itNode)->FastGetSolutionStepValue(ACCELERATION, 0);
            array_1d<double, 3 > & CurrentVelocity            = (itNode)->FastGetSolutionStepValue(VELOCITY,     0);
            array_1d<double, 3 > & CurrentDisplacement        = (itNode)->FastGetSolutionStepValue(DISPLACEMENT, 0);

            // Updating time derivatives ::: Please note that displacements and its time derivatives can not be consistently fixed separately
            noalias(DeltaDisplacement) = CurrentDisplacement - PreviousDisplacement;

            this->UpdateVelocity     (CurrentVelocity,     DeltaDisplacement, PreviousVelocity, PreviousAcceleration);

            this->UpdateAcceleration (CurrentAcceleration, DeltaDisplacement, PreviousVelocity, PreviousAcceleration);
        }

        KRATOS_CATCH( "" );
    }

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

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
    
    //Member variables
    double mNewmark0,mNewmark1,mNewmark2,mNewmark3,mNewmark4,mNewmark5;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    inline void UpdateVelocity(
        array_1d<double, 3 > & CurrentVelocity,
        const array_1d<double, 3 > & DeltaDisplacement,
        const array_1d<double, 3 > & PreviousVelocity,
        const array_1d<double, 3 > & PreviousAcceleration
    )
    {
        noalias(CurrentVelocity) =  (mNewmark1 * DeltaDisplacement - mNewmark4 * PreviousVelocity
                                     - mNewmark5 * PreviousAcceleration);
    }

    //----------------------------------------------------------------------------------------------------------------

    inline void UpdateAcceleration(
        array_1d<double, 3 > & CurrentAcceleration,
        const array_1d<double, 3 > & DeltaDisplacement,
        const array_1d<double, 3 > & PreviousVelocity,
        const array_1d<double, 3 > & PreviousAcceleration
    )
    {
        noalias(CurrentAcceleration) =  (mNewmark0 * DeltaDisplacement - mNewmark2 * PreviousVelocity
                                         -  mNewmark3 * PreviousAcceleration);
    }

}; // Class IncrementalUpdateStaticSmoothingScheme
}  // namespace Kratos

#endif // KRATOS_INCREMENTAL_UPDATE_STATIC_SMOOTHING_SCHEME defined
