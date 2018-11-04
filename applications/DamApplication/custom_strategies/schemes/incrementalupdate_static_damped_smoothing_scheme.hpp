//
//   Project Name:        KratosDamApplication   $
//   Last Modified by:    $Author:Lorenzo Gracia $
//   Date:                $Date:    November 2016$
//   Revision:            $Revision:         1.0 $
//

#if !defined(KRATOS_INCREMENTAL_UPDATE_STATIC_DAMPED_SMOOTHING_SCHEME )
#define  KRATOS_INCREMENTAL_UPDATE_STATIC_DAMPED_SMOOTHING_SCHEME

// Application includes
//~ #include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "custom_strategies/schemes/incrementalupdate_static_smoothing_scheme.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class IncrementalUpdateStaticDampedSmoothingScheme : public IncrementalUpdateStaticSmoothingScheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( IncrementalUpdateStaticDampedSmoothingScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>     BaseType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::DofsArrayType     DofsArrayType;
    typedef ModelPart::NodesContainerType        NodesArrayType;
    using Scheme<TSparseSpace,TDenseSpace>::mSchemeIsInitialized;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    IncrementalUpdateStaticDampedSmoothingScheme(double rayleigh_m, double rayleigh_k)
        : IncrementalUpdateStaticSmoothingScheme<TSparseSpace,TDenseSpace>()

    {
        mRayleighAlpha = rayleigh_m;
        mRayleighBeta = rayleigh_k;

        //Allocate auxiliary memory
        int NumThreads = OpenMPUtils::GetNumThreads();
        mDampingMatrix.resize(NumThreads);
        mVelocityVector.resize(NumThreads);

    }

    //------------------------------------------------------------------------------------

    ///Destructor
    virtual ~IncrementalUpdateStaticDampedSmoothingScheme() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(ModelPart& r_model_part) override
    {
        KRATOS_TRY

        int ierr = Scheme<TSparseSpace,TDenseSpace>::Check(r_model_part);
        if(ierr != 0) return ierr;

        if ( RAYLEIGH_ALPHA.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "RAYLEIGH_ALPHA has Key zero! (check if the application is correctly registered", "" )
        if ( RAYLEIGH_BETA.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "RAYLEIGH_BETA has Key zero! (check if the application is correctly registered", "" )

        // Check rayleigh coefficients
        if( mRayleighAlpha < 0.0 || mRayleighBeta < 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument,"Some of the rayleigh coefficients has an invalid value ", "" )

        return 0;

        KRATOS_CATCH( "" )
    }


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize(ModelPart& r_model_part) override
    {
        KRATOS_TRY

        mDeltaTime = r_model_part.GetProcessInfo()[DELTA_TIME];

        // Pure Newmark scheme
        // Note: this is a quasistatic scheme. The velocity and acceleration are stored without affecting the displacement
        mbeta = 0.25;
        mgamma = 0.5;

        mNewmark0 = ( 1.0 / (mbeta * mDeltaTime * mDeltaTime) );
        mNewmark1 = ( mgamma / (mbeta * mDeltaTime) );
        mNewmark2 = ( 1.0 / (mbeta * mDeltaTime) );
        mNewmark3 = ( 0.5 / (mbeta) - 1.0 );
        mNewmark4 = ( (mgamma / mbeta) - 1.0  );
        mNewmark5 = ( mDeltaTime * 0.5 * ( ( mgamma / mbeta ) - 2.0 ) );

        r_model_part.GetProcessInfo()[RAYLEIGH_ALPHA] = mRayleighAlpha;
        r_model_part.GetProcessInfo()[RAYLEIGH_BETA] = mRayleighBeta;

        mSchemeIsInitialized = true;

        KRATOS_CATCH("")
    }


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
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

    void Update(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b ) override
    {
        KRATOS_TRY;

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector DofSetPartition;
        OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, DofSetPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            typename DofsArrayType::iterator DofsBegin = rDofSet.begin() + DofSetPartition[k];
            typename DofsArrayType::iterator DofsEnd = rDofSet.begin() + DofSetPartition[k+1];

            //Update Displacement and Pressure (DOFs)
            for (typename DofsArrayType::iterator itDof = DofsBegin; itDof != DofsEnd; ++itDof)
            {
                if (itDof->IsFree())
                    itDof->GetSolutionStepValue() += TSparseSpace::GetValue(Dx, itDof->EquationId());
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

    // Note: this is in a parallel loop

    void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        (rCurrentElement) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        (rCurrentElement) -> CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);

        this->AddDampingToLHS (LHS_Contribution, mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDampingToRHS (rCurrentElement, RHS_Contribution, mDampingMatrix[thread], CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void Calculate_RHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        (rCurrentElement) -> CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);

        this->AddDampingToRHS (rCurrentElement, RHS_Contribution, mDampingMatrix[thread], CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void Calculate_LHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        (rCurrentElement) -> CalculateLeftHandSide(LHS_Contribution,CurrentProcessInfo);

        (rCurrentElement) -> CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);

        this->AddDampingToLHS (LHS_Contribution, mDampingMatrix[thread], CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    //Member variables
    double mNewmark0,mNewmark1,mNewmark2,mNewmark3,mNewmark4,mNewmark5;

    double mDeltaTime;
    double mbeta;
    double mgamma;

    double mRayleighAlpha;
    double mRayleighBeta;

    std::vector< Matrix > mDampingMatrix;
    std::vector< Vector > mVelocityVector;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddDampingToLHS(LocalSystemMatrixType& LHS_Contribution,LocalSystemMatrixType& C,ProcessInfo& CurrentProcessInfo)
    {
        // adding damping contribution
        if (C.size1() != 0)
        {
            noalias(LHS_Contribution) += mgamma/(mbeta*mDeltaTime)*C;
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddDampingToRHS(Element::Pointer rCurrentElement,LocalSystemVectorType& RHS_Contribution, LocalSystemMatrixType& C,ProcessInfo& CurrentProcessInfo)
    {
        int thread = OpenMPUtils::ThisThread();

        //adding damping contribution
        if (C.size1() != 0)
        {
            rCurrentElement->GetFirstDerivativesVector(mVelocityVector[thread], 0);

            noalias(RHS_Contribution) -= prod(C, mVelocityVector[thread]);

        }
    }

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

}; // Class IncrementalUpdateStaticDampedSmoothingScheme
}  // namespace Kratos

#endif // KRATOS_INCREMENTAL_UPDATE_STATIC_DAMPED_SMOOTHING_SCHEME defined
