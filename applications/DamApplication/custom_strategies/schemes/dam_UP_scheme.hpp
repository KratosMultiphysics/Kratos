//
//   Project Name:        			KratosDamApplication $
//   Last Modified by:    $Author:    	  Lorenzo Gracia $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_DAM_UP_SCHEME )
#define  KRATOS_DAM_UP_SCHEME

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "includes/element.h"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class DamUPScheme : public Scheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( DamUPScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;
    typedef typename BaseType::DofsArrayType                 DofsArrayType;
    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType         TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    DamUPScheme(double beta, double gamma, double rayleigh_m ,double rayleigh_k  ): Scheme<TSparseSpace,TDenseSpace>()
    {
        mBeta = beta;
        mGamma = gamma;
        mrayleigh_m = rayleigh_m;
        mrayleigh_k = rayleigh_k;


        //Allocate auxiliary memory
        int NumThreads = OpenMPUtils::GetNumThreads();
        mMassMatrix.resize(NumThreads);
        mAccelerationVector.resize(NumThreads);
        mDampingMatrix.resize(NumThreads);
        mVelocityVector.resize(NumThreads);
    }

    //------------------------------------------------------------------------------------

    ///Destructor
    virtual ~DamUPScheme() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(ModelPart& r_model_part) override
    {
        KRATOS_TRY

        //check for variables keys (verify that the variables are correctly initialized)
        if(DISPLACEMENT.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered", "" )
        if(VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"VELOCITY has Key zero! (check if the application is correctly registered", "" )
        if(ACCELERATION.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"ACCELERATION has Key zero! (check if the application is correctly registered", "" )
        if(PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument, "PRESSURE has Key zero! (check if the application is correctly registered", "" )
        if(Dt_PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument, "Dt_PRESSURE has Key zero! (check if the application is correctly registered", "" )
        if(Dt2_PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument, "Dt2_PRESSURE has Key zero! (check if the application is correctly registered", "" )
        if ( VELOCITY_PRESSURE_COEFFICIENT.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "VELOCITY_PRESSURE_COEFFICIENT has Key zero! (check if the application is correctly registered", "" )
        if ( ACCELERATION_PRESSURE_COEFFICIENT.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "ACCELERATION_PRESSURE_COEFFICIENT has Key zero! (check if the application is correctly registered", "" )

        //check that variables are correctly allocated
        for(ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); it++)
        {
			if(it->SolutionStepsDataHas(DISPLACEMENT) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_THROW_ERROR( std::logic_error, "VELOCITY variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(ACCELERATION) == false)
                KRATOS_THROW_ERROR( std::logic_error, "ACCELERATION variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(PRESSURE) == false)
                KRATOS_THROW_ERROR( std::logic_error, "PRESSURE variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(Dt_PRESSURE) == false)
                KRATOS_THROW_ERROR( std::logic_error, "Dt_PRESSURE variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(Dt2_PRESSURE) == false)
                KRATOS_THROW_ERROR( std::logic_error, "Dt2_PRESSURE variable is not allocated for node ", it->Id() )


            if(it->HasDofFor(DISPLACEMENT_X) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_X dof on node ",it->Id() )
            if(it->HasDofFor(DISPLACEMENT_Y) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_Y dof on node ",it->Id() )
            if(it->HasDofFor(DISPLACEMENT_Z) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_Z dof on node ",it->Id() )
            if(it->HasDofFor(PRESSURE) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing PRESSURE dof on node ",it->Id() )
        }

        //check for minimum value of the buffer index.
        if (r_model_part.GetBufferSize() < 2)
            KRATOS_THROW_ERROR( std::logic_error, "insufficient buffer size. Buffer size should be greater than 2. Current size is", r_model_part.GetBufferSize() )

        // Check beta, gamma and theta
        if(mBeta <= 0.0 || mGamma<= 0.0)
            KRATOS_THROW_ERROR( std::invalid_argument,"Some of the scheme variables: beta or  gamma has an invalid value ", "" )

        return 0;

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize(ModelPart& r_model_part) override
    {
        KRATOS_TRY

        mDeltaTime = r_model_part.GetProcessInfo()[DELTA_TIME];
        r_model_part.GetProcessInfo()[VELOCITY_PRESSURE_COEFFICIENT] = mGamma/(mBeta*mDeltaTime);
        r_model_part.GetProcessInfo()[ACCELERATION_PRESSURE_COEFFICIENT] = 1.0/(mBeta*mDeltaTime*mDeltaTime);

        r_model_part.GetProcessInfo()[RAYLEIGH_ALPHA] = mrayleigh_m;
        r_model_part.GetProcessInfo()[RAYLEIGH_BETA] = mrayleigh_k;

        KRATOS_WATCH(mrayleigh_m)
        KRATOS_WATCH(mrayleigh_k)

        BaseType::mSchemeIsInitialized = true;

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        mDeltaTime = r_model_part.GetProcessInfo()[DELTA_TIME];
        r_model_part.GetProcessInfo()[VELOCITY_PRESSURE_COEFFICIENT] = mGamma/(mBeta*mDeltaTime);
        r_model_part.GetProcessInfo()[ACCELERATION_PRESSURE_COEFFICIENT] = 1.0/(mBeta*mDeltaTime*mDeltaTime);

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        int NElems = static_cast<int>(r_model_part.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator itElem = el_begin + i;
            itElem -> InitializeSolutionStep(CurrentProcessInfo);
        }

        int NCons = static_cast<int>(r_model_part.Conditions().size());
        ModelPart::ConditionsContainerType::iterator con_begin = r_model_part.ConditionsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NCons; i++)
        {
            ModelPart::ConditionsContainerType::iterator itCond = con_begin + i;
            itCond -> InitializeSolutionStep(CurrentProcessInfo);
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Predict(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {

      KRATOS_TRY

        // Updating Velocities, Acclerations, DtPressure and Dt2Pressure
        array_1d<double,3> DeltaDisplacement;
        double DeltaPressure;

        const int NNodes = static_cast<int>(r_model_part.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = r_model_part.NodesBegin();

        #pragma omp parallel for private(DeltaDisplacement,DeltaPressure)
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNode = node_begin + i;

            // Termes related to displacement field

            array_1d<double,3>& CurrentDisplacement = itNode->FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double,3>& CurrentAcceleration = itNode->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double,3>& CurrentVelocity = itNode->FastGetSolutionStepValue(VELOCITY);

            const array_1d<double,3>& PreviousDisplacement = itNode->FastGetSolutionStepValue(DISPLACEMENT, 1);
            const array_1d<double,3>& PreviousAcceleration = itNode->FastGetSolutionStepValue(ACCELERATION, 1);
            const array_1d<double,3>& PreviousVelocity = itNode->FastGetSolutionStepValue(VELOCITY, 1);

            if (itNode -> IsFixed(ACCELERATION_X))
            {
                CurrentDisplacement[0] = PreviousDisplacement[0] + mDeltaTime * PreviousVelocity[0] + std::pow(mDeltaTime, 2) * ( ( 0.5 - mBeta) * PreviousAcceleration[0] + mBeta * CurrentAcceleration[0] );
            }
            else if (itNode -> IsFixed(VELOCITY_X))
            {
                CurrentDisplacement[0] = PreviousDisplacement[0] + mDeltaTime*(mBeta/mGamma*(CurrentVelocity[0]-PreviousVelocity[0])+PreviousVelocity[0]);
            }
            else if (itNode -> IsFixed(DISPLACEMENT_X) == false)
            {
                CurrentDisplacement[0] = PreviousDisplacement[0] + mDeltaTime * PreviousVelocity[0] + 0.5 * std::pow(mDeltaTime, 2) * PreviousAcceleration[0];
            }

            if (itNode -> IsFixed(ACCELERATION_Y))
            {
                CurrentDisplacement[1] = PreviousDisplacement[1] + mDeltaTime * PreviousVelocity[1] + std::pow(mDeltaTime, 2) * ( ( 0.5 - mBeta) * PreviousAcceleration[1] + mBeta * CurrentAcceleration[1] );
            }
            else if (itNode -> IsFixed(VELOCITY_Y))
            {
                CurrentDisplacement[1] = PreviousDisplacement[1] + mDeltaTime*(mBeta/mGamma*(CurrentVelocity[1]-PreviousVelocity[1])+PreviousVelocity[1]);
            }
            else if (itNode -> IsFixed(DISPLACEMENT_Y) == false)
            {
                CurrentDisplacement[1] = PreviousDisplacement[1] + mDeltaTime * PreviousVelocity[1] + 0.5 * std::pow(mDeltaTime, 2) * PreviousAcceleration[1];
            }

            // For 3D cases
            if (itNode -> HasDofFor(DISPLACEMENT_Z))
            {
                if (itNode -> IsFixed(ACCELERATION_Z))
                {
                    CurrentDisplacement[2] = PreviousDisplacement[2] + mDeltaTime * PreviousVelocity[2] + std::pow(mDeltaTime, 2) * ( ( 0.5 - mBeta) * PreviousAcceleration[2] + mBeta * CurrentAcceleration[2] );
                }
                else if (itNode -> IsFixed(VELOCITY_Z))
                {
                    CurrentDisplacement[2] = PreviousDisplacement[2] + mDeltaTime*(mBeta/mGamma*(CurrentVelocity[2]-PreviousVelocity[2])+PreviousVelocity[2]);
                }
                else if (itNode -> IsFixed(DISPLACEMENT_Z) == false)
                {
                    CurrentDisplacement[2] = PreviousDisplacement[2] + mDeltaTime * PreviousVelocity[2] + 0.5 * std::pow(mDeltaTime, 2) * PreviousAcceleration[2];
                }
            }

            noalias(DeltaDisplacement) = CurrentDisplacement - PreviousDisplacement;

            noalias(CurrentAcceleration) = 1.0/(mBeta*mDeltaTime*mDeltaTime)*(DeltaDisplacement - mDeltaTime*PreviousVelocity - (0.5-mBeta)*mDeltaTime*mDeltaTime*PreviousAcceleration);
            noalias(CurrentVelocity) = PreviousVelocity + (1.0-mGamma)*mDeltaTime*PreviousAcceleration + mGamma*mDeltaTime*CurrentAcceleration;

            // Terms related to Pressure field

            double& CurrentDt2Pressure = itNode->FastGetSolutionStepValue(Dt2_PRESSURE);
            double& CurrentDtPressure = itNode->FastGetSolutionStepValue(Dt_PRESSURE);
            DeltaPressure = itNode->FastGetSolutionStepValue(PRESSURE) - itNode->FastGetSolutionStepValue(PRESSURE, 1);
            const double& PreviousDt2Pressure = itNode->FastGetSolutionStepValue(Dt2_PRESSURE, 1);
            const double& PreviousDtPressure = itNode->FastGetSolutionStepValue(Dt_PRESSURE, 1);

            CurrentDt2Pressure = 1.0/(mBeta*mDeltaTime*mDeltaTime)*(DeltaPressure - mDeltaTime*PreviousDtPressure - (0.5-mBeta)*mDeltaTime*mDeltaTime*PreviousDt2Pressure);
            CurrentDtPressure = PreviousDtPressure + (1.0-mGamma)*mDeltaTime*PreviousDt2Pressure + mGamma*mDeltaTime*CurrentDt2Pressure;

        }

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeNonLinIteration(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        int NElems = static_cast<int>(r_model_part.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator itElem = el_begin + i;
            itElem -> InitializeNonLinearIteration(CurrentProcessInfo);
        }

        int NCons = static_cast<int>(r_model_part.Conditions().size());
        ModelPart::ConditionsContainerType::iterator con_begin = r_model_part.ConditionsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NCons; i++)
        {
            ModelPart::ConditionsContainerType::iterator itCond = con_begin + i;
            itCond -> InitializeNonLinearIteration(CurrentProcessInfo);
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void FinalizeNonLinIteration(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        int NElems = static_cast<int>(r_model_part.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator itElem = el_begin + i;
            itElem -> FinalizeNonLinearIteration(CurrentProcessInfo);
        }

        int NCons = static_cast<int>(r_model_part.Conditions().size());
        ModelPart::ConditionsContainerType::iterator con_begin = r_model_part.ConditionsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NCons; i++)
        {
            ModelPart::ConditionsContainerType::iterator itCond = con_begin + i;
            itCond -> FinalizeNonLinearIteration(CurrentProcessInfo);
        }

        KRATOS_CATCH("")
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

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        (rCurrentElement) -> CalculateMassMatrix(mMassMatrix[thread],CurrentProcessInfo);

        (rCurrentElement) -> CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);

        AddDynamicsToLHS (LHS_Contribution, mDampingMatrix[thread], mMassMatrix[thread], CurrentProcessInfo);

        AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mDampingMatrix[thread], mMassMatrix[thread], CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void Condition_CalculateSystemContributions(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentCondition) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        (rCurrentCondition) -> EquationIdVector(EquationId,CurrentProcessInfo);

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

        (rCurrentElement) -> CalculateMassMatrix(mMassMatrix[thread], CurrentProcessInfo);

        (rCurrentElement) -> CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mDampingMatrix[thread], mMassMatrix[thread], CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        (rCurrentCondition) -> EquationIdVector(EquationId, CurrentProcessInfo);


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

        (rCurrentElement) -> CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        (rCurrentElement) -> CalculateMassMatrix(mMassMatrix[thread],CurrentProcessInfo);

        (rCurrentElement) -> CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);

        AddDynamicsToLHS (LHS_Contribution, mDampingMatrix[thread], mMassMatrix[thread], CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void Condition_Calculate_LHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentCondition) -> CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        (rCurrentCondition) -> EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

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

        this->UpdateVariablesDerivatives(r_model_part);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    double mBeta;
    double mGamma;
    double mDeltaTime;
    double mrayleigh_m;
    double mrayleigh_k;

    std::vector< Matrix > mMassMatrix;
    std::vector< Vector > mAccelerationVector;
    std::vector< Matrix > mDampingMatrix;
    std::vector< Vector > mVelocityVector;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    inline void UpdateVariablesDerivatives(ModelPart& r_model_part)
    {
        KRATOS_TRY

        // Updating Velocities, Acclerations, DtPressure and Dt2Pressure
        array_1d<double,3> DeltaDisplacement;
        double DeltaPressure;

        const int NNodes = static_cast<int>(r_model_part.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = r_model_part.NodesBegin();

        #pragma omp parallel for private(DeltaDisplacement,DeltaPressure)
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNode = node_begin + i;

            // Termes related to displacement field

            array_1d<double,3>& CurrentDisplacement = itNode->FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double,3>& CurrentAcceleration = itNode->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double,3>& CurrentVelocity = itNode->FastGetSolutionStepValue(VELOCITY);

            const array_1d<double,3>& PreviousDisplacement = itNode->FastGetSolutionStepValue(DISPLACEMENT, 1);
            const array_1d<double,3>& PreviousAcceleration = itNode->FastGetSolutionStepValue(ACCELERATION, 1);
            const array_1d<double,3>& PreviousVelocity = itNode->FastGetSolutionStepValue(VELOCITY, 1);

            noalias(DeltaDisplacement) = CurrentDisplacement - PreviousDisplacement;

            noalias(CurrentAcceleration) = 1.0/(mBeta*mDeltaTime*mDeltaTime)*(DeltaDisplacement - mDeltaTime*PreviousVelocity - (0.5-mBeta)*mDeltaTime*mDeltaTime*PreviousAcceleration);
            noalias(CurrentVelocity) = PreviousVelocity + (1.0-mGamma)*mDeltaTime*PreviousAcceleration + mGamma*mDeltaTime*CurrentAcceleration;

            // Terms related to Pressure field

            double& CurrentDt2Pressure = itNode->FastGetSolutionStepValue(Dt2_PRESSURE);
            double& CurrentDtPressure = itNode->FastGetSolutionStepValue(Dt_PRESSURE);
            DeltaPressure = itNode->FastGetSolutionStepValue(PRESSURE) - itNode->FastGetSolutionStepValue(PRESSURE, 1);
            const double& PreviousDt2Pressure = itNode->FastGetSolutionStepValue(Dt2_PRESSURE, 1);
            const double& PreviousDtPressure = itNode->FastGetSolutionStepValue(Dt_PRESSURE, 1);

            CurrentDt2Pressure = 1.0/(mBeta*mDeltaTime*mDeltaTime)*(DeltaPressure - mDeltaTime*PreviousDtPressure - (0.5-mBeta)*mDeltaTime*mDeltaTime*PreviousDt2Pressure);
            CurrentDtPressure = PreviousDtPressure + (1.0-mGamma)*mDeltaTime*PreviousDt2Pressure + mGamma*mDeltaTime*CurrentDt2Pressure;

        }

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddDynamicsToLHS(LocalSystemMatrixType& LHS_Contribution,LocalSystemMatrixType& M,LocalSystemMatrixType& C,ProcessInfo& CurrentProcessInfo)
    {
        // adding mass contribution
        if (M.size1() != 0)
            noalias(LHS_Contribution) += 1.0/(mBeta*mDeltaTime*mDeltaTime)*M;

        // adding damping contribution
        if (C.size1() != 0)
            noalias(LHS_Contribution) += mGamma/(mBeta*mDeltaTime)*C;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddDynamicsToRHS(Element::Pointer rCurrentElement,LocalSystemVectorType& RHS_Contribution,
                            LocalSystemMatrixType& M,LocalSystemMatrixType& C,ProcessInfo& CurrentProcessInfo)
    {
        int thread = OpenMPUtils::ThisThread();

        //adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentElement->GetSecondDerivativesVector(mAccelerationVector[thread], 0);

            noalias(RHS_Contribution) -= prod(M, mAccelerationVector[thread]);
        }

        //adding damping contribution
        if (C.size1() != 0)
        {
            rCurrentElement->GetFirstDerivativesVector(mVelocityVector[thread], 0);

            noalias(RHS_Contribution) -= prod(C, mVelocityVector[thread]);
        }
    }


}; // Class DamUPScheme
}  // namespace Kratos

#endif // KRATOS_DAM_UP_SCHEME defined

