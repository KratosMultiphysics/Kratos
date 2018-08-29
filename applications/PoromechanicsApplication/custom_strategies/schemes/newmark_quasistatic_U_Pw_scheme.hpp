
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana


#if !defined(KRATOS_NEWMARK_QUASISTATIC_U_PW_SCHEME )
#define  KRATOS_NEWMARK_QUASISTATIC_U_PW_SCHEME

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class NewmarkQuasistaticUPwScheme : public Scheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( NewmarkQuasistaticUPwScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;
    typedef typename BaseType::DofsArrayType                 DofsArrayType;
    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType         TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    NewmarkQuasistaticUPwScheme(double beta, double gamma, double theta) : Scheme<TSparseSpace,TDenseSpace>()
    {
        mBeta = beta;
        mGamma = gamma;
        mTheta = theta;
    }

    //------------------------------------------------------------------------------------

    ///Destructor
    ~NewmarkQuasistaticUPwScheme() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(ModelPart& r_model_part) override
    {
        KRATOS_TRY

        //check for variables keys (verify that the variables are correctly initialized)
        if(DELTA_TIME.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"DELTA_TIME Key is 0. Check if all applications were correctly registered.", "")
        if(DISPLACEMENT.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"DISPLACEMENT Key is 0. Check if all applications were correctly registered.", "" )
        if(VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"VELOCITY Key is 0. Check if all applications were correctly registered.", "" )
        if(ACCELERATION.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"ACCELERATION Key is 0. Check if all applications were correctly registered.", "" )
        if(WATER_PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument, "WATER_PRESSURE Key is 0. Check if all applications were correctly registered.", "" )
        if(DT_WATER_PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument, "DT_WATER_PRESSURE Key is 0. Check if all applications were correctly registered.", "" )
        if ( VELOCITY_COEFFICIENT.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "VELOCITY_COEFFICIENT Key is 0. Check if all applications were correctly registered.", "" )
        if ( DT_PRESSURE_COEFFICIENT.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "DT_PRESSURE_COEFFICIENT Key is 0. Check if all applications were correctly registered.", "" )

        //check that variables are correctly allocated
        for(ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); it++)
        {
            if(it->SolutionStepsDataHas(DISPLACEMENT) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_THROW_ERROR( std::logic_error, "VELOCITY variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(ACCELERATION) == false)
                KRATOS_THROW_ERROR( std::logic_error, "ACCELERATION variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(WATER_PRESSURE) == false)
                KRATOS_THROW_ERROR( std::logic_error, "WATER_PRESSURE variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(DT_WATER_PRESSURE) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DT_WATER_PRESSURE variable is not allocated for node ", it->Id() )

            if(it->HasDofFor(DISPLACEMENT_X) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_X dof on node ",it->Id() )
            if(it->HasDofFor(DISPLACEMENT_Y) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_Y dof on node ",it->Id() )
            if(it->HasDofFor(DISPLACEMENT_Z) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_Z dof on node ",it->Id() )
            if(it->HasDofFor(WATER_PRESSURE) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing WATER_PRESSURE dof on node ",it->Id() )
        }

        //check for minimum value of the buffer index.
        if (r_model_part.GetBufferSize() < 2)
            KRATOS_THROW_ERROR( std::logic_error, "insufficient buffer size. Buffer size should be greater than 2. Current size is", r_model_part.GetBufferSize() )

        // Check beta, gamma and theta
        if(mBeta <= 0.0 || mGamma<= 0.0 || mTheta <= 0.0)
            KRATOS_THROW_ERROR( std::invalid_argument,"Some of the scheme variables: beta, gamma or theta has an invalid value ", "" )

        return 0;

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize(ModelPart& r_model_part) override
    {
        KRATOS_TRY

        mDeltaTime = r_model_part.GetProcessInfo()[DELTA_TIME];
        r_model_part.GetProcessInfo()[VELOCITY_COEFFICIENT] = mGamma/(mBeta*mDeltaTime);
        r_model_part.GetProcessInfo()[DT_PRESSURE_COEFFICIENT] = 1.0/(mTheta*mDeltaTime);

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
        r_model_part.GetProcessInfo()[VELOCITY_COEFFICIENT] = mGamma/(mBeta*mDeltaTime);
        r_model_part.GetProcessInfo()[DT_PRESSURE_COEFFICIENT] = 1.0/(mTheta*mDeltaTime);

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
        this->UpdateVariablesDerivatives(r_model_part);
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

    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        if(rModelPart.GetProcessInfo()[NODAL_SMOOTHING] == true)
        {
            unsigned int Dim = rModelPart.GetProcessInfo()[DOMAIN_SIZE];

            const int NNodes = static_cast<int>(rModelPart.Nodes().size());
            ModelPart::NodesContainerType::iterator node_begin = rModelPart.NodesBegin();

            // Clear nodal variables
            #pragma omp parallel for
            for(int i = 0; i < NNodes; i++)
            {
                ModelPart::NodesContainerType::iterator itNode = node_begin + i;

                itNode->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
                Matrix& rNodalStress = itNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                if(rNodalStress.size1() != Dim)
                    rNodalStress.resize(Dim,Dim,false);
                noalias(rNodalStress) = ZeroMatrix(Dim,Dim);
                itNode->FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) = 0.0;
                itNode->FastGetSolutionStepValue(NODAL_JOINT_AREA) = 0.0;
                itNode->FastGetSolutionStepValue(NODAL_JOINT_WIDTH) = 0.0;
                itNode->FastGetSolutionStepValue(NODAL_JOINT_DAMAGE) = 0.0;
            }

            BaseType::FinalizeSolutionStep(rModelPart,A,Dx,b);

            // Compute smoothed nodal variables
            #pragma omp parallel for
            for(int i = 0; i < NNodes; i++)
            {
                ModelPart::NodesContainerType::iterator itNode = node_begin + i;

                const double& NodalArea = itNode->FastGetSolutionStepValue(NODAL_AREA);
                if (NodalArea>1.0e-20)
                {
                    const double InvNodalArea = 1.0/NodalArea;
                    Matrix& rNodalStress = itNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                    for(unsigned int i = 0; i<Dim; i++)
                    {
                        for(unsigned int j = 0; j<Dim; j++)
                        {
                            rNodalStress(i,j) *= InvNodalArea;
                        }
                    }
                    double& NodalDamage = itNode->FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE);
                    NodalDamage *= InvNodalArea;
                }

                const double& NodalJointArea = itNode->FastGetSolutionStepValue(NODAL_JOINT_AREA);
                if (NodalJointArea>1.0e-20)
                {
                    const double InvNodalJointArea = 1.0/NodalJointArea;
                    double& NodalJointWidth = itNode->FastGetSolutionStepValue(NODAL_JOINT_WIDTH);
                    NodalJointWidth *= InvNodalJointArea;
                    double& NodalJointDamage = itNode->FastGetSolutionStepValue(NODAL_JOINT_DAMAGE);
                    NodalJointDamage *= InvNodalJointArea;
                }
            }
        }
        else
        {
            BaseType::FinalizeSolutionStep(rModelPart,A,Dx,b);
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

        (rCurrentElement) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

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

        (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

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

        (rCurrentElement) -> CalculateLeftHandSide(LHS_Contribution,CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

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
    double mTheta;
    double mDeltaTime;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    inline void UpdateVariablesDerivatives(ModelPart& r_model_part)
    {
        KRATOS_TRY

        //Update Acceleration, Velocity and DtPressure

        array_1d<double,3> DeltaDisplacement;
        double DeltaPressure;

        const int NNodes = static_cast<int>(r_model_part.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = r_model_part.NodesBegin();

        #pragma omp parallel for private(DeltaDisplacement,DeltaPressure)
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNode = node_begin + i;

            array_1d<double,3>& CurrentAcceleration = itNode->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double,3>& CurrentVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
            noalias(DeltaDisplacement) = itNode->FastGetSolutionStepValue(DISPLACEMENT) - itNode->FastGetSolutionStepValue(DISPLACEMENT, 1);
            const array_1d<double,3>& PreviousAcceleration = itNode->FastGetSolutionStepValue(ACCELERATION, 1);
            const array_1d<double,3>& PreviousVelocity = itNode->FastGetSolutionStepValue(VELOCITY, 1);

            noalias(CurrentAcceleration) = 1.0/(mBeta*mDeltaTime*mDeltaTime)*(DeltaDisplacement - mDeltaTime*PreviousVelocity - (0.5-mBeta)*mDeltaTime*mDeltaTime*PreviousAcceleration);
            noalias(CurrentVelocity) = PreviousVelocity + (1.0-mGamma)*mDeltaTime*PreviousAcceleration + mGamma*mDeltaTime*CurrentAcceleration;

            double& CurrentDtPressure = itNode->FastGetSolutionStepValue(DT_WATER_PRESSURE);
            DeltaPressure = itNode->FastGetSolutionStepValue(WATER_PRESSURE) - itNode->FastGetSolutionStepValue(WATER_PRESSURE, 1);
            const double& PreviousDtPressure = itNode->FastGetSolutionStepValue(DT_WATER_PRESSURE, 1);

            CurrentDtPressure = 1.0/(mTheta*mDeltaTime)*(DeltaPressure - (1.0-mTheta)*mDeltaTime*PreviousDtPressure);
        }

        KRATOS_CATCH( "" )
    }

}; // Class NewmarkQuasistaticUPwScheme
}  // namespace Kratos

#endif // KRATOS_NEWMARK_QUASISTATIC_U_PW_SCHEME defined
