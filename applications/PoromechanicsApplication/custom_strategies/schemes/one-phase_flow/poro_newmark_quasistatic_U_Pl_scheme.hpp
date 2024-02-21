
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


#if !defined(KRATOS_PORO_NEWMARK_QUASISTATIC_U_PL_SCHEME )
#define  KRATOS_PORO_NEWMARK_QUASISTATIC_U_PL_SCHEME

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class PoroNewmarkQuasistaticUPlScheme : public Scheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( PoroNewmarkQuasistaticUPlScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;
    typedef typename BaseType::DofsArrayType                 DofsArrayType;
    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType         TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    PoroNewmarkQuasistaticUPlScheme(double theta_u,
                                    double theta_p,
                                    double beta,
                                    double gamma) :
                                    Scheme<TSparseSpace,TDenseSpace>()
    {
        mTheta_u = theta_u;
        mTheta_p = theta_p;
        // Here mBeta and mGamma are not used for the GN11 scheme
        mBeta = 1.0;
        mGamma = 1.0;
    }

    //------------------------------------------------------------------------------------

    ///Destructor
    ~PoroNewmarkQuasistaticUPlScheme() override {}

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
        if(LIQUID_PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument, "LIQUID_PRESSURE Key is 0. Check if all applications were correctly registered.", "" )
        if(DT_LIQUID_PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument, "DT_LIQUID_PRESSURE Key is 0. Check if all applications were correctly registered.", "" )
        if ( VELOCITY_COEFFICIENT.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "VELOCITY_COEFFICIENT Key is 0. Check if all applications were correctly registered.", "" )
        if ( DT_LIQUID_PRESSURE_COEFFICIENT.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "DT_LIQUID_PRESSURE_COEFFICIENT Key is 0. Check if all applications were correctly registered.", "" )

        //check that variables are correctly allocated
        for(ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); it++)
        {
            if(it->SolutionStepsDataHas(DISPLACEMENT) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_THROW_ERROR( std::logic_error, "VELOCITY variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(LIQUID_PRESSURE) == false)
                KRATOS_THROW_ERROR( std::logic_error, "LIQUID_PRESSURE variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(DT_LIQUID_PRESSURE) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DT_LIQUID_PRESSURE variable is not allocated for node ", it->Id() )

            if(it->HasDofFor(DISPLACEMENT_X) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_X dof on node ",it->Id() )
            if(it->HasDofFor(DISPLACEMENT_Y) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_Y dof on node ",it->Id() )
            if(it->HasDofFor(DISPLACEMENT_Z) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_Z dof on node ",it->Id() )
            if(it->HasDofFor(LIQUID_PRESSURE) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing LIQUID_PRESSURE dof on node ",it->Id() )
        }

        //check for minimum value of the buffer index.
        if (r_model_part.GetBufferSize() < 2)
            KRATOS_THROW_ERROR( std::logic_error, "insufficient buffer size. Buffer size should be greater than 2. Current size is", r_model_part.GetBufferSize() )

        // Check theta
        if(mTheta_u <= 0.0 || mTheta_p <= 0.0)
            KRATOS_THROW_ERROR( std::invalid_argument,"Some of the scheme variables: theta_u or theta_p has an invalid value ", "" )

        return 0;

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize(ModelPart& r_model_part) override
    {
        KRATOS_TRY

        mDeltaTime = r_model_part.GetProcessInfo()[DELTA_TIME];
        // Note: VELOCITY_COEFFICIENT is updated according to the GN11 scheme used here
        r_model_part.GetProcessInfo().SetValue(VELOCITY_COEFFICIENT,1.0/(mTheta_u*mDeltaTime));
        r_model_part.GetProcessInfo().SetValue(DT_LIQUID_PRESSURE_COEFFICIENT,1.0/(mTheta_p*mDeltaTime));

        const int NNodes = static_cast<int>(r_model_part.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = r_model_part.NodesBegin();

        // Initialize INITIAL_STRESS_TENSOR
        #pragma omp parallel for
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNode = node_begin + i;

            Matrix& rInitialStress = itNode->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR);
            if(rInitialStress.size1() != 3)
                rInitialStress.resize(3,3,false);
            noalias(rInitialStress) = ZeroMatrix(3,3);
        }

        BaseType::mSchemeIsInitialized = true;

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeElements( ModelPart& rModelPart) override
    {
        KRATOS_TRY

        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        int NElems = static_cast<int>(rModelPart.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();

        // #pragma omp parallel for
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator itElem = el_begin + i;
            itElem -> Initialize(rCurrentProcessInfo);
        }

        this->SetElementsAreInitialized();

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

        const ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();

        int NElems = static_cast<int>(r_model_part.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator itElem = el_begin + i;
            itElem -> InitializeSolutionStep(rCurrentProcessInfo);
        }

        int NCons = static_cast<int>(r_model_part.Conditions().size());
        ModelPart::ConditionsContainerType::iterator con_begin = r_model_part.ConditionsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NCons; i++)
        {
            ModelPart::ConditionsContainerType::iterator itCond = con_begin + i;
            itCond -> InitializeSolutionStep(rCurrentProcessInfo);
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

        const ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();

        int NElems = static_cast<int>(r_model_part.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator itElem = el_begin + i;
            itElem -> InitializeNonLinearIteration(rCurrentProcessInfo);
        }

        int NCons = static_cast<int>(r_model_part.Conditions().size());
        ModelPart::ConditionsContainerType::iterator con_begin = r_model_part.ConditionsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NCons; i++)
        {
            ModelPart::ConditionsContainerType::iterator itCond = con_begin + i;
            itCond -> InitializeNonLinearIteration(rCurrentProcessInfo);
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

        const ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();

        int NElems = static_cast<int>(r_model_part.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator itElem = el_begin + i;
            itElem -> FinalizeNonLinearIteration(rCurrentProcessInfo);
        }

        int NCons = static_cast<int>(r_model_part.Conditions().size());
        ModelPart::ConditionsContainerType::iterator con_begin = r_model_part.ConditionsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NCons; i++)
        {
            ModelPart::ConditionsContainerType::iterator itCond = con_begin + i;
            itCond -> FinalizeNonLinearIteration(rCurrentProcessInfo);
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
            const int NNodes = static_cast<int>(rModelPart.Nodes().size());
            ModelPart::NodesContainerType::iterator node_begin = rModelPart.NodesBegin();

            // Clear nodal variables
            #pragma omp parallel for
            for(int i = 0; i < NNodes; i++)
            {
                ModelPart::NodesContainerType::iterator itNode = node_begin + i;

                itNode->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
                Matrix& rNodalStress = itNode->FastGetSolutionStepValue(NODAL_EFFECTIVE_STRESS_TENSOR);
                if(rNodalStress.size1() != 3)
                    rNodalStress.resize(3,3,false);
                noalias(rNodalStress) = ZeroMatrix(3,3);
                itNode->FastGetSolutionStepValue(NODAL_JOINT_AREA) = 0.0;
                itNode->FastGetSolutionStepValue(NODAL_JOINT_WIDTH) = 0.0;
                itNode->FastGetSolutionStepValue(NODAL_JOINT_DAMAGE) = 0.0;
            }

            BaseType::FinalizeSolutionStep(rModelPart,A,Dx,b);

            // Compute smoothed nodal variables
            #pragma omp parallel for
            for(int n = 0; n < NNodes; n++)
            {
                ModelPart::NodesContainerType::iterator itNode = node_begin + n;

                const double& NodalArea = itNode->FastGetSolutionStepValue(NODAL_AREA);
                if (NodalArea>1.0e-20)
                {
                    const double InvNodalArea = 1.0/NodalArea;
                    Matrix& rNodalStress = itNode->FastGetSolutionStepValue(NODAL_EFFECTIVE_STRESS_TENSOR);
                    for(unsigned int i = 0; i<3; i++)
                    {
                        for(unsigned int j = 0; j<3; j++)
                        {
                            rNodalStress(i,j) *= InvNodalArea;
                        }
                    }
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
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        rCurrentElement.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,rCurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        rCurrentCondition.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);

        rCurrentCondition.EquationIdVector(EquationId,rCurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        rCurrentElement.CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,rCurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        rCurrentCondition.CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);

        rCurrentCondition.EquationIdVector(EquationId, rCurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateLHSContribution(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        rCurrentElement.CalculateLeftHandSide(LHS_Contribution,rCurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,rCurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateLHSContribution(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        rCurrentCondition.CalculateLeftHandSide(LHS_Contribution, rCurrentProcessInfo);

        rCurrentCondition.EquationIdVector(EquationId, rCurrentProcessInfo);

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

        int NumThreads = ParallelUtilities::GetNumThreads();
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

    double mTheta_u;
    double mTheta_p;
    double mBeta;
    double mGamma;
    double mDeltaTime;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    virtual inline void UpdateVariablesDerivatives(ModelPart& r_model_part)
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

            array_1d<double,3>& CurrentVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
            noalias(DeltaDisplacement) = itNode->FastGetSolutionStepValue(DISPLACEMENT) - itNode->FastGetSolutionStepValue(DISPLACEMENT, 1);
            const array_1d<double,3>& PreviousVelocity = itNode->FastGetSolutionStepValue(VELOCITY, 1);

            // Note: Since acceleration is not used in the quasi-static scheme, here we use a GN11 scheme to update velocity. 
            // A GN22 is used in the dynamic scheme
            noalias(CurrentVelocity) = 1.0/(mTheta_u*mDeltaTime)*(DeltaDisplacement - (1.0-mTheta_u)*mDeltaTime*PreviousVelocity);

            double& CurrentDtPressure = itNode->FastGetSolutionStepValue(DT_LIQUID_PRESSURE);
            DeltaPressure = itNode->FastGetSolutionStepValue(LIQUID_PRESSURE) - itNode->FastGetSolutionStepValue(LIQUID_PRESSURE, 1);
            const double& PreviousDtPressure = itNode->FastGetSolutionStepValue(DT_LIQUID_PRESSURE, 1);

            CurrentDtPressure = 1.0/(mTheta_p*mDeltaTime)*(DeltaPressure - (1.0-mTheta_p)*mDeltaTime*PreviousDtPressure);
        }

        KRATOS_CATCH( "" )
    }

}; // Class PoroNewmarkQuasistaticUPlScheme
}  // namespace Kratos

#endif // KRATOS_PORO_NEWMARK_QUASISTATIC_U_PL_SCHEME defined
