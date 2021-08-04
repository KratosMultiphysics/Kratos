// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#if !defined(KRATOS_NEWMARK_QUASISTATIC_U_PW_SCHEME )
#define  KRATOS_NEWMARK_QUASISTATIC_U_PW_SCHEME

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"

// Application includes
#include "geo_mechanics_application_variables.h"

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

    int Check(const ModelPart& r_model_part) const override
    {
        KRATOS_TRY

        BaseType::Check(r_model_part);

        //check that variables are correctly allocated
        for (ModelPart::NodesContainerType::const_iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); it++)
        {
            if(it->SolutionStepsDataHas(DISPLACEMENT) == false)
                KRATOS_ERROR << "DISPLACEMENT variable is not allocated for node "
                             << it->Id()
                             << std::endl;

            if(it->SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_ERROR << "VELOCITY variable is not allocated for node "
                             << it->Id()
                             << std::endl;

            if(it->SolutionStepsDataHas(ACCELERATION) == false)
                KRATOS_ERROR << "ACCELERATION variable is not allocated for node "
                             << it->Id()
                             << std::endl;

            if(it->SolutionStepsDataHas(WATER_PRESSURE) == false)
                KRATOS_ERROR << "WATER_PRESSURE variable is not allocated for node "
                             << it->Id()
                             << std::endl;

            if(it->SolutionStepsDataHas(DT_WATER_PRESSURE) == false)
                KRATOS_ERROR << "DT_WATER_PRESSURE variable is not allocated for node "
                             << it->Id()
                             << std::endl;

            if(it->HasDofFor(DISPLACEMENT_X) == false)
                KRATOS_ERROR << "missing DISPLACEMENT_X dof on node "
                             << it->Id()
                             << std::endl;

            if(it->HasDofFor(DISPLACEMENT_Y) == false)
                KRATOS_ERROR << "missing DISPLACEMENT_Y dof on node "
                             << it->Id()
                             << std::endl;

            if(it->HasDofFor(DISPLACEMENT_Z) == false)
                KRATOS_ERROR << "missing DISPLACEMENT_Z dof on node "
                             << it->Id()
                             << std::endl;

            if(it->HasDofFor(WATER_PRESSURE) == false)
                KRATOS_ERROR << "missing WATER_PRESSURE dof on node "
                             << it->Id()
                             << std::endl;
        }

        //check for minimum value of the buffer index.
        if (r_model_part.GetBufferSize() < 2)
            KRATOS_ERROR << "insufficient buffer size. Buffer size should be greater than 2. Current size is"
                         << r_model_part.GetBufferSize()
                         << std::endl;

        // Check beta, gamma and theta
        if (mBeta <= 0.0 || mGamma<= 0.0 || mTheta <= 0.0)
            KRATOS_ERROR << "Some of the scheme variables: beta, gamma or theta has an invalid value"
                         << std::endl;

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

        const ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        int NElems = static_cast<int>(r_model_part.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator it = el_begin + i;
            const bool entity_is_active = (it->IsDefined(ACTIVE)) ? it->Is(ACTIVE) : true;
            if (entity_is_active) {
                it -> InitializeSolutionStep(CurrentProcessInfo);
            }
        }

        int NCons = static_cast<int>(r_model_part.Conditions().size());
        ModelPart::ConditionsContainerType::iterator con_begin = r_model_part.ConditionsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NCons; i++)
        {
            ModelPart::ConditionsContainerType::iterator it = con_begin + i;
            const bool entity_is_active = (it->IsDefined(ACTIVE)) ? it->Is(ACTIVE) : true;
            if (entity_is_active) {
                it -> InitializeSolutionStep(CurrentProcessInfo);
            }
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

        const ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        int NElems = static_cast<int>(r_model_part.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator it = el_begin + i;
            const bool entity_is_active = (it->IsDefined(ACTIVE)) ? it->Is(ACTIVE) : true;
            if (entity_is_active) {
                it -> InitializeNonLinearIteration(CurrentProcessInfo);
            }
        }

        int NCons = static_cast<int>(r_model_part.Conditions().size());
        ModelPart::ConditionsContainerType::iterator con_begin = r_model_part.ConditionsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NCons; i++)
        {
            ModelPart::ConditionsContainerType::iterator it = con_begin + i;
            const bool entity_is_active = (it->IsDefined(ACTIVE)) ? it->Is(ACTIVE) : true;
            if (entity_is_active) {
                it -> InitializeNonLinearIteration(CurrentProcessInfo);
            }
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

        const ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        int NElems = static_cast<int>(r_model_part.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator it = el_begin + i;
            const bool entity_is_active = (it->IsDefined(ACTIVE)) ? it->Is(ACTIVE) : true;
            if (entity_is_active) {
                it -> FinalizeNonLinearIteration(CurrentProcessInfo);
            }
        }

        int NCons = static_cast<int>(r_model_part.Conditions().size());
        ModelPart::ConditionsContainerType::iterator con_begin = r_model_part.ConditionsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NCons; i++)
        {
            ModelPart::ConditionsContainerType::iterator it = con_begin + i;
            const bool entity_is_active = (it->IsDefined(ACTIVE)) ? it->Is(ACTIVE) : true;
            if (entity_is_active) {
                it -> FinalizeNonLinearIteration(CurrentProcessInfo);
            }
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

        if (rModelPart.GetProcessInfo()[NODAL_SMOOTHING] == true)
        {
            unsigned int Dim = rModelPart.GetProcessInfo()[DOMAIN_SIZE];

            SizeType StressTensorSize = STRESS_TENSOR_SIZE_2D;
            if (Dim == N_DIM_3D) StressTensorSize = STRESS_TENSOR_SIZE_3D; 

            const int NNodes = static_cast<int>(rModelPart.Nodes().size());
            ModelPart::NodesContainerType::iterator node_begin = rModelPart.NodesBegin();

            // Clear nodal variables
            #pragma omp parallel for
            for(int n = 0; n < NNodes; n++)
            {
                ModelPart::NodesContainerType::iterator itNode = node_begin + n;

                itNode->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
                Matrix& rNodalStress = itNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                if (rNodalStress.size1() != StressTensorSize)
                    rNodalStress.resize(StressTensorSize,StressTensorSize,false);
                noalias(rNodalStress) = ZeroMatrix(StressTensorSize, StressTensorSize);
                itNode->FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) = 0.0;
                itNode->FastGetSolutionStepValue(NODAL_JOINT_AREA) = 0.0;
                itNode->FastGetSolutionStepValue(NODAL_JOINT_WIDTH) = 0.0;
                itNode->FastGetSolutionStepValue(NODAL_JOINT_DAMAGE) = 0.0;
            }

            FinalizeSolutionStepActiveEntities(rModelPart,A,Dx,b);

            // Compute smoothed nodal variables
            #pragma omp parallel for
            for(int i = 0; i < NNodes; i++)
            {
                ModelPart::NodesContainerType::iterator itNode = node_begin + i;

                const double& NodalArea = itNode->FastGetSolutionStepValue(NODAL_AREA);
                if (NodalArea > 1.0e-20)
                {
                    const double InvNodalArea = 1.0/NodalArea;
                    Matrix& rNodalStress = itNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                    for(unsigned int i = 0; i < rNodalStress.size1(); i++)
                    {
                        for(unsigned int j = 0; j < rNodalStress.size2(); j++)
                        {
                            rNodalStress(i,j) *= InvNodalArea;
                        }
                    }
                    double& NodalDamage = itNode->FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE);
                    NodalDamage *= InvNodalArea;
                }

                const double& NodalJointArea = itNode->FastGetSolutionStepValue(NODAL_JOINT_AREA);
                if (NodalJointArea > 1.0e-20)
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
            FinalizeSolutionStepActiveEntities(rModelPart,A,Dx,b);
        }

        KRATOS_CATCH("")
    }

    //------------------------------------------------------------------------------------------
    void FinalizeSolutionStepActiveEntities(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        const ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        int NElems = static_cast<int>(r_model_part.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator it = el_begin + i;
            const bool entity_is_active = (it->IsDefined(ACTIVE)) ? it->Is(ACTIVE) : true;
            if (entity_is_active) {
                it -> FinalizeSolutionStep(CurrentProcessInfo);
            }
        }

        int NCons = static_cast<int>(r_model_part.Conditions().size());
        ModelPart::ConditionsContainerType::iterator con_begin = r_model_part.ConditionsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NCons; i++)
        {
            ModelPart::ConditionsContainerType::iterator it = con_begin + i;
            const bool entity_is_active = (it->IsDefined(ACTIVE)) ? it->Is(ACTIVE) : true;
            if (entity_is_active) {
                it -> FinalizeSolutionStep(CurrentProcessInfo);
            }
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
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentElement).CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        (rCurrentElement).EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentCondition).CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        (rCurrentCondition).EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateRHSContribution(
        Element &rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentElement).CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        (rCurrentElement).EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateRHSContribution(
        Condition &rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentCondition).CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        (rCurrentCondition).EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateLHSContribution(
        Element &rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentElement).CalculateLeftHandSide(LHS_Contribution,CurrentProcessInfo);

        (rCurrentElement).EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateLHSContribution(
        Condition &rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentCondition).CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        (rCurrentCondition).EquationIdVector(EquationId, CurrentProcessInfo);

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

    double mBeta;
    double mGamma;
    double mTheta;
    double mDeltaTime;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    NewmarkQuasistaticUPwScheme() : Scheme<TSparseSpace,TDenseSpace>()
    {
        mBeta = 0.25;
        mGamma = 0.5;
        mTheta = 0.5;
    }

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

            array_1d<double,3>& CurrentAcceleration = itNode->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double,3>& CurrentVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
            noalias(DeltaDisplacement) = itNode->FastGetSolutionStepValue(DISPLACEMENT) - itNode->FastGetSolutionStepValue(DISPLACEMENT, 1);
            const array_1d<double,3>& PreviousAcceleration = itNode->FastGetSolutionStepValue(ACCELERATION, 1);
            const array_1d<double,3>& PreviousVelocity = itNode->FastGetSolutionStepValue(VELOCITY, 1);

            noalias(CurrentAcceleration) =  (  DeltaDisplacement
                                             - mDeltaTime*PreviousVelocity 
                                             - (0.5-mBeta)*mDeltaTime*mDeltaTime*PreviousAcceleration )
                                          / (mBeta*mDeltaTime*mDeltaTime);

            noalias(CurrentVelocity) =  PreviousVelocity
                                      + (1.0-mGamma)*mDeltaTime*PreviousAcceleration
                                      + mGamma*mDeltaTime*CurrentAcceleration;

            double& CurrentDtPressure = itNode->FastGetSolutionStepValue(DT_WATER_PRESSURE);
            DeltaPressure = itNode->FastGetSolutionStepValue(WATER_PRESSURE) - itNode->FastGetSolutionStepValue(WATER_PRESSURE, 1);
            const double& PreviousDtPressure = itNode->FastGetSolutionStepValue(DT_WATER_PRESSURE, 1);

            CurrentDtPressure = (DeltaPressure - (1.0-mTheta)*mDeltaTime*PreviousDtPressure) / (mTheta*mDeltaTime);
        }

        KRATOS_CATCH( "" )
    }

}; // Class NewmarkQuasistaticUPwScheme
}  // namespace Kratos

#endif // KRATOS_NEWMARK_QUASISTATIC_U_PW_SCHEME defined
