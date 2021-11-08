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
#include "utilities/parallel_utilities.h"

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

    int Check(const ModelPart& rModelPart) const override
    {
        KRATOS_TRY

        BaseType::Check(rModelPart);

        //check that variables are correctly allocated
        for (const auto& rNode : rModelPart.Nodes())
        {
            if(rNode.SolutionStepsDataHas(DISPLACEMENT) == false)
                KRATOS_ERROR << "DISPLACEMENT variable is not allocated for node "
                             << rNode.Id()
                             << std::endl;

            if(rNode.SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_ERROR << "VELOCITY variable is not allocated for node "
                             << rNode.Id()
                             << std::endl;

            if(rNode.SolutionStepsDataHas(ACCELERATION) == false)
                KRATOS_ERROR << "ACCELERATION variable is not allocated for node "
                             << rNode.Id()
                             << std::endl;

            if(rNode.SolutionStepsDataHas(WATER_PRESSURE) == false)
                KRATOS_ERROR << "WATER_PRESSURE variable is not allocated for node "
                             << rNode.Id()
                             << std::endl;

            if(rNode.SolutionStepsDataHas(DT_WATER_PRESSURE) == false)
                KRATOS_ERROR << "DT_WATER_PRESSURE variable is not allocated for node "
                             << rNode.Id()
                             << std::endl;

            if(rNode.HasDofFor(DISPLACEMENT_X) == false)
                KRATOS_ERROR << "missing DISPLACEMENT_X dof on node "
                             << rNode.Id()
                             << std::endl;

            if(rNode.HasDofFor(DISPLACEMENT_Y) == false)
                KRATOS_ERROR << "missing DISPLACEMENT_Y dof on node "
                             << rNode.Id()
                             << std::endl;

            if(rNode.HasDofFor(DISPLACEMENT_Z) == false)
                KRATOS_ERROR << "missing DISPLACEMENT_Z dof on node "
                             << rNode.Id()
                             << std::endl;

            if(rNode.HasDofFor(WATER_PRESSURE) == false)
                KRATOS_ERROR << "missing WATER_PRESSURE dof on node "
                             << rNode.Id()
                             << std::endl;
        }

        //check for minimum value of the buffer index.
        if (rModelPart.GetBufferSize() < 2)
            KRATOS_ERROR << "insufficient buffer size. Buffer size should be greater than 2. Current size is"
                         << rModelPart.GetBufferSize()
                         << std::endl;

        // Check beta, gamma and theta
        if (mBeta <= 0.0 || mGamma<= 0.0 || mTheta <= 0.0)
            KRATOS_ERROR << "Some of the scheme variables: beta, gamma or theta has an invalid value"
                         << std::endl;

        return 0;

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        mDeltaTime = rModelPart.GetProcessInfo()[DELTA_TIME];
        rModelPart.GetProcessInfo()[VELOCITY_COEFFICIENT] = mGamma/(mBeta*mDeltaTime);
        rModelPart.GetProcessInfo()[DT_PRESSURE_COEFFICIENT] = 1.0/(mTheta*mDeltaTime);

        BaseType::mSchemeIsInitialized = true;

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        mDeltaTime = rModelPart.GetProcessInfo()[DELTA_TIME];
        rModelPart.GetProcessInfo()[VELOCITY_COEFFICIENT] = mGamma/(mBeta*mDeltaTime);
        rModelPart.GetProcessInfo()[DT_PRESSURE_COEFFICIENT] = 1.0/(mTheta*mDeltaTime);

        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        block_for_each(rModelPart.Elements(), [&rCurrentProcessInfo](Element& rElement) {
            const bool isActive = (rElement.IsDefined(ACTIVE)) ? rElement.Is(ACTIVE) : true;
            if (isActive) rElement.InitializeSolutionStep(rCurrentProcessInfo);
        });

        block_for_each(rModelPart.Conditions(), [&rCurrentProcessInfo](Condition& rCondition) {
            const bool isActive = (rCondition.IsDefined(ACTIVE)) ? rCondition.Is(ACTIVE) : true;
            if (isActive) rCondition.InitializeSolutionStep(rCurrentProcessInfo);
        });

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        this->UpdateVariablesDerivatives(rModelPart);
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        block_for_each(rModelPart.Elements(), [&rCurrentProcessInfo](Element& rElement) {
            const bool isActive = (rElement.IsDefined(ACTIVE)) ? rElement.Is(ACTIVE) : true;
            if (isActive) rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
        });

        block_for_each(rModelPart.Conditions(), [&rCurrentProcessInfo](Condition& rCondition) {
            const bool isActive = (rCondition.IsDefined(ACTIVE)) ? rCondition.Is(ACTIVE) : true;
            if (isActive) rCondition.InitializeNonLinearIteration(rCurrentProcessInfo);
        });

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

        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        block_for_each(rModelPart.Elements(), [&rCurrentProcessInfo](Element& rElement) {
            const bool isActive = (rElement.IsDefined(ACTIVE)) ? rElement.Is(ACTIVE) : true;
            if (isActive) rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
        });

        block_for_each(rModelPart.Conditions(), [&rCurrentProcessInfo](Condition& rCondition) {
            const bool isActive = (rCondition.IsDefined(ACTIVE)) ? rCondition.Is(ACTIVE) : true;
            if (isActive) rCondition.FinalizeNonLinearIteration(rCurrentProcessInfo);
        });

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

            // Clear nodal variables
            block_for_each(rModelPart.Nodes(), [&StressTensorSize](Node<3>& rNode) {
                rNode.FastGetSolutionStepValue(NODAL_AREA) = 0.0;
                Matrix& rNodalStress = rNode.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                if (rNodalStress.size1() != StressTensorSize)
                    rNodalStress.resize(StressTensorSize,StressTensorSize,false);
                noalias(rNodalStress) = ZeroMatrix(StressTensorSize, StressTensorSize);
                rNode.FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) = 0.0;
                rNode.FastGetSolutionStepValue(NODAL_JOINT_AREA) = 0.0;
                rNode.FastGetSolutionStepValue(NODAL_JOINT_WIDTH) = 0.0;
                rNode.FastGetSolutionStepValue(NODAL_JOINT_DAMAGE) = 0.0;
            });


            FinalizeSolutionStepActiveEntities(rModelPart,A,Dx,b);


            // Compute smoothed nodal variables
            block_for_each(rModelPart.Nodes(), [&](Node<3>& rNode) {
                const double& NodalArea = rNode.FastGetSolutionStepValue(NODAL_AREA);
                if (NodalArea > 1.0e-20)
                {
                    const double InvNodalArea = 1.0/NodalArea;
                    Matrix& rNodalStress = rNode.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                    for(unsigned int i = 0; i < rNodalStress.size1(); i++)
                    {
                        for(unsigned int j = 0; j < rNodalStress.size2(); j++)
                        {
                            rNodalStress(i,j) *= InvNodalArea;
                        }
                    }
                    rNode.FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) *= InvNodalArea;
                }

                const double& NodalJointArea = rNode.FastGetSolutionStepValue(NODAL_JOINT_AREA);
                if (NodalJointArea > 1.0e-20)
                {
                    const double InvNodalJointArea = 1.0/NodalJointArea;
                    rNode.FastGetSolutionStepValue(NODAL_JOINT_WIDTH)  *= InvNodalJointArea;
                    rNode.FastGetSolutionStepValue(NODAL_JOINT_DAMAGE) *= InvNodalJointArea;
                }
            });

        }
        else
        {
            FinalizeSolutionStepActiveEntities(rModelPart,A,Dx,b);
        }

        KRATOS_CATCH("")
    }

    //------------------------------------------------------------------------------------------
    void FinalizeSolutionStepActiveEntities(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        block_for_each(rModelPart.Elements(), [&rCurrentProcessInfo](Element& rElement) {
            const bool isActive = (rElement.IsDefined(ACTIVE)) ? rElement.Is(ACTIVE) : true;
            if (isActive) rElement.FinalizeSolutionStep(rCurrentProcessInfo);
        });

        block_for_each(rModelPart.Conditions(), [&rCurrentProcessInfo](Condition& rCondition) {
            const bool isActive = (rCondition.IsDefined(ACTIVE)) ? rCondition.Is(ACTIVE) : true;
            if (isActive) rCondition.FinalizeSolutionStep(rCurrentProcessInfo);
        });

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
        ModelPart& rModelPart,
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

        this->UpdateVariablesDerivatives(rModelPart);

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
        mDeltaTime = 0.0;
    }

    virtual inline void UpdateVariablesDerivatives(ModelPart& rModelPart)
    {
        KRATOS_TRY

        //Update Acceleration, Velocity and DtPressure
        block_for_each(rModelPart.Nodes(), [&](Node<3>& rNode) {
            noalias(rNode.FastGetSolutionStepValue(ACCELERATION)) =  ((  rNode.FastGetSolutionStepValue(DISPLACEMENT)
                                                                      - rNode.FastGetSolutionStepValue(DISPLACEMENT, 1))
                                                                   - mDeltaTime * rNode.FastGetSolutionStepValue(VELOCITY, 1) 
                                                                   - (0.5-mBeta) * mDeltaTime * mDeltaTime
                                                                   * rNode.FastGetSolutionStepValue(ACCELERATION, 1) )
                                                                   / (mBeta*mDeltaTime*mDeltaTime);

            noalias(rNode.FastGetSolutionStepValue(VELOCITY)) =  rNode.FastGetSolutionStepValue(VELOCITY, 1)
                                                                + (1.0-mGamma)*mDeltaTime
                                                                * rNode.FastGetSolutionStepValue(ACCELERATION, 1)
                                                                + mGamma * mDeltaTime
                                                                * rNode.FastGetSolutionStepValue(ACCELERATION);

            const double DeltaPressure =  rNode.FastGetSolutionStepValue(WATER_PRESSURE)
                                        - rNode.FastGetSolutionStepValue(WATER_PRESSURE, 1);

            rNode.FastGetSolutionStepValue(DT_WATER_PRESSURE) = ( DeltaPressure
                                                                 - (1.0-mTheta) * mDeltaTime
                                                                  * rNode.FastGetSolutionStepValue(DT_WATER_PRESSURE, 1))
                                                                / (mTheta*mDeltaTime);
        });

        KRATOS_CATCH( "" )
    }

}; // Class NewmarkQuasistaticUPwScheme
}  // namespace Kratos

#endif // KRATOS_NEWMARK_QUASISTATIC_U_PW_SCHEME defined
