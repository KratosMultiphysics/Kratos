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

#if !defined(KRATOS_NEWMARK_DYNAMIC_U_PW_SCHEME )
#define  KRATOS_NEWMARK_DYNAMIC_U_PW_SCHEME

#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_strategies/schemes/newmark_quasistatic_U_Pw_scheme.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class NewmarkDynamicUPwScheme : public NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( NewmarkDynamicUPwScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;
    typedef typename BaseType::DofsArrayType                 DofsArrayType;
    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType         TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mDeltaTime;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mBeta;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mGamma;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mTheta;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    NewmarkDynamicUPwScheme(double beta, double gamma, double theta)
        : NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>(beta, gamma, theta)
    {
        //Allocate auxiliary memory
        int NumThreads = ParallelUtilities::GetNumThreads();
        mMassMatrix.resize(NumThreads);
        mAccelerationVector.resize(NumThreads);
        mDampingMatrix.resize(NumThreads);
        mVelocityVector.resize(NumThreads);
    }

    //------------------------------------------------------------------------------------

    ///Destructor
    ~NewmarkDynamicUPwScheme() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        // Predict Displacements on free nodes and update Acceleration, Velocity and DtPressure
        block_for_each(rModelPart.Nodes(), [&](Node& rNode){
            if (rNode.IsFixed(ACCELERATION_X))
            {
                const double &PreviousDisplacement = rNode.FastGetSolutionStepValue(DISPLACEMENT_X, 1);
                const double &PreviousVelocity     = rNode.FastGetSolutionStepValue(VELOCITY_X, 1);
                const double &PreviousAcceleration = rNode.FastGetSolutionStepValue(ACCELERATION_X, 1);
                const double &CurrentAcceleration  = rNode.FastGetSolutionStepValue(ACCELERATION_X);

                rNode.FastGetSolutionStepValue(DISPLACEMENT_X) =  PreviousDisplacement
                                                                + mDeltaTime * PreviousVelocity
                                                                + mDeltaTime * mDeltaTime
                                                                * ( ( 0.5 - mBeta) * PreviousAcceleration + mBeta * CurrentAcceleration );
            }
            else if (rNode.IsFixed(VELOCITY_X))
            {
                const double &PreviousDisplacement = rNode.FastGetSolutionStepValue(DISPLACEMENT_X, 1);
                const double &PreviousVelocity     = rNode.FastGetSolutionStepValue(VELOCITY_X, 1);
                const double &CurrentVelocity      = rNode.FastGetSolutionStepValue(VELOCITY_X);
                rNode.FastGetSolutionStepValue(DISPLACEMENT_X) =  PreviousDisplacement
                                                                 + mDeltaTime*(mBeta/mGamma*(CurrentVelocity - PreviousVelocity)
                                                                 + PreviousVelocity);
            }
            else if (!rNode.IsFixed(DISPLACEMENT_X))
            {
                const double &PreviousDisplacement = rNode.FastGetSolutionStepValue(DISPLACEMENT_X, 1);
                const double &PreviousVelocity     = rNode.FastGetSolutionStepValue(VELOCITY_X, 1);
                const double &PreviousAcceleration = rNode.FastGetSolutionStepValue(ACCELERATION_X, 1);

                rNode.FastGetSolutionStepValue(DISPLACEMENT_X) =   PreviousDisplacement
                                                                 + mDeltaTime * PreviousVelocity
                                                                 + 0.5 * mDeltaTime * mDeltaTime
                                                                 * PreviousAcceleration;
            }

            if (rNode.IsFixed(ACCELERATION_Y))
            {
                const double &PreviousDisplacement = rNode.FastGetSolutionStepValue(DISPLACEMENT_Y, 1);
                const double &PreviousVelocity     = rNode.FastGetSolutionStepValue(VELOCITY_Y, 1);
                const double &PreviousAcceleration = rNode.FastGetSolutionStepValue(ACCELERATION_Y, 1);
                const double &CurrentAcceleration  = rNode.FastGetSolutionStepValue(ACCELERATION_Y);

                rNode.FastGetSolutionStepValue(DISPLACEMENT_Y) =  PreviousDisplacement
                                                                + mDeltaTime * PreviousVelocity
                                                                + mDeltaTime * mDeltaTime
                                                                * ( ( 0.5 - mBeta) * PreviousAcceleration + mBeta * CurrentAcceleration );
            }
            else if (rNode.IsFixed(VELOCITY_Y))
            {
                const double &PreviousDisplacement = rNode.FastGetSolutionStepValue(DISPLACEMENT_Y, 1);
                const double &PreviousVelocity     = rNode.FastGetSolutionStepValue(VELOCITY_Y, 1);
                const double &CurrentVelocity      = rNode.FastGetSolutionStepValue(VELOCITY_Y);
                rNode.FastGetSolutionStepValue(DISPLACEMENT_Y) =  PreviousDisplacement
                                                                 + mDeltaTime*(mBeta/mGamma*(CurrentVelocity - PreviousVelocity)
                                                                 + PreviousVelocity);
            }
            else if (!rNode.IsFixed(DISPLACEMENT_Y))
            {
                const double &PreviousDisplacement = rNode.FastGetSolutionStepValue(DISPLACEMENT_Y, 1);
                const double &PreviousVelocity     = rNode.FastGetSolutionStepValue(VELOCITY_Y, 1);
                const double &PreviousAcceleration = rNode.FastGetSolutionStepValue(ACCELERATION_Y, 1);

                rNode.FastGetSolutionStepValue(DISPLACEMENT_Y) =   PreviousDisplacement
                                                                 + mDeltaTime * PreviousVelocity
                                                                 + 0.5 * mDeltaTime * mDeltaTime
                                                                 * PreviousAcceleration;
            }

            // For 3D cases
            if (rNode.HasDofFor(DISPLACEMENT_Z))
            {
                if (rNode.IsFixed(ACCELERATION_Z))
                {
                    const double &PreviousDisplacement = rNode.FastGetSolutionStepValue(DISPLACEMENT_Z, 1);
                    const double &PreviousVelocity     = rNode.FastGetSolutionStepValue(VELOCITY_Z, 1);
                    const double &PreviousAcceleration = rNode.FastGetSolutionStepValue(ACCELERATION_Z, 1);
                    const double &CurrentAcceleration  = rNode.FastGetSolutionStepValue(ACCELERATION_Z);

                    rNode.FastGetSolutionStepValue(DISPLACEMENT_Z) =  PreviousDisplacement
                                                                    + mDeltaTime * PreviousVelocity
                                                                    + mDeltaTime * mDeltaTime
                                                                    * ( ( 0.5 - mBeta) * PreviousAcceleration + mBeta * CurrentAcceleration );
                }
                else if (rNode.IsFixed(VELOCITY_Z))
                {
                    const double &PreviousDisplacement = rNode.FastGetSolutionStepValue(DISPLACEMENT_Z, 1);
                    const double &PreviousVelocity     = rNode.FastGetSolutionStepValue(VELOCITY_Z, 1);
                    const double &CurrentVelocity      = rNode.FastGetSolutionStepValue(VELOCITY_Z);
                    rNode.FastGetSolutionStepValue(DISPLACEMENT_Z) =  PreviousDisplacement
                                                                    + mDeltaTime*(mBeta/mGamma*(CurrentVelocity - PreviousVelocity)
                                                                    + PreviousVelocity);
                }
                else if (!rNode.IsFixed(DISPLACEMENT_Z))
                {
                    const double &PreviousDisplacement = rNode.FastGetSolutionStepValue(DISPLACEMENT_Z, 1);
                    const double &PreviousVelocity     = rNode.FastGetSolutionStepValue(VELOCITY_Z, 1);
                    const double &PreviousAcceleration = rNode.FastGetSolutionStepValue(ACCELERATION_Z, 1);

                    rNode.FastGetSolutionStepValue(DISPLACEMENT_Z) =   PreviousDisplacement
                                                                    + mDeltaTime * PreviousVelocity
                                                                    + 0.5 * mDeltaTime * mDeltaTime
                                                                    * PreviousAcceleration;
                }
            }

            noalias(rNode.FastGetSolutionStepValue(ACCELERATION)) =   1.0/(mBeta*mDeltaTime*mDeltaTime)
                                                                    * ( (rNode.FastGetSolutionStepValue(DISPLACEMENT) - rNode.FastGetSolutionStepValue(DISPLACEMENT,1))
                                                                     - mDeltaTime
                                                                     * rNode.FastGetSolutionStepValue(VELOCITY, 1)
                                                                     - (0.5-mBeta) * mDeltaTime * mDeltaTime
                                                                     * rNode.FastGetSolutionStepValue(ACCELERATION,1));

            noalias(rNode.FastGetSolutionStepValue(VELOCITY)) =   rNode.FastGetSolutionStepValue(VELOCITY, 1)
                                                                 + (1.0-mGamma) * mDeltaTime
                                                                 * rNode.FastGetSolutionStepValue(ACCELERATION,1)
                                                                 + mGamma * mDeltaTime
                                                                 * rNode.FastGetSolutionStepValue(ACCELERATION);

            const double DeltaPressure =  rNode.FastGetSolutionStepValue(WATER_PRESSURE)
                                        - rNode.FastGetSolutionStepValue(WATER_PRESSURE, 1);

            rNode.FastGetSolutionStepValue(DT_WATER_PRESSURE) =  1.0/(mTheta*mDeltaTime)
                                                                * ( DeltaPressure
                                                                   - (1.0-mTheta) * mDeltaTime 
                                                                   * rNode.FastGetSolutionStepValue(DT_WATER_PRESSURE,1));
        });


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

        int thread = OpenMPUtils::ThisThread();

        (rCurrentCondition).CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        (rCurrentCondition).CalculateMassMatrix(mMassMatrix[thread], CurrentProcessInfo);

        (rCurrentCondition).CalculateDampingMatrix(mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDynamicsToLHS(LHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);

        (rCurrentCondition).EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

// Note: this is in a parallel loop

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        (rCurrentElement).CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        (rCurrentElement).CalculateMassMatrix(mMassMatrix[thread],CurrentProcessInfo);

        (rCurrentElement).CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);

        this->AddDynamicsToLHS(LHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);

        (rCurrentElement).EquationIdVector(EquationId,CurrentProcessInfo);

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

        int thread = OpenMPUtils::ThisThread();

        (rCurrentElement).CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        (rCurrentElement).CalculateMassMatrix(mMassMatrix[thread], CurrentProcessInfo);

        (rCurrentElement).CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);

        this->AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);

        (rCurrentElement).EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationIds,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentCondition.CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);

        rCurrentCondition.CalculateMassMatrix(mMassMatrix[thread], rCurrentProcessInfo);

        rCurrentCondition.CalculateDampingMatrix(mDampingMatrix[thread], rCurrentProcessInfo);

        this->AddDynamicsToRHS(rCurrentCondition, rRHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], rCurrentProcessInfo);

        rCurrentCondition.EquationIdVector(rEquationIds, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop


    void CalculateLHSContribution(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        (rCurrentCondition).CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        (rCurrentCondition).CalculateMassMatrix(mMassMatrix[thread], CurrentProcessInfo);

        (rCurrentCondition).CalculateDampingMatrix(mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDynamicsToLHS(LHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);

        (rCurrentCondition).EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateLHSContribution(
        Element &rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        (rCurrentElement).CalculateLeftHandSide(LHS_Contribution,CurrentProcessInfo);

        (rCurrentElement).CalculateMassMatrix(mMassMatrix[thread],CurrentProcessInfo);

        (rCurrentElement).CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);

        this->AddDynamicsToLHS(LHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);

        (rCurrentElement).EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    std::vector< Matrix > mMassMatrix;
    std::vector< Vector > mAccelerationVector;
    std::vector< Matrix > mDampingMatrix;
    std::vector< Vector > mVelocityVector;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddDynamicsToLHS(LocalSystemMatrixType& LHS_Contribution,
                          LocalSystemMatrixType& M,
                          LocalSystemMatrixType& C,
                          const ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        // adding mass contribution
        if (M.size1() != 0)
            noalias(LHS_Contribution) += 1.0/(mBeta*mDeltaTime*mDeltaTime)*M;

        // adding damping contribution
        if (C.size1() != 0)
            noalias(LHS_Contribution) += mGamma/(mBeta*mDeltaTime)*C;

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddDynamicsToRHS(Condition& rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& M,
        LocalSystemMatrixType& C,
        const ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

            int thread = OpenMPUtils::ThisThread();

        //adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentCondition.GetSecondDerivativesVector(mAccelerationVector[thread], 0);
            noalias(RHS_Contribution) -= prod(M, mAccelerationVector[thread]);
        }

        //adding damping contribution
        if (C.size1() != 0)
        {
            rCurrentCondition.GetFirstDerivativesVector(mVelocityVector[thread], 0);
            noalias(RHS_Contribution) -= prod(C, mVelocityVector[thread]);
        }

        KRATOS_CATCH("")
    }

    void AddDynamicsToRHS(Element &rCurrentElement,
                          LocalSystemVectorType& RHS_Contribution,
                          LocalSystemMatrixType& M,
                          LocalSystemMatrixType& C,
                          const ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentElement.GetSecondDerivativesVector(mAccelerationVector[thread], 0);
            noalias(RHS_Contribution) -= prod(M, mAccelerationVector[thread]);
        }

        //adding damping contribution
        if (C.size1() != 0)
        {
            rCurrentElement.GetFirstDerivativesVector(mVelocityVector[thread], 0);
            noalias(RHS_Contribution) -= prod(C, mVelocityVector[thread]);
        }

        KRATOS_CATCH( "" )
    }

}; // Class NewmarkDynamicUPwScheme
}  // namespace Kratos

#endif // KRATOS_NEWMARK_DYNAMIC_U_PW_SCHEME defined
