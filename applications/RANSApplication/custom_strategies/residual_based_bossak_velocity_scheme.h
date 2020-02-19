//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Jordi Cotela
//                  Suneth Warnakulasuriya
//

#if !defined(KRATOS_RESIDUAL_BASED_BOSSAK_VELOCITY_SCHEME_H_INCLUDED)
#define KRATOS_RESIDUAL_BASED_BOSSAK_VELOCITY_SCHEME_H_INCLUDED

// System includes
#include <limits>
#include <vector>

// External includes

// Project includes
#include "custom_strategies/relaxed_dof_updater.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/time_discretization.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// A scheme for steady and dynamic equations, using Bossak time integration.
/**
 * It can be used for either first- or second-order time derivatives. Elements
 * and conditions must provide a specialization of SchemeExtension via
 * their data value container, which allows the scheme to operate independently
 * of the variable arrangements in the element or condition.
 */
template <class TSparseSpace, class TDenseSpace>
class ResidualBasedBossakVelocityScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBossakVelocityScheme);

    using BaseType = Scheme<TSparseSpace, TDenseSpace>;

    using SystemMatrixType = typename BaseType::TSystemMatrixType;

    using SystemVectorType = typename BaseType::TSystemVectorType;

    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    using DofsArrayType = typename BaseType::DofsArrayType;

    using NodeType = ModelPart::NodeType;

    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.

    ResidualBasedBossakVelocityScheme(
        const double AlphaBossak,
        const double RelaxationFactor,
        const std::vector<Variable<double> const*> rDisplacementVariables,
        const std::vector<Variable<double> const*> rVelocityVariables,
        const std::vector<Variable<double> const*> rAccelerationVariables,
        const std::vector<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> const*> rDisplacementComponentVariables,
        const std::vector<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> const*> rVelocityComponentVariables,
        const std::vector<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> const*> rAccelerationComponentVariables)
        : mAlphaBossak(AlphaBossak),
          mUpdateAcceleration(rAccelerationVariables.size() > 0 ||
                              rAccelerationComponentVariables.size() > 0),
          mUpdateDisplacement(rDisplacementVariables.size() > 0 ||
                              rDisplacementComponentVariables.size() > 0),
          mRelaxationFactor(RelaxationFactor),
          mDisplacementVariables(rDisplacementVariables),
          mVelocityVariables(rVelocityVariables),
          mAccelerationVariables(rAccelerationVariables),
          mDisplacementComponentVariables(rDisplacementComponentVariables),
          mVelocityComponentVariables(rVelocityComponentVariables),
          mAccelerationComponentVariables(rAccelerationComponentVariables)
    {
        KRATOS_INFO("ResidualBasedBossakVelocityScheme")
            << " Using bossak velocity scheme with alpha_bossak = " << std::scientific
            << mAlphaBossak << " [UpdateAcceleration: " << mUpdateAcceleration
            << ", UpdateDisplacement: " << mUpdateDisplacement << "]\n";

        // Allocate auxiliary memory.
        const int num_threads = OpenMPUtils::GetNumThreads();

        mMassMatrix.resize(num_threads);
        mDampingMatrix.resize(num_threads);
        mValuesVector.resize(num_threads);
        mSecondDerivativeValuesVector.resize(num_threads);
        mSecondDerivativeValuesVectorOld.resize(num_threads);
    }

    /// Destructor.
    ~ResidualBasedBossakVelocityScheme() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void InitializeSolutionStep(ModelPart& rModelPart,
                                SystemMatrixType& rA,
                                SystemVectorType& rDx,
                                SystemVectorType& rb) override
    {
        KRATOS_TRY;

        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];

        KRATOS_ERROR_IF(delta_time < std::numeric_limits<double>::epsilon())
            << "detected delta_time = 0 in the Bossak Scheme ... "
               "check if the time step is created correctly for "
               "the current model part.";

        ResidualBasedBossakVelocityScheme::CalculateBossakConstants(
            mBossak, mAlphaBossak, delta_time);

#pragma omp critical
        {
            rModelPart.GetProcessInfo()[BOSSAK_ALPHA] = mBossak.Alpha;
        }

        KRATOS_CATCH("");
    }

    void Update(ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                SystemMatrixType& rA,
                SystemVectorType& rDx,
                SystemVectorType& rb) override
    {
        KRATOS_TRY;

        mpDofUpdater->UpdateDofs(rDofSet, rDx, mRelaxationFactor);

        this->UpdateTimeSchemeVariables(rModelPart);

        KRATOS_CATCH("");
    }

    void CalculateSystemContributions(Element::Pointer pCurrentElement,
                                      LocalSystemMatrixType& rLHS_Contribution,
                                      LocalSystemVectorType& rRHS_Contribution,
                                      Element::EquationIdVectorType& rEquationId,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        const int k = OpenMPUtils::ThisThread();

        (pCurrentElement)->InitializeNonLinearIteration(rCurrentProcessInfo);
        (pCurrentElement)->CalculateLocalSystem(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
        (pCurrentElement)->CalculateLocalVelocityContribution(mDampingMatrix[k], rRHS_Contribution, rCurrentProcessInfo);

        if (mUpdateAcceleration)
        {
            (pCurrentElement)->CalculateMassMatrix(mMassMatrix[k], rCurrentProcessInfo);
            AddDynamicsToRHS(pCurrentElement, rRHS_Contribution, mDampingMatrix[k],
                             mMassMatrix[k], rCurrentProcessInfo);
        }
        AddDynamicsToLHS(rLHS_Contribution, mDampingMatrix[k], mMassMatrix[k],
                         rCurrentProcessInfo);

        (pCurrentElement)->EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Calculate_RHS_Contribution(Element::Pointer pCurrentElement,
                                    LocalSystemVectorType& rRHS_Contribution,
                                    Element::EquationIdVectorType& rEquationId,
                                    ProcessInfo& rCurrentProcessInfo) override
    {
        const int k = OpenMPUtils::ThisThread();

        // Initializing the non linear iteration for the current element
        (pCurrentElement)->InitializeNonLinearIteration(rCurrentProcessInfo);

        // basic operations for the element considered
        (pCurrentElement)->CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);
        (pCurrentElement)->CalculateLocalVelocityContribution(mDampingMatrix[k], rRHS_Contribution, rCurrentProcessInfo);
        (pCurrentElement)->EquationIdVector(rEquationId, rCurrentProcessInfo);

        // adding the dynamic contributions (static is already included)
        if (mUpdateAcceleration)
        {
            (pCurrentElement)->CalculateMassMatrix(mMassMatrix[k], rCurrentProcessInfo);
            AddDynamicsToRHS(pCurrentElement, rRHS_Contribution, mDampingMatrix[k],
                             mMassMatrix[k], rCurrentProcessInfo);
        }
    }

    void Condition_CalculateSystemContributions(Condition::Pointer pCurrentCondition,
                                                LocalSystemMatrixType& rLHS_Contribution,
                                                LocalSystemVectorType& rRHS_Contribution,
                                                Condition::EquationIdVectorType& rEquationId,
                                                ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY
        const int k = OpenMPUtils::ThisThread();

        (pCurrentCondition)->InitializeNonLinearIteration(rCurrentProcessInfo);
        (pCurrentCondition)->CalculateLocalSystem(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
        (pCurrentCondition)->CalculateLocalVelocityContribution(mDampingMatrix[k], rRHS_Contribution, rCurrentProcessInfo);
        (pCurrentCondition)->EquationIdVector(rEquationId, rCurrentProcessInfo);

        if (mUpdateAcceleration)
        {
            (pCurrentCondition)->CalculateMassMatrix(mMassMatrix[k], rCurrentProcessInfo);
            AddDynamicsToRHS(pCurrentCondition, rRHS_Contribution,
                             mDampingMatrix[k], mMassMatrix[k], rCurrentProcessInfo);
        }
        AddDynamicsToLHS(rLHS_Contribution, mDampingMatrix[k], mMassMatrix[k],
                         rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void Condition_Calculate_RHS_Contribution(Condition::Pointer pCurrentCondition,
                                              LocalSystemVectorType& rRHS_Contribution,
                                              Element::EquationIdVectorType& rEquationId,
                                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        const int k = OpenMPUtils::ThisThread();

        (pCurrentCondition)->InitializeNonLinearIteration(rCurrentProcessInfo);
        (pCurrentCondition)->CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);
        (pCurrentCondition)->CalculateLocalVelocityContribution(mDampingMatrix[k], rRHS_Contribution, rCurrentProcessInfo);
        (pCurrentCondition)->EquationIdVector(rEquationId, rCurrentProcessInfo);

        // adding the dynamic contributions (static is already included)
        if (mUpdateAcceleration)
        {
            (pCurrentCondition)->CalculateMassMatrix(mMassMatrix[k], rCurrentProcessInfo);
            AddDynamicsToRHS(pCurrentCondition, rRHS_Contribution,
                             mDampingMatrix[k], mMassMatrix[k], rCurrentProcessInfo);
        }
        KRATOS_CATCH("");
    }

    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ResidualBasedBossakVelocityScheme";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }
    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    struct BossakConstants
    {
        double Alpha;
        double Gamma;
        double Beta;
        double C0;
        double C1;
        double C2;
        double C3;
        double C4;
        double C5;
        double C6;
    };

    ///@}
    ///@name Protected member Variables
    ///@{

    std::vector<LocalSystemVectorType> mSecondDerivativeValuesVectorOld;
    std::vector<LocalSystemVectorType> mSecondDerivativeValuesVector;
    std::vector<LocalSystemVectorType> mValuesVector;

    std::vector<LocalSystemMatrixType> mMassMatrix;
    std::vector<LocalSystemMatrixType> mDampingMatrix;

    const double mAlphaBossak;
    bool mUpdateAcceleration;
    bool mUpdateDisplacement;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    //****************************************************************************

    /**
    Kdyn = am*M + D + a1*K
     */
    void AddDynamicsToLHS(LocalSystemMatrixType& rLHS_Contribution,
                          LocalSystemMatrixType& rDampingMatrix,
                          LocalSystemMatrixType& rMassMatrix,
                          ProcessInfo& CurrentProcessInfo)
    {
        // multipling time scheme factor
        rLHS_Contribution *= mBossak.C1;

        // adding mass contribution to the dynamic stiffness
        if (rMassMatrix.size1() != 0 && mUpdateAcceleration) // if M matrix declared
        {
            noalias(rLHS_Contribution) += mBossak.C0 * rMassMatrix;
        }

        // adding  damping contribution
        if (rDampingMatrix.size1() != 0) // if M matrix declared
        {
            noalias(rLHS_Contribution) += rDampingMatrix;
        }
    }

    //****************************************************************************

    /// Add Bossak contributions from the inertial term to the RHS vector.
    /** This essentially performs bdyn = b - M*acc for the current element.
     *  Note that viscous/pressure contributions to the RHS are expected to be added by the element itself.
     *  @param[in] rCurrentElement The fluid element we are assembling.
     *  @param[in/out] rRHS_Contribution The right hand side term where the contribution will be added.
     *  @param[in] rD The elemental velocity/pressure LHS matrix.
     *  @param[in] rM The elemental acceleration LHS matrix.
     *  @param[in] rCurrentProcessInfo ProcessInfo instance for the containing ModelPart.
     */
    void AddDynamicsToRHS(Element::Pointer rCurrentElement,
                          LocalSystemVectorType& rRHS_Contribution,
                          LocalSystemMatrixType& rDampingMatrix,
                          LocalSystemMatrixType& rMassMatrix,
                          ProcessInfo& rCurrentProcessInfo)
    {
        // adding inertia contribution
        if (rMassMatrix.size1() != 0)
        {
            const int k = OpenMPUtils::ThisThread();
            rCurrentElement->GetSecondDerivativesVector(
                mSecondDerivativeValuesVector[k], 0);
            (mSecondDerivativeValuesVector[k]) *= (1.00 - mBossak.Alpha);
            rCurrentElement->GetSecondDerivativesVector(
                mSecondDerivativeValuesVectorOld[k], 1);
            noalias(mSecondDerivativeValuesVector[k]) +=
                mBossak.Alpha * mSecondDerivativeValuesVectorOld[k];
            noalias(rRHS_Contribution) -=
                prod(rMassMatrix, mSecondDerivativeValuesVector[k]);
        }
    }

    /// Add Bossak contributions from the inertial term to the RHS vector.
    /** This essentially performs bdyn = b - M*acc for the current condition.
     *  Note that viscous/pressure contributions to the RHS are expected to be added by the element condition.
     *  @param[in] rCurrentCondition The fluid condition we are assembling.
     *  @param[in/out] rRHS_Contribution The right hand side term where the contribution will be added.
     *  @param[in] rD The elemental velocity/pressure LHS matrix.
     *  @param[in] rM The elemental acceleration LHS matrix.
     *  @param[in] rCurrentProcessInfo ProcessInfo instance for the containing ModelPart.
     */
    void AddDynamicsToRHS(Condition::Pointer rCurrentCondition,
                          LocalSystemVectorType& rRHS_Contribution,
                          LocalSystemMatrixType& rDampingMatrix,
                          LocalSystemMatrixType& rMassMatrix,
                          ProcessInfo& rCurrentProcessInfo)
    {
        // adding inertia contribution
        if (rMassMatrix.size1() != 0)
        {
            const int k = OpenMPUtils::ThisThread();
            rCurrentCondition->GetSecondDerivativesVector(
                mSecondDerivativeValuesVector[k], 0);
            (mSecondDerivativeValuesVector[k]) *= (1.00 - mBossak.Alpha);
            rCurrentCondition->GetSecondDerivativesVector(
                mSecondDerivativeValuesVectorOld[k], 1);
            noalias(mSecondDerivativeValuesVector[k]) +=
                mBossak.Alpha * mSecondDerivativeValuesVectorOld[k];

            noalias(rRHS_Contribution) -=
                prod(rMassMatrix, mSecondDerivativeValuesVector[k]);
        }
    }

    void UpdateTimeSchemeVariables(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        UpdateAcceleration<Variable<double>>(rModelPart, mVelocityVariables,
                                             mAccelerationVariables);
        UpdateAcceleration<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>(
            rModelPart, mVelocityComponentVariables, mAccelerationComponentVariables);
        UpdateDisplacement<Variable<double>>(rModelPart, mDisplacementVariables,
                                             mVelocityVariables, mAccelerationVariables);
        UpdateDisplacement<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>(
            rModelPart, mDisplacementComponentVariables,
            mVelocityComponentVariables, mAccelerationComponentVariables);

        KRATOS_CATCH("");
    }

    void UpdateAcceleration(double& rCurrentAcceleration,
                            const double CurrentVelocity,
                            const double OldVelocity,
                            const double OldAcceleration) const
    {
        rCurrentAcceleration = mBossak.C2 * (CurrentVelocity - OldVelocity) -
                               mBossak.C3 * OldAcceleration;
    }

    void UpdateDisplacement(double& rCurrentDisplacement,
                            const double OldDisplacement,
                            const double OldVelocity,
                            const double CurrentAcceleration,
                            const double OldAcceleration) const
    {
        rCurrentDisplacement = OldDisplacement + mBossak.C6 * OldVelocity +
                               mBossak.C4 * OldAcceleration + mBossak.C5 * CurrentAcceleration;
    }

    static void CalculateBossakConstants(BossakConstants& rBossakConstants,
                                         const double Alpha,
                                         const double DeltaTime)
    {
        TimeDiscretization::Bossak bossak(Alpha, 0.25, 0.5);
        rBossakConstants.Alpha = bossak.GetAlphaM();
        rBossakConstants.Gamma = bossak.GetGamma();
        rBossakConstants.Beta = bossak.GetBeta();

        rBossakConstants.C0 =
            (1.0 - rBossakConstants.Alpha) / (rBossakConstants.Gamma * DeltaTime);
        rBossakConstants.C1 =
            DeltaTime / (rBossakConstants.Beta * rBossakConstants.Gamma);
        rBossakConstants.C2 = 1.0 / (rBossakConstants.Gamma * DeltaTime);
        rBossakConstants.C3 = (1.0 - rBossakConstants.Gamma) / rBossakConstants.Gamma;
        rBossakConstants.C4 =
            std::pow(DeltaTime, 2) * (-2.0 * rBossakConstants.Beta + 1.0) / 2.0;
        rBossakConstants.C5 = std::pow(DeltaTime, 2) * rBossakConstants.Beta;
        rBossakConstants.C6 = DeltaTime;
    }

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    using DofUpdaterType = RelaxedDofUpdater<TSparseSpace>;
    using DofUpdaterPointerType = typename DofUpdaterType::UniquePointer;

    DofUpdaterPointerType mpDofUpdater = Kratos::make_unique<DofUpdaterType>();

    double mRelaxationFactor;

    const std::vector<Variable<double> const*> mDisplacementVariables;
    const std::vector<Variable<double> const*> mVelocityVariables;
    const std::vector<Variable<double> const*> mAccelerationVariables;

    const std::vector<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> const*> mDisplacementComponentVariables;
    const std::vector<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> const*> mVelocityComponentVariables;
    const std::vector<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> const*> mAccelerationComponentVariables;

    BossakConstants mBossak;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    // class to hold all the derivatives for updated target variable

    template <class TVariableType>
    void UpdateAcceleration(ModelPart& rModelPart,
                            const std::vector<TVariableType const*>& pVelocityVariables,
                            const std::vector<TVariableType const*>& pAccelerationVariables)
    {
        if (!mUpdateAcceleration)
            return;
        const int number_of_nodes = rModelPart.NumberOfNodes();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(rModelPart.NodesBegin() + i_node);
            for (IndexType i_var = 0; i_var < pAccelerationVariables.size(); ++i_var)
            {
                double& r_current_acceleration =
                    r_node.FastGetSolutionStepValue(*pAccelerationVariables[i_var]);
                const double old_acceleration = r_node.FastGetSolutionStepValue(
                    *pAccelerationVariables[i_var], 1);
                const double current_velocity =
                    r_node.FastGetSolutionStepValue(*pVelocityVariables[i_var]);
                const double old_velocity =
                    r_node.FastGetSolutionStepValue(*pVelocityVariables[i_var], 1);
                UpdateAcceleration(r_current_acceleration, current_velocity,
                                   old_velocity, old_acceleration);
            }
        }
    }

    template <class TVariableType>
    void UpdateDisplacement(ModelPart& rModelPart,
                            const std::vector<TVariableType const*>& pDisplacementVariables,
                            const std::vector<TVariableType const*>& pVelocityVariables,
                            const std::vector<TVariableType const*>& pAccelerationVariables)
    {
        if (!mUpdateDisplacement)
            return;
        const int number_of_nodes = rModelPart.NumberOfNodes();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(rModelPart.NodesBegin() + i_node);
            for (IndexType i_var = 0; i_var < pDisplacementVariables.size(); ++i_var)
            {
                double& r_current_displacement =
                    r_node.FastGetSolutionStepValue(*pDisplacementVariables[i_var]);
                const double old_displacement = r_node.FastGetSolutionStepValue(
                    *pDisplacementVariables[i_var], 1);
                const double current_acceleration =
                    r_node.FastGetSolutionStepValue(*pAccelerationVariables[i_var]);
                const double old_acceleration = r_node.FastGetSolutionStepValue(
                    *pAccelerationVariables[i_var], 1);
                const double old_velocity =
                    r_node.FastGetSolutionStepValue(*pVelocityVariables[i_var], 1);
                UpdateDisplacement(r_current_displacement, old_displacement, old_velocity,
                                   current_acceleration, old_acceleration);
            }
        }
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; /* Class ResidualBasedBossakVelocityScheme */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BOSSAK_VELOCITY_SCHEME_H_INCLUDED defined */
