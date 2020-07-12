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

#if !defined(KRATOS_GENERIC_RESIDUAL_BASED_BOSSAK_VELOCITY_SCALAR_SCHEME_H_INCLUDED)
#define KRATOS_GENERIC_RESIDUAL_BASED_BOSSAK_VELOCITY_SCALAR_SCHEME_H_INCLUDED

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
class GenericResidualBasedBossakVelocityScalarScheme
    : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(GenericResidualBasedBossakVelocityScalarScheme);

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

    GenericResidualBasedBossakVelocityScalarScheme(const double AlphaBossak,
                                                   const double RelaxationFactor,
                                                   const Variable<double>& rScalarVariable,
                                                   const Variable<double>& rScalarRateVariable,
                                                   const Variable<double>& rRelaxedScalarRateVariable)
        : mAlphaBossak(AlphaBossak),
          mUpdateAcceleration(true),
          mRelaxationFactor(RelaxationFactor),
          mrScalarVariable(rScalarVariable),
          mrScalarRateVariable(rScalarRateVariable),
          mrRelaxedScalarRateVariable(rRelaxedScalarRateVariable)
    {
        KRATOS_INFO("GenericResidualBasedBossakVelocityScalarScheme")
            << " Using bossak velocity scheme with alpha_bossak = " << std::scientific
            << mAlphaBossak << " [UpdateAcceleration: " << mUpdateAcceleration << "]\n";

        // Allocate auxiliary memory.
        const int num_threads = OpenMPUtils::GetNumThreads();

        mMassMatrix.resize(num_threads);
        mDampingMatrix.resize(num_threads);
        mValuesVector.resize(num_threads);
        mSecondDerivativeValuesVector.resize(num_threads);
        mSecondDerivativeValuesVectorOld.resize(num_threads);
    }

    /// Destructor.
    ~GenericResidualBasedBossakVelocityScalarScheme() override = default;

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

        GenericResidualBasedBossakVelocityScalarScheme::CalculateBossakConstants(
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

        this->UpdateScalarRateVariables(rModelPart);

        KRATOS_CATCH("");
    }

    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        int value = BaseType::Check(rModelPart);

        KRATOS_ERROR_IF(!rModelPart.HasNodalSolutionStepVariable(mrScalarVariable))
            << mrScalarVariable.Name() << " not in nodal solution step variable list of "
            << rModelPart.Name() << ".\n";
        KRATOS_ERROR_IF(!rModelPart.HasNodalSolutionStepVariable(mrScalarRateVariable))
            << mrScalarRateVariable.Name() << " not in nodal solution step variable list of "
            << rModelPart.Name() << ".\n";
        KRATOS_ERROR_IF(!rModelPart.HasNodalSolutionStepVariable(mrRelaxedScalarRateVariable))
            << mrRelaxedScalarRateVariable.Name() << " not in nodal solution step variable list of "
            << rModelPart.Name() << ".\n";

        return value;
        KRATOS_CATCH("");
    }

    void CalculateSystemContributions(Element& rElement,
                                      LocalSystemMatrixType& rLHS_Contribution,
                                      LocalSystemVectorType& rRHS_Contribution,
                                      Element::EquationIdVectorType& rEquationIdVector,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateDynamicSystem<Element>(rElement, rLHS_Contribution, rRHS_Contribution,
                                        rEquationIdVector, rCurrentProcessInfo);
    }

    void CalculateRHSContribution(Element& rElement,
                                  LocalSystemVectorType& rRHS_Contribution,
                                  Element::EquationIdVectorType& rEquationIdVector,
                                  const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateDynamicRHS<Element>(rElement, rRHS_Contribution,
                                     rEquationIdVector, rCurrentProcessInfo);
    }

    void CalculateSystemContributions(Condition& rCondition,
                                      LocalSystemMatrixType& rLHS_Contribution,
                                      LocalSystemVectorType& rRHS_Contribution,
                                      Condition::EquationIdVectorType& rEquationIdVector,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateDynamicSystem<Condition>(rCondition, rLHS_Contribution, rRHS_Contribution,
                                          rEquationIdVector, rCurrentProcessInfo);
    }

    void CalculateRHSContribution(Condition& rCondition,
                                  LocalSystemVectorType& rRHS_Contribution,
                                  Element::EquationIdVectorType& rEquationIdVector,
                                  const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateDynamicRHS<Condition>(rCondition, rRHS_Contribution,
                                       rEquationIdVector, rCurrentProcessInfo);
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
        return "GenericResidualBasedBossakVelocityScalarScheme";
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

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

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

    const Variable<double>& mrScalarVariable;
    const Variable<double>& mrScalarRateVariable;
    const Variable<double>& mrRelaxedScalarRateVariable;

    BossakConstants mBossak;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    template <class TItemType>
    void CalculateDampingSystem(TItemType& rItem,
                                LocalSystemMatrixType& rLHS_Contribution,
                                LocalSystemVectorType& rRHS_Contribution,
                                typename TItemType::EquationIdVectorType& rEquationIdVector,
                                const ProcessInfo& rCurrentProcessInfo,
                                const int ThreadId)
    {
        KRATOS_TRY;

        rItem.InitializeNonLinearIteration(rCurrentProcessInfo);
        rItem.CalculateLocalSystem(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
        rItem.CalculateLocalVelocityContribution(
            mDampingMatrix[ThreadId], rRHS_Contribution, rCurrentProcessInfo);
        rItem.EquationIdVector(rEquationIdVector, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    template <class TItemType>
    void AddMassMatrixToRHS(TItemType& rItem, LocalSystemVectorType& rRHS_Contribution, const int ThreadId)
    {
        rItem.GetSecondDerivativesVector(mSecondDerivativeValuesVector[ThreadId], 0);
        (mSecondDerivativeValuesVector[ThreadId]) *= (1.00 - mBossak.Alpha);
        rItem.GetSecondDerivativesVector(mSecondDerivativeValuesVectorOld[ThreadId], 1);
        noalias(mSecondDerivativeValuesVector[ThreadId]) +=
            mBossak.Alpha * mSecondDerivativeValuesVectorOld[ThreadId];

        noalias(rRHS_Contribution) -=
            prod(mMassMatrix[ThreadId], mSecondDerivativeValuesVector[ThreadId]);
    }

    template <class TItemType>
    void CalculateDynamicSystem(TItemType& rItem,
                                LocalSystemMatrixType& rLHS_Contribution,
                                LocalSystemVectorType& rRHS_Contribution,
                                typename TItemType::EquationIdVectorType& rEquationIdVector,
                                const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        this->CalculateSteadySystem<TItemType>(rItem, rLHS_Contribution, rRHS_Contribution,
                                               rEquationIdVector, rCurrentProcessInfo);

        const int k = OpenMPUtils::ThisThread();

        rItem.CalculateMassMatrix(mMassMatrix[k], rCurrentProcessInfo);
        // adding mass contribution to the dynamic stiffness
        if (mMassMatrix[k].size1() != 0) // if M matrix declared
        {
            AddMassMatrixToRHS<TItemType>(rItem, rRHS_Contribution, k);

            noalias(rLHS_Contribution) += mBossak.C0 * mMassMatrix[k];
        }

        KRATOS_CATCH("");
    }

    template <class TItemType>
    void CalculateSteadySystem(TItemType& rItem,
                               LocalSystemMatrixType& rLHS_Contribution,
                               LocalSystemVectorType& rRHS_Contribution,
                               typename TItemType::EquationIdVectorType& rEquationIdVector,
                               const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const int k = OpenMPUtils::ThisThread();

        CalculateDampingSystem<TItemType>(rItem, rLHS_Contribution, rRHS_Contribution,
                                          rEquationIdVector, rCurrentProcessInfo, k);

        if (mDampingMatrix[k].size1() != 0)
        {
            noalias(rLHS_Contribution) += mDampingMatrix[k];
        }

        KRATOS_CATCH("");
    }

    template <class TItemType>
    void CalculateSteadyRHS(TItemType& rItem,
                            LocalSystemVectorType& rRHS_Contribution,
                            typename TItemType::EquationIdVectorType& rEquationIdVector,
                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const int k = OpenMPUtils::ThisThread();

        CalculateDampingSystem<TItemType>(rItem, mDampingMatrix[k], rRHS_Contribution,
                                          rEquationIdVector, rCurrentProcessInfo, k);

        KRATOS_CATCH("");
    }

    template <class TItemType>
    void CalculateDynamicRHS(TItemType& rItem,
                             LocalSystemVectorType& rRHS_Contribution,
                             typename TItemType::EquationIdVectorType& rEquationIdVector,
                             const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const int k = OpenMPUtils::ThisThread();

        CalculateDampingSystem<TItemType>(rItem, mDampingMatrix[k], rRHS_Contribution,
                                          rEquationIdVector, rCurrentProcessInfo, k);

        rItem.CalculateMassMatrix(mMassMatrix[k], rCurrentProcessInfo);
        // adding mass contribution to the dynamic stiffness
        if (mMassMatrix[k].size1() != 0) // if M matrix declared
        {
            AddMassMatrixToRHS<TItemType>(rItem, rRHS_Contribution, k);
        }

        KRATOS_CATCH("");
    }

    // class to hold all the derivatives for updated target variable

    void UpdateScalarRateVariables(ModelPart& rModelPart)
    {
        if (!mUpdateAcceleration)
            return;
        const int number_of_nodes = rModelPart.NumberOfNodes();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(rModelPart.NodesBegin() + i_node);
            double& r_current_rate = r_node.FastGetSolutionStepValue(mrScalarRateVariable);
            const double old_rate =
                r_node.FastGetSolutionStepValue(mrScalarRateVariable, 1);
            const double current_value = r_node.FastGetSolutionStepValue(mrScalarVariable);
            const double old_value =
                r_node.FastGetSolutionStepValue(mrScalarVariable, 1);

            // update scalar rate variable
            r_current_rate = mBossak.C2 * (current_value - old_value) - mBossak.C3 * old_rate;

            // update relaxed scalar rate variable
            r_node.FastGetSolutionStepValue(mrRelaxedScalarRateVariable) =
                this->mAlphaBossak * old_rate + (1.0 - this->mAlphaBossak) * r_current_rate;
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

}; /* Class GenericResidualBasedBossakVelocityScalarScheme */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_GENERIC_RESIDUAL_BASED_BOSSAK_VELOCITY_SCALAR_SCHEME_H_INCLUDED defined */
