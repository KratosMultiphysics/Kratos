//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_STEADY_SCALAR_TRANSPORT_SCHEME)
#define KRATOS_STEADY_SCALAR_TRANSPORT_SCHEME

// System includes
#include <iomanip>
#include <sstream>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "custom_strategies/relaxed_dof_updater.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <class TSparseSpace, class TDenseSpace>
class SteadyScalarScheme
    : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SteadyScalarScheme);

    using BaseType = Scheme<TSparseSpace, TDenseSpace>;

    using DofType = typename BaseType::TDofType;

    using DofsArrayType = typename BaseType::DofsArrayType;

    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    using TSystemVectorType = typename BaseType::TSystemVectorType;

    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    SteadyScalarScheme(
        const double RelaxationFactor)
    : mRelaxationFactor(RelaxationFactor)
    {
        const int num_threads = OpenMPUtils::GetNumThreads();
        mDampingMatrix.resize(num_threads);

        mpDofUpdater = Kratos::make_unique<DofUpdaterType>(mRelaxationFactor);
    }

    ~SteadyScalarScheme() override = default;

    ///@}
    ///@name Operators
    ///@{

    void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dv,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        if (rModelPart.GetBufferSize() > 1) {
            block_for_each(rDofSet, [](DofType& pDof) {
                if (pDof.IsFree()) {
                    const double value = pDof.GetSolutionStepValue(1);
                    pDof.GetSolutionStepValue() = value;
                }
            });
        }

        KRATOS_CATCH("");
    }

    void Update(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb) override
    {
        KRATOS_TRY;

        mpDofUpdater->UpdateDofs(rDofSet, rDx);

        KRATOS_CATCH("");
    }

    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

    void CalculateSystemContributions(
        Element& rElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        const int k = OpenMPUtils::ThisThread();

        this->CalculateDampingSystem(rElement, rLHS_Contribution, rRHS_Contribution,
                                     rEquationIdVector, rCurrentProcessInfo, k);

        if (mDampingMatrix[k].size1() != 0) {
            noalias(rLHS_Contribution) += mDampingMatrix[k];
        }

        KRATOS_CATCH("");
    }

    void CalculateSystemContributions(
        Condition& rCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Condition::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        const int k = OpenMPUtils::ThisThread();

        this->CalculateDampingSystem(rCondition, rLHS_Contribution, rRHS_Contribution,
                                     rEquationIdVector, rCurrentProcessInfo, k);

        if (mDampingMatrix[k].size1() != 0) {
            noalias(rLHS_Contribution) += mDampingMatrix[k];
        }

        KRATOS_CATCH("");
    }

    void CalculateRHSContribution(
        Element& rElement,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        const int k = OpenMPUtils::ThisThread();

        this->CalculateDampingSystem(rElement, mDampingMatrix[k], rRHS_Contribution,
                                     rEquationIdVector, rCurrentProcessInfo, k);

        KRATOS_CATCH("");
    }

    void CalculateRHSContribution(
        Condition& rCondition,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        const int k = OpenMPUtils::ThisThread();

        this->CalculateDampingSystem(rCondition, mDampingMatrix[k], rRHS_Contribution,
                                     rEquationIdVector, rCurrentProcessInfo, k);

        KRATOS_CATCH("");
    }

    void InitializeDofUpdater() {
        mpDofUpdater->InitializeAitken();
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.

    std::string Info() const override
    {
        std::stringstream msg;
        msg << "Using generic residual based steady scalar transport scheme "
               "with\n"
            << "     Relaxation factor           : " << this->mRelaxationFactor;
        return msg.str();
    }

    ///@}

protected:
    ///@name Protected Member Variables
    ///@{

    std::vector<LocalSystemMatrixType> mDampingMatrix;
    double mRelaxationFactor;

    ///@}
    ///@name Protected Operators
    ///@{

    template <class TItemType>
    void CalculateDampingSystem(
        TItemType& rItem,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        typename TItemType::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo,
        const int ThreadId)
    {
        KRATOS_TRY;

        rItem.CalculateLocalSystem(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
        rItem.CalculateLocalVelocityContribution(
            mDampingMatrix[ThreadId], rRHS_Contribution, rCurrentProcessInfo);
        rItem.EquationIdVector(rEquationIdVector, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    using DofUpdaterType = RelaxedDofUpdater<TSparseSpace>;
    using DofUpdaterPointerType = typename DofUpdaterType::UniquePointer;

    DofUpdaterPointerType mpDofUpdater;

    ///@}

private:
    ///@name Member Variables
    ///@{



    ///@}
};

///@}

} // namespace Kratos

#endif /* KRATOS_STEADY_SCALAR_TRANSPORT_SCHEME defined */
