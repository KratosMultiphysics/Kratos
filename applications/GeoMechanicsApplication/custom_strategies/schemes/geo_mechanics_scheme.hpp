// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//
#pragma once
#include "solving_strategies/schemes/scheme.h"

namespace Kratos {

template <class TSparseSpace, class TDenseSpace>
class GeoMechanicsScheme : public Scheme<TSparseSpace, TDenseSpace> {
public:
    void GetDofList(const Element& rElement,
                    Element::DofsVectorType& rDofList,
                    const ProcessInfo& rCurrentProcessInfo) override
    {
        const bool isActive = (rElement.IsDefined(ACTIVE)) ? rElement.Is(ACTIVE) : true;
        if (isActive)
            rElement.GetDofList(rDofList, rCurrentProcessInfo);
    }

    void GetDofList(const Condition& rCondition,
                    Element::DofsVectorType& rDofList,
                    const ProcessInfo& rCurrentProcessInfo) override
    {
        const bool isActive = (rCondition.IsDefined(ACTIVE)) ? rCondition.Is(ACTIVE) : true;
        if (isActive)
            rCondition.GetDofList(rDofList, rCurrentProcessInfo);
    }

    void EquationId(const Element& rElement,
                    Element::EquationIdVectorType& rEquationId,
                    const ProcessInfo& rCurrentProcessInfo) override
    {
        const bool isActive = (rElement.IsDefined(ACTIVE)) ? rElement.Is(ACTIVE) : true;
        if (isActive)
            rElement.EquationIdVector(rEquationId, rCurrentProcessInfo);
    }

    void EquationId(const Condition& rCondition,
                    Element::EquationIdVectorType& rEquationId,
                    const ProcessInfo& rCurrentProcessInfo) override
    {
        const bool isActive = (rCondition.IsDefined(ACTIVE)) ? rCondition.Is(ACTIVE) : true;
        if (isActive)
            rCondition.EquationIdVector(rEquationId, rCurrentProcessInfo);
    }

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        SetTimeFactors(rModelPart);

        Scheme<TSparseSpace, TDenseSpace>::mSchemeIsInitialized = true;

        KRATOS_CATCH("")
    }

    void InitializeSolutionStep(ModelPart& rModelPart,
                                TSystemMatrixType& A,
                                TSystemVectorType& Dx,
                                TSystemVectorType& b) override
    {
        KRATOS_TRY

        SetTimeFactors(rModelPart);

        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        block_for_each(rModelPart.Elements(), [&rCurrentProcessInfo](Element& rElement) {
            const bool isActive = (rElement.IsDefined(ACTIVE)) ? rElement.Is(ACTIVE) : true;
            if (isActive)
                rElement.InitializeSolutionStep(rCurrentProcessInfo);
        });

        block_for_each(rModelPart.Conditions(), [&rCurrentProcessInfo](Condition& rCondition) {
            const bool isActive =
                (rCondition.IsDefined(ACTIVE)) ? rCondition.Is(ACTIVE) : true;
            if (isActive)
                rCondition.InitializeSolutionStep(rCurrentProcessInfo);
        });

        KRATOS_CATCH("")
    }

    void Predict(ModelPart& rModelPart,
                 DofsArrayType& rDofSet,
                 TSystemMatrixType& A,
                 TSystemVectorType& Dx,
                 TSystemVectorType& b) override
    {
        this->UpdateVariablesDerivatives(rModelPart);
    }

    void InitializeNonLinIteration(ModelPart& rModelPart,
                                   TSystemMatrixType& A,
                                   TSystemVectorType& Dx,
                                   TSystemVectorType& b) override
    {
        KRATOS_TRY

        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        block_for_each(rModelPart.Elements(), [&rCurrentProcessInfo](Element& rElement) {
            const bool isActive = (rElement.IsDefined(ACTIVE)) ? rElement.Is(ACTIVE) : true;
            if (isActive)
                rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
        });

        block_for_each(rModelPart.Conditions(), [&rCurrentProcessInfo](Condition& rCondition) {
            const bool isActive =
                (rCondition.IsDefined(ACTIVE)) ? rCondition.Is(ACTIVE) : true;
            if (isActive)
                rCondition.InitializeNonLinearIteration(rCurrentProcessInfo);
        });

        KRATOS_CATCH("")
    }

    void FinalizeNonLinIteration(ModelPart& rModelPart,
                                 TSystemMatrixType& A,
                                 TSystemVectorType& Dx,
                                 TSystemVectorType& b) override
    {
        KRATOS_TRY

        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        block_for_each(rModelPart.Elements(), [&rCurrentProcessInfo](Element& rElement) {
            const bool isActive = (rElement.IsDefined(ACTIVE)) ? rElement.Is(ACTIVE) : true;
            if (isActive)
                rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
        });

        block_for_each(rModelPart.Conditions(), [&rCurrentProcessInfo](Condition& rCondition) {
            const bool isActive =
                (rCondition.IsDefined(ACTIVE)) ? rCondition.Is(ACTIVE) : true;
            if (isActive)
                rCondition.FinalizeNonLinearIteration(rCurrentProcessInfo);
        });

        KRATOS_CATCH("")
    }

    void FinalizeSolutionStepActiveEntities(ModelPart& rModelPart,
                                            TSystemMatrixType& A,
                                            TSystemVectorType& Dx,
                                            TSystemVectorType& b)
    {
        KRATOS_TRY

        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        block_for_each(rModelPart.Elements(), [&rCurrentProcessInfo](Element& rElement) {
            const bool isActive = (rElement.IsDefined(ACTIVE)) ? rElement.Is(ACTIVE) : true;
            if (isActive)
                rElement.FinalizeSolutionStep(rCurrentProcessInfo);
        });

        block_for_each(rModelPart.Conditions(), [&rCurrentProcessInfo](Condition& rCondition) {
            const bool isActive =
                (rCondition.IsDefined(ACTIVE)) ? rCondition.Is(ACTIVE) : true;
            if (isActive)
                rCondition.FinalizeSolutionStep(rCurrentProcessInfo);
        });

        KRATOS_CATCH("")
    }

    void CalculateSystemContributions(Element& rCurrentElement,
                                      LocalSystemMatrixType& LHS_Contribution,
                                      LocalSystemVectorType& RHS_Contribution,
                                      Element::EquationIdVectorType& EquationId,
                                      const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        rCurrentElement.CalculateLocalSystem(LHS_Contribution, RHS_Contribution,
                                             CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateSystemContributions(Condition& rCurrentCondition,
                                      LocalSystemMatrixType& LHS_Contribution,
                                      LocalSystemVectorType& RHS_Contribution,
                                      Element::EquationIdVectorType& EquationId,
                                      const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        rCurrentCondition.CalculateLocalSystem(
            LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        rCurrentCondition.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(Element& rCurrentElement,
                                  LocalSystemVectorType& RHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        rCurrentElement.CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(Condition& rCurrentCondition,
                                  LocalSystemVectorType& RHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        rCurrentCondition.CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        rCurrentCondition.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateLHSContribution(Element& rCurrentElement,
                                  LocalSystemMatrixType& LHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        rCurrentElement.CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateLHSContribution(Condition& rCurrentCondition,
                                  LocalSystemMatrixType& LHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        rCurrentCondition.CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        rCurrentCondition.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void Update(ModelPart& rModelPart,
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

            typename DofsArrayType::iterator DofsBegin =
                rDofSet.begin() + DofSetPartition[k];
            typename DofsArrayType::iterator DofsEnd =
                rDofSet.begin() + DofSetPartition[k + 1];

            // Update Displacement and Pressure (DOFs)
            for (typename DofsArrayType::iterator itDof = DofsBegin;
                 itDof != DofsEnd; ++itDof) {
                if (itDof->IsFree())
                    itDof->GetSolutionStepValue() +=
                        TSparseSpace::GetValue(Dx, itDof->EquationId());
            }
        }

        this->UpdateVariablesDerivatives(rModelPart);

        KRATOS_CATCH("")
    }

protected:
    void CheckBufferSize(const ModelPart& rModelPart) const
    {
        KRATOS_ERROR_IF(rModelPart.GetBufferSize() < 2)
            << "insufficient buffer size. Buffer size should be "
               "greater or equal to 2. Current size is "
            << rModelPart.GetBufferSize() << std::endl;
    }

    virtual inline void SetTimeFactors(ModelPart& rModelPart) = 0;
    virtual inline void UpdateVariablesDerivatives(ModelPart& rModelPart) = 0;

    double GetDeltaTime() const
    {
        return mDeltaTime;
    }

    void SetDeltaTime(double DeltaTime)
    {
        mDeltaTime = DeltaTime;
    }

private:
    double mDeltaTime = 0.0;
};

} // namespace Kratos
