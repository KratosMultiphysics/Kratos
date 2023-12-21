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

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class GeoMechanicsTimeIntegrationScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    using BaseType = Scheme<TSparseSpace, TDenseSpace>;
    using DofsArrayType = typename BaseType::DofsArrayType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    GeoMechanicsTimeIntegrationScheme(const Variable<double>& rVariable,
                                      const Variable<double>& rDeltaTimeVariable,
                                      const Variable<double>& rDeltaTimeVariableCoefficient)
        : mVariable(rVariable),
          mDeltaTimeVariable(rDeltaTimeVariable),
          mDeltaTimeVariableCoefficient(rDeltaTimeVariableCoefficient)
    {
    }

    int Check(const ModelPart& rModelPart) const final
    {
        KRATOS_TRY

        Scheme<TSparseSpace, TDenseSpace>::Check(rModelPart);
        CheckAllocatedVariables(rModelPart);
        CheckBufferSize(rModelPart);

        return 0;

        KRATOS_CATCH("")
    }

    void GetDofList(const Element& rElement,
                    Element::DofsVectorType& rDofList,
                    const ProcessInfo& rCurrentProcessInfo) override
    {
        GetDofListImpl(rElement, rDofList, rCurrentProcessInfo);
    }

    void GetDofList(const Condition& rCondition,
                    Condition::DofsVectorType& rDofList,
                    const ProcessInfo& rCurrentProcessInfo) override
    {
        GetDofListImpl(rCondition, rDofList, rCurrentProcessInfo);
    }

    template <typename T>
    void GetDofListImpl(const T& rElementOrCondition,
                        typename T::DofsVectorType& rDofList,
                        const ProcessInfo& rCurrentProcessInfo)
    {
        if (IsActive(rElementOrCondition))
            rElementOrCondition.GetDofList(rDofList, rCurrentProcessInfo);
    }

    void EquationId(const Element& rElement,
                    Element::EquationIdVectorType& rEquationId,
                    const ProcessInfo& rCurrentProcessInfo) override
    {
        EquationIdImpl(rElement, rEquationId, rCurrentProcessInfo);
    }

    void EquationId(const Condition& rCondition,
                    Condition::EquationIdVectorType& rEquationId,
                    const ProcessInfo& rCurrentProcessInfo) override
    {
        EquationIdImpl(rCondition, rEquationId, rCurrentProcessInfo);
    }

    template <typename T>
    void EquationIdImpl(const T& rElementOrCondition,
                        typename T::EquationIdVectorType& rEquationId,
                        const ProcessInfo& rCurrentProcessInfo)
    {
        if (IsActive(rElementOrCondition))
            rElementOrCondition.EquationIdVector(rEquationId, rCurrentProcessInfo);
    }

    void Initialize(ModelPart& rModelPart) override
    {
        Scheme<TSparseSpace, TDenseSpace>::Initialize(rModelPart);

        KRATOS_TRY
        SetTimeFactors(rModelPart);
        KRATOS_CATCH("")
    }

    void Predict(ModelPart& rModelPart,
                 DofsArrayType&,
                 TSystemMatrixType&,
                 TSystemVectorType&,
                 TSystemVectorType&) override
    {
        this->UpdateVariablesDerivatives(rModelPart);
    }

    void InitializeSolutionStep(ModelPart& rModelPart,
                                TSystemMatrixType&,
                                TSystemVectorType&,
                                TSystemVectorType&) override
    {
        KRATOS_TRY

        SetTimeFactors(rModelPart);

        BlockForEachActiveElement(rModelPart, &Element::InitializeSolutionStep);
        BlockForEachActiveCondition(rModelPart, &Condition::InitializeSolutionStep);

        KRATOS_CATCH("")
    }

    void InitializeNonLinIteration(ModelPart& rModelPart,
                                   TSystemMatrixType&,
                                   TSystemVectorType&,
                                   TSystemVectorType&) override
    {
        KRATOS_TRY

        BlockForEachActiveElement(rModelPart, &Element::InitializeNonLinearIteration);
        BlockForEachActiveCondition(rModelPart, &Condition::InitializeNonLinearIteration);

        KRATOS_CATCH("")
    }

    void FinalizeNonLinIteration(ModelPart& rModelPart,
                                 TSystemMatrixType&,
                                 TSystemVectorType&,
                                 TSystemVectorType&) override
    {
        KRATOS_TRY

        BlockForEachActiveElement(rModelPart, &Element::FinalizeNonLinearIteration);
        BlockForEachActiveCondition(rModelPart, &Condition::FinalizeNonLinearIteration);

        KRATOS_CATCH("")
    }

    void FinalizeSolutionStepActiveEntities(ModelPart& rModelPart,
                                            TSystemMatrixType&,
                                            TSystemVectorType&,
                                            TSystemVectorType&)
    {
        KRATOS_TRY

        BlockForEachActiveElement(rModelPart, &Element::FinalizeSolutionStep);
        BlockForEachActiveCondition(rModelPart, &Condition::FinalizeSolutionStep);

        KRATOS_CATCH("")
    }

    template <typename MemFuncPtr>
    void BlockForEachActiveElement(ModelPart& rModelPart, MemFuncPtr pMemberFunction)
    {
        const auto& r_current_process_info = rModelPart.GetProcessInfo();
        block_for_each(rModelPart.Elements(),
                       [&r_current_process_info, pMemberFunction](auto& rElement)
        {
            if (IsActive(rElement))
            {
                (rElement.*pMemberFunction)(r_current_process_info);
            }
        });
    }

    template <typename MemFuncPtr>
    void BlockForEachActiveCondition(ModelPart& rModelPart, MemFuncPtr pMemberFunction)
    {
        const auto& r_current_process_info = rModelPart.GetProcessInfo();
        block_for_each(rModelPart.Conditions(),
                       [&r_current_process_info, pMemberFunction](Condition& rCondition)
        {
            if (IsActive(rCondition))
                (rCondition.*pMemberFunction)(r_current_process_info);
        });
    }

    template <class T>
    static bool IsActive(const T& rComponent)
    {
        return !(rComponent.IsDefined(ACTIVE)) || rComponent.Is(ACTIVE);
    }

    void FinalizeSolutionStep(ModelPart& rModelPart,
                              TSystemMatrixType& A,
                              TSystemVectorType& Dx,
                              TSystemVectorType& b) override
    {
        FinalizeSolutionStepActiveEntities(rModelPart, A, Dx, b);
    }

    void CalculateSystemContributions(Element& rCurrentElement,
                                      LocalSystemMatrixType& LHS_Contribution,
                                      LocalSystemVectorType& RHS_Contribution,
                                      Element::EquationIdVectorType& EquationId,
                                      const ProcessInfo& CurrentProcessInfo) override
    {
        CalculateSystemContributionsImpl(rCurrentElement, LHS_Contribution, RHS_Contribution,
                                         EquationId, CurrentProcessInfo);
    }

    void CalculateSystemContributions(Condition& rCurrentCondition,
                                      LocalSystemMatrixType& LHS_Contribution,
                                      LocalSystemVectorType& RHS_Contribution,
                                      Element::EquationIdVectorType& EquationId,
                                      const ProcessInfo& CurrentProcessInfo) override
    {
        CalculateSystemContributionsImpl(rCurrentCondition, LHS_Contribution, RHS_Contribution,
                                         EquationId, CurrentProcessInfo);
    }

    template <typename T>
    void CalculateSystemContributionsImpl(T& rCurrentComponent,
                                          LocalSystemMatrixType& LHS_Contribution,
                                          LocalSystemVectorType& RHS_Contribution,
                                          typename T::EquationIdVectorType& EquationId,
                                          const ProcessInfo& CurrentProcessInfo)

    {
        KRATOS_TRY

        rCurrentComponent.CalculateLocalSystem(
            LHS_Contribution, RHS_Contribution, CurrentProcessInfo);
        rCurrentComponent.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(Element& rCurrentElement,
                                  LocalSystemVectorType& RHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo& CurrentProcessInfo) override
    {
        CalculateRHSContributionImpl(rCurrentElement, RHS_Contribution,
                                     EquationId, CurrentProcessInfo);
    }

    void CalculateRHSContribution(Condition& rCurrentCondition,
                                  LocalSystemVectorType& RHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo& CurrentProcessInfo) override
    {
        CalculateRHSContributionImpl(rCurrentCondition, RHS_Contribution,
                                     EquationId, CurrentProcessInfo);
    }

    template <typename T>
    void CalculateRHSContributionImpl(T& rCurrentComponent,
                                      LocalSystemVectorType& RHS_Contribution,
                                      typename T::EquationIdVectorType& EquationId,
                                      const ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        rCurrentComponent.CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);
        rCurrentComponent.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateLHSContribution(Element& rCurrentElement,
                                  LocalSystemMatrixType& LHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo& CurrentProcessInfo) override
    {
        CalculateLHSContributionImpl(rCurrentElement, LHS_Contribution,
                                     EquationId, CurrentProcessInfo);
    }

    void CalculateLHSContribution(Condition& rCurrentCondition,
                                  LocalSystemMatrixType& LHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo& CurrentProcessInfo) override
    {
        CalculateLHSContributionImpl(rCurrentCondition, LHS_Contribution,
                                     EquationId, CurrentProcessInfo);
    }

    template <typename T>
    void CalculateLHSContributionImpl(T& rCurrentComponent,
                                      LocalSystemMatrixType& LHS_Contribution,
                                      typename T::EquationIdVectorType& EquationId,
                                      const ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        rCurrentComponent.CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);
        rCurrentComponent.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void Update(ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                TSystemMatrixType&,
                TSystemVectorType& Dx,
                TSystemVectorType&) override
    {
        KRATOS_TRY

        block_for_each(rDofSet, [&Dx](auto& dof)
        {
            if (dof.IsFree())
            {
                dof.GetSolutionStepValue() +=
                    TSparseSpace::GetValue(Dx, dof.EquationId());
            }
        });

        this->UpdateVariablesDerivatives(rModelPart);

        KRATOS_CATCH("")
    }

protected:
    void CheckBufferSize(const ModelPart& rModelPart) const
    {
        constexpr auto minimum_buffer_size = ModelPart::IndexType{2};
        KRATOS_ERROR_IF(rModelPart.GetBufferSize() < minimum_buffer_size)
            << "insufficient buffer size. Buffer size should be "
               "greater than or equal to "
            << minimum_buffer_size << ". Current size is "
            << rModelPart.GetBufferSize() << std::endl;
    }

    template <class T>
    void CheckSolutionStepsData(const Node& rNode, const Variable<T>& rVariable) const
    {
        KRATOS_ERROR_IF_NOT(rNode.SolutionStepsDataHas(rVariable))
            << rVariable.Name() << " variable is not allocated for node "
            << rNode.Id() << std::endl;
    }

    template <class T>
    void CheckDof(const Node& rNode, const Variable<T>& rVariable) const
    {
        KRATOS_ERROR_IF_NOT(rNode.HasDofFor(rVariable))
            << "missing " << rVariable.Name() << " dof on node " << rNode.Id()
            << std::endl;
    }

    virtual void CheckAllocatedVariables(const ModelPart& rModelPart) const
    {
        for (const auto& r_node : rModelPart.Nodes())
        {
            this->CheckSolutionStepsData(r_node, mVariable);
            this->CheckSolutionStepsData(r_node, mDeltaTimeVariable);
            this->CheckDof(r_node, mVariable);
        }
    }

    virtual inline void SetTimeFactors(ModelPart& rModelPart) = 0;
    virtual inline void UpdateVariablesDerivatives(ModelPart& rModelPart)
    {
        KRATOS_TRY

        block_for_each(rModelPart.Nodes(),
                       [this](Node& rNode) {
            UpdateScalarTimeDerivative(rNode, mVariable, mDeltaTimeVariable);
        });

        KRATOS_CATCH("")
    }


    virtual void UpdateScalarTimeDerivative(Node& rNode,
                                            const Variable<double>& rVariable,
                                            const Variable<double>& rDtVariable) const = 0;

    [[nodiscard]] double GetDeltaTime() const { return mDeltaTime; }

    void SetDeltaTime(double DeltaTime) { mDeltaTime = DeltaTime; }

    Variable<double> mVariable;
    Variable<double> mDeltaTimeVariable;
    Variable<double> mDeltaTimeVariableCoefficient;

private:
    double mDeltaTime = 1.0;
};

} // namespace Kratos
