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
class GeoMechanicsTimeIntegrationScheme : public Scheme<TSparseSpace, TDenseSpace> {
public:
    using BaseType = Scheme<TSparseSpace, TDenseSpace>;
    using DofsArrayType = typename BaseType::DofsArrayType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

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
                       [&r_current_process_info, pMemberFunction](auto& rElement) {
                           if (IsActive(rElement)) {
                               (rElement.*pMemberFunction)(r_current_process_info);
                           }
                       });
    }

    template <typename MemFuncPtr>
    void BlockForEachActiveCondition(ModelPart& rModelPart, MemFuncPtr pMemberFunction)
    {
        const auto& r_current_process_info = rModelPart.GetProcessInfo();
        block_for_each(rModelPart.Conditions(),
                       [&r_current_process_info, pMemberFunction](Condition& rCondition) {
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

        int num_threads = ParallelUtilities::GetNumThreads();
        OpenMPUtils::PartitionVector dof_set_partition;
        OpenMPUtils::DivideInPartitions(static_cast<int>(rDofSet.size()),
                                        num_threads, dof_set_partition);

#pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            typename DofsArrayType::iterator dofs_begin =
                rDofSet.begin() + dof_set_partition[k];
            typename DofsArrayType::iterator dofs_end =
                rDofSet.begin() + dof_set_partition[k + 1];

            // Update Displacement and Pressure (DOFs)
            for (typename DofsArrayType::iterator it_dof = dofs_begin;
                 it_dof != dofs_end; ++it_dof) {
                if (it_dof->IsFree())
                    it_dof->GetSolutionStepValue() +=
                        TSparseSpace::GetValue(Dx, it_dof->EquationId());
            }
        }

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

    void CheckSolutionStepsData(const Node& r_node, const Variable<double>& variable) const
    {
        KRATOS_ERROR_IF_NOT(r_node.SolutionStepsDataHas(variable))
            << variable.Name() << " variable is not allocated for node "
            << r_node.Id() << std::endl;
    }

    void CheckDof(const Node& r_node, const Variable<double>& variable) const
    {
        KRATOS_ERROR_IF_NOT(r_node.HasDofFor(variable))
            << "missing " << variable.Name() << " dof on node " << r_node.Id()
            << std::endl;
    }

    virtual inline void SetTimeFactors(ModelPart& rModelPart)
    {
        // intentionally empty
    }
    virtual inline void UpdateVariablesDerivatives(ModelPart& rModelPart)
    {
        // intentionally empty
    }

    [[nodiscard]] double GetDeltaTime() const
    {
        return mDeltaTime;
    }

    void SetDeltaTime(double DeltaTime)
    {
        mDeltaTime = DeltaTime;
    }

private:
    double mDeltaTime = 1.0;
};

} // namespace Kratos
