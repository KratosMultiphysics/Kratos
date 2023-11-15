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
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    int Check(const ModelPart& rModelPart) const override
    {
        KRATOS_TRY

        Scheme<TSparseSpace, TDenseSpace>::Check(rModelPart);
        CheckAllocatedVariables(rModelPart);
        this->CheckBufferSize(rModelPart);

        return 0;

        KRATOS_CATCH("")
    }

    void GetDofList(const Element& rElement,
                    Element::DofsVectorType& rDofList,
                    const ProcessInfo& rCurrentProcessInfo) override
    {
        if (IsActive(rElement))
            rElement.GetDofList(rDofList, rCurrentProcessInfo);
    }

    void GetDofList(const Condition& rCondition,
                    Condition::DofsVectorType& rDofList,
                    const ProcessInfo& rCurrentProcessInfo) override
    {
        if (IsActive(rCondition))
            rCondition.GetDofList(rDofList, rCurrentProcessInfo);
    }

    void EquationId(const Element& rElement,
                    Element::EquationIdVectorType& rEquationId,
                    const ProcessInfo& rCurrentProcessInfo) override
    {
        if (IsActive(rElement))
            rElement.EquationIdVector(rEquationId, rCurrentProcessInfo);
    }

    void EquationId(const Condition& rCondition,
                    Condition::EquationIdVectorType& rEquationId,
                    const ProcessInfo& rCurrentProcessInfo) override
    {
        if (IsActive(rCondition))
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

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        block_for_each(rModelPart.Elements(), [&r_current_process_info, this](Element& rElement) {
            if (IsActive(rElement))
                rElement.InitializeSolutionStep(r_current_process_info);
        });

        block_for_each(rModelPart.Conditions(), [&r_current_process_info,
                                                 this](Condition& rCondition) {
            if (IsActive(rCondition))
                rCondition.InitializeSolutionStep(r_current_process_info);
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

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        block_for_each(rModelPart.Elements(), [&r_current_process_info, this](Element& rElement) {
            if (IsActive(rElement))
                rElement.InitializeNonLinearIteration(r_current_process_info);
        });

        block_for_each(rModelPart.Conditions(), [&r_current_process_info,
                                                 this](Condition& rCondition) {
            if (IsActive(rCondition))
                rCondition.InitializeNonLinearIteration(r_current_process_info);
        });

        KRATOS_CATCH("")
    }

    void FinalizeNonLinIteration(ModelPart& rModelPart,
                                 TSystemMatrixType& A,
                                 TSystemVectorType& Dx,
                                 TSystemVectorType& b) override
    {
        KRATOS_TRY

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        block_for_each(rModelPart.Elements(), [&r_current_process_info, this](Element& rElement) {
            if (IsActive(rElement))
                rElement.FinalizeNonLinearIteration(r_current_process_info);
        });

        block_for_each(rModelPart.Conditions(), [&r_current_process_info,
                                                 this](Condition& rCondition) {
            if (IsActive(rCondition))
                rCondition.FinalizeNonLinearIteration(r_current_process_info);
        });

        KRATOS_CATCH("")
    }

    void FinalizeSolutionStepActiveEntities(ModelPart& rModelPart,
                                            TSystemMatrixType& A,
                                            TSystemVectorType& Dx,
                                            TSystemVectorType& b)
    {
        KRATOS_TRY

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        block_for_each(rModelPart.Elements(), [&r_current_process_info, this](Element& rElement) {
            if (IsActive(rElement))
                rElement.FinalizeSolutionStep(r_current_process_info);
        });

        block_for_each(rModelPart.Conditions(),
                       [&r_current_process_info, this](Condition& rCondition) {
                           if (IsActive(rCondition))
                               rCondition.FinalizeSolutionStep(r_current_process_info);
                       });

        KRATOS_CATCH("")
    }

    template <class T>
    bool IsActive(const T& rComponent) const
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
            for (typename DofsArrayType::iterator itDof = dofs_begin;
                 itDof != dofs_end; ++itDof) {
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
        constexpr int minimum_buffer_size = 2;
        KRATOS_ERROR_IF(rModelPart.GetBufferSize() < minimum_buffer_size)
            << "insufficient buffer size. Buffer size should be "
               "greater or equal to "
            << minimum_buffer_size << ". Current size is "
            << rModelPart.GetBufferSize() << std::endl;
    }

    template <class T>
    void CheckSolutionStepsData(const Node& r_node, const Variable<T>& variable) const
    {
        KRATOS_ERROR_IF(!r_node.SolutionStepsDataHas(variable))
            << variable.Name() << " variable is not allocated for node "
            << r_node.Id() << std::endl;
    }

    template <class T>
    void CheckDof(const Node& r_node, const Variable<T>& variable) const
    {
        KRATOS_ERROR_IF(!r_node.HasDofFor(variable))
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

    virtual void CheckAllocatedVariables(const ModelPart& rModelPart) const
    {
        // intentionally empty
    }

    virtual void UpdateScalarTimeDerivative(Node& rNode,
                                            const Variable<double>& variable,
                                            const Variable<double>& dt_variable) const
    {
        // intentionally empty
    }
;
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
