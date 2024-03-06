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

#include "geo_mechanics_application_variables.h"
#include "solving_strategies/schemes/scheme.h"

namespace Kratos
{

struct FirstOrderScalarVariable {
    Variable<double> instance;
    Variable<double> first_time_derivative;
    Variable<double> delta_time_coefficient;

    FirstOrderScalarVariable(const Variable<double>& instance,
                             const Variable<double>& first_time_derivative,
                             const Variable<double>& delta_time_coefficient)
        : instance(instance),
          first_time_derivative(first_time_derivative),
          delta_time_coefficient(delta_time_coefficient)
    {
    }
};

struct SecondOrderVectorVariable {
    Variable<array_1d<double, 3>> instance;
    Variable<array_1d<double, 3>> first_time_derivative;
    Variable<array_1d<double, 3>> second_time_derivative;

    explicit SecondOrderVectorVariable(const Variable<array_1d<double, 3>>& instance)
        : instance(instance),
          first_time_derivative(instance.GetTimeDerivative()),
          second_time_derivative(first_time_derivative.GetTimeDerivative())
    {
    }
};

template <class TSparseSpace, class TDenseSpace>
class GeoMechanicsTimeIntegrationScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    using BaseType              = Scheme<TSparseSpace, TDenseSpace>;
    using DofsArrayType         = typename BaseType::DofsArrayType;
    using TSystemVectorType     = typename BaseType::TSystemVectorType;
    using TSystemMatrixType     = typename BaseType::TSystemMatrixType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    GeoMechanicsTimeIntegrationScheme(const std::vector<FirstOrderScalarVariable>& rFirstOrderScalarVariables,
                                      const std::vector<SecondOrderVectorVariable>& rSecondOrderVectorVariables)
        : mFirstOrderScalarVariables(rFirstOrderScalarVariables),
          mSecondOrderVectorVariables(rSecondOrderVectorVariables)
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

    void GetDofList(const Element& rElement, Element::DofsVectorType& rDofList, const ProcessInfo& rCurrentProcessInfo) override
    {
        GetDofListImpl(rElement, rDofList, rCurrentProcessInfo);
    }

    void GetDofList(const Condition& rCondition, Condition::DofsVectorType& rDofList, const ProcessInfo& rCurrentProcessInfo) override
    {
        GetDofListImpl(rCondition, rDofList, rCurrentProcessInfo);
    }

    template <typename T>
    void GetDofListImpl(const T& rElementOrCondition, typename T::DofsVectorType& rDofList, const ProcessInfo& rCurrentProcessInfo)
    {
        if (IsActive(rElementOrCondition))
            rElementOrCondition.GetDofList(rDofList, rCurrentProcessInfo);
    }

    void EquationId(const Element&                 rElement,
                    Element::EquationIdVectorType& rEquationId,
                    const ProcessInfo&             rCurrentProcessInfo) override
    {
        EquationIdImpl(rElement, rEquationId, rCurrentProcessInfo);
    }

    void EquationId(const Condition&                 rCondition,
                    Condition::EquationIdVectorType& rEquationId,
                    const ProcessInfo&               rCurrentProcessInfo) override
    {
        EquationIdImpl(rCondition, rEquationId, rCurrentProcessInfo);
    }

    template <typename T>
    void EquationIdImpl(const T&                          rElementOrCondition,
                        typename T::EquationIdVectorType& rEquationId,
                        const ProcessInfo&                rCurrentProcessInfo)
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

    void Predict(ModelPart& rModelPart, DofsArrayType&, TSystemMatrixType&, TSystemVectorType&, TSystemVectorType&) override
    {
        this->UpdateVariablesDerivatives(rModelPart);
    }

    void InitializeSolutionStep(ModelPart& rModelPart, TSystemMatrixType&, TSystemVectorType&, TSystemVectorType&) override
    {
        KRATOS_TRY

        SetTimeFactors(rModelPart);

        BlockForEachActiveElement(rModelPart, &Element::InitializeSolutionStep);
        BlockForEachActiveCondition(rModelPart, &Condition::InitializeSolutionStep);

        KRATOS_CATCH("")
    }

    void InitializeNonLinIteration(ModelPart& rModelPart, TSystemMatrixType&, TSystemVectorType&, TSystemVectorType&) override
    {
        KRATOS_TRY

        BlockForEachActiveElement(rModelPart, &Element::InitializeNonLinearIteration);
        BlockForEachActiveCondition(rModelPart, &Condition::InitializeNonLinearIteration);

        KRATOS_CATCH("")
    }

    void FinalizeNonLinIteration(ModelPart& rModelPart, TSystemMatrixType&, TSystemVectorType&, TSystemVectorType&) override
    {
        KRATOS_TRY

        BlockForEachActiveElement(rModelPart, &Element::FinalizeNonLinearIteration);
        BlockForEachActiveCondition(rModelPart, &Condition::FinalizeNonLinearIteration);

        KRATOS_CATCH("")
    }

    void FinalizeSolutionStepActiveEntities(ModelPart& rModelPart, TSystemMatrixType&, TSystemVectorType&, TSystemVectorType&)
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
        block_for_each(rModelPart.Elements(), [&r_current_process_info, pMemberFunction](auto& rElement) {
            if (IsActive(rElement)) {
                (rElement.*pMemberFunction)(r_current_process_info);
            }
        });
    }

    template <typename MemFuncPtr>
    void BlockForEachActiveCondition(ModelPart& rModelPart, MemFuncPtr pMemberFunction)
    {
        const auto& r_current_process_info = rModelPart.GetProcessInfo();
        block_for_each(rModelPart.Conditions(), [&r_current_process_info, pMemberFunction](Condition& rCondition) {
            if (IsActive(rCondition)) (rCondition.*pMemberFunction)(r_current_process_info);
        });
    }

    template <class T>
    static bool IsActive(const T& rComponent)
    {
        return !(rComponent.IsDefined(ACTIVE)) || rComponent.Is(ACTIVE);
    }

    void FinalizeSolutionStep(ModelPart& rModelPart, TSystemMatrixType& A, TSystemVectorType& Dx, TSystemVectorType& b) override
    {
        FinalizeSolutionStepActiveEntities(rModelPart, A, Dx, b);
    }

    void CalculateSystemContributions(Element&                       rCurrentElement,
                                      LocalSystemMatrixType&         LHS_Contribution,
                                      LocalSystemVectorType&         RHS_Contribution,
                                      Element::EquationIdVectorType& EquationId,
                                      const ProcessInfo&             CurrentProcessInfo) override
    {
        CalculateSystemContributionsImpl(rCurrentElement, LHS_Contribution, RHS_Contribution,
                                         EquationId, CurrentProcessInfo);
    }

    void CalculateSystemContributions(Condition&                     rCurrentCondition,
                                      LocalSystemMatrixType&         LHS_Contribution,
                                      LocalSystemVectorType&         RHS_Contribution,
                                      Element::EquationIdVectorType& EquationId,
                                      const ProcessInfo&             CurrentProcessInfo) override
    {
        CalculateSystemContributionsImpl(rCurrentCondition, LHS_Contribution, RHS_Contribution,
                                         EquationId, CurrentProcessInfo);
    }

    template <typename T>
    void CalculateSystemContributionsImpl(T&                                rCurrentComponent,
                                          LocalSystemMatrixType&            LHS_Contribution,
                                          LocalSystemVectorType&            RHS_Contribution,
                                          typename T::EquationIdVectorType& EquationId,
                                          const ProcessInfo&                CurrentProcessInfo)

    {
        KRATOS_TRY

        rCurrentComponent.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);
        rCurrentComponent.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(Element&                       rCurrentElement,
                                  LocalSystemVectorType&         RHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo&             CurrentProcessInfo) override
    {
        CalculateRHSContributionImpl(rCurrentElement, RHS_Contribution, EquationId, CurrentProcessInfo);
    }

    void CalculateRHSContribution(Condition&                     rCurrentCondition,
                                  LocalSystemVectorType&         RHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo&             CurrentProcessInfo) override
    {
        CalculateRHSContributionImpl(rCurrentCondition, RHS_Contribution, EquationId, CurrentProcessInfo);
    }

    template <typename T>
    void CalculateRHSContributionImpl(T&                                rCurrentComponent,
                                      LocalSystemVectorType&            RHS_Contribution,
                                      typename T::EquationIdVectorType& EquationId,
                                      const ProcessInfo&                CurrentProcessInfo)
    {
        KRATOS_TRY

        rCurrentComponent.CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);
        rCurrentComponent.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateLHSContribution(Element&                       rCurrentElement,
                                  LocalSystemMatrixType&         LHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo&             CurrentProcessInfo) override
    {
        CalculateLHSContributionImpl(rCurrentElement, LHS_Contribution, EquationId, CurrentProcessInfo);
    }

    void CalculateLHSContribution(Condition&                     rCurrentCondition,
                                  LocalSystemMatrixType&         LHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo&             CurrentProcessInfo) override
    {
        CalculateLHSContributionImpl(rCurrentCondition, LHS_Contribution, EquationId, CurrentProcessInfo);
    }

    template <typename T>
    void CalculateLHSContributionImpl(T&                                rCurrentComponent,
                                      LocalSystemMatrixType&            LHS_Contribution,
                                      typename T::EquationIdVectorType& EquationId,
                                      const ProcessInfo&                CurrentProcessInfo)
    {
        KRATOS_TRY

        rCurrentComponent.CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);
        rCurrentComponent.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void Update(ModelPart& rModelPart, DofsArrayType& rDofSet, TSystemMatrixType&, TSystemVectorType& Dx, TSystemVectorType&) override
    {
        KRATOS_TRY

        block_for_each(rDofSet, [&Dx](auto& dof) {
            if (dof.IsFree()) {
                dof.GetSolutionStepValue() += TSparseSpace::GetValue(Dx, dof.EquationId());
            }
        });

        UpdateVariablesDerivatives(rModelPart);

        KRATOS_CATCH("")
    }

protected:
    const Variable<double>& GetComponentFromVectorVariable(const Variable<array_1d<double, 3>>& rSource,
                                                           const std::string& rComponent) const
    {
        return KratosComponents<Variable<double>>::Get(rSource.Name() + "_" + rComponent);
    }

    virtual inline void SetTimeFactors(ModelPart& rModelPart)
    {
        mDeltaTime = rModelPart.GetProcessInfo()[DELTA_TIME];
    }

    virtual inline void UpdateVariablesDerivatives(ModelPart& rModelPart) = 0;

    [[nodiscard]] double GetDeltaTime() const { return mDeltaTime; }

    [[nodiscard]] const std::vector<SecondOrderVectorVariable>& GetSecondOrderVectorVariables() const
    {
        return mSecondOrderVectorVariables;
    }

    [[nodiscard]] const std::vector<FirstOrderScalarVariable>& GetFirstOrderScalarVariables() const
    {
        return mFirstOrderScalarVariables;
    }

private:
    void CheckAllocatedVariables(const ModelPart& rModelPart) const
    {
        for (const auto& r_node : rModelPart.Nodes()) {
            for (const auto& r_first_order_scalar_variable : mFirstOrderScalarVariables) {
                this->CheckSolutionStepsData(r_node, r_first_order_scalar_variable.instance);
                this->CheckSolutionStepsData(r_node, r_first_order_scalar_variable.first_time_derivative);
                this->CheckDof(r_node, r_first_order_scalar_variable.instance);
            }
        }

        for (const auto& r_node : rModelPart.Nodes()) {
            for (const auto& r_second_order_vector_variable : this->mSecondOrderVectorVariables) {
                if (!rModelPart.HasNodalSolutionStepVariable(r_second_order_vector_variable.instance))
                    continue;

                this->CheckSolutionStepsData(r_node, r_second_order_vector_variable.instance);
                this->CheckSolutionStepsData(r_node, r_second_order_vector_variable.first_time_derivative);
                this->CheckSolutionStepsData(r_node, r_second_order_vector_variable.second_time_derivative);

                // We don't check for "Z", since it is optional (in case of a 2D problem)
                std::vector<std::string> components{"X", "Y"};
                for (const auto& component : components) {
                    const auto& variable_component =
                        GetComponentFromVectorVariable(r_second_order_vector_variable.instance, component);
                    this->CheckDof(r_node, variable_component);
                }
            }
        }
    }

    void CheckBufferSize(const ModelPart& rModelPart) const
    {
        constexpr auto minimum_buffer_size = ModelPart::IndexType{2};
        KRATOS_ERROR_IF(rModelPart.GetBufferSize() < minimum_buffer_size)
            << "insufficient buffer size. Buffer size should be "
               "greater than or equal to "
            << minimum_buffer_size << ". Current size is " << rModelPart.GetBufferSize() << std::endl;
    }

    template <class T>
    void CheckSolutionStepsData(const Node& rNode, const Variable<T>& rVariable) const
    {
        KRATOS_ERROR_IF_NOT(rNode.SolutionStepsDataHas(rVariable))
            << rVariable.Name() << " variable is not allocated for node " << rNode.Id() << std::endl;
    }

    template <class T>
    void CheckDof(const Node& rNode, const Variable<T>& rVariable) const
    {
        KRATOS_ERROR_IF_NOT(rNode.HasDofFor(rVariable))
            << "missing " << rVariable.Name() << " dof on node " << rNode.Id() << std::endl;
    }

    double                                 mDeltaTime = 1.0;
    std::vector<FirstOrderScalarVariable>  mFirstOrderScalarVariables;
    std::vector<SecondOrderVectorVariable> mSecondOrderVectorVariables;
};

} // namespace Kratos
