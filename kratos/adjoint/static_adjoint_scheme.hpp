//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// --- Kratos Core Includes ---
#include "adjoint/adjoint_scheme.hpp"


namespace Kratos {


template <class TSparse, class TDense>
class KRATOS_API(KRATOS_CORE) StaticAdjointScheme : public AdjointScheme<TSparse,TDense> {
public:
    using Base = AdjointScheme<TSparse,TDense>;

    using Base::Base;

    KRATOS_CLASS_POINTER_DEFINITION(StaticAdjointScheme);

private:
    /// @details Computes @f$K@f$ and @f$-\frac{\partial j}{\partial u}@f$
    template <class T>
    requires (std::is_same_v<T,Element> || std::is_same_v<T,Condition>)
    void CalculateSystemContributionsImpl(
        T& rElement,
        typename Base::LocalSystemMatrixType& rLhsContribution,
        typename Base::LocalSystemVectorType& rRhsContribution,
        Element::EquationIdVectorType& rEquationIds,
        const ProcessInfo& rProcessInfo) {
            KRATOS_TRY
                // Collect DoFs and convert them to IAdjoint::DynamicVariable.
                Element::DofsVectorType dofs;
                rElement.GetDofList(
                    dofs,
                    rProcessInfo);

                const Geometry<Node>& r_geometry = rElement.GetGeometry();
                std::vector<IAdjoint::DynamicVariable> state_variables;
                state_variables.reserve(dofs.size());
                std::transform(
                    dofs.begin(),
                    dofs.end(),
                    std::back_inserter(state_variables),
                    [&r_geometry] (const Dof<double>* p_dof) -> IAdjoint::DynamicVariable {
                        // Find the node's local index within the geometry.
                        const auto it_node = std::find_if(
                            r_geometry.begin(),
                            r_geometry.end(),
                            [target_id = p_dof->Id()] (const Node& r_node) {
                                return r_node == target_id;
                            });
                        const auto i_node = std::distance(
                            r_geometry.begin(),
                            it_node);

                        // Construct a variable representing the DoF.
                        return IAdjoint::DynamicVariable(
                            p_dof->GetVariable(),
                            i_node);
                    });

                // Convert DoFs to equation Ids.
                rEquationIds.resize(dofs.size());
                std::transform(
                    dofs.begin(),
                    dofs.end(),
                    rEquationIds.begin(),
                    [] (const Dof<double>* p_dof) -> Dof<double>::EquationIdType {
                        return p_dof->EquationId();
                    });

                // Compute K.
                rElement.CalculateLeftHandSide(
                    rLhsContribution,
                    rProcessInfo);

                // Compute -\frac{\partial j}{\partial u}.
                this->GetResponseFunction().ComputeDerivative(
                    rRhsContribution,
                    rElement,
                    state_variables,
                    rProcessInfo,
                    0);
                rRhsContribution *= -1.0;
            KRATOS_CATCH("while processing element " + std::to_string(rElement.Id()))
    }

    /// @details Computes @f$K@f$.
    template <class T>
    requires (std::is_same_v<T,Element> || std::is_same_v<T,Condition>)
    void CalculateLHSContributionsImpl(
        T& rElement,
        typename Base::LocalSystemMatrixType& rLhsContribution,
        Element::EquationIdVectorType& rEquationIds,
        const ProcessInfo& rProcessInfo) {
            KRATOS_TRY
                rElement.CalculateLeftHandSide(
                    rLhsContribution,
                    rProcessInfo);
                rElement.EquationIdVector(
                    rEquationIds,
                    rProcessInfo);
            KRATOS_CATCH("while processing element " + std::to_string(rElement.Id()))
    }

    /// @details Computes @f$-\frac{\partial j}{\partial u}@f$
    template <class T>
    requires (std::is_same_v<T,Element> || std::is_same_v<T,Condition>)
    void CalculateRHSContributionsImpl(
        T& rElement,
        typename Base::LocalSystemVectorType& rRhsContribution,
        Element::EquationIdVectorType& rEquationIds,
        const ProcessInfo& rProcessInfo) {
            KRATOS_TRY
                // Collect state variables from the response function.
                std::vector<IAdjoint::DynamicVariable> state_variables;
                this->GetResponseFunction().GetStateVariables(
                    state_variables,
                    rElement,
                    rProcessInfo);

                // Compute the response function's derivatives
                // with respect to the collected state variables.
                this->GetResponseFunction().ComputeDerivative(
                    rRhsContribution,
                    rElement,
                    state_variables,
                    rProcessInfo,
                    0);

                // Get equation IDs for the collected state variables.
                this->template GetEquationIds<T>(
                    rEquationIds,
                    state_variables,
                    rElement);

                rRhsContribution *= -1.0;
            KRATOS_CATCH("while processing element " + std::to_string(rElement.Id()))
    }

public:
    /// @copydoc Base::CalculateSystemContributions
    void CalculateSystemContributions(
        Element& rElement,
        typename Base::LocalSystemMatrixType& rLhsContribution,
        typename Base::LocalSystemVectorType& rRhsContribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rProcessInfo) override {
            return this->template CalculateSystemContributionsImpl<Element>(
                rElement,
                rLhsContribution,
                rRhsContribution,
                rEquationId,
                rProcessInfo);
    }

    /// @copydoc Base::CalculateLHSContribution
    void CalculateLHSContribution(
        Element& rElement,
        typename Base::LocalSystemMatrixType& rLhsContribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rProcessInfo) override {
            return this->template CalculateLHSContributionsImpl<Element>(
                rElement,
                rLhsContribution,
                rEquationId,
                rProcessInfo);
    }

    /// @copydoc Base::CalculateRHSContribution
    void CalculateRHSContribution(
        Element& rElement,
        typename Base::LocalSystemVectorType& rRhsContribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rProcessInfo) override {
            return this->template CalculateRHSContributionsImpl<Element>(
                rElement,
                rRhsContribution,
                rEquationId,
                rProcessInfo);
    }

    /// @copydoc Base::CalculateSystemContributions
    void CalculateSystemContributions(
        Condition& rCondition,
        typename Base::LocalSystemMatrixType& rLhsContribution,
        typename Base::LocalSystemVectorType& rRhsContribution,
        Condition::EquationIdVectorType& rEquationId,
        const ProcessInfo& rProcessInfo) override {
            return this->template CalculateSystemContributionsImpl<Condition>(
                rCondition,
                rLhsContribution,
                rRhsContribution,
                rEquationId,
                rProcessInfo);
    }

    /// @copydoc Base::CalculateLHSContribution
    void CalculateLHSContribution(
        Condition& rCondition,
        typename Base::LocalSystemMatrixType& rLhsContribution,
        Condition::EquationIdVectorType& rEquationId,
        const ProcessInfo& rProcessInfo) override {
            return this->template CalculateLHSContributionsImpl<Condition>(
                rCondition,
                rLhsContribution,
                rEquationId,
                rProcessInfo);
    }

    /// @copydoc Base::CalculateRHSContribution
    void CalculateRHSContribution(
        Condition& rCondition,
        typename Base::LocalSystemVectorType& rRhsContribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rProcessInfo) override {
            return this->template CalculateRHSContributionsImpl<Condition>(
                rCondition,
                rRhsContribution,
                rEquationId,
                rProcessInfo);
    }


}; // class StaticAdjointScheme


} // namespace Kratos
