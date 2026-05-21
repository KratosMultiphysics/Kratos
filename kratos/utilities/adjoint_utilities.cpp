
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

// --- Kratos Includes ---
#include "utilities/adjoint_utilities.hpp"


namespace Kratos {


template <IAdjoint::ResidualTerm Term, class TEntity>
requires (std::is_same_v<TEntity,Element> || std::is_same_v<TEntity,Condition>)
template <Globals::DataLocation TContext>
void AdjointFiniteDifferenceUtility<Term,TEntity>::Perturb(
    TEntity& rEntity,
    const IAdjoint::DynamicVariable& rVariable,
    double PerturbationMagnitude,
    std::size_t iBuffer) const {
        if constexpr (TContext == Globals::DataLocation::NodeHistorical) {
            KRATOS_TRY
                const std::size_t i_node = rVariable.GetDynamicIndex();
                KRATOS_ERROR_IF_NOT(i_node < rEntity.GetGeometry().size());
                Node& r_node = rEntity.GetGeometry()[i_node];
                KRATOS_ERROR_IF_NOT(r_node.SolutionStepsDataHas(rVariable));
                r_node.SolutionStepData().GetValue<double>(rVariable, iBuffer) += PerturbationMagnitude;
            KRATOS_CATCH("")
        } /*if node historical*/ else if constexpr (TContext == Globals::DataLocation::NodeNonHistorical) {
            const std::size_t i_node = rVariable.GetDynamicIndex();
            KRATOS_ERROR_IF_NOT(i_node < rEntity.GetGeometry().size());
            Node& r_node = rEntity.GetGeometry()[i_node];
            switch (rVariable.Key()) {
                case SHAPE_X.Key():
                    r_node.X0() += PerturbationMagnitude;
                    break;
                case SHAPE_Y.Key():
                    r_node.Y0() += PerturbationMagnitude;
                    break;
                case SHAPE_Z.Key():
                    r_node.Z0() += PerturbationMagnitude;
                    break;
                default:
                    r_node.GetData().GetValue<double>(rVariable) += PerturbationMagnitude;
            }
        } /*if node non-historical*/ else if constexpr (TContext == Globals::DataLocation::Element) {
            // Assume the element refers to its properties.
            rEntity.GetProperties().Data().template GetValue<double>(rVariable) += PerturbationMagnitude;
        } /*if element*/ else {
            KRATOS_ERROR << "unsupported data location to perturb";
        }
}


template <IAdjoint::ResidualTerm Term, class TEntity>
requires (std::is_same_v<TEntity,Element> || std::is_same_v<TEntity,Condition>)
void AdjointFiniteDifferenceUtility<Term,TEntity>::ComputeFiniteDifferences(
    Value& rOutput,
    const Value& rReferenceState,
    const Value& rPerturbedState,
    double PerturbationMagnitude,
    std::size_t OutputBegin) const {
        KRATOS_ERROR_IF_NOT(PerturbationMagnitude);
        const double scale = 1 / PerturbationMagnitude;
        if constexpr (Term == IAdjoint::ResidualTerm::Load) {
            KRATOS_ERROR_IF_NOT(rReferenceState.size() == rPerturbedState.size());
            KRATOS_ERROR_IF_NOT(OutputBegin + rPerturbedState.size() <= rOutput.size());
            for (std::size_t i_component=0ul; i_component<rReferenceState.size(); ++i_component) {
                rOutput[i_component + OutputBegin] = scale * (rPerturbedState[i_component] - rReferenceState[i_component]);
            } // for i_component
        } else {
            KRATOS_ERROR_IF_NOT(rReferenceState.size1() == rPerturbedState.size1());
            KRATOS_ERROR_IF_NOT(rReferenceState.size2() == rPerturbedState.size2());
            KRATOS_ERROR_IF_NOT(rOutput.size2() == rPerturbedState.size2());
            KRATOS_ERROR_IF_NOT(OutputBegin + rPerturbedState.size1() <= rOutput.size1());
            for (std::size_t i_row=0ul; i_row<rReferenceState.size1(); ++i_row) {
                row(rOutput, i_row + OutputBegin) = scale * (row(rPerturbedState, i_row) - row(rReferenceState, i_row));
            } // for i_row
        }
}


template <IAdjoint::ResidualTerm Term, class TEntity>
requires (std::is_same_v<TEntity,Element> || std::is_same_v<TEntity,Condition>)
void AdjointFiniteDifferenceUtility<Term,TEntity>::ResizeOutput(
    Matrix& rOutput,
    std::span<const std::size_t> SingleVariateDerivativeShape,
    std::size_t VariableCount) {
        KRATOS_ERROR_IF(SingleVariateDerivativeShape.empty());
        if constexpr (std::is_same_v<Value,Vector>) {
            KRATOS_ERROR_IF_NOT(SingleVariateDerivativeShape.size() == 1);
            rOutput.resize(
                SingleVariateDerivativeShape.front(),
                VariableCount,
                false);
        } else if constexpr (std::is_same_v<Value,Matrix>) {
            KRATOS_ERROR_IF_NOT(SingleVariateDerivativeShape.size() == 2);
            rOutput.resize(
                SingleVariateDerivativeShape.front(),
                VariableCount,
                false);
        } else static_assert(std::is_same_v<Value,void>, "invalid derivative value type");
}


template <IAdjoint::ResidualTerm Term, class TEntity>
requires (std::is_same_v<TEntity,Element> || std::is_same_v<TEntity,Condition>)
void AdjointFiniteDifferenceUtility<Term,TEntity>::EvaluateTerm(
    TEntity& rEntity,
    const ProcessInfo& rProcessInfo,
    Value& rOutput) {
        KRATOS_TRY
            if constexpr (Term == IAdjoint::ResidualTerm::Load) {
                rEntity.CalculateRightHandSide(rOutput, rProcessInfo);
            } else if constexpr (Term == IAdjoint::ResidualTerm::Stiffness) {
                rEntity.CalculateLeftHandSide(rOutput, rProcessInfo);
            } else if constexpr (Term == IAdjoint::ResidualTerm::Damping) {
                rEntity.CalculateDampingMatrix(rOutput, rProcessInfo);
            } else if constexpr (Term == IAdjoint::ResidualTerm::Mass) {
                rEntity.CalculateMassMatrix(rOutput, rProcessInfo);
            } else static_assert(Term == IAdjoint::ResidualTerm::Load, "invalid residual term");
        KRATOS_CATCH("")
}


template <IAdjoint::ResidualTerm Term, class TEntity>
requires (std::is_same_v<TEntity,Element> || std::is_same_v<TEntity,Condition>)
void AdjointFiniteDifferenceUtility<Term,TEntity>::FiniteDifferenceDerivative(
    const TEntity& rEntity,
    const Vector& rValues,
    std::span<const Perturbation> Perturbations,
    Matrix& rOutput,
    std::size_t iBuffer,
    const ProcessInfo& rProcessInfo) const
requires (Term != IAdjoint::ResidualTerm::Load) {
        this->FiniteDifferenceDerivative(
            rEntity,
            &rValues,
            Perturbations,
            rOutput,
            iBuffer,
            rProcessInfo);
}


template <IAdjoint::ResidualTerm Term, class TEntity>
requires (std::is_same_v<TEntity,Element> || std::is_same_v<TEntity,Condition>)
void AdjointFiniteDifferenceUtility<Term,TEntity>::FiniteDifferenceDerivative(
    const TEntity& rEntity,
    std::span<const Perturbation> Perturbations,
    Matrix& rOutput,
    std::size_t iBuffer,
    const ProcessInfo& rProcessInfo) const
requires (Term == IAdjoint::ResidualTerm::Load) {
        this->FiniteDifferenceDerivative(
            rEntity,
            nullptr,
            Perturbations,
            rOutput,
            iBuffer,
            rProcessInfo);
}


template <IAdjoint::ResidualTerm Term, class TEntity>
requires (std::is_same_v<TEntity,Element> || std::is_same_v<TEntity,Condition>)
void AdjointFiniteDifferenceUtility<Term,TEntity>::FiniteDifferenceDerivative(
    const TEntity& rEntity,
    const Vector* pValues,
    std::span<const Perturbation> Perturbations,
    Matrix& rOutput,
    std::size_t iBuffer,
    const ProcessInfo& rProcessInfo) const {
        KRATOS_TRY
            // This wonderful piece of software does not care about petty concepts
            // such as const-correctness, so I have to do a disgusting const cast here.
            TEntity& r_mutable_entity = *const_cast<TEntity*>(&rEntity);

            // Make a copy of the entity, its geometry, nodes, properties and constitutive law.
            // Since the naming and design in Kratos is deliberately misleading, TEntity::Clone
            // does not actually perform a deep copy of the entity, so the above mentioned
            // components have to be copied separately.
            typename TEntity::Pointer p_entity_clone;

            {
                // Copy nodes.
                typename TEntity::NodesArrayType nodes;
                for (Node& r_node : r_mutable_entity.GetGeometry())
                    nodes.push_back(r_node.Clone(r_node.Id()));

                // Clone the entity.
                p_entity_clone = rEntity.Clone(
                    rEntity.Id(),
                    nodes);

                // Copy the DataValueContainer attached to the geometry.
                p_entity_clone->GetGeometry().GetData() = rEntity.GetGeometry().GetData();

                // Copy the entity's properties.
                Properties::Pointer p_properties_clone(new Properties(rEntity.GetProperties()));
            }

            // Compute the unperturbed reference state
            // and initialize the output container.
            rOutput.clear();
            Value reference_state, perturbed_state, finite_difference;
            this->EvaluateTerm(*p_entity_clone, rProcessInfo, reference_state);

            if constexpr (Term == IAdjoint::ResidualTerm::Load) {
                const std::array<std::size_t,1> reference_state_shape {reference_state.size()};
                this->ResizeOutput(
                    rOutput,
                    reference_state_shape,
                    Perturbations.size());
            } else {
                const std::array<std::size_t,2> reference_state_shape {
                    reference_state.size1(),
                    reference_state.size2()};
                this->ResizeOutput(
                    rOutput,
                    reference_state_shape,
                    Perturbations.size());
            }

            // Loop through the perturbations and compute finite differences.
            for (std::size_t i_perturbation=0ul; i_perturbation<Perturbations.size(); ++i_perturbation) {
                const auto& [r_variable, r_context, r_magnitude] = Perturbations[i_perturbation];

                // Perturb, compute the perturbed state, then undo the perturbation.
                switch (r_context) {
                    case Globals::DataLocation::NodeHistorical: {
                        this->Perturb<Globals::DataLocation::NodeHistorical>(
                            *p_entity_clone,
                            r_variable,
                            r_magnitude,
                            iBuffer);
                        this->EvaluateTerm(
                            *p_entity_clone,
                            rProcessInfo,
                            perturbed_state);
                        this->Perturb<Globals::DataLocation::NodeHistorical>(
                            *p_entity_clone,
                            r_variable,
                            -r_magnitude,
                            iBuffer);
                        break;
                    } // case r_context == Globals::DataLocation::NodeHistorical
                    case Globals::DataLocation::NodeNonHistorical: {
                        this->Perturb<Globals::DataLocation::NodeNonHistorical>(
                            *p_entity_clone,
                            r_variable,
                            r_magnitude,
                            iBuffer);
                        this->EvaluateTerm(
                            *p_entity_clone,
                            rProcessInfo,
                            perturbed_state);
                        this->Perturb<Globals::DataLocation::NodeNonHistorical>(
                            *p_entity_clone,
                            r_variable,
                            -r_magnitude,
                            iBuffer);
                        break;
                    } // case r_context == Globals::DataLocation::NodeNonHistorical
                    case Globals::DataLocation::Element: {
                        this->Perturb<Globals::DataLocation::Element>(
                            *p_entity_clone,
                            r_variable,
                            r_magnitude,
                            iBuffer);
                        this->EvaluateTerm(
                            *p_entity_clone,
                            rProcessInfo,
                            perturbed_state);
                        this->Perturb<Globals::DataLocation::Element>(
                            *p_entity_clone,
                            r_variable,
                            -r_magnitude,
                            iBuffer);
                        break;
                    } // case r_context == Globals::DataLocation::Element
                    default: KRATOS_ERROR
                        << "unsupported context " << (int)r_context;
                } // switch r_context

                // Compute finite differences.
                const double scale = 1.0 / r_magnitude;
                if constexpr (Term == IAdjoint::ResidualTerm::Load) {
                    column(rOutput, i_perturbation) = scale * (perturbed_state - reference_state);
                } else {
                    // There's an extra multiplication by -1 because the element
                    // computes the load vector assuming it's on the right hand
                    // side of the equations, while the rest of the terms are on
                    // the left hand side.
                    // Since we're interested in individual terms of the residual,
                    // all terms are on the same side:
                    // R = F - (Ma + Dv + Ku)
                    perturbed_state = -scale * (perturbed_state - reference_state);
                    column(rOutput, i_perturbation) = prod(perturbed_state, *pValues);
                }
            } // for r_variable, r_location in variables
        KRATOS_CATCH("")
}


#define KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY(Term, TEntity)         \
    template class AdjointFiniteDifferenceUtility<Term,TEntity>;

KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY(IAdjoint::ResidualTerm::Load, Element)
KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY(IAdjoint::ResidualTerm::Mass, Element)
KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY(IAdjoint::ResidualTerm::Damping, Element)
KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY(IAdjoint::ResidualTerm::Stiffness, Element)
KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY(IAdjoint::ResidualTerm::Load, Condition)
KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY(IAdjoint::ResidualTerm::Mass, Condition)
KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY(IAdjoint::ResidualTerm::Damping, Condition)
KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY(IAdjoint::ResidualTerm::Stiffness, Condition)


#undef KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY


} // namespace Kratos
