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

// --- External Includes ---
#include <pybind11/stl.h>

// --- Core Includes ---
#include "adjoint/adjoint_interface.hpp"

// --- STL Includes ---
#include <memory>


namespace Kratos::Python {


void AddAdjointInterfaceToPython(pybind11::module_& rModule) {
    auto adjoint_interface_bindings = pybind11::class_<IAdjoint,std::shared_ptr<IAdjoint>>(
        rModule,
        "IAdjoint");
    adjoint_interface_bindings
        .def(pybind11::init<>())
        .def_static(
            "TermName",
            &IAdjoint::TermName)
        ;

    pybind11::class_<
        IAdjoint::DynamicVariable,
        VariableData
    >(adjoint_interface_bindings, "DynamicVariable")
        .def(pybind11::init<>())
        .def(pybind11::init<const VariableData&>())
        .def(pybind11::init<const VariableData&,std::size_t>())
        .def(
            "GetDynamicIndex",
            &IAdjoint::DynamicVariable::GetDynamicIndex)
        .def(
            "SetDynamicIndex",
            &IAdjoint::DynamicVariable::SetDynamicIndex)
        ;

    pybind11::enum_<IAdjoint::ResidualTerm>(
        adjoint_interface_bindings,
        "ResidualTerm")
            .value("Load", IAdjoint::ResidualTerm::Load)
            .value("Stiffness", IAdjoint::ResidualTerm::Stiffness)
            .value("Damping", IAdjoint::ResidualTerm::Damping)
            .value("Mass", IAdjoint::ResidualTerm::Mass)
            ;

    pybind11::class_<
        IAdjointElement,
        IAdjoint,
        std::shared_ptr<IAdjointElement>
    >(rModule, "IAdjointElement")
        .def(
            "GetInfluencingVariables",
            [] (
                const IAdjointElement& rThis,
                IAdjoint::ResidualTerm Term,
                const ProcessInfo& rProcessInfo) {
                    std::vector<IAdjoint::DynamicVariable> output;
                    switch (Term) {
                        case IAdjoint::ResidualTerm::Load:
                            rThis.GetInfluencingVariables<IAdjoint::ResidualTerm::Load>(
                                output,
                                rProcessInfo);
                            break;
                        case IAdjoint::ResidualTerm::Stiffness:
                            rThis.GetInfluencingVariables<IAdjoint::ResidualTerm::Stiffness>(
                                output,
                                rProcessInfo);
                            break;
                        case IAdjoint::ResidualTerm::Damping:
                            rThis.GetInfluencingVariables<IAdjoint::ResidualTerm::Damping>(
                                output,
                                rProcessInfo);
                            break;
                        case IAdjoint::ResidualTerm::Mass:
                            rThis.GetInfluencingVariables<IAdjoint::ResidualTerm::Mass>(
                                output,
                                rProcessInfo);
                            break;
                        default:
                            KRATOS_ERROR << "invalid residual term '" << IAdjoint::TermName(Term) << "'";
                    }
                return output;})
        .def(
            "ComputeDerivative",
            [] (
                const IAdjointElement& rThis,
                IAdjoint::ResidualTerm Term,
                const std::vector<IAdjoint::DynamicVariable> Variables,
                const ProcessInfo& rProcessInfo,
                int iBuffer) {
                    Matrix output;
                    switch (Term) {
                        case IAdjoint::ResidualTerm::Stiffness:
                            rThis.ComputeDerivative<IAdjoint::ResidualTerm::Stiffness>(
                                output,
                                Variables,
                                rProcessInfo,
                                iBuffer);
                            break;
                        case IAdjoint::ResidualTerm::Damping:
                            rThis.ComputeDerivative<IAdjoint::ResidualTerm::Damping>(
                                output,
                                Variables,
                                rProcessInfo,
                                iBuffer);
                            break;
                        case IAdjoint::ResidualTerm::Mass:
                            rThis.ComputeDerivative<IAdjoint::ResidualTerm::Mass>(
                                output,
                                Variables,
                                rProcessInfo,
                                iBuffer);
                            break;
                        case IAdjoint::ResidualTerm::Load:
                            rThis.ComputeDerivative<IAdjoint::ResidualTerm::Load>(
                                output,
                                Variables,
                                rProcessInfo,
                                iBuffer);
                            break;
                        default:
                            KRATOS_ERROR << "invalid residual term '" << IAdjoint::TermName(Term) << "'";
                    }
                    return output;})
        ;
} // AddAdjointInterfaceToPython


} // namespace Kratos::Python
