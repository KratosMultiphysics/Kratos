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
#include <pybind11/operators.h>

// --- Core Includes ---
#include "add_adjoint_interface_to_python.hpp"
#include "adjoint/adjoint_interface.hpp"
#include "adjoint/response_function.hpp"
#include "adjoint/sensor_response.hpp"

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
        .def(pybind11::self == pybind11::self)
        .def(pybind11::self == VariableData())
        .def(VariableData() == pybind11::self)
        .def(pybind11::self < pybind11::self)
        .def(pybind11::self < VariableData())
        .def(VariableData() < pybind11::self)
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
        std::shared_ptr<IAdjointElement>,
        IAdjoint
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
                return output;},
                pybind11::arg("ResidualTerm"),
                pybind11::arg("rProcessInfo"))
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
                    return output;},
                pybind11::arg("ResidualTerm"),
                pybind11::arg("Variables"),
                pybind11::arg("rProcessInfo"),
                pybind11::arg("iBuffer"))
        ;

    pybind11::class_<
        ResponseFunction,
        ResponseFunction::Pointer,
        IAdjoint
    >(rModule, "ResponseFunction")
        .def(pybind11::init<>())
        .def(pybind11::init([] (std::vector<const Variable<double>*> DesignVariableTypes) {
            return ResponseFunction::Pointer(new ResponseFunction(DesignVariableTypes));
        }))
        .def(
            "ComputeCache",
            &ResponseFunction::ComputeCache,
            pybind11::arg("rModelPart"))
        .def(
            "ClearCache",
            &ResponseFunction::ClearCache)
        .def(
            "ComputeValue",
            pybind11::overload_cast<const ModelPart&,int>(&ResponseFunction::ComputeValue, pybind11::const_),
            pybind11::arg("rModelPart"),
            pybind11::arg("iBuffer"))
        .def(
            "GetStateVariables",
            [] (const ResponseFunction& rThis, const Element& rElement, const ProcessInfo& rProcessInfo) -> std::vector<IAdjoint::DynamicVariable> {
                std::vector<IAdjoint::DynamicVariable> output;
                rThis.GetStateVariables(
                    output,
                    rElement,
                    rProcessInfo);
                return output;
            },
            pybind11::arg("rElement"),
            pybind11::arg("rProcessInfo"))
        .def(
            "GetStateVariables",
            [] (const ResponseFunction& rThis, const Condition& rCondition, const ProcessInfo& rProcessInfo) -> std::vector<IAdjoint::DynamicVariable> {
                std::vector<IAdjoint::DynamicVariable> output;
                rThis.GetStateVariables(
                    output,
                    rCondition,
                    rProcessInfo);
                return output;
            },
            pybind11::arg("rElement"),
            pybind11::arg("rProcessInfo"))
        .def(
            "GetDesignVariables",
            [] (const ResponseFunction& rThis, const Element& rElement, const ProcessInfo& rProcessInfo) -> std::vector<IAdjoint::DynamicVariable> {
                std::vector<IAdjoint::DynamicVariable> output;
                rThis.GetDesignVariables(
                    output,
                    rElement,
                    rProcessInfo);
                return output;
            },
            pybind11::arg("rElement"),
            pybind11::arg("rProcessInfo"))
        .def(
            "GetDesignVariables",
            [] (const ResponseFunction& rThis, const Condition& rCondition, const ProcessInfo& rProcessInfo) -> std::vector<IAdjoint::DynamicVariable> {
                std::vector<IAdjoint::DynamicVariable> output;
                rThis.GetDesignVariables(
                    output,
                    rCondition,
                    rProcessInfo);
                return output;
            },
            pybind11::arg("rElement"),
            pybind11::arg("rProcessInfo"))
        .def(
            "ComputeDerivative",
            [] (
                const ResponseFunction& rThis,
                const Element& rElement,
                const std::vector<IAdjoint::DynamicVariable>& rVariables,
                const ProcessInfo& rProcessInfo,
                int iBuffer) -> Vector {
                    Vector output;
                    rThis.ComputeDerivative(
                        output,
                        rElement,
                        rVariables,
                        rProcessInfo,
                        iBuffer);
                    return output;
                },
            pybind11::arg("rElement"),
            pybind11::arg("rVariables"),
            pybind11::arg("rProcessInfo"),
            pybind11::arg("iBuffer"))
        .def(
            "ComputeDerivative",
            [] (
                const ResponseFunction& rThis,
                const Condition& rCondition,
                const std::vector<IAdjoint::DynamicVariable>& rVariables,
                const ProcessInfo& rProcessInfo,
                int iBuffer) -> Vector {
                    Vector output;
                    rThis.ComputeDerivative(
                        output,
                        rCondition,
                        rVariables,
                        rProcessInfo,
                        iBuffer);
                    return output;
                },
            pybind11::arg("rCondition"),
            pybind11::arg("rVariables"),
            pybind11::arg("rProcessInfo"),
            pybind11::arg("iBuffer"))
        .def(
            "SetDesignVariableTypes",
            [] (ResponseFunction& rThis, const std::vector<const Variable<double>*>& rVariableTypes) -> void {
                rThis.SetDesignVariableTypes(rVariableTypes);
            },
            pybind11::arg("rVariableTypes"))
        .def(
            "GetDesignVariableTypes",
            [] (const ResponseFunction& rThis) -> std::vector<const Variable<double>*> {
                const auto variables = rThis.GetDesignVariableTypes();
                return std::vector<const Variable<double>*>(variables.begin(), variables.end());
            })
        ;

    pybind11::class_<
        SensorResponse,
        SensorResponse::Pointer,
        ResponseFunction
    >(rModule, "SensorResponse")
        .def(pybind11::init<>())
        .def(pybind11::init([] (
            const std::vector<Variable<double>*> rDesignVariableTypes,
            const std::string& rName,
            Node::Pointer pNode,
            Element::Pointer pElement) -> SensorResponse::Pointer {
                std::span<const Variable<double>* const> design_variable_types(
                    rDesignVariableTypes.data(),
                    rDesignVariableTypes.size());
                return std::make_shared<SensorResponse>(
                    design_variable_types,
                    rName,
                    pNode,
                    pElement);
            }),
            pybind11::arg("rDesignVariableTypes"),
            pybind11::arg("rName"),
            pybind11::arg("pNode"),
            pybind11::arg("pElement"))
        .def(
            "Name",
            &SensorResponse::Name)
        .def(
            "GetNode",
            [] (SensorResponse& rThis) -> Node& {
                return *rThis.GetNode();
            },
            pybind11::return_value_policy::reference)
        .def(
            "AddTensorAdaptor",
            &SensorResponse::AddTensorAdaptor,
            pybind11::arg("name"),
            pybind11::arg("adaptor"))
        .def(
            "GetTensorAdaptor",
            &SensorResponse::GetTensorAdaptor,
            pybind11::arg("name"))
        .def(
            "ClearTensorAdaptors",
            &SensorResponse::ClearTensorAdaptors)
        .def(
            "GetContainingElement",
            &SensorResponse::GetContainingElement)
        ;
} // AddAdjointInterfaceToPython


} // namespace Kratos::Python
