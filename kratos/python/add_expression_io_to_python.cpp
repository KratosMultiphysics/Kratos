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

// System includes
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

// Project includes
#include "add_expression_io_to_python.h"
#include "containers/container_expression/expressions/io/expression_io.h"
#include "containers/container_expression/expressions/io/variable_expression_io.h"
#include "containers/container_expression/expressions/io/c_array_copy_expression_io.h"
#include "containers/container_expression/expressions/io/c_array_move_expression_input.h"
#include <bits/utility.h>


namespace Kratos::Python {


namespace Detail {
class ExpressionInputTrampoline final : public ExpressionInput
{
public:
    Expression::Pointer Execute() const override
    {
        PYBIND11_OVERRIDE_PURE(
            Expression::Pointer, /*return type*/
            ExpressionInput,     /*class name*/
            Execute              /*function name*/
        );
    }
}; // class ExpressionInputTrampoline

class ExpressionOutputTrampoline final : public ExpressionOutput
{
public:
    void Execute(const Expression& rExpression)
    {
        PYBIND11_OVERRIDE_PURE(
            void,               /*return type*/
            ExpressionOutput,   /*class name*/
            Execute             /*function name*/
        );
    }
}; // class ExpressionOutputTrampoline

template <class TVariant>
struct DerefVariant {};

template <class ...TPointers>
struct DerefVariant<std::variant<TPointers...>>
{
    using type = std::variant<typename std::pointer_traits<TPointers>::element_type*...>;
};
} // namespace Detail


void AddExpressionIOToPython(pybind11::module& rModule)
{
    pybind11::class_<ExpressionInput, Detail::ExpressionInputTrampoline>(rModule, "ExpressionInput")
        .def("Execute", &ExpressionInput::Execute)
        ;

    pybind11::class_<ExpressionOutput, Detail::ExpressionOutputTrampoline>(rModule, "ExpressionOutput")
        .def("Execute", &ExpressionOutput::Execute)
        ;

    auto variable_expression_io = rModule.def_submodule("VariableExpressionIO");

    pybind11::enum_<VariableExpressionIO::ContainerType>(variable_expression_io, "ContainerType")
        .value("NodalHistorical", VariableExpressionIO::NodalHistorical)
        .value("NodalNonHistorical", VariableExpressionIO::NodalNonHistorical)
        .value("ElementNonHistorical", VariableExpressionIO::ElementNonHistorical)
        .value("ConditionNonHistorical", VariableExpressionIO::ConditionNonHistorical)
        ;

    pybind11::enum_<VariableExpressionIO::MeshType>(variable_expression_io, "MeshType")
        .value("Local", VariableExpressionIO::Local)
        .value("Interface", VariableExpressionIO::Interface)
        .value("Ghost", VariableExpressionIO::Ghost)
        ;

    pybind11::class_<VariableExpressionIO::VariableExpressionInput, ExpressionInput>(variable_expression_io, "Input")
        //.def(pybind11::init([](const ModelPart& rModelPart,
        //                       const VariableExpressionIO::VariableType& rVariable,
        //                       const VariableExpressionIO::ContainerType& rContainerType,
        //                       const VariableExpressionIO::MeshType& rMeshType) {
        //                        return std::visit([&](const auto& rVar){return VariableExpressionIO::VariableExpressionInput(
        //                            rModelPart,
        //                            *rVar,
        //                            rContainerType,
        //                            rMeshType
        //                        );}, rVariable);
        //                    }),
        //     pybind11::arg("model_part"),
        //     pybind11::arg("variable"),
        //     pybind11::arg("container_type"),
        //     pybind11::arg("mesh_type") = VariableExpressionIO::Local)
        .def(pybind11::init<const ModelPart&,
                            const VariableExpressionIO::VariableType&,
                            const VariableExpressionIO::ContainerType&,
                            const VariableExpressionIO::MeshType&>(),
             pybind11::arg("model_part"),
             pybind11::arg("variable"),
             pybind11::arg("container_type"),
             pybind11::arg("mesh_type") = VariableExpressionIO::Local)
        ;
}


} // namespace Kratos::Python
