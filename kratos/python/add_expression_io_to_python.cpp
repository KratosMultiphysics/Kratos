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
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

// Project includes
#include "add_expression_io_to_python.h"
#include "includes/define_python.h"
#include "containers/container_expression/expressions/expression.h"
#include "containers/container_expression/expressions/io/expression_io.h"
#include "containers/container_expression/expressions/io/variable_expression_io.h"
#include "containers/container_expression/expressions/io/c_array_copy_expression_io.h"
#include "containers/container_expression/expressions/io/c_array_move_expression_input.h"
#include "containers/container_expression/expressions/literal/literal_expression.h"
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"
#include "containers/container_expression/expressions/arithmetic_operators.h"
#include "containers/container_expression/expressions/view_operators.h"


namespace Kratos::Python {


namespace Detail {


class ExpressionTrampoline final : public Expression
{
public:
    using Expression::Expression;

    using Expression::IndexType;

    double Evaluate(const IndexType EntityIndex,
                    const IndexType EntityDataBeginIndex,
                    const IndexType ComponentIndex) const override
    {
        PYBIND11_OVERRIDE_PURE(
            double,                 /*return type*/
            Expression,             /*base type*/
            Evaluate,               /*function name*/
            EntityIndex,
            EntityDataBeginIndex,
            ComponentIndex
        );
    }

    const std::vector<IndexType> GetItemShape() const override
    {
        PYBIND11_OVERRIDE_PURE(
            const std::vector<IndexType>,   /*return type*/
            Expression,                     /*base type*/
            GetItemShape                    /*function name*/
        );
    }

    std::string Info() const override
    {
        PYBIND11_OVERRIDE_PURE(
            std::string,    /*return type*/
            Expression,     /*base type*/
            Info            /*function name*/
        );
    }
}; // class ExpressionTrampoline


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
    void Execute(const Expression& rExpression) override
    {
        PYBIND11_OVERRIDE_PURE(
            void,               /*return type*/
            ExpressionOutput,   /*class name*/
            Execute             /*function name*/
        );
    }
}; // class ExpressionOutputTrampoline

} // namespace Detail


void AddExpressionIOToPython(pybind11::module& rModule)
{
    pybind11::class_<Expression, Expression::Pointer, Detail::ExpressionTrampoline>(rModule, "Expression")
        .def("GetItemShape", &Expression::GetItemShape)
        .def("NumberOfEntities", &Expression::NumberOfEntities)
        .def("GetItemComponentCount", &Expression::GetItemComponentCount)
        .def("__add__", [](Expression::Pointer pLeft, double Right) {return pLeft + Right;})
        //.def("__add__", [](double Left, Expression::Pointer pRight) {return Left + pRight;})
        .def("__add__", [](Expression::Pointer pLeft, Expression::Pointer pRight) {return pLeft + pRight;})
        .def("__sub__", [](Expression::Pointer pLeft, double Right) {return pLeft - Right;})
        //.def("__sub__", [](double Left, Expression::Pointer pRight) {return Left - pRight;})
        .def("__sub__", [](Expression::Pointer pLeft, Expression::Pointer pRight) {return pLeft - pRight;})
        .def("__mul__", [](Expression::Pointer pLeft, double Right) {return pLeft * Right;})
        //.def("__mul__", [](double Left, Expression::Pointer pRight) {retmatekelemenurn Left * pRight;})
        .def("__mul__", [](Expression::Pointer pLeft, Expression::Pointer pRight) {return pLeft * pRight;})
        .def("__truediv__", [](Expression::Pointer pLeft, double Right) {return pLeft / Right;})
        //.def("__truediv__", [](double Left, Expression::Pointer pRight) {return Left / pRight;})
        .def("__truediv__", [](Expression::Pointer pLeft, Expression::Pointer pRight) {return pLeft / pRight;})
        .def("__pow__", [](Expression::Pointer pLeft, double Right) {return Power(pLeft, Right);})
        .def("__pow__", [](Expression::Pointer pLeft, Expression::Pointer pRight) {return Power(pLeft, pRight);})
        .def("__neg__", [](Expression::Pointer pOperand) {return -1.0 * pOperand;})
        .def("__str__", &Expression::Info)
        ;

    pybind11::class_<LiteralExpression<double>, Kratos::intrusive_ptr<LiteralExpression<double>>, Expression>(rModule, "LiteralExpression")
        .def_static("Create",
                    LiteralExpression<double>::Create,
                    pybind11::arg("value"),
                    pybind11::arg("number_of_entities"))
        ;

    pybind11::class_<LiteralFlatExpression<double>, Kratos::intrusive_ptr<LiteralFlatExpression<double>>, Expression>(rModule, "LiteralFlatExpression")
        .def_static("Create", [](pybind11::array_t<double>& rArray,
                                 std::size_t NumberOfItems,
                                 const std::vector<std::size_t>& rShape) {
                                    return LiteralFlatExpression<double>::Create(
                                        rArray.mutable_data(),
                                        NumberOfItems,
                                        rShape
                                    );
                                 },
                    pybind11::arg("array"),
                    pybind11::arg("number_of_items_in _array"),
                    pybind11::arg("item_shape"))
        ;

    auto variable_expression_io = rModule.def_submodule("VariableExpressionIO");

    pybind11::class_<ExpressionInput, Detail::ExpressionInputTrampoline, ExpressionInput::Pointer>(variable_expression_io, "ExpressionInput")
        .def("Execute", &ExpressionInput::Execute)
        ;

    pybind11::class_<ExpressionOutput, Detail::ExpressionOutputTrampoline, ExpressionOutput::Pointer>(variable_expression_io, "ExpressionOutput")
        .def("Execute", &ExpressionOutput::Execute)
        ;

    pybind11::enum_<VariableExpressionIO::ContainerType>(variable_expression_io, "ContainerType")
        .value("NodalHistorical", VariableExpressionIO::NodalHistorical)
        .value("NodalNonHistorical", VariableExpressionIO::NodalNonHistorical)
        .value("ElementNonHistorical", VariableExpressionIO::ElementNonHistorical)
        .value("ConditionNonHistorical", VariableExpressionIO::ConditionNonHistorical)
        ;

    pybind11::enum_<MeshType>(variable_expression_io, "MeshType")
        .value("Local", MeshType::Local)
        .value("Interface", MeshType::Interface)
        .value("Ghost", MeshType::Ghost)
        ;

    pybind11::class_<VariableExpressionIO::VariableExpressionInput, VariableExpressionIO::VariableExpressionInput::Pointer, ExpressionInput>(variable_expression_io, "Input")
        .def(pybind11::init<const ModelPart&,
                            const VariableExpressionIO::VariableType&,
                            const VariableExpressionIO::ContainerType&>(),
             pybind11::arg("model_part"),
             pybind11::arg("variable"),
             pybind11::arg("container_type"))
        .def(pybind11::init<const ContainerExpression<ModelPart::NodesContainerType>&,
                            const VariableExpressionIO::VariableType&,
                            const bool>(),
             pybind11::arg("nodal_container_expression"),
             pybind11::arg("variable"),
             pybind11::arg("is_historical"))
        .def(pybind11::init<const ContainerExpression<ModelPart::ConditionsContainerType>&,
                            const VariableExpressionIO::VariableType&>(),
             pybind11::arg("condition_container_expression"),
             pybind11::arg("variable"))
        .def(pybind11::init<const ContainerExpression<ModelPart::ElementsContainerType>&,
                            const VariableExpressionIO::VariableType&>(),
             pybind11::arg("element_container_expression"),
             pybind11::arg("variable"))
        ;


    pybind11::class_<VariableExpressionIO::VariableExpressionOutput, VariableExpressionIO::VariableExpressionOutput::Pointer, ExpressionOutput>(variable_expression_io, "Output")
        .def(pybind11::init<ModelPart&,
                            const VariableExpressionIO::VariableType&,
                            const VariableExpressionIO::ContainerType&>(),
             pybind11::arg("model_part"),
             pybind11::arg("variable"),
             pybind11::arg("container_type"))
        .def(pybind11::init<ContainerExpression<ModelPart::NodesContainerType>&,
                            const VariableExpressionIO::VariableType&,
                            const bool>(),
             pybind11::arg("nodal_container_expression"),
             pybind11::arg("variable"),
             pybind11::arg("is_historical"))
        .def(pybind11::init<ContainerExpression<ModelPart::ConditionsContainerType>&,
                            const VariableExpressionIO::VariableType&>(),
             pybind11::arg("condition_container_expression"),
             pybind11::arg("variable"))
        .def(pybind11::init<ContainerExpression<ModelPart::ElementsContainerType>&,
                            const VariableExpressionIO::VariableType&>(),
             pybind11::arg("element_container_expression"),
             pybind11::arg("variable"))
        ;

    pybind11::class_<CArrayExpressionInput, CArrayExpressionInput::Pointer, ExpressionInput>(rModule, "CArrayExpressionInput")
        .def(pybind11::init([](const pybind11::array_t<double>& rArray,
                               int NumberOfEntities,
                               const std::vector<int>& rShape){
                                return CArrayExpressionInput(rArray.data(),
                                                             NumberOfEntities,
                                                             rShape.data(),
                                                             rShape.size());
                            }),
             pybind11::arg("array").noconvert(),
             pybind11::arg("number_of_items"),
             pybind11::arg("shape"))
        ;

    pybind11::class_<CArrayExpressionOutput, CArrayExpressionOutput::Pointer, ExpressionOutput>(rModule, "CArrayExpressionOutput")
        .def(pybind11::init([](pybind11::array_t<double>& rArray) {
                                return CArrayExpressionOutput(rArray.mutable_data(), rArray.size());
                               }),
             pybind11::arg("target_array").noconvert());
        ;
}


} // namespace Kratos::Python
