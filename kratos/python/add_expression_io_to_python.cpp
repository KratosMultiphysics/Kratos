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
#include "expression/expression.h"
#include "expression/expression_io.h"
#include "expression/variable_expression_io.h"
#include "expression/c_array_expression_io.h"
#include "expression/literal_expression_input.h"
#include "expression/literal_expression.h"
#include "expression/literal_flat_expression.h"
#include "expression/arithmetic_operators.h"
#include "expression/view_operators.h"


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

template<class TContainerType, class TRawDataType>
void AddCArrayExpressionIOMethods(pybind11::module& rModule)
{
    std::string expression_name;
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        expression_name = "nodal_expression";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        expression_name = "condition_expression";
    } else {
        expression_name = "element_expression";
    }

    rModule.def(
        "Read", [](
            ContainerExpression<TContainerType>& rContainerExpression,
            const pybind11::array_t<TRawDataType>& rArray) {

            KRATOS_ERROR_IF(rArray.ndim() == 0) << "Passed data is not compatible.\n";

            // dimension of the numpy array is always one dimension greater than the kratos stored dimension for each
            // entity. That is because, first dimension of the numpy array shows how many entities are there
            // in the numpy array to be read in. If the numpy array dimension is [45, 3, 4] then it shows
            // there are 45 entities each having matrices of shape [3, 4].
            std::vector<int> shape(rArray.ndim() - 1);
            std::copy(rArray.shape() + 1, rArray.shape() + rArray.ndim(), shape.begin());

            rContainerExpression.SetExpression(
                CArrayExpressionIO::Input(
                    rArray.data(), rArray.shape()[0], shape.data(), shape.size())
                    .Execute());

        },
        pybind11::arg(expression_name.c_str()),
        pybind11::arg("array").noconvert());

    rModule.def(
        "Move", [](
            ContainerExpression<TContainerType>& rContainerExpression,
            pybind11::array_t<TRawDataType>& rArray) {

            KRATOS_ERROR_IF(rArray.ndim() == 0) << "Passed data is not compatible.\n";

            // dimension of the numpy array is always one dimension greater than the kratos stored dimension for each
            // entity. That is because, first dimension of the numpy array shows how many entities are there
            // in the numpy array to be read in. If the numpy array dimension is [45, 3, 4] then it shows
            // there are 45 entities each having matrices of shape [3, 4].
            std::vector<int> shape(rArray.ndim() - 1);
            std::copy(rArray.shape() + 1, rArray.shape() + rArray.ndim(), shape.begin());

            rContainerExpression.SetExpression(
                CArrayExpressionIO::MoveInput(
                    rArray.mutable_data(), rArray.shape()[0], shape.data(), shape.size())
                    .Execute());
        },
        pybind11::arg(expression_name.c_str()),
        pybind11::arg("array").noconvert());

    rModule.def(
        "Write",
        [](const ContainerExpression<TContainerType>& rContainerExpression,
           pybind11::array_t<TRawDataType>& rArray) {
            CArrayExpressionIO::Output(rArray.mutable_data(),
                                                       rArray.size())
                .Execute(rContainerExpression.GetExpression());
        },
        pybind11::arg(expression_name.c_str()),
        pybind11::arg("target_array").noconvert());
}


template<class TContainerType, class TRawDataType>
void AddLiteralExpressionIOMethods(pybind11::module& rModule)
{
    std::string expression_name;
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        expression_name = "nodal_expression";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        expression_name = "condition_expression";
    } else {
        expression_name = "element_expression";
    }

    rModule.def("SetData", &LiteralExpressionIO::SetData<TContainerType>,
                    pybind11::arg(expression_name.c_str()),
                    pybind11::arg("value"));
    rModule.def("SetDataToZero", &LiteralExpressionIO::SetDataToZero<TContainerType>,
                    pybind11::arg(expression_name.c_str()),
                    pybind11::arg("variable"));
}

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
    variable_expression_io.def("Read",
                               &VariableExpressionIO::Read<MeshType::Local>,
                               pybind11::arg("nodal_container_expression"),
                               pybind11::arg("variable"),
                               pybind11::arg("is_historical"));
    variable_expression_io.def("Read",
                               &VariableExpressionIO::Read<ModelPart::ConditionsContainerType, MeshType::Local>,
                               pybind11::arg("condition_container_expression"),
                               pybind11::arg("variable"));
    variable_expression_io.def("Read",
                               &VariableExpressionIO::Read<ModelPart::ElementsContainerType, MeshType::Local>,
                               pybind11::arg("element_container_expression"),
                               pybind11::arg("variable"));
    variable_expression_io.def("Write",
                               &VariableExpressionIO::Write<MeshType::Local>,
                               pybind11::arg("nodal_container_expression"),
                               pybind11::arg("variable"),
                               pybind11::arg("is_historical"));
    variable_expression_io.def("Write",
                               &VariableExpressionIO::Write<ModelPart::ConditionsContainerType, MeshType::Local>,
                               pybind11::arg("condition_container_expression"),
                               pybind11::arg("variable"));
    variable_expression_io.def("Write",
                               &VariableExpressionIO::Write<ModelPart::ElementsContainerType, MeshType::Local>,
                               pybind11::arg("element_container_expression"),
                               pybind11::arg("variable"));

    pybind11::class_<ExpressionInput, Detail::ExpressionInputTrampoline, ExpressionInput::Pointer>(variable_expression_io, "ExpressionInput")
        .def("Execute", &ExpressionInput::Execute)
        ;

    pybind11::class_<ExpressionOutput, Detail::ExpressionOutputTrampoline, ExpressionOutput::Pointer>(variable_expression_io, "ExpressionOutput")
        .def("Execute", &ExpressionOutput::Execute)
        ;

    pybind11::enum_<ContainerType>(rModule, "ContainerType")
        .value("NodalHistorical", ContainerType::NodalHistorical)
        .value("NodalNonHistorical", ContainerType::NodalNonHistorical)
        .value("ElementNonHistorical", ContainerType::ElementNonHistorical)
        .value("ConditionNonHistorical", ContainerType::ConditionNonHistorical)
        ;

    pybind11::enum_<MeshType>(rModule, "MeshType")
        .value("Local", MeshType::Local)
        .value("Interface", MeshType::Interface)
        .value("Ghost", MeshType::Ghost)
        ;

    pybind11::class_<VariableExpressionIO::VariableExpressionInput, VariableExpressionIO::VariableExpressionInput::Pointer, ExpressionInput>(variable_expression_io, "Input")
        .def(pybind11::init<const ModelPart&,
                            const VariableExpressionIO::VariableType&,
                            const ContainerType&>(),
             pybind11::arg("model_part"),
             pybind11::arg("variable"),
             pybind11::arg("container_type"))
        ;

    pybind11::class_<VariableExpressionIO::VariableExpressionOutput, VariableExpressionIO::VariableExpressionOutput::Pointer, ExpressionOutput>(variable_expression_io, "Output")
        .def(pybind11::init<ModelPart&,
                            const VariableExpressionIO::VariableType&,
                            const ContainerType&>(),
             pybind11::arg("model_part"),
             pybind11::arg("variable"),
             pybind11::arg("container_type"))
        ;

    auto carray_expression_io = rModule.def_submodule("CArrayExpressionIO");
    Detail::AddCArrayExpressionIOMethods<ModelPart::NodesContainerType, int>(carray_expression_io);
    Detail::AddCArrayExpressionIOMethods<ModelPart::NodesContainerType, double>(carray_expression_io);
    Detail::AddCArrayExpressionIOMethods<ModelPart::ConditionsContainerType, int>(carray_expression_io);
    Detail::AddCArrayExpressionIOMethods<ModelPart::ConditionsContainerType, double>(carray_expression_io);
    Detail::AddCArrayExpressionIOMethods<ModelPart::ElementsContainerType, int>(carray_expression_io);
    Detail::AddCArrayExpressionIOMethods<ModelPart::ElementsContainerType, double>(carray_expression_io);

    pybind11::class_<CArrayExpressionIO::Input, CArrayExpressionIO::Input::Pointer, ExpressionInput>(
        carray_expression_io, "Input")
        .def(pybind11::init([](const pybind11::array_t<double>& rArray,
                               int NumberOfEntities,
                               const std::vector<int>& rShape) {
                 return CArrayExpressionIO::Input(
                     rArray.data(), NumberOfEntities, rShape.data(), rShape.size());
             }),
             pybind11::arg("array").noconvert(),
             pybind11::arg("number_of_items"),
             pybind11::arg("shape"));

    pybind11::class_<CArrayExpressionIO::Output, CArrayExpressionIO::Output::Pointer, ExpressionOutput>(
        carray_expression_io, "Output")
        .def(pybind11::init([](pybind11::array_t<double>& rArray) {
                 return CArrayExpressionIO::Output(
                     rArray.mutable_data(), rArray.size());
             }),
             pybind11::arg("target_array").noconvert());
    ;

    auto data_expression_io = rModule.def_submodule("LiteralExpressionIO");
    Detail::AddLiteralExpressionIOMethods<ModelPart::NodesContainerType, int>(data_expression_io);
    Detail::AddLiteralExpressionIOMethods<ModelPart::NodesContainerType, double>(data_expression_io);
    Detail::AddLiteralExpressionIOMethods<ModelPart::ConditionsContainerType, int>(data_expression_io);
    Detail::AddLiteralExpressionIOMethods<ModelPart::ConditionsContainerType, double>(data_expression_io);
    Detail::AddLiteralExpressionIOMethods<ModelPart::ElementsContainerType, int>(data_expression_io);
    Detail::AddLiteralExpressionIOMethods<ModelPart::ElementsContainerType, double>(data_expression_io);

    pybind11::class_<LiteralExpressionIO::LiteralExpressionInput, LiteralExpressionIO::LiteralExpressionInput::Pointer, ExpressionInput>(
        data_expression_io, "Input")
        .def(pybind11::init<const ModelPart&,
                            const LiteralExpressionIO::DataType&,
                            const ContainerType&>(),
             pybind11::arg("model_part"),
             pybind11::arg("variable"),
             pybind11::arg("container_type"))
        ;
}


} // namespace Kratos::Python
