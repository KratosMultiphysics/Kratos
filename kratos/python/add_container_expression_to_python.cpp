//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <numeric>

// External includes
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

// Project includes
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"
#include "includes/define_python.h"
#include "containers/container_expression/container_data_io.h"
#include "containers/container_expression/container_expression.h"
#include "containers/container_expression/specialized_container_expression.h"
#include "add_container_expression_to_python_utils.h"

// Include base h
#include "add_container_expression_to_python.h"

namespace Kratos::Python
{

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
} // namespace Detail


void  AddContainerExpressionToPython(pybind11::module& m)
{
    auto sub_module = m.def_submodule("ContainerExpression");

    pybind11::class_<Expression, Kratos::intrusive_ptr<Expression>, Detail::ExpressionTrampoline>(sub_module, "Expression")
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
        //.def("__mul__", [](double Left, Expression::Pointer pRight) {return Left * pRight;})
        .def("__mul__", [](Expression::Pointer pLeft, Expression::Pointer pRight) {return pLeft * pRight;})
        .def("__truediv__", [](Expression::Pointer pLeft, double Right) {return pLeft / Right;})
        //.def("__truediv__", [](double Left, Expression::Pointer pRight) {return Left / pRight;})
        .def("__truediv__", [](Expression::Pointer pLeft, Expression::Pointer pRight) {return pLeft / pRight;})
        .def("__pow__", [](Expression::Pointer pLeft, double Right) {return Pow(pLeft, Right);})
        .def("__pow__", [](Expression::Pointer pLeft, Expression::Pointer pRight) {return Pow(pLeft, pRight);})
        .def("__neg__", [](Expression::Pointer pOperand) {return -1.0 * pOperand;})
        .def("__str__", &Expression::Info)
        ;

    pybind11::class_<LiteralExpression<double>, Kratos::intrusive_ptr<LiteralExpression<double>>, Expression>(sub_module, "LiteralExpression")
        .def_static("Create",
                    LiteralExpression<double>::Create,
                    pybind11::arg("value"),
                    pybind11::arg("number_of_entities"))
        ;

    pybind11::class_<LiteralFlatExpression<double>, Kratos::intrusive_ptr<LiteralFlatExpression<double>>, Expression>(sub_module, "LiteralFlatExpression")
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

    AddContainerExpressionToPython<ModelPart::NodesContainerType>(sub_module, "NodalExpression");
    AddContainerExpressionToPython<ModelPart::ConditionsContainerType>(sub_module, "ConditionExpression");
    AddContainerExpressionToPython<ModelPart::ElementsContainerType>(sub_module, "ElementExpression");

    AddSpecializedContainerExpressionToPython<ModelPart::NodesContainerType, ContainerDataIOTags::Historical>(sub_module, "HistoricalExpression");
    AddSpecializedContainerExpressionToPython<ModelPart::NodesContainerType, ContainerDataIOTags::NonHistorical>(sub_module, "NodalNonHistoricalExpression");
    AddSpecializedContainerExpressionToPython<ModelPart::ConditionsContainerType, ContainerDataIOTags::NonHistorical>(sub_module, "ConditionNonHistoricalExpression");
    AddSpecializedContainerExpressionToPython<ModelPart::ElementsContainerType, ContainerDataIOTags::NonHistorical>(sub_module, "ElementNonHistoricalExpression");
}

} // namespace Kratos::Python
