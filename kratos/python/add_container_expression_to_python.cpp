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

// Project includes
#include "add_expression_io_to_python.h"
#include "expression/arithmetic_operators.h"
#include "expression/c_array_expression_io.h"
#include "expression/container_expression.h"
#include "expression/expression_utils.h"
#include "includes/define_python.h"
#include "numpy_utils.h"

// Include base h
#include "add_container_expression_to_python.h"

namespace Kratos::Python
{

template<class TContainerType>
void AddContainerExpressionToPython(pybind11::module& m, const std::string& rName)
{
    namespace py = pybind11;

    using container_expression_holder_base = ContainerExpression<TContainerType>;
    py::class_<container_expression_holder_base, typename container_expression_holder_base::Pointer>(m, rName.c_str())
        .def(py::init<ModelPart&>(), py::arg("model_part"))
        .def("CopyFrom", &container_expression_holder_base::CopyFrom, py::arg("origin_container_expression"))
        .def("SetExpression", &container_expression_holder_base::SetExpression)
        .def("HasExpression", &container_expression_holder_base::HasExpression)
        .def("GetExpression", &container_expression_holder_base::pGetExpression)
        .def("GetModelPart", py::overload_cast<>(&container_expression_holder_base::GetModelPart), py::return_value_policy::reference)
        .def("GetContainer", py::overload_cast<>(&container_expression_holder_base::GetContainer), py::return_value_policy::reference)
        .def("GetItemShape", &container_expression_holder_base::GetItemShape)
        .def("GetItemComponentCount", &container_expression_holder_base::GetItemComponentCount)
        .def("GetMaxDepth", &container_expression_holder_base::GetMaxDepth)
        .def("Evaluate", [](const container_expression_holder_base& rSelf){
            const auto& r_shape = rSelf.GetItemShape();
            auto array = AllocateNumpyArray<double>(rSelf.GetContainer().size(), r_shape);
            CArrayExpressionIO::Write(rSelf, array.mutable_data(), array.size());
            return array;
        })
        .def("Clone", &container_expression_holder_base::Clone)
        .def("__add__", [](const container_expression_holder_base& rSelf, const container_expression_holder_base& rOther) { return rSelf + rOther; })
        .def("__iadd__", [](container_expression_holder_base& rSelf, const container_expression_holder_base& rOther) { rSelf = rSelf + rOther; return rSelf; })
        .def("__add__", [](const container_expression_holder_base& rSelf, const double Value) { return rSelf + Value; })
        .def("__iadd__", [](container_expression_holder_base& rSelf, const double Value) { rSelf = rSelf + Value; return rSelf; })
        .def("__sub__", [](const container_expression_holder_base& rSelf, const container_expression_holder_base& rOther) { return rSelf - rOther; })
        .def("__isub__", [](container_expression_holder_base& rSelf, const container_expression_holder_base& rOther) { rSelf = rSelf - rOther; return rSelf; })
        .def("__sub__", [](const container_expression_holder_base& rSelf, const double Value) { return rSelf - Value; })
        .def("__isub__", [](container_expression_holder_base& rSelf, const double Value) { rSelf = rSelf - Value; return rSelf; })
        .def("__mul__", [](const container_expression_holder_base& rSelf, const container_expression_holder_base& rOther) { return rSelf * rOther; })
        .def("__imul__", [](container_expression_holder_base& rSelf, const container_expression_holder_base& rOther) { rSelf = rSelf * rOther; return rSelf; })
        .def("__mul__", [](const container_expression_holder_base& rSelf, const double Value) { return rSelf * Value; })
        .def("__imul__", [](container_expression_holder_base& rSelf, const double Value) { rSelf = rSelf * Value; return rSelf; })
        .def("__truediv__", [](const container_expression_holder_base& rSelf, const container_expression_holder_base& rOther) { return rSelf / rOther; })
        .def("__itruediv__", [](container_expression_holder_base& rSelf, const container_expression_holder_base& rOther) { rSelf = rSelf / rOther; return rSelf; })
        .def("__truediv__", [](const container_expression_holder_base& rSelf, const double Value) { return rSelf / Value; })
        .def("__itruediv__", [](container_expression_holder_base& rSelf, const double Value) { rSelf = rSelf / Value; return rSelf; })
        .def("__pow__", [](container_expression_holder_base& rSelf, const container_expression_holder_base& rInput) { container_expression_holder_base result(rSelf.GetModelPart()); result = ExpressionUtils::Pow(rSelf, rInput); return result; })
        .def("__ipow__", [](container_expression_holder_base& rSelf, const container_expression_holder_base& rInput) { rSelf = ExpressionUtils::Pow(rSelf, rInput); return rSelf; })
        .def("__pow__", [](container_expression_holder_base& rSelf, const double Value) { container_expression_holder_base result(rSelf.GetModelPart()); result = ExpressionUtils::Pow(rSelf, Value); return result; })
        .def("__ipow__", [](container_expression_holder_base& rSelf, const double Value) { rSelf = ExpressionUtils::Pow(rSelf, Value); return rSelf; })
        .def("__neg__", [](container_expression_holder_base& rSelf) { return rSelf * -1.0; })
        .def("PrintData", &container_expression_holder_base::PrintData)
        .def("__str__", &container_expression_holder_base::Info)
        ;
}

template<class TContainerType>
void AddContainerExpressionUtilsToPython(pybind11::module& m, const std::string& rName)
{
     namespace py = pybind11;

     m.def("Collapse", &ExpressionUtils::Collapse<TContainerType>, py::arg(rName.c_str()));
     m.def("Abs", &ExpressionUtils::Abs<TContainerType>, py::arg(rName.c_str()));
     m.def("EntityMin", &ExpressionUtils::EntityMin<TContainerType>, py::arg(rName.c_str()));
     m.def("EntityMax", &ExpressionUtils::EntityMax<TContainerType>, py::arg(rName.c_str()));
     m.def("EntitySum", &ExpressionUtils::EntitySum<TContainerType>, py::arg(rName.c_str()));
     m.def("Sum", &ExpressionUtils::Sum<TContainerType>, py::arg(rName.c_str()));
     m.def("NormInf", &ExpressionUtils::NormInf<TContainerType>, py::arg(rName.c_str()));
     m.def("NormL2", &ExpressionUtils::NormL2<TContainerType>, py::arg(rName.c_str()));
     m.def("NormP", &ExpressionUtils::NormP<TContainerType>, py::arg(rName.c_str()), py::arg("p_value"));
     m.def("Scale", py::overload_cast<const ContainerExpression<TContainerType>&, const double>(&ExpressionUtils::Scale<TContainerType>), py::arg(rName.c_str()), py::arg("scaling_coeff"));
     m.def("Scale", py::overload_cast<const ContainerExpression<TContainerType>&, const ContainerExpression<TContainerType>&>(&ExpressionUtils::Scale<TContainerType>), py::arg(rName.c_str()), py::arg(("scaling_" + rName).c_str()));
     m.def("Pow", py::overload_cast<const ContainerExpression<TContainerType>&, const double>(&ExpressionUtils::Pow<TContainerType>), py::arg(rName.c_str()), py::arg("power"));
     m.def("Pow", py::overload_cast<const ContainerExpression<TContainerType>&, const ContainerExpression<TContainerType>&>(&ExpressionUtils::Pow<TContainerType>), py::arg(rName.c_str()), py::arg(("power_" + rName).c_str()));
     m.def("Slice", &ExpressionUtils::Slice<TContainerType>, py::arg(rName.c_str()), py::arg("offset"), py::arg("stride"));
     m.def("Reshape", py::overload_cast<const ContainerExpression<TContainerType>&, const std::vector<std::size_t>&>(&ExpressionUtils::Reshape<TContainerType>), py::arg(rName.c_str()), py::arg("new_shape"));
     m.def("Comb", py::overload_cast<const std::vector<typename ContainerExpression<TContainerType>::Pointer>&>(&ExpressionUtils::Comb<TContainerType>), py::arg(("other_" + rName + "s").c_str()));
     m.def("InnerProduct", &ExpressionUtils::InnerProduct<TContainerType>, py::arg((rName + "_1").c_str()), py::arg((rName + "_2").c_str()));
}

void  AddContainerExpressionToPython(pybind11::module& m)
{
    auto container_exp_sub_module = m.def_submodule("Expression");

    AddContainerExpressionToPython<ModelPart::NodesContainerType>(container_exp_sub_module, "NodalExpression");
    AddContainerExpressionToPython<ModelPart::ConditionsContainerType>(container_exp_sub_module, "ConditionExpression");
    AddContainerExpressionToPython<ModelPart::ElementsContainerType>(container_exp_sub_module, "ElementExpression");

    AddExpressionIOToPython(container_exp_sub_module);

    auto utils = container_exp_sub_module.def_submodule("Utils");
    AddContainerExpressionUtilsToPython<ModelPart::NodesContainerType>(utils, "nodal_expression");
    AddContainerExpressionUtilsToPython<ModelPart::ConditionsContainerType>(utils, "condition_expression");
    AddContainerExpressionUtilsToPython<ModelPart::ElementsContainerType>(utils, "element_expression");
}

} // namespace Kratos::Python
