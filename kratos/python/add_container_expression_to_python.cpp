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
        .def("Slice",
             &container_expression_holder_base::Slice,
             py::arg("offset"),
             py::arg("stride"))
        .def("Reshape",
             &container_expression_holder_base::Reshape,
             py::arg("new_shape"))
        .def("Comb",
             [](container_expression_holder_base& rSelf,
                const container_expression_holder_base& rOther)
                {return rSelf.Comb(rOther);},
             py::arg("other"))
        .def("Comb",
             [](container_expression_holder_base& rSelf,
                const std::vector<typename container_expression_holder_base::Pointer>& rOthers)
                {return rSelf.Comb(rOthers);},
             py::arg("others"))
        .def("Evaluate", [](const container_expression_holder_base& rSelf){
            const auto& r_shape = rSelf.GetItemShape();
            auto array = AllocateNumpyArray<double>(rSelf.GetContainer().size(), r_shape);
            CArrayExpressionIO::Write(rSelf, array.mutable_data(), array.size());
            return array;
        })
        .def("Clone", &container_expression_holder_base::Clone)
        .def("Scale", [](const container_expression_holder_base& rSelf, const container_expression_holder_base& rOther){auto copy = rSelf; copy.SetExpression(Scale(rSelf.pGetExpression(), rOther.pGetExpression())); return copy;})
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
        .def("__pow__", [](container_expression_holder_base& rSelf, const container_expression_holder_base& rInput) { container_expression_holder_base result(rSelf.GetModelPart()); result = Power(rSelf, rInput); return result; })
        .def("__ipow__", [](container_expression_holder_base& rSelf, const container_expression_holder_base& rInput) { rSelf = Power(rSelf, rInput); return rSelf; })
        .def("__pow__", [](container_expression_holder_base& rSelf, const double Value) { container_expression_holder_base result(rSelf.GetModelPart()); result = Power(rSelf, Value); return result; })
        .def("__ipow__", [](container_expression_holder_base& rSelf, const double Value) { rSelf = Power(rSelf, Value); return rSelf; })
        .def("__neg__", [](container_expression_holder_base& rSelf) { return rSelf *= -1.0; })
        .def("PrintData", &container_expression_holder_base::PrintData)
        .def("__str__", &container_expression_holder_base::Info)
        ;
}

void  AddContainerExpressionToPython(pybind11::module& m)
{
    auto container_exp_sub_module = m.def_submodule("Expression");

    AddContainerExpressionToPython<ModelPart::NodesContainerType>(container_exp_sub_module, "NodalExpression");
    AddContainerExpressionToPython<ModelPart::ConditionsContainerType>(container_exp_sub_module, "ConditionExpression");
    AddContainerExpressionToPython<ModelPart::ElementsContainerType>(container_exp_sub_module, "ElementExpression");

    AddExpressionIOToPython(container_exp_sub_module);
}

} // namespace Kratos::Python
