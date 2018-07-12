//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


// System includes

// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "nurbs_brep_application_variables.h"
#include "custom_utilities/NodeCurveGeometry3D.h"

namespace py = pybind11;

namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(py::module& m)
{
    using namespace pybind11::literals;

    using Vector3D = Kratos::array_1d<double, 3>;

    // register Interval
    {
        using Type = ANurbs::Interval<double>;

        py::class_<Type>(m, "Interval")
            .def(py::init<double, double>())
            .def_property_readonly("T0", &Type::T0)
            .def_property_readonly("T1", &Type::T1)
        ;
    }

    // register NodeCurveGeometry
    {
        using Type = NodeCurveGeometry3D;
        using Pointer = std::shared_ptr<Type>;

        pybind11::class_<Type, Pointer>(m, "NodeCurveGeometry3D")
            .def(pybind11::init<int, int>(),
                "Degree"_a,
                "NumberOfNodes"_a)
            .def("Knot", &Type::Knot,
                "Index"_a)
            .def("SetKnot", &Type::SetKnot,
                "Index"_a,
                "Value"_a)
            .def("Knots", &Type::Knots)
            .def("Node", &Type::Node,
                "Index"_a)
            .def("SetNode", &Type::SetNode,
                "Index"_a,
                "Value"_a)
            .def("Pole", &Type::Pole,
                "Index"_a)
            .def("SetPole", &Type::SetPole,
                "Index"_a,
                "Value"_a)
            .def("Weight", &Type::Weight,
                "Index"_a)
            .def("SetWeight", &Type::SetWeight,
                "Index"_a,
                "Value"_a)
            .def("Spans", &Type::Spans)
            .def("PointAt", &Type::PointAt,
                "T"_a)
            .def("DerivativesAt", &Type::DerivativesAt,
                "T"_a,
                "Order"_a)
            .def("ValueAt", (double (Type::*)(const Variable<double>&,
                const double)) &Type::ValueAt<double>,
                "Variable"_a,
                "T"_a)
            .def("ValueAt", (std::vector<double> (Type::*)(
                const Variable<double>&, const double, const int))
                &Type::ValueAt<double>,
                "Variable"_a,
                "T"_a,
                "Order"_a)
            .def("ValueAt", (Vector3D (Type::*)(const Variable<Vector3D>&,
                const double)) &Type::ValueAt<Vector3D>,
                "Variable"_a,
                "T"_a)
            .def("ValueAt", (std::vector<Vector3D> (Type::*)(
                const Variable<Vector3D>&, const double, const int))
                &Type::ValueAt<Vector3D>,
                "Variable"_a,
                "T"_a,
                "Order"_a)
        ;
    }
}

} // namespace Python.
} // Namespace Kratos
