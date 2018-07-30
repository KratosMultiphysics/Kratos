/*
//  KRATOS .___  ________    _____
//         |   |/  _____/   /  _  \
//         |   /   \  ___  /  /_\  \
//         |   \    \_\  \/    |    \
//         |___|\______  /\____|__  /
//                     \/         \/  Application
//
//  License: BSD License
//           Kratos default license: kratos/license.txt
*/

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "iga_application_variables.h"
#include "custom_utilities/node_curve_geometry_3d.h"
#include "custom_utilities/node_surface_geometry_3d.h"

namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Array3D = Kratos::array_1d<double, 3>;

    // register Interval
    {
        using Type = ANurbs::Interval<double>;

        py::class_<Type>(m, "Interval")
            .def(py::init<double, double>())
            .def_property_readonly("T0", &Type::T0)
            .def_property_readonly("T1", &Type::T1)
            .def_property_readonly("Min", &Type::Min)
            .def_property_readonly("Max", &Type::Max)
            .def_property_readonly("Delta", &Type::Delta)
            .def_property_readonly("Length", &Type::Length)
            .def("NormalizedAt", &Type::NormalizedAt,
                "T"_a)
            .def("ParameterAtNormalized", &Type::ParameterAtNormalized,
                "T"_a)
            .def("Clamp", (double (Type::*)(const double) const)  &Type::Clamp,
                "T"_a)
        ;
    }

    // register NodeCurveGeometry
    {
        using Type = NodeCurveGeometry3D;
        using Pointer = std::shared_ptr<Type>;

        using VariableComponent = Kratos::VariableComponent<
            Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3>>>;

        pybind11::class_<Type, Pointer>(m, "NodeCurveGeometry3D")
            .def(pybind11::init<int, int>(),
                "Degree"_a,
                "NumberOfNodes"_a)
            .def_property_readonly("Degree", &Type::Degree)
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
                const double) const) &Type::ValueAt<double>,
                "Variable"_a,
                "T"_a)
            .def("ValueAt", (std::vector<double> (Type::*)(
                const Variable<double>&, const double, const int) const)
                &Type::ValueAt<double>,
                "Variable"_a,
                "T"_a,
                "Order"_a)
            .def("ValueAt", (Array3D (Type::*)(const Variable<Array3D>&,
                const double) const) &Type::ValueAt<Array3D>,
                "Variable"_a,
                "T"_a)
            .def("ValueAt", (std::vector<Array3D> (Type::*)(
                const Variable<Array3D>&, const double, const int) const)
                &Type::ValueAt<Array3D>,
                "Variable"_a,
                "T"_a,
                "Order"_a)
            .def("ValueAt", (double (Type::*)(const VariableComponent&,
                const double) const) &Type::ValueAt<double, VariableComponent>,
                "Variable"_a,
                "T"_a)
            .def("ValueAt", (std::vector<double> (Type::*)(
                const VariableComponent&, const double,
                const int) const) &Type::ValueAt<double, VariableComponent>,
                "Variable"_a,
                "T"_a,
                "Order"_a)
        ;
    }

    // register NodeSurfaceGeometry
    {
        using Type = NodeSurfaceGeometry3D;
        using Pointer = std::shared_ptr<Type>;

        using VariableComponent = Kratos::VariableComponent<
            Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3>>>;

        pybind11::class_<Type, Pointer>(m, "NodeSurfaceGeometry3D")
            .def(pybind11::init<int, int, int, int>(),
                "DegreeU"_a,
                "DegreeV"_a,
                "NumberOfNodesU"_a,
                "NumberOfNodesV"_a)
            .def_property_readonly("DegreeU", &Type::DegreeU)
            .def_property_readonly("DegreeV", &Type::DegreeV)
            .def("SetKnotU", &Type::SetKnotU,
                "Index"_a,
                "Value"_a)
            .def("SetKnotV", &Type::SetKnotV,
                "Index"_a,
                "Value"_a)
            .def("Node", &Type::Node,
                "IndexU"_a,
                "IndexV"_a)
            .def("SetNode", &Type::SetNode,
                "IndexU"_a,
                "IndexV"_a,
                "Value"_a)
            .def("Pole", &Type::Pole,
                "IndexU"_a,
                "IndexV"_a)
            .def("SetPole", &Type::SetPole,
                "IndexU"_a,
                "IndexV"_a,
                "Value"_a)
            .def("Weight", &Type::Weight,
                "IndexU"_a,
                "IndexV"_a)
            .def("SetWeight", &Type::SetWeight,
                "IndexU"_a,
                "IndexV"_a,
                "Value"_a)
            .def("PointAt", &Type::PointAt,
                "U"_a,
                "V"_a)
            .def("DerivativesAt", &Type::DerivativesAt,
                "U"_a,
                "V"_a,
                "Order"_a)
            .def("ValueAt", (double (Type::*)(const Variable<double>&,
                const double, const double) const) &Type::ValueAt<double>,
                "Variable"_a,
                "U"_a,
                "V"_a)
            .def("ValueAt", (std::vector<double> (Type::*)(
                const Variable<double>&, const double, const double, const int)
                const)
                &Type::ValueAt<double>,
                "Variable"_a,
                "U"_a,
                "V"_a,
                "Order"_a)
            .def("ValueAt", (Array3D (Type::*)(const Variable<Array3D>&,
                const double, const double) const) &Type::ValueAt<Array3D>,
                "Variable"_a,
                "U"_a,
                "V"_a)
            .def("ValueAt", (std::vector<Array3D> (Type::*)(
                const Variable<Array3D>&, const double, const double,
                const int) const) &Type::ValueAt<Array3D>,
                "Variable"_a,
                "U"_a,
                "V"_a,
                "Order"_a)
            .def("ValueAt", (double (Type::*)(const VariableComponent&,
                const double, const double) const) &Type::ValueAt<double,
                VariableComponent>,
                "Variable"_a,
                "U"_a,
                "V"_a)
            .def("ValueAt", (std::vector<double> (Type::*)(
                const VariableComponent&, const double, const double, const int)
                const) &Type::ValueAt<double, VariableComponent>,
                "Variable"_a,
                "U"_a,
                "V"_a,
                "Order"_a)
        ;
    }
}

} // namespace Python
} // Namespace Kratos
