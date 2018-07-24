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

namespace py = pybind11;

namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(py::module& m)
{
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COORDINATES)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TANGENTS)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CROSS_AREA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRESTRESS_CAUCHY)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHAPE_FUNCTION_VALUES)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHAPE_FUNCTION_LOCAL_DERIVATIVES)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RAYLEIGH_ALPHA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RAYLEIGH_BETA)

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

    // register CurveShapeEvaluator
    {
        using Type = ANurbs::CurveShapeEvaluator<double>;

        pybind11::class_<Type>(m, "CurveShapeEvaluator")
            .def(pybind11::init<int, int>(),
                "Degree"_a,
                "Order"_a)
            .def("Resize", &Type::Resize,
                "Degree"_a,
                "Order"_a)
            .def_property_readonly("Degree", &Type::Degree)
            .def_property_readonly("Order", &Type::Order)
            .def_property_readonly("NumberOfNonzeroPoles",
                &Type::NbNonzeroPoles)
            .def_property_readonly("FirstNonzeroPole",
                &Type::FirstNonzeroPole)
            .def_property_readonly("LastNonzeroPole",
                &Type::LastNonzeroPole)
            .def_property_readonly("NumberOfShapes", &Type::NbShapes)
            .def("__call__", &Type::Value,
                "Order"_a,
                "Pole"_a)
            .def("Compute", &Type::Compute<std::vector<double>>,
                "Knots"_a,
                "T"_a)
            .def("Compute", &Type::Compute<std::vector<double>,
                std::vector<double>>,
                "Knots"_a,
                "Weights"_a,
                "T"_a)
        ;
    }

    // register CurveShapeEvaluator
    {
        using Type = ANurbs::SurfaceShapeEvaluator<double>;

        pybind11::class_<Type>(m, "SurfaceShapeEvaluator")
            .def(pybind11::init<int, int, int>(),
                "DegreeU"_a,
                "DegreeV"_a,
                "Order"_a)
            .def("Resize", &Type::Resize,
                "DegreeU"_a,
                "DegreeV"_a,
                "Order"_a)
            .def_property_readonly("DegreeU", &Type::DegreeU)
            .def_property_readonly("DegreeV", &Type::DegreeV)
            .def_property_readonly("Order", &Type::Order)
            .def_property_readonly("NumberOfShapes", (int (Type::*)(void) const)
                &Type::NbShapes)
            .def_property_readonly("NumberOfNonzeroPolesU",
                &Type::NbNonzeroPolesU)
            .def_property_readonly("NumberOfNonzeroPolesV",
                &Type::NbNonzeroPolesV)
            .def_property_readonly("FirstNonzeroPoleU",
                &Type::FirstNonzeroPoleU)
            .def_property_readonly("FirstNonzeroPoleV",
                &Type::FirstNonzeroPoleV)
            .def_property_readonly("LastNonzeroPoleU",
                &Type::LastNonzeroPoleU)
            .def_property_readonly("LastNonzeroPoleV",
                &Type::LastNonzeroPoleV)
            .def("__call__", (double (Type::*)(const int, const int, const int)
                const) &Type::operator(),
                "Derivative"_a,
                "PoleU"_a,
                "PoleV"_a)
            .def("Compute", &Type::Compute<std::vector<double>>,
                "KnotsU"_a,
                "KnotsV"_a,
                "U"_a,
                "V"_a)
            .def("Compute", &Type::Compute<std::vector<double>,
                ANurbs::Grid<double>>,
                "KnotsU"_a,
                "KnotsV"_a,
                "Weights"_a,
                "U"_a,
                "V"_a)
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

    // register IntegrationPoint1
    {
        using Type = ANurbs::IntegrationPoint1<double>;

        pybind11::class_<Type>(m, "IntegrationPoint1")
            .def("__iter__", 
                [](const Type &self) {
                    return pybind11::make_iterator(&self.t, &self.t + 2);
                }, pybind11::keep_alive<0, 1>())
            .def_readwrite("t", &Type::t)
            .def_readwrite("weight", &Type::weight)
        ;
    }

    // register IntegrationPoint2
    {
        using Type = ANurbs::IntegrationPoint2<double>;

        pybind11::class_<Type>(m, "IntegrationPoint2")
            .def_readwrite("u", &Type::u)
            .def_readwrite("v", &Type::v)
            .def_readwrite("weight", &Type::weight)
        ;
    }

    // register IntegrationPoints
    {
        using Type = ANurbs::IntegrationPoints<double>;

        pybind11::class_<Type>(m, "IntegrationPoints")
            .def_static("Points1", &Type::Points1,
                "Degree"_a,
                "Domain"_a)
            .def_static("Points2", &Type::Points2,
                "DegreeU"_a,
                "DegreeV"_a,
                "DomainU"_a,
                "DomainV"_a)
        ;
    }
}

} // namespace Python
} // Namespace Kratos
