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

template <int TDimension>
void RegisterCurveGeometryBase(
    pybind11::module& m,
    const std::string& name)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    using VectorType = ANurbs::Point<double, TDimension>;

    using Type = ANurbs::CurveGeometryBase<double, VectorType>;
    using Holder = ANurbs::Pointer<Type>;

    pybind11::class_<Type, Holder>(m, name.c_str());
}

template <int TDimension>
void RegisterCurveGeometry(
    pybind11::module& m,
    const std::string& name)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    using VectorType = ANurbs::Point<double, TDimension>;

    using Type = ANurbs::CurveGeometry<double, TDimension>;
    using Holder = ANurbs::Pointer<Type>;
    using Base = ANurbs::CurveGeometryBase<double, VectorType>;

    pybind11::class_<Type, Base, Holder>(m, name.c_str())
        .def(pybind11::init<int, int, bool>(),
            "Degree"_a,
            "NumberOfPoles"_a,
            "IsRational"_a)
        .def_property_readonly("Degree", &Type::Degree)
        .def_property_readonly("Domain", &Type::Domain)
        .def("Knot", &Type::Knot,
            "Index"_a)
        .def("SetKnot", &Type::SetKnot,
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
    ;
}

template <int TDimension>
void RegisterSurfaceGeometryBase(
    pybind11::module& m,
    const std::string& name)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    using VectorType = ANurbs::Point<double, TDimension>;

    using Type = ANurbs::SurfaceGeometryBase<double, VectorType>;
    using Holder = ANurbs::Pointer<Type>;

    pybind11::class_<Type, Holder>(m, name.c_str());
}

template <int TDimension>
void RegisterSurfaceGeometry(
    pybind11::module& m,
    const std::string& name)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    using VectorType = ANurbs::Point<double, TDimension>;

    using Type = ANurbs::SurfaceGeometry<double, TDimension>;
    using Holder = ANurbs::Pointer<Type>;
    using Base = ANurbs::SurfaceGeometryBase<double, VectorType>;

    pybind11::class_<Type, Base, Holder>(m, name.c_str())
        .def(pybind11::init<int, int, int, int, bool>(),
            "DegreeU"_a,
            "DegreeV"_a,
            "NumberOfPolesU"_a,
            "NumberOfPolesV"_a,
            "IsRational"_a)
        .def_property_readonly("DegreeU", &Type::DegreeU)
        .def_property_readonly("DegreeV", &Type::DegreeV)
        .def_property_readonly("DomainU", &Type::DomainU)
        .def_property_readonly("DomainV", &Type::DomainV)
        .def_property_readonly("NbKnotsU", &Type::NbKnotsU)
        .def_property_readonly("NbKnotsV", &Type::NbKnotsV)
        .def("KnotU", &Type::KnotU,
            "Index"_a)
        .def("KnotV", &Type::KnotV,
            "Index"_a)
        .def("SetKnotU", &Type::SetKnotU,
            "Index"_a,
            "Value"_a)
        .def("SetKnotV", &Type::SetKnotV,
            "Index"_a,
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
    ;
}

template <int TDimension>
void RegisterCurveBase(
    pybind11::module& m,
    const std::string& name)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    using VectorType = ANurbs::Point<double, TDimension>;

    using Type = ANurbs::CurveBase<double, VectorType>;
    using Holder = ANurbs::Pointer<Type>;

    py::class_<Type, Holder>(m, name.c_str())
        .def_property_readonly("Domain", &Type::Domain)
        .def("PointAt", &Type::PointAt,
            "T"_a)
        .def("DerivativesAt", &Type::DerivativesAt,
            "T"_a,
            "Order"_a)
        .def("Spans", &Type::Spans)
    ;
}

template <int TDimension>
void RegisterCurve(
    pybind11::module& m,
    const std::string& name)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    using VectorType = ANurbs::Point<double, TDimension>;
    using GeometryType = ANurbs::CurveGeometry<double, TDimension>;

    using Type = ANurbs::Curve<GeometryType>;
    using Holder = ANurbs::Pointer<Type>;
    using Base = ANurbs::CurveBase<double, VectorType>;

    pybind11::class_<Type, Base, Holder>(m, name.c_str())
        .def(pybind11::init<ANurbs::Pointer<GeometryType>,
            ANurbs::Interval<double>>(),
            "CurveGeometry"_a,
            "Domain"_a)
        .def_property_readonly("Domain", &Type::Domain)
        .def("PointAt", &Type::PointAt,
            "T"_a)
        .def("DerivativesAt", &Type::DerivativesAt,
            "T"_a,
            "Order"_a)
        .def("Spans", &Type::Spans)
    ;
}

template <int TDimension>
void RegisterSurfaceBase(
    pybind11::module& m,
    const std::string& name)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    using VectorType = ANurbs::Point<double, TDimension>;

    using Type = ANurbs::SurfaceBase<double, VectorType>;
    using Holder = ANurbs::Pointer<Type>;

    py::class_<Type, Holder>(m, name.c_str())
        .def_property_readonly("DomainU", &Type::DomainU)
        .def_property_readonly("DomainV", &Type::DomainV)
        .def("PointAt", &Type::PointAt,
            "U"_a,
            "V"_a)
        .def("DerivativesAt", &Type::DerivativesAt,
            "U"_a,
            "V"_a,
            "Order"_a)
        .def("SpansU", &Type::SpansU)
        .def("SpansV", &Type::SpansV)
    ;
}

template <int TDimension>
void RegisterSurface(
    pybind11::module& m,
    const std::string& name)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    using VectorType = ANurbs::Point<double, TDimension>;
    using GeometryType = ANurbs::SurfaceGeometry<double, TDimension>;

    using Type = ANurbs::Surface<GeometryType>;
    using Holder = ANurbs::Pointer<Type>;
    using Base = ANurbs::SurfaceBase<double, VectorType>;

    pybind11::class_<Type, Base, Holder>(m, name.c_str())
        .def(pybind11::init<ANurbs::Pointer<GeometryType>,
            ANurbs::Interval<double>, ANurbs::Interval<double>>(),
            "SurfaceGeometry"_a,
            "DomainU"_a,
            "DomainV"_a)
        .def_property_readonly("DomainU", &Type::DomainU)
        .def_property_readonly("DomainV", &Type::DomainV)
        .def("PointAt", &Type::PointAt,
            "U"_a,
            "V"_a)
        .def("DerivativesAt", &Type::DerivativesAt,
            "U"_a,
            "V"_a,
            "Order"_a)
        .def("SpansU", &Type::SpansU)
        .def("SpansV", &Type::SpansV)
    ;
}

template <int TDimension>
void RegisterCurveOnSurface(
    pybind11::module& m,
    const std::string& name)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Vector2Type = ANurbs::Point<double, 2>;
    using VectorType = ANurbs::Point<double, TDimension>;

    using CurveGeometryBaseType = ANurbs::CurveGeometryBase<double, Vector2Type>;
    using SurfaceGeometryBaseType = ANurbs::SurfaceGeometryBase<double, VectorType>;

    using Type = ANurbs::CurveOnSurface<CurveGeometryBaseType, SurfaceGeometryBaseType>;
    using Pointer = ANurbs::Pointer<Type>;
    using Base = ANurbs::CurveBase<double, VectorType>;

    pybind11::class_<Type, Base, Pointer>(m, name.c_str())
        .def(pybind11::init<ANurbs::Pointer<CurveGeometryBaseType>,
            ANurbs::Pointer<SurfaceGeometryBaseType>,
            ANurbs::Interval<double>>(),
            "CurveGeometry"_a,
            "SurfaceGeometry"_a,
            "Domain"_a)
        .def_property_readonly("Domain", &Type::Domain)
        .def("PointAt", &Type::PointAt,
            "T"_a)
        .def("DerivativesAt", &Type::DerivativesAt,
            "T"_a,
            "Order"_a)
        .def("Spans", &Type::Spans)
    ;
}

template <int TDimension>
void RegisterPointOnCurveProjection(
    pybind11::module& m,
    const std::string& name)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    using VectorType = ANurbs::Point<double, TDimension>;
    using CurveBaseType = ANurbs::CurveBase<double, VectorType>;

    using Type = ANurbs::PointOnCurveProjection<double, VectorType>;
    using Pointer = ANurbs::Pointer<Type>;

    py::class_<Type, Pointer>(m, name.c_str())
        .def(py::init<ANurbs::Pointer<CurveBaseType>, double>(),
            "Curve"_a,
            "Tolerance"_a)
        .def("Compute",
            &Type::Compute)
        .def_property_readonly("Parameter",
            &Type::Parameter)
        .def_property_readonly("Point",
            &Type::Point)
    ;
}

template <int TDimension>
void RegisterCurveTessellation(
    pybind11::module& m,
    const std::string& name)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    using VectorType = ANurbs::Point<double, TDimension>;

    using Type = ANurbs::CurveTessellation<double, VectorType>;
    using Pointer = ANurbs::Pointer<Type>;

    py::class_<Type, Pointer>(m, name.c_str())
        .def(py::init<>())
        .def("Compute", &Type::Compute,
            "Curve"_a,
            "Tolerance"_a)
        .def_property_readonly("NbPoints", &Type::NbPoints)
        .def("Parameter", &Type::Parameter,
            "index"_a)
        .def("Point", &Type::Point,
            "index"_a)
    ;
}

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
            .def("ParameterAtNormalized", (double (Type::*)(const double) const)
                &Type::ParameterAtNormalized,
                "T"_a)
                "T"_a)
            .def("Clamp", (double (Type::*)(const double) const)  &Type::Clamp,
                "T"_a)
        ;
    }

    // register NodeCurveGeometry
    {
        using Type = NodeCurveGeometry3D;
        using Pointer = ANurbs::Pointer<Type>;

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
        using Pointer = ANurbs::Pointer<Type>;

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

    // register Point1D
    {
        using Type = ANurbs::Point1D;

        pybind11::class_<Type>(m, "Point1D")
            .def(pybind11::init<double>(),
                "X"_a)
            .def_property_readonly("X", &Type::X)
        ;
    }

    // register Point2D
    {
        using Type = ANurbs::Point2D;

        pybind11::class_<Type>(m, "Point2D")
            .def(pybind11::init<double, double>(),
                "X"_a,
                "Y"_a)
            .def_property_readonly("X", &Type::X)
            .def_property_readonly("Y", &Type::Y)
        ;
    }

    // register Point3D
    {
        using Type = ANurbs::Point3D;

        pybind11::class_<Type>(m, "Point3D")
            .def(pybind11::init<double, double, double>(),
                "X"_a,
                "Y"_a,
                "Z"_a)
            .def_property_readonly("X", &Type::X)
            .def_property_readonly("Y", &Type::Y)
            .def_property_readonly("Z", &Type::Z)
        ;
    }

    RegisterCurveGeometryBase<1>(m, "CurveGeometryBase1D");
    RegisterCurveGeometryBase<2>(m, "CurveGeometryBase2D");
    RegisterCurveGeometryBase<3>(m, "CurveGeometryBase3D");

    RegisterCurveGeometry<1>(m, "CurveGeometry1D");
    RegisterCurveGeometry<2>(m, "CurveGeometry2D");
    RegisterCurveGeometry<3>(m, "CurveGeometry3D");
    
    RegisterSurfaceGeometryBase<1>(m, "SurfaceGeometryBase1D");
    RegisterSurfaceGeometryBase<2>(m, "SurfaceGeometryBase2D");
    RegisterSurfaceGeometryBase<3>(m, "SurfaceGeometryBase3D");

    RegisterSurfaceGeometry<1>(m, "SurfaceGeometry1D");
    RegisterSurfaceGeometry<2>(m, "SurfaceGeometry2D");
    RegisterSurfaceGeometry<3>(m, "SurfaceGeometry3D");

    RegisterCurveBase<1>(m, "CurveBase1D");
    RegisterCurveBase<2>(m, "CurveBase2D");
    RegisterCurveBase<3>(m, "CurveBase3D");

    RegisterCurve<1>(m, "Curve1D");
    RegisterCurve<2>(m, "Curve2D");
    RegisterCurve<3>(m, "Curve3D");

    RegisterSurfaceBase<1>(m, "SurfaceBase1D");
    RegisterSurfaceBase<2>(m, "SurfaceBase2D");
    RegisterSurfaceBase<3>(m, "SurfaceBase3D");

    RegisterSurface<1>(m, "Surface1D");
    RegisterSurface<2>(m, "Surface2D");
    RegisterSurface<3>(m, "Surface3D");

    RegisterCurveOnSurface<1>(m, "CurveOnSurface1D");
    RegisterCurveOnSurface<2>(m, "CurveOnSurface2D");
    RegisterCurveOnSurface<3>(m, "CurveOnSurface3D");

    RegisterPointOnCurveProjection<2>(m, "PointOnCurveProjection2D");
    RegisterPointOnCurveProjection<3>(m, "PointOnCurveProjection3D");

    RegisterCurveTessellation<2>(m, "CurveTessellation2D");
    RegisterCurveTessellation<3>(m, "CurveTessellation3D");
}

} // namespace Python
} // Namespace Kratos
