/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "iga_application_variables.h"
#include "custom_utilities/anurbs.h"
#include "custom_utilities/node_curve_geometry_3d.h"
#include "custom_utilities/node_surface_geometry_3d.h"

#include "custom_utilities/brep_json_io.h"
#include "custom_utilities/nurbs_brep_modeler.h"


namespace Kratos {
namespace Python {

void RegisterPoint1D(
    pybind11::module& m,
    const std::string& name)
{

    using namespace pybind11::literals;

    using Type = Kratos::array_1d<double, 1>;

    pybind11::class_<Type>(m, name.c_str())
        .def(pybind11::init<>([](const double x) -> Type {
                Type vector;
                vector[0] = x;
                return vector;
            }),
            "X"_a)
        .def_property_readonly("X", [](const Type& self) -> double {
                return self[0];
            })
    ;
}

void RegisterPoint2D(
    pybind11::module& m,
    const std::string& name)
{
    using namespace pybind11::literals;

    using Type = Kratos::array_1d<double, 2>;

    pybind11::class_<Type>(m, name.c_str())
        .def(pybind11::init<>([](const double x, const double y) -> Type {
                Type vector;
                vector[0] = x;
                vector[1] = y;
                return vector;
            }),
            "X"_a,
            "Y"_a)
        .def_property_readonly("X", [](const Type& self) -> double {
                return self[0];
            })
        .def_property_readonly("Y", [](const Type& self) -> double {
                return self[1];
            })
    ;
}

void RegisterInterval(
    pybind11::module& m,
    const std::string& name)
{
    using namespace pybind11::literals;

    using Type = ANurbs::Interval<double>;

    pybind11::class_<Type>(m, name.c_str())
        .def(pybind11::init<double, double>())
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
        .def("Clamp", (double (Type::*)(const double) const) &Type::Clamp,
            "T"_a)
    ;
}

void RegisterCurveShapeEvaluator(
    pybind11::module& m,
    const std::string& name)
{
    using namespace pybind11::literals;

    using Type = ANurbs::CurveShapeEvaluator<double>;

    pybind11::class_<Type>(m, name.c_str())
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

void RegisterSurfaceShapeEvaluator(
    pybind11::module& m,
    const std::string& name)
{
    using namespace pybind11::literals;

    using Type = ANurbs::SurfaceShapeEvaluator<double>;

    pybind11::class_<Type>(m, name.c_str())
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
        .def_property_readonly("NumberOfNonzeroPoles",
            &Type::NbNonzeroPoles)
        .def_property_readonly("FirstNonzeroPoleU",
            &Type::FirstNonzeroPoleU)
        .def_property_readonly("FirstNonzeroPoleV",
            &Type::FirstNonzeroPoleV)
        .def_property_readonly("LastNonzeroPoleU",
            &Type::LastNonzeroPoleU)
        .def_property_readonly("LastNonzeroPoleV",
            &Type::LastNonzeroPoleV)
        .def_property_readonly("NonzeroPoleIndices", &Type::NonzeroPoleIndices)
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

template <int TDimension>
void RegisterCurveGeometryBase(
    pybind11::module& m,
    const std::string& name)
{
    using namespace pybind11::literals;

    using Type = CurveGeometryBase<TDimension>;
    using Holder = shared_ptr<Type>;

    pybind11::class_<Type, Holder>(m, name.c_str())
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
void RegisterCurveGeometry(
    pybind11::module& m,
    const std::string& name)
{
    using namespace pybind11::literals;

    using Type = CurveGeometry<TDimension>;
    using Holder = shared_ptr<Type>;
    using Base = CurveGeometryBase<TDimension>;

    pybind11::class_<Type, Base, Holder>(m, name.c_str())
        .def(pybind11::init<int, int, bool>(),
            "Degree"_a,
            "NumberOfPoles"_a,
            "IsRational"_a)
    ;
}

template <int TDimension>
void RegisterSurfaceGeometryBase(
    pybind11::module& m,
    const std::string& name)
{
    using namespace pybind11::literals;

    using VectorType = Kratos::array_1d<double, TDimension>;

    using Type = SurfaceGeometryBase<TDimension>;
    using Holder = shared_ptr<Type>;

    pybind11::class_<Type, Holder>(m, name.c_str())
        .def_property_readonly("DegreeU", &Type::DegreeU)
        .def_property_readonly("DegreeV", &Type::DegreeV)
        .def_property_readonly("DomainU", &Type::DomainU)
        .def_property_readonly("DomainV", &Type::DomainV)
        .def_property_readonly("NumberOfKnotsU", &Type::NbKnotsU)
        .def_property_readonly("NumberOfKnotsV", &Type::NbKnotsV)
        .def_property_readonly("NumberOfPolesU", &Type::NbPolesU)
        .def_property_readonly("NumberOfPolesV", &Type::NbPolesV)
        .def("KnotsU", &Type::KnotsU)
        .def("KnotsV", &Type::KnotsV)
        .def("SpansU", &Type::SpansU)
        .def("SpansV", &Type::SpansV)
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
        .def("Pole", (VectorType (Type::*)(const int, const int) const)
            &Type::Pole,
            "IndexU"_a,
            "IndexV"_a)
        .def("Pole", (VectorType (Type::*)(const int) const) &Type::Pole,
            "Index"_a)
        .def("SetPole", (void (Type::*)(const int, const VectorType&))
            &Type::SetPole,
            "Index"_a,
            "Value"_a)
        .def("SetPole", (void (Type::*)(const int, const int,
            const VectorType&)) &Type::SetPole,
            "IndexU"_a,
            "IndexV"_a,
            "Value"_a)
        .def("Weight", &Type::Weight,
            "IndexU"_a,
            "IndexV"_a)
        .def("SetWeight", (void (Type::*)(const int, const double))
            &Type::SetWeight,
            "Index"_a,
            "Value"_a)
        .def("SetWeight", (void (Type::*)(const int, const int,
            const double)) &Type::SetWeight,
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
void RegisterSurfaceGeometry(
    pybind11::module& m,
    const std::string& name)
{
    using namespace pybind11::literals;

    using Type = SurfaceGeometry<TDimension>;
    using Holder = shared_ptr<Type>;
    using Base = SurfaceGeometryBase<TDimension>;

    pybind11::class_<Type, Base, Holder>(m, name.c_str())
        .def(pybind11::init<int, int, int, int, bool>(),
            "DegreeU"_a,
            "DegreeV"_a,
            "NumberOfPolesU"_a,
            "NumberOfPolesV"_a,
            "IsRational"_a);
    ;
}

void RegisterNodeCurveGeometry(
    pybind11::module& m,
    const std::string& name)
{
    using namespace pybind11::literals;

    using VectorType = Kratos::array_1d<double, 3>;

    using Type = NodeCurveGeometry3D;
    using Holder = shared_ptr<Type>;
    using Base = CurveGeometryBase<3>;

    using VariableComponent = Kratos::VariableComponent<
        Kratos::VectorComponentAdaptor<VectorType>>;

    pybind11::class_<Type, Base, Holder>(m, name.c_str())
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
        .def("ValueAt", (VectorType (Type::*)(const Variable<VectorType>&,
            const double) const) &Type::ValueAt<VectorType>,
            "Variable"_a,
            "T"_a)
        .def("ValueAt", (std::vector<VectorType> (Type::*)(
            const Variable<VectorType>&, const double, const int) const)
            &Type::ValueAt<VectorType>,
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

void RegisterNodeSurfaceGeometry(
    pybind11::module& m,
    const std::string& name)
{
    using namespace pybind11::literals;

    using VectorType = Kratos::array_1d<double, 3>;

    using Type = NodeSurfaceGeometry3D;
    using Holder = shared_ptr<Type>;
    using Base = ANurbs::SurfaceGeometryBase<VectorType>;

    using VariableComponent = Kratos::VariableComponent<
        Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3>>>;

    pybind11::class_<Type, Base, Holder>(m, name.c_str())
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
        .def("ValueAt", (VectorType (Type::*)(const Variable<VectorType>&,
            const double, const double) const) &Type::ValueAt<VectorType>,
            "Variable"_a,
            "U"_a,
            "V"_a)
        .def("ValueAt", (std::vector<VectorType> (Type::*)(
            const Variable<VectorType>&, const double, const double,
            const int) const) &Type::ValueAt<VectorType>,
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

template <int TDimension>
void RegisterCurveBase(
    pybind11::module& m,
    const std::string& name)
{
    using namespace pybind11::literals;

    using Type = CurveBase<TDimension>;
    using Holder = shared_ptr<Type>;

    pybind11::class_<Type, Holder>(m, name.c_str())
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
    using namespace pybind11::literals;

    using GeometryType = CurveGeometry<TDimension>;

    using Type = Curve<TDimension>;
    using Holder = shared_ptr<Type>;
    using Base = CurveBase<TDimension>;

    pybind11::class_<Type, Base, Holder>(m, name.c_str())
        .def(pybind11::init<shared_ptr<GeometryType>,
            ANurbs::Interval<double>>(),
            "CurveGeometry"_a,
            "Domain"_a)
    ;
}

template <int TDimension>
void RegisterSurfaceBase(
    pybind11::module& m,
    const std::string& name)
{
    using namespace pybind11::literals;

    using Type = SurfaceBase<TDimension>;
    using Holder = shared_ptr<Type>;

    pybind11::class_<Type, Holder>(m, name.c_str())
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
    using namespace pybind11::literals;

    using GeometryType = SurfaceGeometry<TDimension>;

    using Type = Surface<TDimension>;
    using Holder = shared_ptr<Type>;
    using Base = SurfaceBase<TDimension>;

    pybind11::class_<Type, Base, Holder>(m, name.c_str())
        .def(pybind11::init<shared_ptr<GeometryType>,
            ANurbs::Interval<double>, ANurbs::Interval<double>>(),
            "SurfaceGeometry"_a,
            "DomainU"_a,
            "DomainV"_a)
    ;
}

template <int TDimension>
void RegisterCurveOnSurface(
    pybind11::module& m,
    const std::string& name)
{
    using namespace pybind11::literals;

    using Type = CurveOnSurface<TDimension>;
    using Holder = shared_ptr<Type>;
    using Base = CurveBase<TDimension>;

    pybind11::class_<Type, Base, Holder>(m, name.c_str())
        .def(pybind11::init<shared_ptr<CurveGeometryBase<2>>,
            shared_ptr<SurfaceGeometryBase<TDimension>>,
            ANurbs::Interval<double>>(),
            "CurveGeometry"_a,
            "SurfaceGeometry"_a,
            "Domain"_a)
    ;
}

template <int TDimension>
void RegisterPointOnCurveProjection(
    pybind11::module& m,
    const std::string& name)
{
    using namespace pybind11::literals;

    using VectorType = Kratos::array_1d<double, TDimension>;
    using CurveBaseType = CurveBase<TDimension>;

    using Type = ANurbs::PointOnCurveProjection<VectorType>;
    using Holder = shared_ptr<Type>;

    pybind11::class_<Type, Holder>(m, name.c_str())
        .def(pybind11::init<shared_ptr<CurveBaseType>, double>(),
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
    using namespace pybind11::literals;

    using VectorType = Kratos::array_1d<double, TDimension>;

    using Type = ANurbs::CurveTessellation<VectorType>;
    using Holder = shared_ptr<Type>;

    pybind11::class_<Type, Holder>(m, name.c_str())
        .def(pybind11::init<>())
        .def("Compute", &Type::Compute,
            "Curve"_a,
            "Tolerance"_a)
        .def_property_readonly("NumberOfPoints", &Type::NbPoints)
        .def("Parameter", &Type::Parameter,
            "index"_a)
        .def("Point", &Type::Point,
            "index"_a)
    ;
}

void RegisterIntegrationPoint1(
    pybind11::module& m,
    const std::string& name)
{
    using namespace pybind11::literals;

    using Type = ANurbs::IntegrationPoint1<double>;

    pybind11::class_<Type>(m, name.c_str())
        .def("__iter__",
            [](const Type &self) {
                return pybind11::make_iterator(&self.t, &self.t + 2);
            }, pybind11::keep_alive<0, 1>())
        .def_readwrite("t", &Type::t)
        .def_readwrite("weight", &Type::weight)
    ;
}

void RegisterIntegrationPoint2(
    pybind11::module& m,
    const std::string& name)
{
    using namespace pybind11::literals;

    using Type = ANurbs::IntegrationPoint2<double>;

    pybind11::class_<Type>(m, name.c_str())
        .def("__iter__",
            [](const Type &self) {
                return pybind11::make_iterator(&self.u, &self.u + 3);
            }, pybind11::keep_alive<0, 1>())
        .def_readwrite("u", &Type::u)
        .def_readwrite("v", &Type::v)
        .def_readwrite("weight", &Type::weight)
    ;
}

void RegisterIntegrationPoints(
    pybind11::module& m,
    const std::string& name)
{
    using namespace pybind11::literals;

    using Type = ANurbs::IntegrationPoints<double>;

    pybind11::class_<Type>(m, name.c_str())
        .def_static("Points1D", &Type::Points1,
            "Degree"_a,
            "Domain"_a)
        .def_static("Points2D", (std::vector<ANurbs::IntegrationPoint2<double>>
            (*)(const size_t, const size_t, const ANurbs::Interval<double>&,
            const ANurbs::Interval<double>&)) &Type::Points2,
            "DegreeU"_a,
            "DegreeV"_a,
            "DomainU"_a,
            "DomainV"_a)
        .def_static("Points2D", (std::vector<ANurbs::IntegrationPoint2<double>>
            (*)(const size_t, const size_t,
            const std::vector<ANurbs::Interval<double>>&,
            const std::vector<ANurbs::Interval<double>>&)) &Type::Points2,
            "DegreeU"_a,
            "DegreeV"_a,
            "DomainsU"_a,
            "DomainsV"_a)
    ;
}

void AddCustomUtilitiesToPython(
    pybind11::module& m)
{
    using namespace pybind11::literals;

    RegisterInterval(m, "Interval");

    RegisterCurveShapeEvaluator(m, "CurveShapeEvaluator");
    RegisterSurfaceShapeEvaluator(m, "SurfaceShapeEvaluator");

    RegisterPoint1D(m, "Point1D");
    RegisterPoint2D(m, "Point2D");

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

    RegisterNodeCurveGeometry(m, "NodeCurveGeometry3D");
    RegisterNodeSurfaceGeometry(m, "NodeSurfaceGeometry3D");

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

    RegisterIntegrationPoint1(m, "IntegrationPoint1");
    RegisterIntegrationPoint2(m, "IntegrationPoint2");
    RegisterIntegrationPoints(m, "IntegrationPoints");

    py::class_<BrepJsonIO, typename BrepJsonIO::Pointer>(m, "BrepJsonIO")
        .def(py::init<>())
        ;

    py::class_<NurbsBrepModeler, typename NurbsBrepModeler::Pointer>(m, "NurbsBrepModeler")
        .def(py::init<ModelPart&>())
        .def("ImportGeometry", &NurbsBrepModeler::ImportGeometry)
        .def("ImportModelPart", &NurbsBrepModeler::ImportModelPart)
        ;
}

} // namespace Python
} // Namespace Kratos
