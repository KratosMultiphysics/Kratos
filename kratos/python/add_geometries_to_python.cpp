//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "geometries/point.h"
#include "includes/node.h"
#include "python/containers_interface.h"
#include "python/add_geometries_to_python.h"
#include "geometries/geometry.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/prism_3d_6.h"
// #include "geometries/prism_3d_15.h"
#include "geometries/hexahedra_3d_8.h"
// #include "geometries/hexahedra_3d_20.h"
// #include "geometries/hexahedra_3d_27.h"
// Nurbs Geometries
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_curve_geometry.h"

namespace Kratos
{

namespace Python
{
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef Node<3> NodeType;
    typedef PointerVector<NodeType> NodeContainerType;
    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::PointsArrayType PointsArrayType;
    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef typename Point::CoordinatesArrayType CoordinatesArrayType;

    const PointerVector< Node<3> >& ConstGetPoints( GeometryType& geom ) { return geom.Points(); }
    PointerVector< Node<3> >& GetPoints( GeometryType& geom ) { return geom.Points(); }

    // Id utilities
    void SetId1(
        GeometryType& dummy, IndexType geometry_id)
    {
        return(dummy.SetId(geometry_id));
    }

    void SetId2(
        GeometryType& dummy, const std::string& geometry_name)
    {
        return(dummy.SetId(geometry_name));
    }

    bool IsIdGeneratedFromString1(GeometryType& dummy)
    {
        return(dummy.IsIdGeneratedFromString());
    }

    bool IsIdSelfAssigned1(GeometryType& dummy)
    {
        return(dummy.IsIdSelfAssigned());
    }

    array_1d<double,3> GetNormal(
        GeometryType& dummy,
        CoordinatesArrayType& LocalCoords
        )
    {
        return( dummy.Normal(LocalCoords) );
    }

    array_1d<double,3> FastGetNormal(GeometryType& dummy)
    {
        CoordinatesArrayType LocalCoords;
        LocalCoords.clear();
        return( dummy.Normal(LocalCoords) );
    }

    array_1d<double,3> GetUnitNormal(
        GeometryType& dummy,
        CoordinatesArrayType& LocalCoords
        )
    {
        return( dummy.UnitNormal(LocalCoords) );
    }

    array_1d<double,3> FastGetUnitNormal(GeometryType& dummy)
    {
        CoordinatesArrayType LocalCoords;
        LocalCoords.clear();
        return( dummy.UnitNormal(LocalCoords) );
    }

void  AddGeometriesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef Node<3> NodeType;
    typedef NodeType::Pointer pNodeType;
    typedef Geometry<NodeType > GeometryType;

    py::class_<GeometryType, GeometryType::Pointer >(m,"Geometry")
    .def(py::init<>())
    .def(py::init< IndexType >())
    .def(py::init< std::string >())
    .def(py::init< GeometryType::PointsArrayType& >())
    .def(py::init< IndexType, GeometryType::PointsArrayType& >())
    .def(py::init< std::string, GeometryType::PointsArrayType& >())
    // Id functions
    .def_property("Id", &GeometryType::Id, SetId1)
    .def("SetId", SetId1)
    .def("SetId", SetId2)
    .def("IsIdGeneratedFromString", IsIdGeneratedFromString1)
    .def("IsIdSelfAssigned", IsIdSelfAssigned1)
    .def_static("GenerateId", &GeometryType::GenerateId)
    // Dimension access
    .def("WorkingSpaceDimension",&GeometryType::WorkingSpaceDimension)
    .def("LocalSpaceDimension",&GeometryType::LocalSpaceDimension)
    .def("Dimension", &GeometryType::Dimension)
    .def("DomainSize",&GeometryType::DomainSize)
    .def("PointsNumber",&GeometryType::PointsNumber)
    .def("Normal",GetNormal)
    .def("Normal",FastGetNormal)
    .def("UnitNormal",GetUnitNormal)
    .def("UnitNormal",FastGetUnitNormal)
    .def("Center",&GeometryType::Center)
    .def("Length",&GeometryType::Length)
    .def("Area",&GeometryType::Area)
    .def("Volume",&GeometryType::Volume)
    .def("__str__", PrintObject<GeometryType>)
    .def("__getitem__", [](GeometryType& self, unsigned int i){return self(i);} )
    .def("__iter__",    [](GeometryType& self){return py::make_iterator(self.begin(), self.end());},  py::keep_alive<0,1>())
    .def("__len__",     [](GeometryType& self){return self.PointsNumber();} )
//     .def("Points", &GeometryType::ConstGetPoints)
//     .def("Points", &GeometryType::GetPoints)
    ;

    // 2D
    py::class_<Line2D2<NodeType>, Line2D2<NodeType>::Pointer,  GeometryType  >(m,"Line2D2").def(py::init<pNodeType, pNodeType>())
    ;
    py::class_<Line2D3<NodeType>, Line2D3<NodeType>::Pointer,  GeometryType  >(m,"Line2D3").def(py::init<pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Triangle2D3<NodeType>, Triangle2D3<NodeType>::Pointer,  GeometryType  >(m,"Triangle2D3").def(py::init<pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Triangle2D6<NodeType>, Triangle2D6<NodeType>::Pointer,  GeometryType  >(m,"Triangle2D6").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Quadrilateral2D4<NodeType>, Quadrilateral2D4<NodeType>::Pointer,  GeometryType  >(m,"Quadrilateral2D4").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Quadrilateral2D8<NodeType>, Quadrilateral2D8<NodeType>::Pointer,  GeometryType  >(m,"Quadrilateral2D8").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Quadrilateral2D9<NodeType>, Quadrilateral2D9<NodeType>::Pointer,  GeometryType  >(m,"Quadrilateral2D9").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;

    // 3D
    py::class_<Line3D2<NodeType>, Line3D2<NodeType>::Pointer,  GeometryType  >(m,"Line3D2").def(py::init<pNodeType, pNodeType>())
    ;
    py::class_<Line3D3<NodeType>, Line3D3<NodeType>::Pointer,  GeometryType  >(m,"Line3D3").def(py::init<pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Triangle3D3<NodeType>, Triangle3D3<NodeType>::Pointer,  GeometryType  >(m,"Triangle3D3").def(py::init<pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Triangle3D6<NodeType>, Triangle3D6<NodeType>::Pointer,  GeometryType  >(m,"Triangle3D6").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Quadrilateral3D4<NodeType>, Quadrilateral3D4<NodeType>::Pointer,  GeometryType  >(m,"Quadrilateral3D4").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Quadrilateral3D8<NodeType>, Quadrilateral3D8<NodeType>::Pointer,  GeometryType  >(m,"Quadrilateral3D8").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Quadrilateral3D9<NodeType>, Quadrilateral3D9<NodeType>::Pointer,  GeometryType  >(m,"Quadrilateral3D9").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Tetrahedra3D4<NodeType>, Tetrahedra3D4<NodeType>::Pointer,  GeometryType  >(m,"Tetrahedra3D4").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Tetrahedra3D10<NodeType>, Tetrahedra3D10<NodeType>::Pointer,  GeometryType  >(m,"Tetrahedra3D10").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Prism3D6<NodeType>, Prism3D6<NodeType>::Pointer,  GeometryType  >(m,"Prism3D6").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
//     py::class_<Prism3D15<NodeType>, Prism3D15<NodeType>::Pointer,  GeometryType  >(m,"Prism3D15").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
//     ;
    py::class_<Hexahedra3D8<NodeType>, Hexahedra3D8<NodeType>::Pointer,  GeometryType  >(m,"Hexahedra3D8").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
//     py::class_<Hexahedra3D20<NodeType>, Hexahedra3D20<NodeType>::Pointer,  GeometryType  >(m,"Hexahedra3D20").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
//     ;
//     py::class_<Hexahedra3D27<NodeType>, Hexahedra3D27<NodeType>::Pointer,  GeometryType  >(m,"Hexahedra3D27").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
//     ;

    // Adding PointsArrayType to interface
    PointerVectorPythonInterface<PointsArrayType>().CreateInterface(m,"NodesVector");

    /// Nurbs Geometries
    // NurbsSurfaceGeometry3D
    py::class_<NurbsSurfaceGeometry<3, NodeContainerType>, NurbsSurfaceGeometry<3, NodeContainerType>::Pointer, GeometryType >(m, "NurbsSurfaceGeometry3D")
        .def(py::init<const PointsArrayType&, const SizeType, const SizeType, const Vector&, const Vector&>())
        .def(py::init<const PointsArrayType&, const SizeType, const SizeType, const Vector&, const Vector&, const Vector&>())
        .def("PolynomialDegreeU", &NurbsSurfaceGeometry<3, NodeContainerType>::PolynomialDegreeU)
        .def("PolynomialDegreeV", &NurbsSurfaceGeometry<3, NodeContainerType>::PolynomialDegreeV)
        .def("KnotsU", &NurbsSurfaceGeometry<3, NodeContainerType>::KnotsU)
        .def("KnotsV", &NurbsSurfaceGeometry<3, NodeContainerType>::KnotsV)
        .def("NumberOfKnotsU", &NurbsSurfaceGeometry<3, NodeContainerType>::NumberOfKnotsU)
        .def("NumberOfKnotsV", &NurbsSurfaceGeometry<3, NodeContainerType>::NumberOfKnotsV)
        .def("IsRational", &NurbsSurfaceGeometry<3, NodeContainerType>::IsRational)
        .def("Weights", &NurbsSurfaceGeometry<3, NodeContainerType>::Weights)
        .def("NumberOfControlPointsU", &NurbsSurfaceGeometry<3, NodeContainerType>::NumberOfControlPointsU)
        .def("NumberOfControlPointsV", &NurbsSurfaceGeometry<3, NodeContainerType>::NumberOfControlPointsV)
        ;

    // NurbsCurveGeometry3D
    py::class_<NurbsCurveGeometry<3, NodeContainerType>, NurbsCurveGeometry<3, NodeContainerType>::Pointer, GeometryType >(m, "NurbsCurveGeometry3D")
        .def(py::init<const PointsArrayType&, const SizeType, const Vector&>())
        .def(py::init<const PointsArrayType&, const SizeType, const Vector&, const Vector&>())
        .def("PolynomialDegree", &NurbsCurveGeometry<3, NodeContainerType>::PolynomialDegree)
        .def("Knots", &NurbsCurveGeometry<3, NodeContainerType>::Knots)
        .def("NumberOfKnots", &NurbsCurveGeometry<3, NodeContainerType>::NumberOfKnots)
        .def("NumberOfControlPoints", &NurbsCurveGeometry<3, NodeContainerType>::NumberOfNonzeroControlPoints)
        .def("IsRational", &NurbsCurveGeometry<3, NodeContainerType>::IsRational)
        .def("Weights", &NurbsCurveGeometry<3, NodeContainerType>::Weights)
        ;

    // NurbsCurveGeometry2D
    py::class_<NurbsCurveGeometry<2, NodeContainerType>, NurbsCurveGeometry<2, NodeContainerType>::Pointer, GeometryType >(m, "NurbsCurveGeometry2D")
        .def(py::init<const PointsArrayType&, const SizeType, const Vector>())
        .def(py::init<const PointsArrayType&, const SizeType, const Vector, const Vector>())
        .def("PolynomialDegree", &NurbsCurveGeometry<2, NodeContainerType>::PolynomialDegree)
        .def("Knots", &NurbsCurveGeometry<2, NodeContainerType>::Knots)
        .def("NumberOfKnots", &NurbsCurveGeometry<2, NodeContainerType>::NumberOfKnots)
        .def("NumberOfControlPoints", &NurbsCurveGeometry<2, NodeContainerType>::NumberOfNonzeroControlPoints)
        .def("IsRational", &NurbsCurveGeometry<2, NodeContainerType>::IsRational)
        .def("Weights", &NurbsCurveGeometry<2, NodeContainerType>::Weights)
        ;

}

}  // namespace Python.

} // Namespace Kratos

