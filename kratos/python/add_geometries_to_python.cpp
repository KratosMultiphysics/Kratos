//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
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
#include "geometries/pyramid_3d_5.h"
#include "geometries/pyramid_3d_13.h"
#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"
// Nurbs Geometries
#include "geometries/nurbs_volume_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_curve_geometry.h"
#include "geometries/surface_in_nurbs_volume_geometry.h"

namespace Kratos::Python
{
    using IndexType = std::size_t;
    using SizeType = std::size_t;
    using NodeType = Node;
    using NodeContainerType = PointerVector<NodeType>;
    using GeometryType = Geometry<NodeType>;
    using PointsArrayType = typename GeometryType::PointsArrayType;
    using IntegrationPointsArrayType = typename GeometryType::IntegrationPointsArrayType;
    using GeometriesArrayType = typename GeometryType::GeometriesArrayType;
    using CoordinatesArrayType = typename Point::CoordinatesArrayType;

    const PointerVector< Node >& ConstGetPoints( GeometryType& geom ) { return geom.Points(); }
    PointerVector< Node >& GetPoints( GeometryType& geom ) { return geom.Points(); }

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

    ///@}
    ///@name Set and Calculate
    ///@{

    template< class TDataType >
    void Assign(
        GeometryType& dummy, const Variable<TDataType>& rVariable, TDataType Value)
    {
        dummy.Assign(rVariable, Value);
    }

    template< class TDataType >
    TDataType Calculate(
        GeometryType& dummy, const Variable<TDataType>& rVariable)
    {
        TDataType Output;
        dummy.Calculate(rVariable, Output);
        return Output;
    }

    ///@}

void  AddGeometriesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef Node NodeType;
    typedef NodeType::Pointer pNodeType;
    typedef Geometry<NodeType > GeometryType;

    py::class_<GeometryType, GeometryType::Pointer >(m, "Geometry")
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
    .def("WorkingSpaceDimension", &GeometryType::WorkingSpaceDimension)
    .def("LocalSpaceDimension", &GeometryType::LocalSpaceDimension)
    .def("DomainSize", &GeometryType::DomainSize)
    .def("EdgesNumber", &GeometryType::EdgesNumber)
    .def("PointsNumber", &GeometryType::PointsNumber)
    .def("PointsNumberInDirection", &GeometryType::PointsNumberInDirection)
    .def("PolynomialDegree", &GeometryType::PolynomialDegree)
    // Geometry data
    .def("GetDefaultIntegrationMethod", &GeometryType::GetDefaultIntegrationMethod)
    .def("GetGeometryFamily", &GeometryType::GetGeometryFamily)
    .def("GetGeometryType", &GeometryType::GetGeometryType)
    // Geometry Parts
    .def("GetGeometryPart", [](GeometryType& self, IndexType Index)
        { return(self.GetGeometryPart(Index)); })
    .def_property_readonly_static("BACKGROUND_GEOMETRY_INDEX", [](py::object)
        { return GeometryType::BACKGROUND_GEOMETRY_INDEX; })
    // Integration
    .def("IntegrationPointsNumber", [](GeometryType& self)
        { return(self.IntegrationPointsNumber()); })
    // Quadrature points
    .def("CreateQuadraturePointGeometries", [](GeometryType& self,
        GeometriesArrayType& rResultGeometries, IndexType NumberOfShapeFunctionDerivatives)
        {
            IntegrationInfo integration_info = self.GetDefaultIntegrationInfo();
            return(self.CreateQuadraturePointGeometries(rResultGeometries, NumberOfShapeFunctionDerivatives, integration_info));
        })
    .def("CreateQuadraturePointGeometries", [](GeometryType& self,
        GeometriesArrayType& rResultGeometries, IndexType NumberOfShapeFunctionDerivatives, std::vector<std::array<double,4>>& rIntegrationPoints)
        {
            IntegrationPointsArrayType integration_points(rIntegrationPoints.size());
            for( IndexType i = 0; i < rIntegrationPoints.size(); ++i){
                IntegrationPoint<3> point_tmp(rIntegrationPoints[i][0],rIntegrationPoints[i][1],rIntegrationPoints[i][2],rIntegrationPoints[i][3]);
                integration_points[i] = point_tmp;
            }
            IntegrationInfo integration_info = self.GetDefaultIntegrationInfo();
            return(self.CreateQuadraturePointGeometries(rResultGeometries, NumberOfShapeFunctionDerivatives, integration_points, integration_info)); })
    // Normal
    .def("Normal", [](GeometryType& self)
        { const auto& r_type = self.GetGeometryType();
          KRATOS_WARNING_IF("Geometry", !(r_type == GeometryData::KratosGeometryType::Kratos_Line2D2 || r_type == GeometryData::KratosGeometryType::Kratos_Triangle3D3))
              << "WARNING:: Your geometry is not linear, please provide local coordinates or a GP index. Normal will be computed at zero local coordinates." << std::endl;
          CoordinatesArrayType LocalCoords;
          LocalCoords.clear();
          return(self.Normal(LocalCoords)); })
    .def("Normal", [](GeometryType& self, CoordinatesArrayType& LocalCoords)
        { return(self.Normal(LocalCoords)); })
    .def("Normal", [](GeometryType& self, IndexType IntegrationPointIndex)
        { return(self.Normal(IntegrationPointIndex)); })
    .def("UnitNormal", [](GeometryType& self)
        { const auto& r_type = self.GetGeometryType();
          KRATOS_WARNING_IF("Geometry", !(r_type == GeometryData::KratosGeometryType::Kratos_Line2D2 || r_type == GeometryData::KratosGeometryType::Kratos_Triangle3D3))
              << "Your geometry is not linear, please provide local coordinates or a GP index. Normal will be computed at zero local coordinates." << std::endl;
          CoordinatesArrayType LocalCoords;
          LocalCoords.clear();
          return(self.UnitNormal(LocalCoords)); })
    .def("UnitNormal", [](GeometryType& self, CoordinatesArrayType& LocalCoords)
        { return(self.UnitNormal(LocalCoords)); })
    .def("UnitNormal", [](GeometryType& self, IndexType IntegrationPointIndex)
        { return(self.UnitNormal(IntegrationPointIndex)); })
     // Jacobian
    .def("Jacobian", [](GeometryType& self, IndexType IntegrationPointIndex)
        { Matrix results; return(self.Jacobian(results, IntegrationPointIndex)); })
    .def("DeterminantOfJacobian", [](GeometryType& self)
        { Vector results; return(self.DeterminantOfJacobian(results)); })
    .def("DeterminantOfJacobian", [](GeometryType& self, IndexType IntegrationPointIndex)
        { return(self.DeterminantOfJacobian(IntegrationPointIndex)); })
    // ShapeFunctionsValues
    .def("ShapeFunctionsValues", [](GeometryType& self)
        { return(self.ShapeFunctionsValues()); })
    .def("ShapeFunctionDerivatives", [](GeometryType& self, IndexType DerivativeOrderIndex,
        IndexType IntegrationPointIndex)
        { return(self.ShapeFunctionDerivatives(DerivativeOrderIndex, IntegrationPointIndex, self.GetDefaultIntegrationMethod())); })
    // Mapping
    .def("GlobalCoordinates", [](GeometryType& self, CoordinatesArrayType& LocalCoordinates)
        {
        CoordinatesArrayType result = ZeroVector( 3 );
        return(self.GlobalCoordinates(result, LocalCoordinates)); })
    // Geometrical
    .def("Center",&GeometryType::Center)
    .def("Length",&GeometryType::Length)
    .def("Area",&GeometryType::Area)
    .def("Volume",&GeometryType::Volume)
    // Assign
    .def("Assign", Assign<bool>)
    .def("Assign", Assign<int>)
    .def("Assign", Assign<double>)
    .def("Assign", Assign<array_1d<double, 2>>)
    .def("Assign", Assign<array_1d<double, 3>>)
    .def("Assign", Assign<array_1d<double, 6>>)
    .def("Assign", Assign<Vector>)
    .def("Assign", Assign<Matrix>)
    // Calculate
    .def("Calculate", Calculate<bool>)
    .def("Calculate", Calculate<int>)
    .def("Calculate", Calculate<double>)
    .def("Calculate", Calculate<array_1d<double, 2>>)
    .def("Calculate", Calculate<array_1d<double, 3>>)
    .def("Calculate", Calculate<array_1d<double, 6>>)
    .def("Calculate", Calculate<Vector>)
    .def("Calculate", Calculate<Matrix>)
    // Info
    .def("Info",&GeometryType::Info)
    // Print
    .def("__str__", PrintObject<GeometryType>)
    // Access to nodes
    .def("__getitem__", [](GeometryType& self, unsigned int i){return self(i);} )
    .def("__iter__",    [](GeometryType& self){return py::make_iterator(self.begin(), self.end());},  py::keep_alive<0,1>())
    .def("__len__",     [](GeometryType& self){return self.PointsNumber();} )
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
    py::class_<Pyramid3D5<NodeType>, Pyramid3D5<NodeType>::Pointer,  GeometryType  >(m,"Pyramid3D5").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Pyramid3D13<NodeType>, Pyramid3D13<NodeType>::Pointer,  GeometryType  >(m,"Pyramid3D13").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Prism3D6<NodeType>, Prism3D6<NodeType>::Pointer,  GeometryType  >(m,"Prism3D6").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Prism3D15<NodeType>, Prism3D15<NodeType>::Pointer,  GeometryType  >(m,"Prism3D15").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Hexahedra3D8<NodeType>, Hexahedra3D8<NodeType>::Pointer,  GeometryType  >(m,"Hexahedra3D8").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Hexahedra3D20<NodeType>, Hexahedra3D20<NodeType>::Pointer,  GeometryType  >(m,"Hexahedra3D20").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Hexahedra3D27<NodeType>, Hexahedra3D27<NodeType>::Pointer,  GeometryType  >(m,"Hexahedra3D27").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;

    // Adding PointsArrayType to interface
    PointerVectorPythonInterface<PointsArrayType>().CreateInterface(m,"NodesVector");

    // Adding GeometriesArrayType to interface
    PointerVectorPythonInterface<GeometriesArrayType>().CreateInterface(m, "GeometriesVector");

    /// Nurbs Geometries
    // NurbsVolumeGeometry
    py::class_<NurbsVolumeGeometry<NodeContainerType>, NurbsVolumeGeometry<NodeContainerType>::Pointer, GeometryType >(m, "NurbsVolumeGeometry")
        .def(py::init<const PointsArrayType&, const SizeType, const SizeType, const SizeType, const Vector&, const Vector&, const Vector&>())
        .def("PolynomialDegreeU", &NurbsVolumeGeometry<NodeContainerType>::PolynomialDegreeU)
        .def("PolynomialDegreeV", &NurbsVolumeGeometry<NodeContainerType>::PolynomialDegreeV)
        .def("PolynomialDegreeW", &NurbsVolumeGeometry<NodeContainerType>::PolynomialDegreeW)
        .def("KnotsU", &NurbsVolumeGeometry<NodeContainerType>::KnotsU)
        .def("KnotsV", &NurbsVolumeGeometry<NodeContainerType>::KnotsV)
        .def("KnotsW", &NurbsVolumeGeometry<NodeContainerType>::KnotsW)
        .def("NumberOfKnotsU", &NurbsVolumeGeometry<NodeContainerType>::NumberOfKnotsU)
        .def("NumberOfKnotsV", &NurbsVolumeGeometry<NodeContainerType>::NumberOfKnotsV)
        .def("NumberOfKnotsW", &NurbsVolumeGeometry<NodeContainerType>::NumberOfKnotsW)
        .def("IsRational", &NurbsVolumeGeometry<NodeContainerType>::IsRational)
        .def("NumberOfControlPointsU", &NurbsVolumeGeometry<NodeContainerType>::NumberOfControlPointsU)
        .def("NumberOfControlPointsV", &NurbsVolumeGeometry<NodeContainerType>::NumberOfControlPointsV)
        .def("NumberOfControlPointsW", &NurbsVolumeGeometry<NodeContainerType>::NumberOfControlPointsW)
        ;

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

    py::class_<SurfaceInNurbsVolumeGeometry<3, NodeContainerType>, SurfaceInNurbsVolumeGeometry<3, NodeContainerType>::Pointer, GeometryType>(m, "SurfaceInNurbsVolumeGeometry")
        .def(py::init<NurbsVolumeGeometry<NodeContainerType>::Pointer, GeometryType::Pointer>())
        ;

}

}  // namespace Kratos::Python.
