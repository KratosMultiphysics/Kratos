//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


// System includes


// External includes


// Project includes
#include "includes/define.h"

#include "mor_application.h"
#include "mor_application_variables.h"

#include "geometries/line_2d_2.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/point_3d.h"
#include "geometries/triangle_3d_3.h"


namespace Kratos {

  typedef Node<3> NodeType;

KratosMORApplication::KratosMORApplication():
    KratosApplication("MORApplication"),
      mAcousticElement2D2N(0, Element::GeometryType::Pointer(new Line2D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
      mAcousticElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mAcousticElement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mAcousticElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mAcousticElement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mAcousticPMLElement2D2N(0, Element::GeometryType::Pointer(new Line2D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
      mAcousticPMLElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mAcousticPMLElement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mAcousticPMLElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mAcousticPMLElement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),

      // conditions
      mAcousticLoadConcition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mAcousticLoadConcition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mAcousticLoadConcition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Condition::GeometryType::PointsArrayType(4)))),
      mAcousticRobinConcition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mAcousticRobinConcition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mAcousticRobinConcition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Condition::GeometryType::PointsArrayType(4)))),
      mAcousticStructureCouplingCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mAcousticStructureCouplingCondition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Condition::GeometryType::PointsArrayType(4)))),
      mAcousticStructureCouplingCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mAcousticStructureMappingCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mAcousticStructureMappingCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mAcousticStructureMappingCondition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Condition::GeometryType::PointsArrayType(4)))),
      mDisplacementOutputCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<NodeType >(Condition::GeometryType::PointsArrayType(1)))),
      mPressureOutputCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<NodeType >(Condition::GeometryType::PointsArrayType(1))))
    {}

void KratosMORApplication::Register()
{
     // calling base class register to register Kratos components
     KratosApplication::Register();
     KRATOS_INFO("") << "Initializing KratosMORApplication..." << std::endl;

  KRATOS_REGISTER_VARIABLE( FREQUENCY )
  KRATOS_REGISTER_VARIABLE( BUILD_LEVEL )
  KRATOS_REGISTER_VARIABLE( SCALAR_OUTPUT )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( COMPONENT_OUTPUT )



  // complex dof variables
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( REAL_DISPLACEMENT )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( IMAG_DISPLACEMENT )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( PML_IMAG_DISTANCE )

  KRATOS_REGISTER_VARIABLE( REAL_PRESSURE )
  KRATOS_REGISTER_VARIABLE( IMAG_PRESSURE )

  // elements
  KRATOS_REGISTER_ELEMENT("AcousticElement2D2N", mAcousticElement2D2N)
  KRATOS_REGISTER_ELEMENT("AcousticElement2D3N", mAcousticElement2D3N)
  KRATOS_REGISTER_ELEMENT("AcousticElement2D4N", mAcousticElement2D4N)
  KRATOS_REGISTER_ELEMENT("AcousticElement3D4N", mAcousticElement3D4N)
  KRATOS_REGISTER_ELEMENT("AcousticElement3D8N", mAcousticElement3D8N)
  KRATOS_REGISTER_ELEMENT("AcousticPMLElement2D2N", mAcousticPMLElement2D2N)
  KRATOS_REGISTER_ELEMENT("AcousticPMLElement2D3N", mAcousticPMLElement2D3N)
  KRATOS_REGISTER_ELEMENT("AcousticPMLElement2D4N", mAcousticPMLElement2D4N)
  KRATOS_REGISTER_ELEMENT("AcousticPMLElement3D4N", mAcousticPMLElement3D4N)
  KRATOS_REGISTER_ELEMENT("AcousticPMLElement3D8N", mAcousticPMLElement3D8N)

  // conditions
  KRATOS_REGISTER_CONDITION("AcousticLoadCondition2D2N", mAcousticLoadConcition2D2N)
  KRATOS_REGISTER_CONDITION("AcousticLoadCondition3D3N", mAcousticLoadConcition3D3N)
  KRATOS_REGISTER_CONDITION("AcousticLoadCondition3D4N", mAcousticLoadConcition3D4N)
  KRATOS_REGISTER_CONDITION("AcousticRobinCondition2D2N", mAcousticRobinConcition2D2N)
  KRATOS_REGISTER_CONDITION("AcousticRobinCondition3D3N", mAcousticRobinConcition3D3N)
  KRATOS_REGISTER_CONDITION("AcousticRobinCondition3D4N", mAcousticRobinConcition3D4N)
  KRATOS_REGISTER_CONDITION("AcousticStructureCouplingCondition2D2N", mAcousticStructureCouplingCondition2D2N)
  KRATOS_REGISTER_CONDITION("AcousticStructureCouplingCondition3D4N", mAcousticStructureCouplingCondition3D4N)
  KRATOS_REGISTER_CONDITION("AcousticStructureCouplingCondition3D3N", mAcousticStructureCouplingCondition3D3N)
  KRATOS_REGISTER_CONDITION("AcousticStructureMappingCondition2D2N", mAcousticStructureMappingCondition2D2N)
  KRATOS_REGISTER_CONDITION("AcousticStructureMappingCondition3D3N", mAcousticStructureMappingCondition3D3N)
  KRATOS_REGISTER_CONDITION("AcousticStructureMappingCondition3D4N", mAcousticStructureMappingCondition3D4N)

  KRATOS_REGISTER_CONDITION("DisplacementOutputCondition3D1N", mDisplacementOutputCondition3D1N)
  KRATOS_REGISTER_CONDITION("PressureOutputCondition3D1N", mPressureOutputCondition3D1N)

}
}  // namespace Kratos.
