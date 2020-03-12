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
#include "includes/variables.h"
//#include "includes/constitutive_law.h"

#include "geometries/quadrilateral_2d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"

#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

namespace Kratos {

  typedef Node<3> NodeType;

KratosMORApplication::KratosMORApplication():
    KratosApplication("MORApplication"),
      mAcousticElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mAcousticElement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mAcousticElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mAcousticElement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),

      /* CONDITIONS */
      // Adding point load conditions
      mPointDisplacementCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<NodeType >(Condition::GeometryType::PointsArrayType(1))))
      
    {}

void KratosMORApplication::Register()
{
     // calling base class register to register Kratos components
     KratosApplication::Register();
     KRATOS_INFO("") << "Initializing KratosMORApplication..." << std::endl;

  KRATOS_REGISTER_VARIABLE( FREQUENCY )
  KRATOS_REGISTER_VARIABLE( ACOUSTIC_PRESSURE )
  KRATOS_REGISTER_VARIABLE( PRESSURE_GRADIENT )
  KRATOS_REGISTER_VARIABLE( ACOUSTIC_VELOCITY )
  KRATOS_REGISTER_VARIABLE( ACOUSTIC_DISPLACEMENT )

  KRATOS_REGISTER_ELEMENT("AcousticElement2D3N", mAcousticElement2D3N)
  KRATOS_REGISTER_ELEMENT("AcousticElement2D4N", mAcousticElement2D4N)
  KRATOS_REGISTER_ELEMENT("AcousticElement3D4N", mAcousticElement3D4N)
  KRATOS_REGISTER_ELEMENT("AcousticElement3D8N", mAcousticElement3D8N)
  

// Register linear elastics laws
  KRATOS_REGISTER_CONSTITUTIVE_LAW("AcousticMaterial", mAcousticMaterial);


  // Register the conditions
  // Point loads
  KRATOS_REGISTER_CONDITION("PointDisplacementCondition3D1N", mPointDisplacementCondition3D1N)

}
}  // namespace Kratos.
