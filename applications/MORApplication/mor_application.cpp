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

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/point_3d.h"


namespace Kratos {

  typedef Node<3> NodeType;

KratosMORApplication::KratosMORApplication():
    KratosApplication("MORApplication"),
      mAcousticElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mAcousticElement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      
      // conditions
      mDisplacementOutputCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<NodeType >(Condition::GeometryType::PointsArrayType(1))))
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
  KRATOS_REGISTER_VARIABLE( REAL_PRESSURE )
  KRATOS_REGISTER_VARIABLE( IMAG_PRESSURE )

  KRATOS_REGISTER_ELEMENT("AcousticElement3D4N", mAcousticElement3D4N)
  KRATOS_REGISTER_ELEMENT("AcousticElement3D8N", mAcousticElement3D8N)


  KRATOS_REGISTER_CONDITION("DisplacementOutputCondition3D1N", mDisplacementOutputCondition3D1N)


}
}  // namespace Kratos.
