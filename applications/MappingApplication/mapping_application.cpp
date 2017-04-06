//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"


// System includes


// External includes


// Project includes
#include "mapping_application.h"
#include "mapping_application_variables.h"

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/prism_3d_6.h"
#include "geometries/hexahedra_3d_8.h"


namespace Kratos {

KratosMappingApplication::KratosMappingApplication():
  mVolumeCondition3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> >( Element::GeometryType::PointsArrayType( 4  ) ) ) ),
  mVolumeCondition3D6N( 0, Element::GeometryType::Pointer( new Prism3D6<Node<3> >(Element::GeometryType::PointsArrayType(6)))),
  mVolumeCondition3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8))))
 { }

void KratosMappingApplication::Register() {
 	// calling base class register to register Kratos components
 	KratosApplication::Register();

  std::cout << "    KRATOS ______  ___                      _____  "                          << std::endl;
  std::cout << "           ___   |/  /_____ ___________________(_)_____________ _  "          << std::endl;
  std::cout << "           __  /|_/ /_  __ `/__  __ \\__  __ \\_  /__  __ \\_  __ `/  "       << std::endl;
  std::cout << "           _  /  / / / /_/ /__  /_/ /_  /_/ /  / _  / / /  /_/ /  "           << std::endl;
  std::cout << "           /_/  /_/  \\__,_/ _  .___/_  .___//_/  /_/ /_/_\\__, /  "          << std::endl;
  std::cout << "                            /_/     /_/                 /____/ Application"   << std::endl;

 	std::cout << "Initializing KratosMappingApplication... " << std::endl;
  
  // Needed to exchange Information abt the found neighbors (i.e. only for debugging)
  KRATOS_REGISTER_VARIABLE( NEIGHBOR_RANK )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NEIGHBOR_COORDINATES )

  // Needed for Volume Mapping
  KRATOS_REGISTER_CONDITION( "VolumeCondition3D4N", mVolumeCondition3D4N );
  KRATOS_REGISTER_CONDITION( "VolumeCondition3D6N", mVolumeCondition3D6N );
  KRATOS_REGISTER_CONDITION( "VolumeCondition3D8N", mVolumeCondition3D8N );

}
}  // namespace Kratos.