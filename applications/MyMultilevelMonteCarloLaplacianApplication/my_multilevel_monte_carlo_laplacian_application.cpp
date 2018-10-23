//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand@{KRATOS_APP_AUTHOR}
//


// System includes


// External includes


// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/line_2d.h"
#include "geometries/point_2d.h"
#include "my_multilevel_monte_carlo_laplacian_application.h"
#include "my_multilevel_monte_carlo_laplacian_application_variables.h"


namespace Kratos {

KratosMyMultilevelMonteCarloLaplacianApplication::KratosMyMultilevelMonteCarloLaplacianApplication():
    KratosApplication("MyMultilevelMonteCarloLaplacianApplication"),
    mMyLaplacianElement( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>>( Element::GeometryType::PointsArrayType (3) ) ) )
    {}

void KratosMyMultilevelMonteCarloLaplacianApplication::Register() {
 	// calling base class register to register Kratos components
 	KratosApplication::Register();
 	std::cout << "Initializing KratosMyMultilevelMonteCarloLaplacianApplication... " << std::endl;

  KRATOS_REGISTER_VARIABLE( SOLUTION )
  KRATOS_REGISTER_VARIABLE( FORCING )
  KRATOS_REGISTER_ELEMENT( "MyLaplacianElement", mMyLaplacianElement )

}
}  // namespace Kratos.
