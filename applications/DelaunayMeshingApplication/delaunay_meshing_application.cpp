//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

// System includes

// External includes

// Project includes
#include "delaunay_meshing_application.h"
#include "delaunay_meshing_application_variables.h"


namespace Kratos {

KratosDelaunayMeshingApplication::KratosDelaunayMeshingApplication():
    KratosApplication("DelaunayMeshingApplication")
{}

void KratosDelaunayMeshingApplication::Register() {
  // calling base class register to register Kratos components
  KratosApplication::Register();
  std::cout << "Initializing KratosDelaunayMeshingApplication... " << std::endl;


}
}  // namespace Kratos.
