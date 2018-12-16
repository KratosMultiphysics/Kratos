//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//


// System includes


// External includes


// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/line_2d.h"
#include "geometries/point_2d.h"
#include "multilevel_monte_carlo_application.h"
#include "multilevel_monte_carlo_application_variables.h"


namespace Kratos {

KratosMultilevelMonteCarloApplication::KratosMultilevelMonteCarloApplication():
    KratosApplication("MultilevelMonteCarloApplication")
    {}

void KratosMultilevelMonteCarloApplication::Register() {
 	// calling base class register to register Kratos components
 	KratosApplication::Register();
 	std::cout << "Initializing KratosMultilevelMonteCarloApplication... " << std::endl;


}
}  // namespace Kratos.
