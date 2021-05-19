//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Guillermo Casas, gcasas@cimne.upc.edu $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.3 $
//
//



// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "plasma_dynamics_application.h"
#include "includes/kernel.h"
#include "includes/kratos_flags.h"

// #include "plasma_dynamics_application_variables.h"

namespace Kratos
{

KratosPlasmaDynamicsApplication::KratosPlasmaDynamicsApplication():KratosApplication("PlasmaDynamicsApplication"){}

void KratosPlasmaDynamicsApplication::Register()
{
  KRATOS_INFO("") << "\n"
                  << "Importing PLASMA DYNAMICS APPLICATION... "<<std::endl;

  // std::cout << "Initializing KratosPlasmaDynamicsApplication... " << std::endl;

  /* Define In Global variables.cpp */

}

}  // namespace Kratos.


