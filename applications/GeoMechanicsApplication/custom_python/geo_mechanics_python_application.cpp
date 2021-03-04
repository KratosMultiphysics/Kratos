// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

// System includes

#if defined(KRATOS_PYTHON)

// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"


// Application includes
#include "geo_mechanics_application.h"
#include "geo_mechanics_application_variables.h"


namespace Kratos {
namespace Python {

  using namespace pybind11;

  PYBIND11_MODULE(KratosGeoMechanicsApplication,m)
  {

  class_<KratosGeoMechanicsApplication,
         KratosGeoMechanicsApplication::Pointer,
         KratosApplication>(m, "KratosGeoMechanicsApplication")
        .def(init<>());



    //Registering variables in python


    /* Reset displacement "flag" needed for GeoMechanicalApplication*/



  }


} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
