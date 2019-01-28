//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"

// Utilities
#include "custom_utilities/energy_utilities.h"

namespace Kratos
{

namespace Python
{

void  AddCustomUtilitiesToPython(pybind11::module& m)
{

  namespace py = pybind11;

  py::class_<EnergyUtilities>(m,"EnergyUtilities")
      .def(py::init<>())
      .def("GetTotalKinematicEnergy",&EnergyUtilities::GetTotalKinematicEnergy)
      .def("CalculateNodalMass",&EnergyUtilities::CalculateNodalMass)
      .def("GetTotalStrainEnergy",&EnergyUtilities::GetTotalStrainEnergy)
      .def("GetGravitationalEnergy",&EnergyUtilities::GetGravitationalEnergy)
      .def("GetExternallyAppliedEnergy",&EnergyUtilities::GetExternallyAppliedEnergy)
      ;

}

}  // namespace Python.

} // Namespace Kratos
