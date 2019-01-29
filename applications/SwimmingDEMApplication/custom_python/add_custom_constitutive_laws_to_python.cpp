//
// Author: Guillermo Casas (gcasas@cimne.upc.edu)
//

// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"

#include "../custom_constitutive/hydrodynamic_interaction_law.h"


namespace Kratos {
namespace Python {

namespace py = pybind11;

void AddCustomConstitutiveLawsToPython(pybind11::module& m) {

    py::class_<HydrodynamicInteractionLaw, HydrodynamicInteractionLaw::Pointer>(m, "HydrodynamicInteractionLaw")
        .def(py::init<>())
        .def("Clone", &HydrodynamicInteractionLaw::Clone)
        .def("SetHydrodynamicInteractionLawInProperties", &HydrodynamicInteractionLaw::SetHydrodynamicInteractionLawInProperties)
        .def("GetTypeOfLaw", &HydrodynamicInteractionLaw::GetTypeOfLaw)
        ;
  }

} // namespace Python.
} // namespace Kratos.
