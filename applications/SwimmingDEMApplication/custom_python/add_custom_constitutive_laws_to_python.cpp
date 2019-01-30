//
// Author: Guillermo Casas (gcasas@cimne.upc.edu)
//

// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"

#include "../custom_constitutive/hydrodynamic_interaction_law.h"
#include "../custom_constitutive/drag_law.h"


namespace Kratos {
namespace Python {

namespace py = pybind11;

typedef DragLaw BaseDragLawType;

void AddCustomConstitutiveLawsToPython(pybind11::module& m) {

    py::class_<HydrodynamicInteractionLaw, HydrodynamicInteractionLaw::Pointer>(m, "HydrodynamicInteractionLaw")
        .def(py::init<>())
        .def(py::init<DragLaw&>())
        .def("Clone", &HydrodynamicInteractionLaw::Clone)
        .def("SetHydrodynamicInteractionLawInProperties", &HydrodynamicInteractionLaw::SetHydrodynamicInteractionLawInProperties)
        .def("GetTypeOfLaw", &HydrodynamicInteractionLaw::GetTypeOfLaw)
        ;

    py::class_<DragLaw, DragLaw::Pointer>(m, "DragLaw")
        .def(py::init<>())
        .def("Clone", &DragLaw::Clone)
        .def("SetDragLawInProperties", &DragLaw::SetDragLawInProperties)
        .def("GetTypeOfLaw", &DragLaw::GetTypeOfLaw)
        ;

    py::class_<StokesDragLaw, StokesDragLaw::Pointer, BaseDragLawType>(m, "StokesDragLaw")
        .def(py::init<>())
        ;
  }

} // namespace Python.
} // namespace Kratos.
