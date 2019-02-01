//
// Author: Guillermo Casas (gcasas@cimne.upc.edu)
//

// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"

#include "../custom_constitutive/hydrodynamic_interaction_law.h"
#include "../custom_constitutive/power_law_hydrodynamic_interaction_law.h"
#include "../custom_constitutive/drag_laws/drag_law.h"
#include "../custom_constitutive/drag_laws/stokes_drag_law.h"
#include "../custom_constitutive/drag_laws/beetstra_drag_law.h"
#include "../custom_constitutive/drag_laws/schiller_and_naumann_drag_law.h"
#include "../custom_constitutive/drag_laws/haider_and_levenspiel_drag_law.h"
#include "../custom_constitutive/drag_laws/ganser_drag_law.h"
#include "../custom_constitutive/drag_laws/shah_drag_law.h"

namespace Kratos {
namespace Python {

namespace py = pybind11;

typedef DragLaw BaseDragLawType;

void AddCustomConstitutiveLawsToPython(pybind11::module& m) {

    py::class_<HydrodynamicInteractionLaw, HydrodynamicInteractionLaw::Pointer>(m, "HydrodynamicInteractionLaw")
        .def(py::init<>())
        .def(py::init<const DragLaw&>())
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

    py::class_<BeetstraDragLaw, BeetstraDragLaw::Pointer, BaseDragLawType>(m, "BeetstraDragLaw")
        .def(py::init<>())
        ;

    py::class_<SchillerAndNaumannDragLaw, SchillerAndNaumannDragLaw::Pointer, BaseDragLawType>(m, "SchillerAndNaumannDragLaw")
        .def(py::init<>())
        ;

    py::class_<HaiderAndLevenspielDragLaw, HaiderAndLevenspielDragLaw::Pointer, BaseDragLawType>(m, "HaiderAndLevenspielDragLaw")
        .def(py::init<>())
        ;

    py::class_<GanserDragLaw, GanserDragLaw::Pointer, BaseDragLawType>(m, "GanserDragLaw")
        .def(py::init<>())
        ;

    py::class_<ShahDragLaw, ShahDragLaw::Pointer, BaseDragLawType>(m, "ShahDragLaw")
        .def(py::init<>())
        ;
  }

} // namespace Python.
} // namespace Kratos.
