//
// Author: Guillermo Casas (gcasas@cimne.upc.edu)
//

// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"

// Hydrodynamic Laws
#include "../custom_constitutive/hydrodynamic_interaction_law.h"
#include "../custom_constitutive/power_law_hydrodynamic_interaction_law.h"

// Buoyancy Laws
#include "../custom_constitutive/buoyancy_laws/buoyancy_law.h"
#include "../custom_constitutive/buoyancy_laws/archimedes_buoyancy_law.h"

// Drag laws
#include "../custom_constitutive/drag_laws/drag_law.h"
#include "../custom_constitutive/drag_laws/stokes_drag_law.h"
#include "../custom_constitutive/drag_laws/beetstra_drag_law.h"
#include "../custom_constitutive/drag_laws/schiller_and_naumann_drag_law.h"
#include "../custom_constitutive/drag_laws/haider_and_levenspiel_drag_law.h"
#include "../custom_constitutive/drag_laws/ganser_drag_law.h"
#include "../custom_constitutive/drag_laws/shah_drag_law.h"
#include "../custom_constitutive/drag_laws/newton_drag_law.h"

// Inviscid force laws
#include "../custom_constitutive/inviscid_force_laws/inviscid_force_law.h"
#include "../custom_constitutive/inviscid_force_laws/auton_hunt_prudhomme_inviscid_force_law.h"

// history force laws
#include "../custom_constitutive/history_force_laws/history_force_law.h"
#include "../custom_constitutive/history_force_laws/boussinesq_basset_history_force_law.h"

namespace Kratos {
namespace Python {

namespace py = pybind11;

typedef BuoyancyLaw BaseBuoyancyLawType;
typedef DragLaw BaseDragLawType;
typedef InviscidForceLaw BaseInviscidForceLawType;
typedef HistoryForceLaw BaseHistoryForceLaw;

void AddCustomConstitutiveLawsToPython(pybind11::module& m) {

    py::class_<HydrodynamicInteractionLaw, HydrodynamicInteractionLaw::Pointer>(m, "HydrodynamicInteractionLaw")
        .def(py::init<>())
        .def(py::init<Properties::Pointer, Parameters&>())
        .def("Clone", &HydrodynamicInteractionLaw::Clone)
        .def("SetHydrodynamicInteractionLawInProperties", &HydrodynamicInteractionLaw::SetHydrodynamicInteractionLawInProperties)
        .def("GetTypeOfLaw", &HydrodynamicInteractionLaw::GetTypeOfLaw)
        .def("SetBuoyancyLaw", &HydrodynamicInteractionLaw::SetBuoyancyLaw)
        .def("SetDragLaw", &HydrodynamicInteractionLaw::SetDragLaw)
        .def("SetInviscidForceLaw", &HydrodynamicInteractionLaw::SetInviscidForceLaw)
        .def("SetHistoryForceLaw", &HydrodynamicInteractionLaw::SetHistoryForceLaw)
        ;

    // Buoyancy laws
    py::class_<BuoyancyLaw, BuoyancyLaw::Pointer>(m, "BuoyancyLaw")
        .def(py::init<Parameters&>())
        .def(py::init<>())
        .def("Clone", &BuoyancyLaw::Clone)
        .def("SetBuoyancyLawInProperties", &BuoyancyLaw::SetBuoyancyLawInProperties)
        .def("GetTypeOfLaw", &BuoyancyLaw::GetTypeOfLaw)
        ;

    py::class_<ArchimedesBuoyancyLaw, ArchimedesBuoyancyLaw::Pointer, BaseBuoyancyLawType>(m, "ArchimedesBuoyancyLaw")
        .def(py::init<>())
        .def(py::init<Parameters&>())
        ;

    // Drag laws

    py::class_<DragLaw, DragLaw::Pointer>(m, "DragLaw")
        .def(py::init<>())
        .def(py::init<Parameters&>())
        .def("Clone", &DragLaw::Clone)
        .def("SetDragLawInProperties", &DragLaw::SetDragLawInProperties)
        .def("GetTypeOfLaw", &DragLaw::GetTypeOfLaw)
        ;

    py::class_<StokesDragLaw, StokesDragLaw::Pointer, BaseDragLawType>(m, "StokesDragLaw")
        .def(py::init<>())
        .def(py::init<Parameters&>())
        ;

    py::class_<BeetstraDragLaw, BeetstraDragLaw::Pointer, BaseDragLawType>(m, "BeetstraDragLaw")
        .def(py::init<>())
        .def(py::init<Parameters&>())
        ;

    py::class_<SchillerAndNaumannDragLaw, SchillerAndNaumannDragLaw::Pointer, BaseDragLawType>(m, "SchillerAndNaumannDragLaw")
        .def(py::init<>())
        .def(py::init<Parameters&>())
        ;

    py::class_<HaiderAndLevenspielDragLaw, HaiderAndLevenspielDragLaw::Pointer, BaseDragLawType>(m, "HaiderAndLevenspielDragLaw")
        .def(py::init<>())
        .def(py::init<Parameters&>())
        ;

    py::class_<GanserDragLaw, GanserDragLaw::Pointer, BaseDragLawType>(m, "GanserDragLaw")
        .def(py::init<>())
        .def(py::init<Parameters&>())
        ;

    py::class_<ShahDragLaw, ShahDragLaw::Pointer, BaseDragLawType>(m, "ShahDragLaw")
        .def(py::init<>())
        .def(py::init<Parameters&>())
        ;

    py::class_<NewtonDragLaw, NewtonDragLaw::Pointer, BaseDragLawType>(m, "NewtonDragLaw")
        .def(py::init<>())
        .def(py::init<Parameters&>())
        ;

    // Inviscid force laws
    py::class_<InviscidForceLaw, InviscidForceLaw::Pointer>(m, "InviscidForceLaw")
        .def(py::init<>())
        .def(py::init<Parameters&>())
        .def("Clone", &InviscidForceLaw::Clone)
        .def("SetDragLawInProperties", &InviscidForceLaw::SetInviscidForceLawInProperties)
        .def("GetTypeOfLaw", &InviscidForceLaw::GetTypeOfLaw)
        ;

    py::class_<AutonHuntPrudhommeInviscidForceLaw, AutonHuntPrudhommeInviscidForceLaw::Pointer, BaseInviscidForceLawType>(m, "AutonHuntPrudhommeInviscidForceLaw")
        .def(py::init<>())
        .def(py::init<Parameters&>())
        ;

    // History force laws
    py::class_<HistoryForceLaw, HistoryForceLaw::Pointer>(m, "HistoryForceLaw")
        .def(py::init<>())
        .def(py::init<Parameters&>())
        .def("Clone", &HistoryForceLaw::Clone)
        .def("SetDragLawInProperties", &HistoryForceLaw::SetHistoryForceLawInProperties)
        .def("GetTypeOfLaw", &HistoryForceLaw::GetTypeOfLaw)
        ;

    py::class_<BoussinesqBassetHistoryForceLaw, BoussinesqBassetHistoryForceLaw::Pointer, BaseHistoryForceLaw>(m, "BoussinesqBassetHistoryForceLaw")
        .def(py::init<>())
        .def(py::init<Parameters&>())
        ;
  }

} // namespace Python.
} // namespace Kratos.
