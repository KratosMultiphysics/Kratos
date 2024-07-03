//
// Author: Guillermo Casas (gcasas@cimne.upc.edu)
//

// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"

// Hydrodynamic Laws
#include "custom_constitutive/hydrodynamic_interaction_law.h"
#include "custom_constitutive/power_law_hydrodynamic_interaction_law.h"

// Buoyancy Laws
#include "custom_constitutive/buoyancy_laws/buoyancy_law.h"
#include "custom_constitutive/buoyancy_laws/archimedes_buoyancy_law.h"

// Drag laws
#include "../custom_constitutive/drag_laws/drag_law.h"
#include "../custom_constitutive/drag_laws/stokes_drag_law.h"
#include "../custom_constitutive/drag_laws/beetstra_drag_law.h"
#include "../custom_constitutive/drag_laws/schiller_and_naumann_drag_law.h"
#include "../custom_constitutive/drag_laws/haider_and_levenspiel_drag_law.h"
#include "../custom_constitutive/drag_laws/ganser_drag_law.h"
#include "../custom_constitutive/drag_laws/shah_drag_law.h"
#include "../custom_constitutive/drag_laws/newton_drag_law.h"
#include "../custom_constitutive/drag_laws/chien_drag_law.h"
#include "../custom_constitutive/drag_laws/dallavalle_drag_law.h"
#include "../custom_constitutive/drag_laws/difelice_drag_law.h"
#include "../custom_constitutive/drag_laws/rong_drag_law.h"

// Virtual mass force laws
#include "custom_constitutive/virtual_mass_force_laws/virtual_mass_force_law.h"
#include "custom_constitutive/virtual_mass_force_laws/auton_hunt_prudhomme_virtual_mass_force_law.h"
#include "custom_constitutive/virtual_mass_force_laws/zuber_virtual_mass_force_law.h"

// Undisturbed flow force laws
#include "custom_constitutive/undisturbed_force_laws/undisturbed_force_law.h"
#include "custom_constitutive/undisturbed_force_laws/pressure_gradient_undisturbed_force_law.h"

// History force laws
#include "custom_constitutive/history_force_laws/history_force_law.h"
#include "custom_constitutive/history_force_laws/boussinesq_basset_history_force_law.h"

// Vorticity-induced lift laws
#include "custom_constitutive/vorticity_induced_lift_laws/vorticity_induced_lift_law.h"
#include "custom_constitutive/vorticity_induced_lift_laws/el_samni_lift_law.h"
#include "custom_constitutive/vorticity_induced_lift_laws/saffman_lift_law.h"
#include "custom_constitutive/vorticity_induced_lift_laws/mei_lift_law.h"

// Rotation-induced lift laws
#include "custom_constitutive/rotation_induced_lift_laws/rotation_induced_lift_law.h"
#include "custom_constitutive/rotation_induced_lift_laws/rubinow_and_keller_lift_law.h"
#include "custom_constitutive/rotation_induced_lift_laws/oesterle_dinh_lift_law.h"
#include "custom_constitutive/rotation_induced_lift_laws/loth_rotation_induced_lift_law.h"

// Steady viscous torque laws
#include "custom_constitutive/steady_viscous_torque_laws/steady_viscous_torque_law.h"
#include "custom_constitutive/steady_viscous_torque_laws/rubinow_and_keller_torque_law.h"
#include "custom_constitutive/steady_viscous_torque_laws/loth_steady_viscous_torque_law.h"

namespace Kratos {
namespace Python {

namespace py = pybind11;

typedef HydrodynamicInteractionLaw BaseHydrodynamicInteractionLawType;
typedef BuoyancyLaw BaseBuoyancyLawType;
typedef DragLaw BaseDragLawType;
typedef VirtualMassForceLaw BaseVirtualMassForceLawType;
typedef UndisturbedForceLaw BaseUndisturbedForceLawType;
typedef HistoryForceLaw BaseHistoryForceLaw;
typedef VorticityInducedLiftLaw BaseVorticityInducedLiftLawType;
typedef RotationInducedLiftLaw BaseRotationInducedLiftLawType;
typedef SteadyViscousTorqueLaw BaseSteadyViscousTorqueLaw;

void AddCustomConstitutiveLawsToPython(pybind11::module& m) {

    py::class_<HydrodynamicInteractionLaw, HydrodynamicInteractionLaw::Pointer>(m, "HydrodynamicInteractionLaw")
        .def(py::init<>())
        .def(py::init<Properties::Pointer, Parameters&>())
        .def("Clone", &HydrodynamicInteractionLaw::Clone)
        .def("SetHydrodynamicInteractionLawInProperties", &HydrodynamicInteractionLaw::SetHydrodynamicInteractionLawInProperties)
        .def("GetTypeOfLaw", &HydrodynamicInteractionLaw::GetTypeOfLaw)
        .def("SetBuoyancyLaw", &HydrodynamicInteractionLaw::SetBuoyancyLaw)
        .def("SetDragLaw", &HydrodynamicInteractionLaw::SetDragLaw)
        .def("SetVirtualMassForceLaw", &HydrodynamicInteractionLaw::SetVirtualMassForceLaw)
        .def("SetUndisturbedForceLaw", &HydrodynamicInteractionLaw::SetUndisturbedForceLaw)
        .def("SetHistoryForceLaw", &HydrodynamicInteractionLaw::SetHistoryForceLaw)
        .def("SetVorticityInducedLiftLaw", &HydrodynamicInteractionLaw::SetVorticityInducedLiftLaw)
        .def("SetRotationInducedLiftLaw", &HydrodynamicInteractionLaw::SetRotationInducedLiftLaw)
        .def("SetSteadyViscousTorqueLaw", &HydrodynamicInteractionLaw::SetSteadyViscousTorqueLaw)
        ;

     py::class_<PowerLawFluidHydrodynamicInteractionLaw,
                PowerLawFluidHydrodynamicInteractionLaw::Pointer,
                BaseHydrodynamicInteractionLawType>(m, "PowerLawFluidHydrodynamicInteractionLaw")
        .def(py::init<>())
        .def(py::init<Properties::Pointer, Parameters&>())
        ;

    // Buoyancy laws
    py::class_<BuoyancyLaw, BuoyancyLaw::Pointer>(m, "BuoyancyLaw")
        .def(py::init<Parameters>())
        .def(py::init<>())
        .def("Clone", &BuoyancyLaw::Clone)
        .def("SetBuoyancyLawInProperties", &BuoyancyLaw::SetBuoyancyLawInProperties)
        .def("GetTypeOfLaw", &BuoyancyLaw::GetTypeOfLaw)
        ;

    py::class_<ArchimedesBuoyancyLaw, ArchimedesBuoyancyLaw::Pointer, BaseBuoyancyLawType>(m, "ArchimedesBuoyancyLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;

    // Drag laws

    py::class_<DragLaw, DragLaw::Pointer>(m, "DragLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        .def("Clone", &DragLaw::Clone)
        .def("SetDragLawInProperties", &DragLaw::SetDragLawInProperties)
        .def("GetTypeOfLaw", &DragLaw::GetTypeOfLaw)
        ;

    py::class_<StokesDragLaw, StokesDragLaw::Pointer, BaseDragLawType>(m, "StokesDragLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;

    py::class_<BeetstraDragLaw, BeetstraDragLaw::Pointer, BaseDragLawType>(m, "BeetstraDragLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;

    py::class_<SchillerAndNaumannDragLaw, SchillerAndNaumannDragLaw::Pointer, BaseDragLawType>(m, "SchillerAndNaumannDragLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;

    py::class_<HaiderAndLevenspielDragLaw, HaiderAndLevenspielDragLaw::Pointer, BaseDragLawType>(m, "HaiderAndLevenspielDragLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;

    py::class_<ChienDragLaw, ChienDragLaw::Pointer, BaseDragLawType>(m, "ChienDragLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;

    py::class_<GanserDragLaw, GanserDragLaw::Pointer, BaseDragLawType>(m, "GanserDragLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;

    py::class_<ShahDragLaw, ShahDragLaw::Pointer, BaseDragLawType>(m, "ShahDragLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;

    py::class_<NewtonDragLaw, NewtonDragLaw::Pointer, BaseDragLawType>(m, "NewtonDragLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;

    py::class_<DallavalleDragLaw, DallavalleDragLaw::Pointer, BaseDragLawType>(m, "DallavalleDragLaw")
        .def(py::init<>())
        .def(py::init<Parameters&>())
        ;

    py::class_<DiFeliceDragLaw, DiFeliceDragLaw::Pointer, BaseDragLawType>(m, "DiFeliceDragLaw")
        .def(py::init<>())
        .def(py::init<Parameters&>())
        ;

    py::class_<RongDragLaw, RongDragLaw::Pointer, BaseDragLawType>(m, "RongDragLaw")
        .def(py::init<>())
        .def(py::init<Parameters&>())
        ;

    // Virtual mass force laws
    py::class_<VirtualMassForceLaw, VirtualMassForceLaw::Pointer>(m, "VirtualMassForceLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        .def("Clone", &VirtualMassForceLaw::Clone)
        .def("SetVirtualMassForceLawInProperties", &VirtualMassForceLaw::SetVirtualMassForceLawInProperties)
        .def("GetTypeOfLaw", &VirtualMassForceLaw::GetTypeOfLaw)
        ;

    py::class_<AutonHuntPrudhommeVirtualMassForceLaw, AutonHuntPrudhommeVirtualMassForceLaw::Pointer, BaseVirtualMassForceLawType>(m, "AutonHuntPrudhommeVirtualMassForceLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;


    py::class_<ZuberVirtualMassForceLaw, ZuberVirtualMassForceLaw::Pointer, BaseVirtualMassForceLawType>(m, "ZuberVirtualMassForceLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;

    // Undisturbed flow force laws
    py::class_<UndisturbedForceLaw, UndisturbedForceLaw::Pointer>(m, "UndisturbedForceLaw")
            .def(py::init<>())
            .def(py::init<Parameters>())
            ;


    py::class_<PressureGradientUndisturbedForceLaw, PressureGradientUndisturbedForceLaw::Pointer, BaseUndisturbedForceLawType>(m, "PressureGradientUndisturbedForceLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;


    // History force laws
    py::class_<HistoryForceLaw, HistoryForceLaw::Pointer>(m, "HistoryForceLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        .def("Clone", &HistoryForceLaw::Clone)
        .def("SetHistoryForceLawInProperties", &HistoryForceLaw::SetHistoryForceLawInProperties)
        .def("GetTypeOfLaw", &HistoryForceLaw::GetTypeOfLaw)
        ;

    py::class_<BoussinesqBassetHistoryForceLaw, BoussinesqBassetHistoryForceLaw::Pointer, BaseHistoryForceLaw>(m, "BoussinesqBassetHistoryForceLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;

    // Vorticity-induced lift laws
    py::class_<VorticityInducedLiftLaw, VorticityInducedLiftLaw::Pointer>(m, "VorticityInducedLiftLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        .def("Clone", &VorticityInducedLiftLaw::Clone)
        .def("SetVorticityInducedLiftLawInProperties", &VorticityInducedLiftLaw::SetVorticityInducedLiftLawInProperties)
        .def("GetTypeOfLaw", &VorticityInducedLiftLaw::GetTypeOfLaw)
        ;

    py::class_<ElSamniLiftLaw, ElSamniLiftLaw::Pointer, BaseVorticityInducedLiftLawType>(m, "ElSamniLiftLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;

    py::class_<SaffmanLiftLaw, SaffmanLiftLaw::Pointer, BaseVorticityInducedLiftLawType>(m, "SaffmanLiftLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;

    py::class_<MeiLiftLaw, MeiLiftLaw::Pointer, BaseVorticityInducedLiftLawType>(m, "MeiLiftLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;

    // Rotation-induced lift laws
    py::class_<RotationInducedLiftLaw, RotationInducedLiftLaw::Pointer>(m, "RotationInducedLiftLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        .def("Clone", &RotationInducedLiftLaw::Clone)
        .def("SetRotationInducedLiftLawInProperties", &RotationInducedLiftLaw::SetRotationInducedLiftLawInProperties)
        .def("GetTypeOfLaw", &RotationInducedLiftLaw::GetTypeOfLaw)
        ;

    py::class_<RubinowAndKellerLiftLaw, RubinowAndKellerLiftLaw::Pointer, BaseRotationInducedLiftLawType>(m, "RubinowAndKellerLiftLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;

    py::class_<OesterleAndDinhLiftLaw, OesterleAndDinhLiftLaw::Pointer, BaseRotationInducedLiftLawType>(m, "OesterleAndDinhLiftLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;

    py::class_<LothRotationInducedLiftLaw, LothRotationInducedLiftLaw::Pointer, BaseRotationInducedLiftLawType>(m, "LothRotationInducedLiftLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;

    // Steady viscous torque laws
    py::class_<SteadyViscousTorqueLaw, SteadyViscousTorqueLaw::Pointer>(m, "SteadyViscousTorqueLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        .def("Clone", &SteadyViscousTorqueLaw::Clone)
        .def("SetSteadyViscousTorqueLawInProperties", &SteadyViscousTorqueLaw::SetSteadyViscousTorqueLawInProperties)
        .def("GetTypeOfLaw", &SteadyViscousTorqueLaw::GetTypeOfLaw)
        ;

    py::class_<RubinowAndKellerTorqueLaw, RubinowAndKellerTorqueLaw::Pointer, BaseSteadyViscousTorqueLaw>(m, "RubinowAndKellerTorqueLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;

    py::class_<LothSteadyViscousTorqueLaw, LothSteadyViscousTorqueLaw::Pointer, BaseSteadyViscousTorqueLaw>(m, "LothSteadyViscousTorqueLaw")
        .def(py::init<>())
        .def(py::init<Parameters>())
        ;
  }

} // namespace Python.
} // namespace Kratos.
