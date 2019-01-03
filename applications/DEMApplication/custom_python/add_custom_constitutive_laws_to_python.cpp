//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"

#include "../custom_constitutive/DEM_discontinuum_constitutive_law.h"
#include "../custom_constitutive/DEM_continuum_constitutive_law.h"
#include "../custom_constitutive/DEM_compound_constitutive_law.h"

#include "../custom_constitutive/DEM_D_Linear_viscous_Coulomb_CL.h"
#include "../custom_constitutive/DEM_D_Hertz_viscous_Coulomb_CL.h"
#include "../custom_constitutive/DEM_D_Hertz_viscous_Coulomb_Nestle_CL.h"
#include "../custom_constitutive/DEM_D_Bentonite_Colloid_CL.h"
#include "../custom_constitutive/DEM_D_Linear_viscous_Coulomb_2D_CL.h"
#include "../custom_constitutive/DEM_D_Hertz_viscous_Coulomb_2D_CL.h"
#include "../custom_constitutive/DEM_D_JKR_cohesive_law.h"
#include "../custom_constitutive/DEM_D_DMT_cohesive_law.h"

#include "../custom_constitutive/DEM_D_Hertz_confined_CL.h"
#include "../custom_constitutive/DEM_D_Linear_confined_CL.h"
#include "../custom_constitutive/DEM_D_Linear_HighStiffness_CL.h"

#include "../custom_constitutive/DEM_Dempack_CL.h"
#include "../custom_constitutive/DEM_Dempack_2D_CL.h"
#include "../custom_constitutive/DEM_KDEM_CL.h"
#include "../custom_constitutive/DEM_KDEM_Rankine_CL.h"
#include "../custom_constitutive/DEM_KDEM_Mohr_Coulomb_CL.h"
#include "../custom_constitutive/DEM_KDEM_CamClay_CL.h"
#include "../custom_constitutive/dem_kdem_fissured_rock_cl.h"
#include "../custom_constitutive/DEM_sintering_continuum_CL.h"
#include "../custom_constitutive/DEM_KDEM_fabric_CL.h"
#include "../custom_constitutive/DEM_ExponentialHC_CL.h"
#include "../custom_constitutive/DEM_Dempack_torque_CL.h"
#include "../custom_constitutive/DEM_Dempack_dev_CL.h"
#include "../custom_constitutive/DEM_Dempack_2D_dev_CL.h"
#include "../custom_constitutive/dem_d_linear_custom_constants_cl.h"
#include "../custom_constitutive/DEM_D_Hertz_dependent_friction_CL.h"
#include "../custom_constitutive/dem_kdem_2d_cl.h"
#include "../custom_constitutive/dem_kdem_fabric_2d_cl.h"

namespace Kratos {
namespace Python {

namespace py = pybind11;

void AddCustomConstitutiveLawsToPython(pybind11::module& m) {

    // DEM Discontinuum Constitutive Laws:

    py::class_<DEMDiscontinuumConstitutiveLaw, DEMDiscontinuumConstitutiveLaw::Pointer>(m, "DEMDiscontinuumConstitutiveLaw")
        .def(py::init<>())
        .def("Clone", &DEMDiscontinuumConstitutiveLaw::Clone)
        .def("SetConstitutiveLawInProperties", &DEMDiscontinuumConstitutiveLaw::SetConstitutiveLawInProperties)
        .def("GetTypeOfLaw", &DEMDiscontinuumConstitutiveLaw::GetTypeOfLaw)
        ;

    py::class_<Variable<DEMDiscontinuumConstitutiveLaw::Pointer>, Variable<DEMDiscontinuumConstitutiveLaw::Pointer>::Pointer>(m, "DEMDiscontinuumConstitutiveLawPointerVariable")
        .def("__str__", PrintObject<Variable<DEMDiscontinuumConstitutiveLaw::Pointer>>)
        ;

    py::class_<DEM_D_Linear_viscous_Coulomb, DEM_D_Linear_viscous_Coulomb::Pointer, DEMDiscontinuumConstitutiveLaw>(m, "DEM_D_Linear_viscous_Coulomb")
        .def(py::init<>())
        ;

    py::class_<DEM_D_Bentonite_Colloid, DEM_D_Bentonite_Colloid::Pointer, DEMDiscontinuumConstitutiveLaw>(m, "DEM_D_Bentonite_Colloid")
        .def(py::init<>())
        ;

    py::class_<DEM_D_Linear_viscous_Coulomb2D, DEM_D_Linear_viscous_Coulomb2D::Pointer, DEM_D_Linear_viscous_Coulomb>(m, "DEM_D_Linear_viscous_Coulomb2D")
        .def(py::init<>())
        ;

    py::class_<DEM_D_Hertz_viscous_Coulomb, DEM_D_Hertz_viscous_Coulomb::Pointer, DEMDiscontinuumConstitutiveLaw>(m, "DEM_D_Hertz_viscous_Coulomb")
        .def(py::init<>())
        ;

    py::class_<DEM_D_Hertz_viscous_Coulomb2D, DEM_D_Hertz_viscous_Coulomb2D::Pointer, DEM_D_Hertz_viscous_Coulomb>(m, "DEM_D_Hertz_viscous_Coulomb2D")
        .def(py::init<>())
        ;

    py::class_<DEM_compound_constitutive_law<DEM_D_Hertz_viscous_Coulomb, DEM_D_JKR_Cohesive_Law>, DEM_compound_constitutive_law<DEM_D_Hertz_viscous_Coulomb, DEM_D_JKR_Cohesive_Law>::Pointer, DEM_D_Hertz_viscous_Coulomb>(m, "DEM_D_Hertz_viscous_Coulomb_JKR")
        .def(py::init<>())
        ;

    py::class_<DEM_compound_constitutive_law<DEM_D_Hertz_viscous_Coulomb, DEM_D_DMT_Cohesive_Law>, DEM_compound_constitutive_law<DEM_D_Hertz_viscous_Coulomb, DEM_D_DMT_Cohesive_Law>::Pointer, DEM_D_Hertz_viscous_Coulomb>(m, "DEM_D_Hertz_viscous_Coulomb_DMT")
        .def(py::init<>())
        ;

    py::class_<DEM_compound_constitutive_law<DEM_D_Linear_viscous_Coulomb, DEM_D_JKR_Cohesive_Law>, DEM_compound_constitutive_law<DEM_D_Linear_viscous_Coulomb, DEM_D_JKR_Cohesive_Law>::Pointer, DEM_D_Linear_viscous_Coulomb>(m, "DEM_D_Linear_viscous_Coulomb_JKR")
        .def(py::init<>())
        ;

    py::class_<DEM_compound_constitutive_law<DEM_D_Linear_viscous_Coulomb, DEM_D_DMT_Cohesive_Law>, DEM_compound_constitutive_law<DEM_D_Linear_viscous_Coulomb, DEM_D_DMT_Cohesive_Law>::Pointer, DEM_D_Linear_viscous_Coulomb>(m, "DEM_D_Linear_viscous_Coulomb_DMT")
        .def(py::init<>())
        ;

    py::class_<DEM_D_Linear_Custom_Constants, DEM_D_Linear_Custom_Constants::Pointer, DEM_D_Linear_viscous_Coulomb>(m, "DEM_D_Linear_Custom_Constants")
        .def(py::init<>())
        ;

    py::class_<DEM_D_Hertz_dependent_friction, DEM_D_Hertz_dependent_friction::Pointer, DEMDiscontinuumConstitutiveLaw>(m, "DEM_D_Hertz_dependent_friction")
        .def(py::init<>())
        ;

    py::class_<DEM_D_Hertz_confined, DEM_D_Hertz_confined::Pointer, DEM_D_Hertz_viscous_Coulomb>(m, "DEM_D_Hertz_confined")
        .def(py::init<>())
        ;

    py::class_<DEM_D_Linear_confined, DEM_D_Linear_confined::Pointer, DEM_D_Linear_viscous_Coulomb>(m, "DEM_D_Linear_confined")
        .def(py::init<>())
        ;

    py::class_<DEM_D_Hertz_viscous_Coulomb_Nestle, DEM_D_Hertz_viscous_Coulomb_Nestle::Pointer, DEM_D_Hertz_viscous_Coulomb>(m, "DEM_D_Hertz_viscous_Coulomb_Nestle")
        .def(py::init<>())
        ;

    py::class_<DEM_D_Linear_HighStiffness, DEM_D_Linear_HighStiffness::Pointer, DEMDiscontinuumConstitutiveLaw>(m, "DEM_D_Linear_HighStiffness")
        .def(py::init<>())
        ;
    // DEM Continuum Constitutive Laws:

    // DEM Continuum Constitutive Laws:

    py::class_<DEMContinuumConstitutiveLaw, DEMContinuumConstitutiveLaw::Pointer>(m, "DEMContinuumConstitutiveLaw")
        .def(py::init<>())
        .def("Clone", &DEMContinuumConstitutiveLaw::Clone)
        .def("SetConstitutiveLawInProperties", &DEMContinuumConstitutiveLaw::SetConstitutiveLawInProperties)
        .def("GetTypeOfLaw", &DEMContinuumConstitutiveLaw::GetTypeOfLaw)
        .def("CheckRequirementsOfStressTensor", &DEMContinuumConstitutiveLaw::CheckRequirementsOfStressTensor)
        ;

    py::class_<Variable<DEMContinuumConstitutiveLaw::Pointer>, Variable<DEMContinuumConstitutiveLaw::Pointer>::Pointer>(m, "DEMContinuumConstitutiveLawPointerVariable")
        .def("__str__", PrintObject<Variable<DEMContinuumConstitutiveLaw::Pointer>>)
        ;

    py::class_<DEM_Dempack, DEM_Dempack::Pointer, DEMContinuumConstitutiveLaw>(m, "DEM_Dempack")
        .def(py::init<>())
        ;

    py::class_<DEM_Dempack2D, DEM_Dempack2D::Pointer, DEM_Dempack>(m, "DEM_Dempack2D")
        .def(py::init<>())
        ;

    py::class_<DEM_Dempack_torque, DEM_Dempack_torque::Pointer, DEM_Dempack>(m, "DEM_Dempack_torque")
        .def(py::init<>())
        ;

    py::class_<DEM_Dempack_dev, DEM_Dempack_dev::Pointer, DEM_Dempack>(m, "DEM_Dempack_dev")
        .def(py::init<>())
        ;

    py::class_<DEM_Dempack2D_dev, DEM_Dempack2D_dev::Pointer, DEM_Dempack_dev>(m, "DEM_Dempack2D_dev")
        .def(py::init<>())
        ;

    py::class_<DEM_KDEM, DEM_KDEM::Pointer, DEMContinuumConstitutiveLaw>(m, "DEM_KDEM")
        .def(py::init<>())
        ;

    py::class_<DEM_sintering_continuum, DEM_sintering_continuum::Pointer, DEM_KDEM>(m, "DEM_sintering_continuum")
        .def(py::init<>())
        ;

    py::class_<DEM_KDEMFabric, DEM_KDEMFabric::Pointer, DEM_KDEM>(m, "DEM_KDEMFabric")
        .def(py::init<>())
        ;

    py::class_<DEM_KDEM_Rankine, DEM_KDEM_Rankine::Pointer, DEM_KDEM>(m, "DEM_KDEM_Rankine")
        .def(py::init<>())
        ;

    py::class_<DEM_KDEM_Mohr_Coulomb, DEM_KDEM_Mohr_Coulomb::Pointer, DEM_KDEM_Rankine>(m, "DEM_KDEM_Mohr_Coulomb")
        .def(py::init<>())
        ;
    
    py::class_<DEM_KDEM_CamClay, DEM_KDEM_CamClay::Pointer, DEM_KDEM_Rankine>(m, "DEM_KDEM_CamClay")
        .def(py::init<>())
        ;

    py::class_<DEM_KDEM_Fissured_Rock_CL, DEM_KDEM_Fissured_Rock_CL::Pointer, DEM_KDEM_Rankine>(m, "DEM_KDEM_Fissured_Rock")
        .def(py::init<>())
        ;

    py::class_<DEM_KDEM2D, DEM_KDEM2D::Pointer, DEM_KDEM>(m, "DEM_KDEM2D")
        .def(py::init<>())
        ;

    py::class_<DEM_KDEMFabric2D, DEM_KDEMFabric2D::Pointer, DEM_KDEM2D>(m, "DEM_KDEMFabric2D")
        .def(py::init<>())
        ;

    py::class_<DEM_ExponentialHC, DEM_ExponentialHC::Pointer, DEMContinuumConstitutiveLaw>(m, "DEM_ExponentialHC")
        .def(py::init<>())
        ;
}

} // namespace Python.
} // namespace Kratos.
