//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"

#include "custom_constitutive/DEM_discontinuum_constitutive_law.h"
#include "custom_constitutive/DEM_continuum_constitutive_law.h"

#include "custom_constitutive/DEM_compound_constitutive_law.h"
#include "custom_constitutive/DEM_compound_constitutive_law_for_PBM.h"

#include "custom_constitutive/DEM_D_Linear_viscous_Coulomb_CL.h"
#include "custom_constitutive/DEM_D_Hertz_viscous_Coulomb_CL.h"
#include "custom_constitutive/DEM_D_Hertz_viscous_Coulomb_Nestle_CL.h"
#include "custom_constitutive/DEM_D_Bentonite_Colloid_CL.h"
#include "custom_constitutive/DEM_D_Linear_viscous_Coulomb_2D_CL.h"
#include "custom_constitutive/DEM_D_Hertz_viscous_Coulomb_2D_CL.h"
#include "custom_constitutive/DEM_D_JKR_cohesive_law.h"
#include "custom_constitutive/DEM_D_DMT_cohesive_law.h"
#include "custom_constitutive/DEM_D_Stress_dependent_cohesive_CL.h"
#include "custom_constitutive/DEM_D_Quadratic_CL.h"
#include "custom_constitutive/DEM_D_void_CL.h"

#include "custom_constitutive/DEM_rolling_friction_model.h"

#include "custom_constitutive/DEM_D_Hertz_confined_CL.h"
#include "custom_constitutive/DEM_D_Linear_confined_CL.h"
#include "custom_constitutive/DEM_D_Linear_HighStiffness_CL.h"
#include "custom_constitutive/DEM_D_Linear_HighStiffness_2D_CL.h"
#include "custom_constitutive/DEM_D_Linear_classic_CL.h"

#include "custom_constitutive/DEM_Dempack_CL.h"
#include "custom_constitutive/DEM_Dempack_2D_CL.h"
#include "custom_constitutive/DEM_KDEM_CL.h"
#include "custom_constitutive/DEM_KDEM_soft_torque_CL.h"
#include "custom_constitutive/DEM_KDEM_soft_torque_with_noise_CL.h"
#include "custom_constitutive/DEM_KDEM_with_damage_CL.h"
#include "custom_constitutive/DEM_KDEM_with_damage_parallel_bond_CL.h"
#include "custom_constitutive/DEM_KDEM_with_damage_parallel_bond_capped_CL.h"
#include "custom_constitutive/DEM_KDEM_with_damage_parallel_bond_2D_CL.h"
#include "custom_constitutive/DEM_KDEM_with_damage_parallel_bond_Hertz_CL.h"
#include "custom_constitutive/DEM_KDEM_with_damage_parallel_bond_Hertz_2D_CL.h"
#include "custom_constitutive/DEM_KDEM_Rankine_CL.h"
#include "custom_constitutive/DEM_KDEM_Mohr_Coulomb_CL.h"
#include "custom_constitutive/DEM_KDEM_CamClay_CL.h"
#include "custom_constitutive/dem_kdem_fissured_rock_cl.h"
#include "custom_constitutive/DEM_KDEM_fabric_CL.h"
#include "custom_constitutive/DEM_beam_constitutive_law.h"
#include "custom_constitutive/DEM_ExponentialHC_CL.h"
#include "custom_constitutive/DEM_Dempack_torque_CL.h"
#include "custom_constitutive/DEM_Dempack_dev_CL.h"
#include "custom_constitutive/DEM_Dempack_2D_dev_CL.h"
#include "custom_constitutive/dem_d_linear_custom_constants_cl.h"
#include "custom_constitutive/DEM_D_Conical_damage_CL.h"
#include "custom_constitutive/dem_kdem_2d_cl.h"
#include "custom_constitutive/dem_kdem_fabric_2d_cl.h"
#include "custom_constitutive/DEM_parallel_bond_CL.h"
#include "custom_constitutive/DEM_smooth_joint_CL.h"
#include "custom_constitutive/DEM_parallel_bond_for_membrane_CL.h"
#include "custom_constitutive/DEM_rolling_friction_model_constant_torque.h"
#include "custom_constitutive/DEM_rolling_friction_model_viscous_torque.h"
#include "custom_constitutive/DEM_rolling_friction_model_bounded.h"


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

    py::class_<DEM_D_Stress_Dependent_Cohesive, DEM_D_Stress_Dependent_Cohesive::Pointer, DEMDiscontinuumConstitutiveLaw>(m, "DEM_D_Stress_Dependent_Cohesive")
        .def(py::init<>())
        ;

    py::class_<DEM_D_Conical_damage, DEM_D_Conical_damage::Pointer, DEMDiscontinuumConstitutiveLaw>(m, "DEM_D_Conical_damage")
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

    py::class_<DEM_D_Linear_HighStiffness_2D, DEM_D_Linear_HighStiffness_2D::Pointer, DEMDiscontinuumConstitutiveLaw>(m, "DEM_D_Linear_HighStiffness_2D")
        .def(py::init<>())
        ;

    py::class_<DEM_D_Quadratic, DEM_D_Quadratic::Pointer, DEMDiscontinuumConstitutiveLaw>(m, "DEM_D_Quadratic")
        .def(py::init<>())
        ;

    py::class_<DEM_D_Linear_classic, DEM_D_Linear_classic::Pointer, DEMDiscontinuumConstitutiveLaw>(m, "DEM_D_Linear_classic")
        .def(py::init<>())
        ;

    py::class_<DEM_D_void, DEM_D_void::Pointer, DEMDiscontinuumConstitutiveLaw>(m, "DEM_D_void")
        .def(py::init<>())
        ;

    // DEM Continuum Constitutive Laws:

    py::class_<DEMContinuumConstitutiveLaw, DEMContinuumConstitutiveLaw::Pointer>(m, "DEMContinuumConstitutiveLaw")
        .def(py::init<>())
        .def("Clone", &DEMContinuumConstitutiveLaw::Clone)
        .def("SetConstitutiveLawInProperties", &DEMContinuumConstitutiveLaw::SetConstitutiveLawInProperties)
        .def("SetConstitutiveLawInPropertiesWithParameters", &DEMContinuumConstitutiveLaw::SetConstitutiveLawInPropertiesWithParameters)
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

    py::class_<DEM_KDEM_soft_torque, DEM_KDEM_soft_torque::Pointer, DEM_KDEM>(m, "DEM_KDEM_soft_torque")
        .def(py::init<>())
        ;

    py::class_<DEM_KDEM_soft_torque_with_noise, DEM_KDEM_soft_torque_with_noise::Pointer, DEM_KDEM_soft_torque>(m, "DEM_KDEM_soft_torque_with_noise")
        .def(py::init<>())
        ;

    py::class_<DEM_KDEM_with_damage, DEM_KDEM_with_damage::Pointer, DEM_KDEM_soft_torque>(m, "DEM_KDEM_with_damage")
        .def(py::init<>())
        ;

    py::class_<DEM_KDEM_with_damage_parallel_bond, DEM_KDEM_with_damage_parallel_bond::Pointer, DEM_KDEM_with_damage>(m, "DEM_KDEM_with_damage_parallel_bond")
        .def(py::init<>())
        ;

    py::class_<DEM_KDEM_with_damage_parallel_bond_capped, DEM_KDEM_with_damage_parallel_bond_capped::Pointer, DEM_KDEM_with_damage_parallel_bond>(m, "DEM_KDEM_with_damage_parallel_bond_capped")
        .def(py::init<>())
        ;

    py::class_<DEM_KDEM_with_damage_parallel_bond_2D, DEM_KDEM_with_damage_parallel_bond_2D::Pointer, DEM_KDEM_with_damage_parallel_bond>(m, "DEM_KDEM_with_damage_parallel_bond_2D")
        .def(py::init<>())
        ;

    py::class_<DEM_KDEM_with_damage_parallel_bond_Hertz, DEM_KDEM_with_damage_parallel_bond_Hertz::Pointer, DEM_KDEM_with_damage_parallel_bond>(m, "DEM_KDEM_with_damage_parallel_bond_Hertz")
        .def(py::init<>())
        ;
    
    py::class_<DEM_KDEM_with_damage_parallel_bond_Hertz_2D, DEM_KDEM_with_damage_parallel_bond_Hertz_2D::Pointer, DEM_KDEM_with_damage_parallel_bond_Hertz>(m, "DEM_KDEM_with_damage_parallel_bond_Hertz_2D")
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

    py::class_<DEM_parallel_bond, DEM_parallel_bond::Pointer, DEMContinuumConstitutiveLaw>(m, "DEM_parallel_bond")
        .def(py::init<>())
        ;

    py::class_<DEM_smooth_joint, DEM_smooth_joint::Pointer, DEMContinuumConstitutiveLaw>(m, "DEM_smooth_joint")
        .def(py::init<>())
        ;

    py::class_<DEM_parallel_bond_for_membrane, DEM_parallel_bond_for_membrane::Pointer, DEMContinuumConstitutiveLaw>(m, "DEM_parallel_bond_for_membrane")
        .def(py::init<>())
        ;

    py::class_<DEM_compound_constitutive_law_for_PBM<DEM_parallel_bond, DEM_D_Hertz_viscous_Coulomb>, DEM_compound_constitutive_law_for_PBM<DEM_parallel_bond, DEM_D_Hertz_viscous_Coulomb>::Pointer, DEM_parallel_bond>(m, "DEM_parallel_bond_Hertz")
        .def(py::init<>())
        ;

    py::class_<DEM_compound_constitutive_law_for_PBM<DEM_parallel_bond, DEM_D_Linear_classic>, DEM_compound_constitutive_law_for_PBM<DEM_parallel_bond, DEM_D_Linear_classic>::Pointer, DEM_parallel_bond>(m, "DEM_parallel_bond_Linear")
        .def(py::init<>())
        ;

    py::class_<DEM_compound_constitutive_law_for_PBM<DEM_parallel_bond, DEM_D_Quadratic>, DEM_compound_constitutive_law_for_PBM<DEM_parallel_bond, DEM_D_Quadratic>::Pointer, DEM_parallel_bond>(m, "DEM_parallel_bond_Quadratic")
        .def(py::init<>())
        ;


    //for compound constitutive law

    // DEM Beam Constitutive Laws:

    py::class_<DEMBeamConstitutiveLaw, DEMBeamConstitutiveLaw::Pointer>(m, "DEMBeamConstitutiveLaw")
        .def(py::init<>())
        .def("Clone", &DEMBeamConstitutiveLaw::Clone)
        .def("SetConstitutiveLawInProperties", &DEMBeamConstitutiveLaw::SetConstitutiveLawInProperties)
        .def("GetTypeOfLaw", &DEMBeamConstitutiveLaw::GetTypeOfLaw)
        .def("CheckRequirementsOfStressTensor", &DEMBeamConstitutiveLaw::CheckRequirementsOfStressTensor)
        ;

    py::class_<Variable<DEMBeamConstitutiveLaw::Pointer>, Variable<DEMBeamConstitutiveLaw::Pointer>::Pointer>(m, "DEMBeamConstitutiveLawPointerVariable")
        .def("__str__", PrintObject<Variable<DEMBeamConstitutiveLaw::Pointer>>)
        ;

    // DEM Rolling Friction Model:

    py::class_<DEMRollingFrictionModel, DEMRollingFrictionModel::Pointer>(m, "DEMRollingFrictionModel")
        .def(py::init<>())
        .def("Clone", &DEMRollingFrictionModel::Clone)
        .def("SetAPrototypeOfThisInProperties", &DEMRollingFrictionModel::SetAPrototypeOfThisInProperties)
        //.def("GetTypeOfLaw", &DEMRollingFrictionModel::GetTypeOfLaw)
        ;

    py::class_<Variable<DEMRollingFrictionModel::Pointer>, Variable<DEMRollingFrictionModel::Pointer>::Pointer>(m, "DEMRollingFrictionModelPointerVariable")
        .def("__str__", PrintObject<Variable<DEMRollingFrictionModel::Pointer>>)
        ;

    py::class_<DEMRollingFrictionModelConstantTorque, DEMRollingFrictionModelConstantTorque::Pointer, DEMRollingFrictionModel>(m, "DEMRollingFrictionModelConstantTorque")
    .def(py::init<>())
    ;

    py::class_<DEMRollingFrictionModelViscousTorque, DEMRollingFrictionModelViscousTorque::Pointer, DEMRollingFrictionModel>(m, "DEMRollingFrictionModelViscousTorque")
    .def(py::init<>())
    ;

    py::class_<DEMRollingFrictionModelBounded, DEMRollingFrictionModelBounded::Pointer, DEMRollingFrictionModel>(m, "DEMRollingFrictionModelBounded")
    .def(py::init<>())
    ;
}

} // namespace Python.
} // namespace Kratos.
