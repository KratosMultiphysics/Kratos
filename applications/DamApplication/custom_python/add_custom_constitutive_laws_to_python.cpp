//
//   Project Name:         $
//   Last modified by:    $Author:             $
//   Date:                $Date:               $
//   Revision:            $Revision:           $
//

// System includes

// Project includes
#include "includes/constitutive_law.h"

//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

//constitutive laws
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_stress.hpp"

#include "custom_utilities/nodal_young_modulus_accessor.h"

#include "custom_constitutive/linear_elastic_2D_plane_strain_nodal.hpp"
#include "custom_constitutive/linear_elastic_2D_plane_stress_nodal.hpp"


#include "custom_constitutive/thermal_simo_ju_local_damage_3D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_local_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_local_damage_plane_stress_2D_law.hpp"

#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_plane_stress_2D_law.hpp"

#include "custom_constitutive/thermal_modified_mises_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/thermal_modified_mises_nonlocal_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/thermal_modified_mises_nonlocal_damage_plane_stress_2D_law.hpp"

namespace Kratos
{

    namespace Python
    {

        namespace py = pybind11;

        void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
        {
            py::class_< ThermalLinearElastic3DLaw, ThermalLinearElastic3DLaw::Pointer, ConstitutiveLaw >
            (m, "ThermalLinearElastic3DLaw")
            .def(py::init<>());
            py::class_< ThermalLinearElastic2DPlaneStrain, ThermalLinearElastic2DPlaneStrain::Pointer, ConstitutiveLaw >
            (m, "ThermalLinearElastic2DPlaneStrain")
            .def(py::init<>());
            py::class_< ThermalLinearElastic2DPlaneStress, ThermalLinearElastic2DPlaneStress::Pointer, ConstitutiveLaw >
            (m, "ThermalLinearElastic2DPlaneStress")
            .def(py::init<>());

            py::class_< LinearElastic2DPlaneStrainNodal, LinearElastic2DPlaneStrainNodal::Pointer, ConstitutiveLaw >
            (m, "LinearElastic2DPlaneStrainNodal")
            .def(py::init<>());
            py::class_< LinearElastic2DPlaneStressNodal, LinearElastic2DPlaneStressNodal::Pointer, ConstitutiveLaw >
            (m, "LinearElastic2DPlaneStressNodal")
            .def(py::init<>());

            // The 3D mechanical elastic law is now the standard, accessor-aware
            // StructuralMechanicsApplication flexible law; the thermal nodal laws are
            // the existing Dam thermal adapters (which retrieve E through the
            // accessor-aware Properties::GetValue). The historical nodal names are
            // kept as aliases so existing material files keep working. The spatially
            // varying Young's modulus is supplied by the NodalYoungModulusAccessor
            // installed by the nodal Young's modulus processes.
            m.attr("LinearElastic3DLawNodal") = py::module_::import(
                "KratosMultiphysics.StructuralMechanicsApplication").attr("FlexibleLinearElastic3DLaw");
            m.attr("ThermalLinearElastic3DLawNodal") = m.attr("ThermalLinearElastic3DLaw");
            m.attr("ThermalLinearElastic2DPlaneStrainNodal") = m.attr("ThermalLinearElastic2DPlaneStrain");
            m.attr("ThermalLinearElastic2DPlaneStressNodal") = m.attr("ThermalLinearElastic2DPlaneStress");

            // Installs the NodalYoungModulusAccessor (which exposes the historical
            // NODAL_YOUNG_MODULUS field as YOUNG_MODULUS through the standard
            // Properties::GetValue mechanism) on every property of the model part.
            // The nodal Young's modulus processes already do this; the helper is
            // provided for workflows that prescribe the nodal field by other means.
            m.def("AddNodalYoungModulusAccessor", [](ModelPart& rModelPart) {
                for (auto& r_properties : rModelPart.GetMesh(0).Properties()) {
                    if (!r_properties.HasAccessor(YOUNG_MODULUS)) {
                        r_properties.SetAccessor(
                            YOUNG_MODULUS,
                            Kratos::make_unique<NodalYoungModulusAccessor>());
                    }
                }
            }, py::arg("model_part"));

            py::class_< ThermalSimoJuLocalDamage3DLaw, ThermalSimoJuLocalDamage3DLaw::Pointer, ConstitutiveLaw >
            (m, "ThermalSimoJuLocalDamage3DLaw")
            .def(py::init<>());
            py::class_< ThermalSimoJuLocalDamagePlaneStrain2DLaw, ThermalSimoJuLocalDamagePlaneStrain2DLaw::Pointer, ConstitutiveLaw >
            (m, "ThermalSimoJuLocalDamagePlaneStrain2DLaw")
            .def(py::init<>());
            py::class_< ThermalSimoJuLocalDamagePlaneStress2DLaw, ThermalSimoJuLocalDamagePlaneStress2DLaw::Pointer, ConstitutiveLaw >
            (m, "ThermalSimoJuLocalDamagePlaneStress2DLaw")
            .def(py::init<>());

            py::class_< ThermalSimoJuNonlocalDamage3DLaw, ThermalSimoJuNonlocalDamage3DLaw::Pointer, ConstitutiveLaw >
            (m, "ThermalSimoJuNonlocalDamage3DLaw")
            .def(py::init<>());
            py::class_< ThermalSimoJuNonlocalDamagePlaneStrain2DLaw, ThermalSimoJuNonlocalDamagePlaneStrain2DLaw::Pointer, ConstitutiveLaw >
            (m, "ThermalSimoJuNonlocalDamagePlaneStrain2DLaw")
            .def(py::init<>());
            py::class_< ThermalSimoJuNonlocalDamagePlaneStress2DLaw, ThermalSimoJuNonlocalDamagePlaneStress2DLaw::Pointer, ConstitutiveLaw >
            (m, "ThermalSimoJuNonlocalDamagePlaneStress2DLaw")
            .def(py::init<>());

            py::class_< ThermalModifiedMisesNonlocalDamage3DLaw, ThermalModifiedMisesNonlocalDamage3DLaw::Pointer, ConstitutiveLaw >
            (m, "ThermalModifiedMisesNonlocalDamage3DLaw")
            .def(py::init<>());
            py::class_< ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw, ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw::Pointer, ConstitutiveLaw >
            (m, "ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw")
            .def(py::init<>());
            py::class_< ThermalModifiedMisesNonlocalDamagePlaneStress2DLaw, ThermalModifiedMisesNonlocalDamagePlaneStress2DLaw::Pointer, ConstitutiveLaw >
            (m, "ThermalModifiedMisesNonlocalDamagePlaneStress2DLaw")
            .def(py::init<>());
        }

    }  // namespace Python.
}  // namespace Kratos.
