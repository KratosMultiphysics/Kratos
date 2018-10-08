//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: iga_structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//


// System includes
#include <pybind11/stl.h>

// External includes
#include "includes/define_python.h"
#include "includes/constitutive_law.h"

// Project includes
#include "custom_python/add_constitutive_laws_to_python.h"
// #include "constitutive_laws/damage_tc_plane_stress_2d_law.h"
#include "constitutive_laws/plane_stress_tc_damage_law.h"
#include "custom_constitutive/plane_stress_2d_tc_damage_law.h"
#include "custom_constitutive/plane_stress_2d_kinematically_enriched_law.h"


namespace Kratos
{
namespace Python
{


    using namespace pybind11;

void AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
	// class_< DamageTCPlaneStress2DLaw, typename DamageTCPlaneStress2DLaw::Pointer, ConstitutiveLaw >
	// 	(m, "DamageTCPlaneStress2DLaw").def(init<>())
	// 	;

	//class_< PlaneStressTCDamageLaw, typename PlaneStressTCDamageLaw::Pointer, ConstitutiveLaw >
	//	(m, "PlaneStressTCDamageLaw").def(init<>())
	//	;
    class_< PlaneStress2dTCDamageLaw, typename PlaneStress2dTCDamageLaw::Pointer, ConstitutiveLaw >
        (m, "PlaneStress2dTCDamageLaw").def(init<>())
        ;

    class_< PlaneStress2dKinematicallyEnrichedLaw, typename PlaneStress2dKinematicallyEnrichedLaw::Pointer, ConstitutiveLaw >
        (m, "PlaneStress2dKinematicallyEnrichedLaw").def(init<>())
        ;
    //class_< PlaneStress2dTCDamageLaw, typename PlaneStress2dTCDamageLaw::Pointer, ConstitutiveLaw >
    //    (m, "PlaneStress2dTCDamageLaw").def(init<>())
    //    ;

}

}  // namespace Python.

}  // namespace Kratos.
