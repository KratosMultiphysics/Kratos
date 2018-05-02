/*
==============================================================================
KratosMultiScaleApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:37:00 $
//   Revision:            $Revision: 1.00 $
//
//

// System includes

// External includes
#include <boost/python.hpp>
#include "includes/constitutive_law.h"
#include "includes/properties.h"

// Project includes
#include "add_constitutive_laws_to_python.h"
#include "constitutive_laws/linear_elastic_no_tension_plane_stress_2d_law.h"
#include "constitutive_laws/linear_elastic_thick_shell_law.h"
#include "constitutive_laws/conv_diff_constitutive_law_3d.h"
#include "constitutive_laws/conv_diff_anisotropic_3d_law.h"
#include "constitutive_laws/conv_diff_plane_stress_2d_law.h"
#include "constitutive_laws/conv_diff_anisotropic_2d_law.h"
#include "constitutive_laws/j2_constitutive_law_3d.h"
#include "constitutive_laws/damage_iso_plane_stress_2d_law.h"
#include "constitutive_laws/damage_tc_plane_stress_2d_law.h"
#include "constitutive_laws/damage_tc_3d_law.h"
#include "constitutive_laws/scalar_damage_interface_2d_law.h"
#include "constitutive_laws/scalar_damage_interface_3d_law.h"
#include "constitutive_laws/conv_diff_interface_2d_law.h"
#include "constitutive_laws/conv_diff_interface_3d_law.h"
#include "constitutive_laws/plastic_damage_interface_2d_law.h"
#include "constitutive_laws/shell_from_3d_constitutive_law_adapter.h"
#include "constitutive_laws/planestress_from_3d_constitutive_law_adapter.h"
#include "constitutive_laws/planestrain_from_3d_constitutive_law_adapter.h"
#include "constitutive_laws/interpolated_constitutive_law_2d.h"
#include "multiscale_application_variables.h"

namespace Kratos
{
namespace Python
{

using namespace boost::python;


void AddConstitutiveLawsToPython()
{
	class_< InterpolatedConstitutiveLaw2D, bases< ConstitutiveLaw >, boost::noncopyable >(
		"InterpolatedConstitutiveLaw2D",
		init<const RveMaterialDatabase::Pointer&>())
		;

	class_< LinearElasticNoTensionPlaneStress2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >(
		"LinearElasticNoTensionPlaneStress2DLaw",
		init<>())
		;

	class_< ConvDiffConstitutiveLaw3D, bases< ConstitutiveLaw >, boost::noncopyable >(
		"ConvDiffConstitutiveLaw3D",
		init<>())
		;

	class_< ConvDiffAnisotropic3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >(
		"ConvDiffAnisotropic3DLaw",
		init<>())
		;

	class_< ConvDiffPlaneStress2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >(
		"ConvDiffPlaneStress2DLaw",
		init<>())
		;

	class_< ConvDiffAnisotropic2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >(
		"ConvDiffAnisotropic2DLaw",
		init<>())
		;


	class_< LinearElasticThickShellLaw, bases< ConstitutiveLaw >, boost::noncopyable >(
		"LinearElasticThickShellLaw",
		init<>())
		;

	class_< J2ConstitutiveLaw3D, bases< ConstitutiveLaw >, boost::noncopyable >(
		"J2ConstitutiveLaw3D",
		init<>())
		;

	class_< DamageIsoPlaneStress2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >(
		"DamageIsoPlaneStress2DLaw",
		init<>())
		;

	class_< DamageTCPlaneStress2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >(
		"DamageTCPlaneStress2DLaw",
		init<>())
		;

	class_< DamageTC3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >(
		"DamageTC3DLaw",
		init<>())
		;

	class_< ScalarDamageInterface2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >(
		"ScalarDamageInterface2DLaw",
		init<>())
		;

	class_< ScalarDamageInterface3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >(
		"ScalarDamageInterface3DLaw",
		init<>())
		;

	class_< ConvDiffInterface2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >(
		"ConvDiffInterface2DLaw",
		init<>())
		;

	class_< ConvDiffInterface3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >(
		"ConvDiffInterface3DLaw",
		init<>())
		;

	class_< PlasticDamageInterface2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >(
		"PlasticDamageInterface2DLaw",
		init<>())
		;

	typedef ShellFrom3DConstitutiveLawAdapter<ConstitutiveLaw> ShellFrom3DConstitutiveLawAdapterType;
	class_< ShellFrom3DConstitutiveLawAdapterType, ShellFrom3DConstitutiveLawAdapterType::Pointer,
		    bases< ConstitutiveLaw >, boost::noncopyable >(
			"ShellFrom3DConstitutiveLawAdapter",
			init<const ConstitutiveLaw::Pointer&>())
			;

	typedef PlaneStressFrom3DConstitutiveLawAdapter<ConstitutiveLaw> PlaneStressFrom3DConstitutiveLawAdapterType;
	class_< PlaneStressFrom3DConstitutiveLawAdapterType, PlaneStressFrom3DConstitutiveLawAdapterType::Pointer,
		    bases< ConstitutiveLaw >, boost::noncopyable >(
			"PlaneStressFrom3DConstitutiveLawAdapter",
			init<const ConstitutiveLaw::Pointer&>())
			;

	typedef PlaneStrainFrom3DConstitutiveLawAdapter<ConstitutiveLaw> PlaneStrainFrom3DConstitutiveLawAdapterType;
	class_< PlaneStrainFrom3DConstitutiveLawAdapterType, PlaneStrainFrom3DConstitutiveLawAdapterType::Pointer,
		    bases< ConstitutiveLaw >, boost::noncopyable >(
			"PlaneStrainFrom3DConstitutiveLawAdapter",
			init<const ConstitutiveLaw::Pointer&>())
			;

}

}

}
