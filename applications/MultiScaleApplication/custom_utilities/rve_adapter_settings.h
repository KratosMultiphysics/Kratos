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
//   Date:                $Date: 2013-06-06 10:37:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(RVE_ADAPTER_SETTINGS_H_INCLUDED)
#define RVE_ADAPTER_SETTINGS_H_INCLUDED

#include "includes/constitutive_law.h"

namespace Kratos
{

struct RveAdapterSettings_Thermal_2D
{
	typedef ConstitutiveLaw::SizeType SizeType;

	static SizeType GetStrainSize() { return 2; }
	static SizeType WorkingSpaceDimension() { return 2; }
	static void GetLawFeatures(ConstitutiveLaw::Features rFeatures)
	{
		rFeatures.mOptions.Set(ConstitutiveLaw::PLANE_STRESS_LAW);
		rFeatures.mOptions.Set(ConstitutiveLaw::INFINITESIMAL_STRAINS);

		rFeatures.mStrainMeasures.push_back(ConstitutiveLaw::StrainMeasure_Infinitesimal);

		rFeatures.mStrainSize = GetStrainSize();

		rFeatures.mSpaceDimension = WorkingSpaceDimension();
	}
};

struct RveAdapterSettings_Thermal_3D
{
	typedef ConstitutiveLaw::SizeType SizeType;

	static SizeType GetStrainSize() { return 3; }
	static SizeType WorkingSpaceDimension() { return 3; }
	static void GetLawFeatures(ConstitutiveLaw::Features rFeatures)
	{
		rFeatures.mOptions.Set(ConstitutiveLaw::THREE_DIMENSIONAL_LAW);
		rFeatures.mOptions.Set(ConstitutiveLaw::INFINITESIMAL_STRAINS);

		rFeatures.mStrainMeasures.push_back(ConstitutiveLaw::StrainMeasure_Infinitesimal);

		rFeatures.mStrainSize = GetStrainSize();

		rFeatures.mSpaceDimension = WorkingSpaceDimension();
	}
};

struct RveAdapterSettings_PlaneStress
{
	typedef ConstitutiveLaw::SizeType SizeType;

	static SizeType GetStrainSize() { return 3; }
	static SizeType WorkingSpaceDimension() { return 2; }
	static ConstitutiveLaw::StrainMeasure GetStrainMeasure() { return ConstitutiveLaw::StrainMeasure_Infinitesimal; }
	static ConstitutiveLaw::StressMeasure GetStressMeasure() { return ConstitutiveLaw::StressMeasure_Cauchy; }
	static void GetLawFeatures(ConstitutiveLaw::Features rFeatures)
	{
		rFeatures.mOptions.Set( ConstitutiveLaw::PLANE_STRESS_LAW );
		rFeatures.mOptions.Set( ConstitutiveLaw::INFINITESIMAL_STRAINS );

		rFeatures.mStrainMeasures.push_back(ConstitutiveLaw::StrainMeasure_Infinitesimal);

		rFeatures.mStrainSize = GetStrainSize();

		rFeatures.mSpaceDimension = WorkingSpaceDimension();
	}
};

struct RveAdapterSettings_PlaneStrain
{
	typedef ConstitutiveLaw::SizeType SizeType;

	static SizeType GetStrainSize() { return 4; }
	static SizeType WorkingSpaceDimension() { return 2; }
	static ConstitutiveLaw::StrainMeasure GetStrainMeasure() { return ConstitutiveLaw::StrainMeasure_Infinitesimal; }
	static ConstitutiveLaw::StressMeasure GetStressMeasure() { return ConstitutiveLaw::StressMeasure_Cauchy; }
	static void GetLawFeatures(ConstitutiveLaw::Features rFeatures)
	{
		rFeatures.mOptions.Set( ConstitutiveLaw::PLANE_STRAIN_LAW );
		rFeatures.mOptions.Set( ConstitutiveLaw::INFINITESIMAL_STRAINS );

		rFeatures.mStrainMeasures.push_back(ConstitutiveLaw::StrainMeasure_Infinitesimal);

		rFeatures.mStrainSize = GetStrainSize();

		rFeatures.mSpaceDimension = WorkingSpaceDimension();
	}
};

struct RveAdapterSettings_3D
{
	typedef ConstitutiveLaw::SizeType SizeType;

	static SizeType GetStrainSize() { return 6; }
	static SizeType WorkingSpaceDimension() { return 3; }
	static ConstitutiveLaw::StrainMeasure GetStrainMeasure() { return ConstitutiveLaw::StrainMeasure_Infinitesimal; }
	static ConstitutiveLaw::StressMeasure GetStressMeasure() { return ConstitutiveLaw::StressMeasure_Cauchy; }
	static void GetLawFeatures(ConstitutiveLaw::Features rFeatures)
	{
		rFeatures.mOptions.Set( ConstitutiveLaw::THREE_DIMENSIONAL_LAW );
		rFeatures.mOptions.Set( ConstitutiveLaw::INFINITESIMAL_STRAINS );

		rFeatures.mStrainMeasures.push_back(ConstitutiveLaw::StrainMeasure_Infinitesimal);

		rFeatures.mStrainSize = GetStrainSize();

		rFeatures.mSpaceDimension = WorkingSpaceDimension();
	}
};

struct RveAdapterSettings_ThickShell
{
	typedef ConstitutiveLaw::SizeType SizeType;

	static SizeType GetStrainSize() { return 8; }
	static SizeType WorkingSpaceDimension() { return 3; }
	static ConstitutiveLaw::StrainMeasure GetStrainMeasure() { return ConstitutiveLaw::StrainMeasure_Infinitesimal; }
	static ConstitutiveLaw::StressMeasure GetStressMeasure() { return ConstitutiveLaw::StressMeasure_Cauchy; }
	static void GetLawFeatures(ConstitutiveLaw::Features rFeatures)
	{
		rFeatures.mOptions.Set( ConstitutiveLaw::THREE_DIMENSIONAL_LAW );
		rFeatures.mOptions.Set( ConstitutiveLaw::INFINITESIMAL_STRAINS );

		rFeatures.mStrainMeasures.push_back(ConstitutiveLaw::StrainMeasure_Infinitesimal);

		rFeatures.mStrainSize = GetStrainSize();

		rFeatures.mSpaceDimension = WorkingSpaceDimension();
	}
};

struct RveAdapterSettings_ThinShell
{
	typedef ConstitutiveLaw::SizeType SizeType;

	static SizeType GetStrainSize() { return 6; }
	static SizeType WorkingSpaceDimension() { return 3; }
	static ConstitutiveLaw::StrainMeasure GetStrainMeasure() { return ConstitutiveLaw::StrainMeasure_Infinitesimal; }
	static ConstitutiveLaw::StressMeasure GetStressMeasure() { return ConstitutiveLaw::StressMeasure_Cauchy; }
	static void GetLawFeatures(ConstitutiveLaw::Features rFeatures)
	{
		rFeatures.mOptions.Set( ConstitutiveLaw::PLANE_STRESS_LAW );
		rFeatures.mOptions.Set( ConstitutiveLaw::INFINITESIMAL_STRAINS );

		rFeatures.mStrainMeasures.push_back(ConstitutiveLaw::StrainMeasure_Infinitesimal);

		rFeatures.mStrainSize = GetStrainSize();

		rFeatures.mSpaceDimension = WorkingSpaceDimension();
	}
};

}


#endif // RVE_ADAPTER_SETTINGS_H_INCLUDED
