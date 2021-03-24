//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    JMCarbonell
//					 (adapted to Particle Mechanics by Peter Wilson)
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/johnson_cook_thermal_plastic_plane_strain_2D_law.hpp"
#include "particle_mechanics_application_variables.h"

namespace Kratos
{
    JohnsonCookThermalPlastic2DPlaneStrainLaw::JohnsonCookThermalPlastic2DPlaneStrainLaw()
  : JohnsonCookThermalPlastic3DLaw()
  { }

    JohnsonCookThermalPlastic2DPlaneStrainLaw::JohnsonCookThermalPlastic2DPlaneStrainLaw(const JohnsonCookThermalPlastic2DPlaneStrainLaw& rOther)
  : JohnsonCookThermalPlastic3DLaw(rOther)
  { }

  ConstitutiveLaw::Pointer JohnsonCookThermalPlastic2DPlaneStrainLaw::Clone() const
  {
    return Kratos::make_shared<JohnsonCookThermalPlastic2DPlaneStrainLaw>(*this);
  }

  JohnsonCookThermalPlastic2DPlaneStrainLaw::~JohnsonCookThermalPlastic2DPlaneStrainLaw()
  { }

  void JohnsonCookThermalPlastic2DPlaneStrainLaw::MakeStrainStressVectorFromMatrix(const Matrix& rInput, Vector& rOutput)
  {
	  if (rOutput.size() != GetStrainSize()) rOutput.resize(GetStrainSize(), false);

	  // 2D stress arrangement
	  rOutput[0] = rInput(0, 0);
	  rOutput[1] = rInput(1, 1);
	  rOutput[2] = 2.0 * rInput(0, 1); //xy
  }


  void JohnsonCookThermalPlastic2DPlaneStrainLaw::MakeStrainStressMatrixFromVector(const Vector& rInput, Matrix& rOutput)
  {
	  if (rOutput.size1() != 2 || rOutput.size2() != 2)rOutput.resize(2, 2, false);

	  // 3D stress arrangement
      rOutput.clear();

      // Normal components
	  rOutput(0, 0) = rInput[0];
	  rOutput(1, 1) = rInput[1];

      // Shear components
	  rOutput(0, 1) = 0.5 * rInput[2]; //xy

      // Fill symmetry
	  rOutput(1, 0) = rOutput(0, 1);
  }

  void JohnsonCookThermalPlastic2DPlaneStrainLaw::ComputeCharacteristicLength(
	  const GeometryType& geom,
	  const Properties& rMaterialProperties,
	  double& rCharacteristicLength)
  {

	  // Updated for MPM - we take the material point volume
	  rCharacteristicLength = 0.0;
	  double area = geom.GetValue(MP_VOLUME);

	  if (rMaterialProperties.Has(THICKNESS)) area /= rMaterialProperties[THICKNESS];
	  else KRATOS_ERROR << "2D ANALYSIS SHOULD HAVE THICKNESS IN MATERIAL PROPERTIES!\n";
	  rCharacteristicLength = std::sqrt(area);

	  KRATOS_ERROR_IF(rCharacteristicLength == 0.0) << "Characteristic length not set properly!\n"
		  << "Geom MP_VOLUME = " << geom.GetValue(MP_VOLUME) << "\n";
  }

} // Namespace Kratos
