//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/johnson_cook_thermal_plastic_axisym_2D_law.hpp"

namespace Kratos
{
    JohnsonCookThermalPlastic2DAxisymLaw::JohnsonCookThermalPlastic2DAxisymLaw()
  : JohnsonCookThermalPlastic2DPlaneStrainLaw()
  { }

    JohnsonCookThermalPlastic2DAxisymLaw::JohnsonCookThermalPlastic2DAxisymLaw(const JohnsonCookThermalPlastic2DAxisymLaw& rOther)
  : JohnsonCookThermalPlastic2DPlaneStrainLaw(rOther)
  { }

  ConstitutiveLaw::Pointer JohnsonCookThermalPlastic2DAxisymLaw::Clone() const
  {
    return Kratos::make_shared<JohnsonCookThermalPlastic2DAxisymLaw>(*this);
  }

  JohnsonCookThermalPlastic2DAxisymLaw::~JohnsonCookThermalPlastic2DAxisymLaw()
  { }

  void JohnsonCookThermalPlastic2DAxisymLaw::MakeStrainStressVectorFromMatrix(const Matrix& rInput, Vector& rOutput)
  {
	  if (rOutput.size() != GetStrainSize()) rOutput.resize(GetStrainSize(), false);

	  // 2D stress arrangement
	  rOutput[0] = rInput(0, 0);
	  rOutput[1] = rInput(1, 1);
	  rOutput[2] = rInput(2, 2);
	  rOutput[3] = 2.0 * rInput(0, 1); //xy
  }


  void JohnsonCookThermalPlastic2DAxisymLaw::MakeStrainStressMatrixFromVector(const Vector& rInput, Matrix& rOutput)
  {
	  if (rOutput.size1() != 3 || rOutput.size2() != 3)rOutput.resize(3, 3, false);

	  // 3D stress arrangement
      rOutput.clear();

      // Normal components
	  rOutput(0, 0) = rInput[0];
	  rOutput(1, 1) = rInput[1];
	  rOutput(2, 2) = rInput[2];

      // Shear components
	  rOutput(0, 1) = 0.5 * rInput[3]; //xy

      // Fill symmetry
	  rOutput(1, 0) = rOutput(0, 1);
  }

} // Namespace Kratos
