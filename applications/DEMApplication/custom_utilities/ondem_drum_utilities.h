#pragma once

#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_elements/rigid_body_element.h"

namespace Kratos
{
  class KRATOS_API(DEM_APPLICATION) ONDEMDrumUtilities
  {
    public:
      KRATOS_CLASS_POINTER_DEFINITION(ONDEMDrumUtilities);
      ONDEMDrumUtilities() {}
      ~ONDEMDrumUtilities() {}
      void ExecuteInitialize(ModelPart& particlesMP, ModelPart& wallsMP);
      void Calculate(ModelPart& particlesMP, ModelPart& wallsMP, bool force_execute);
      void ExecuteCalculations(ModelPart& particlesMP, ModelPart& wallsMP);
      void ExecuteFinalize(ModelPart& particlesMP, ModelPart& wallsMP);
      void InitializeCellsMaps_0100(void);
      void InitializeCellsMaps_0125(void);
      void ResetCellsMap(std::vector<std::vector<std::vector<int>>>& map);

      std::ofstream mFile_Contact;
      std::ofstream mFile_Velocity;
      std::ofstream mFile_Energy;
      int mNumParticlesAll;
      int mNumParticlesLarge;
      int mNumParticlesSmall;
      int mNumWalls;
      int mCellsX_0100;
      int mCellsY_0100;
      int mCellsZ_0100;
      int mCellsX_0125;
      int mCellsY_0125;
      int mCellsZ_0125;
      double mRadiusMean;
      std::vector<std::vector<std::vector<int>>> mCountParticlesLarge_0100;
      std::vector<std::vector<std::vector<int>>> mCountParticlesSmall_0100;
      std::vector<std::vector<std::vector<int>>> mCountParticlesLarge_0125;
      std::vector<std::vector<std::vector<int>>> mCountParticlesSmall_0125;
  };
}