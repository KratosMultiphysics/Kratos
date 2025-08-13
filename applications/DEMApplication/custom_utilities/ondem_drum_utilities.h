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
      void Calculate(ModelPart& particlesMP, ModelPart& wallsMP);
      void ExecuteFinalize(ModelPart& particlesMP, ModelPart& wallsMP);
      void InitializeCellsMap(std::vector<std::vector<std::vector<int>>>& map);
      void ResetCellsMap(std::vector<std::vector<std::vector<int>>>& map);

      std::ofstream mFile_Contact;
      std::ofstream mFile_Velocity;
      std::ofstream mFile_Energy;
      int mNumParticlesAll;
      int mNumParticlesLarge;
      int mNumParticlesSmall;
      int mNumWalls;
      int mCellsX;
      int mCellsY;
      int mCellsZ;
      double mRadiusMean;
      std::vector<std::vector<std::vector<int>>> mCountParticlesLarge;
      std::vector<std::vector<std::vector<int>>> mCountParticlesSmall;
  };
}