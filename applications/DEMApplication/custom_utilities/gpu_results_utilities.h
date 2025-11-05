#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_elements/rigid_body_element.h"

namespace Kratos
{
  class KRATOS_API(DEM_APPLICATION) GPUResultsUtilities
  {
    public:
      KRATOS_CLASS_POINTER_DEFINITION(GPUResultsUtilities);
      GPUResultsUtilities() {}
      ~GPUResultsUtilities() {}
      void ExecuteInitialize(ModelPart& particlesMP, ModelPart& wallsMP);
      void Calculate(ModelPart& particlesMP, ModelPart& wallsMP);
      void ExecuteFinalize(ModelPart& particlesMP, ModelPart& wallsMP);

      std::ofstream mFile_Global;
      int mNumParticles;
      int mNumWalls;
      double mTimeToWrite;
      double mRadiusMean;
  };
}