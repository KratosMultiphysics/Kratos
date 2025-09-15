//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#pragma once

// System includes

// External includes
#include "includes/model_part.h"

// Project includes
#include "custom_elements/thermal_spheric_particle.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) HeatMapUtilities
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(HeatMapUtilities);

      // Constructor / destructor methods
      HeatMapUtilities();
      ~HeatMapUtilities();

      // Public methods
      void ExecuteInitialize           (ModelPart& rModelPart);
      void ExecuteFinalizeSolutionStep (ModelPart& rModelPart);
      void ExecuteFinalize             (ModelPart& rModelPart);

      // Public variables
      int mDimX, mDimY, mDimZ;  // Dimensions of heat map matrices
      std::vector<std::vector<std::vector<double>>> mGlobalHeatMapGenerationDampingPP;  // Global heat map matrix for heat generaion by damping between particle-particle
      std::vector<std::vector<std::vector<double>>> mGlobalHeatMapGenerationDampingPW;  // Global heat map matrix for heat generaion by damping between particle-wall
      std::vector<std::vector<std::vector<double>>> mGlobalHeatMapGenerationSlidingPP;  // Global heat map matrix for heat generaion by sliding between particle-particle
      std::vector<std::vector<std::vector<double>>> mGlobalHeatMapGenerationSlidingPW;  // Global heat map matrix for heat generaion by sliding between particle-wall
      std::vector<std::vector<std::vector<double>>> mGlobalHeatMapGenerationRollingPP;  // Global heat map matrix for heat generaion by rolling between particle-particle
      std::vector<std::vector<std::vector<double>>> mGlobalHeatMapGenerationRollingPW;  // Global heat map matrix for heat generaion by rolling between particle-wall

    private:

      // Private methods
      void ResetMap(std::vector<std::vector<std::vector<double>>>& map);

      // Assignment operator
      HeatMapUtilities& operator=(HeatMapUtilities const& rOther);

  }; // Class HeatMapUtilities
} // namespace Kratos
