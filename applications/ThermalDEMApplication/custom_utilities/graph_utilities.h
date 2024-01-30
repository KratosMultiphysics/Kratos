//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#ifndef GRAPH_UTILITIES_H_INCLUDED
#define	GRAPH_UTILITIES_H_INCLUDED

// System includes
#include <iostream>

// External includes
#include "includes/model_part.h"

// Project includes
#include "custom_elements/thermal_spheric_particle.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) GraphUtilities
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(GraphUtilities);

      // Constructor / destructor methods
      GraphUtilities();
      ~GraphUtilities();

      // Public methods
      void ExecuteInitialize(bool ParticleTempMin,
                             bool ParticleTempMax,
                             bool ParticleTempAvg,
                             bool ParticleTempDev,
                             bool ModelTempAvg,
                             bool ParticleHeatFluxContributions,
                             bool ParticleHeatGenContributions,
                             bool ParticleEnergyContributions);
      void ExecuteFinalizeSolutionStep(ModelPart& rModelPart);
      void ExecuteFinalize(void);

    protected:

      // Protected attributes
      bool mGraph_ParticleTempMin;
      bool mGraph_ParticleTempMax;
      bool mGraph_ParticleTempAvg;
      bool mGraph_ParticleTempDev;
      bool mGraph_ModelTempAvg;
      bool mGraph_ParticleHeatFluxContributions;
      bool mGraph_ParticleHeatGenContributions;
      bool mGraph_ParticleEnergyContributions;

      std::ofstream mFile_ParticleTempMin;
      std::ofstream mFile_ParticleTempMax;
      std::ofstream mFile_ParticleTempAvg;
      std::ofstream mFile_ParticleTempDev;
      std::ofstream mFile_ModelTempAvg;
      std::ofstream mFile_ParticleHeatFluxContributions;
      std::ofstream mFile_ParticleHeatGenContributions;
      std::ofstream mFile_ParticleEnergyContributions;

    private:

      // Assignment operator
      GraphUtilities& operator=(GraphUtilities const& rOther);

  }; // Class GraphUtilities
} // namespace Kratos

#endif // GRAPH_UTILITIES_H_INCLUDED
