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
      void ExecuteInitialize(bool ParticleTempAll,
                             bool ParticleTempMin,
                             bool ParticleTempMax,
                             bool ParticleTempAvg,
                             bool ParticleTempAvgVol,
                             bool ParticleTempDev,
                             bool EnergyMechanical,
                             bool EnergyDissipated,
                             bool EnergyThermal,
                             bool HeatFluxContributions,
                             bool HeatGenValues,
                             bool HeatGenContributions);
      void ExecuteFinalizeSolutionStep(ModelPart& rModelPart, const int write_all_temp_freq);
      void ExecuteFinalize(void);

    protected:

      // Protected attributes
      bool mGraph_ParticleTempAll;
      bool mGraph_ParticleTempMin;
      bool mGraph_ParticleTempMax;
      bool mGraph_ParticleTempAvg;
      bool mGraph_ParticleTempAvgVol;
      bool mGraph_ParticleTempDev;
      bool mGraph_EnergyMechanical;
      bool mGraph_EnergyDissipated;
      bool mGraph_EnergyThermal;
      bool mGraph_HeatFluxContributions;
      bool mGraph_HeatGenValues;
      bool mGraph_HeatGenContributions;

      std::ofstream mFile_ParticleTempAll;
      std::ofstream mFile_ParticleTempMin;        // Minimum particle temperature
      std::ofstream mFile_ParticleTempMax;        // Maximum particle temperature
      std::ofstream mFile_ParticleTempAvg;        // Average particle temperature
      std::ofstream mFile_ParticleTempAvgVol;     // Volume-average particle temperature
      std::ofstream mFile_ParticleTempDev;        // Std dev particle temperature
      std::ofstream mFile_EnergyMechanical;       // Mechanical energy components (current)
      std::ofstream mFile_EnergyDissipated;       // Dissipated energy components (accumulated)
      std::ofstream mFile_EnergyThermal;          // Accumulated thermal energy generation components (U[J]=Q[W]*t[s])
      std::ofstream mFile_HeatFluxContributions;  // Relative contributions of heat flux components
      std::ofstream mFile_HeatGenValues;          // Current values of heat generation components (Q[W])
      std::ofstream mFile_HeatGenContributions;   // Relative contributions of current values heat generation components (%)

    private:

      // Private methods
      void OpenFiles(void);
      void CloseFiles(void);
      void WriteGraphs(ModelPart& rModelPart);

      // Assignment operator
      GraphUtilities& operator=(GraphUtilities const& rOther);

  }; // Class GraphUtilities
} // namespace Kratos

#endif // GRAPH_UTILITIES_H_INCLUDED
