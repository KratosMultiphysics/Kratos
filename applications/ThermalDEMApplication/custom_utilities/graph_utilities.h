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
                             bool MechanicalEnergy,
                             bool ThermalEnergy,
                             bool HeatGenValues,
                             bool HeatGenContrib);
      void ExecuteFinalizeSolutionStep(ModelPart& rModelPart);
      void ExecuteFinalize(void);

    protected:

      // Protected attributes
      bool mGraph_ParticleTempMin;
      bool mGraph_ParticleTempMax;
      bool mGraph_ParticleTempAvg;
      bool mGraph_ParticleTempDev;
      bool mGraph_ModelTempAvg;
      bool mGraph_MechanicalEnergy;
      bool mGraph_ThermalEnergy;
      bool mGraph_HeatGenValues;
      bool mGraph_HeatGenContrib;

      std::ofstream mFile_ParticleTempMin;  // Minimum particle temperature
      std::ofstream mFile_ParticleTempMax;  // Maximum particle temperature
      std::ofstream mFile_ParticleTempAvg;  // Average particle temperature
      std::ofstream mFile_ParticleTempDev;  // Std dev particle temperature
      std::ofstream mFile_ModelTempAvg;     // Volume-average particle temperature
      std::ofstream mFile_MechanicalEnergy; // Mechanical energy components (current energies and accumulated dissipations)
      std::ofstream mFile_ThermalEnergy;    // Accumulated thermal energy generation components (U[J]=Q[W]*t[s])
      std::ofstream mFile_HeatGenValues;    // Current values of heat generation components (Q[W])
      std::ofstream mFile_HeatGenContrib;   // Relative contributions of current values heat generation components (%)
      std::ofstream mFile_MassInSilo;       // Relative contributions of current values heat generation components (%)

    private:

      // Assignment operator
      GraphUtilities& operator=(GraphUtilities const& rOther);

  }; // Class GraphUtilities
} // namespace Kratos

#endif // GRAPH_UTILITIES_H_INCLUDED
