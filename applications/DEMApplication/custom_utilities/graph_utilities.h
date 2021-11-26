//
// Author:  Rafael Rangel, rrangel@cimne.upc.edu
// Date:    November 2021
//

#ifndef GRAPH_UTILITIES_H
#define	GRAPH_UTILITIES_H

// System includes
#include <iostream>

// Project includes
#include "custom_elements/thermal_spheric_particle.h"

// External includes

namespace Kratos {

  class KRATOS_API(DEM_APPLICATION) GraphUtilities {

  public:

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
                           bool ParticleHeatFluxContributions);
    void ExecuteFinalizeSolutionStep(ModelPart& r_modelpart);
    void ExecuteFinalize();

  protected:
    // Protected attributes
    bool mGraph_ParticleTempMin;
    bool mGraph_ParticleTempMax;
    bool mGraph_ParticleTempAvg;
    bool mGraph_ParticleTempDev;
    bool mGraph_ModelTempAvg;
    bool mGraph_ParticleHeatFluxContributions;

    std::ofstream mFile_ParticleTempMin;
    std::ofstream mFile_ParticleTempMax;
    std::ofstream mFile_ParticleTempAvg;
    std::ofstream mFile_ParticleTempDev;
    std::ofstream mFile_ModelTempAvg;
    std::ofstream mFile_ParticleHeatFluxContributions;

  private:
    // Assignment operator
    GraphUtilities& operator=(GraphUtilities const& rOther);
  };

} // namespace Kratos

#endif  // GRAPH_UTILITIES_H