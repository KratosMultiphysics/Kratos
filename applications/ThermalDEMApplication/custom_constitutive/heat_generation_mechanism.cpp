//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

// System includes

// External includes

// Project includes
#include "heat_generation_mechanism.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  HeatGenerationMechanism::HeatGenerationMechanism() {}
  HeatGenerationMechanism::~HeatGenerationMechanism() {}

  //------------------------------------------------------------------------------------------------------------
  void HeatGenerationMechanism::SetHeatGenerationMechanismInProperties(Properties::Pointer pProp, bool verbose) const {}

  //------------------------------------------------------------------------------------------------------------
  double HeatGenerationMechanism::ComputeHeatGeneration(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double HeatGenerationMechanism::FinalizeHeatGeneration(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return 0.0;
  }

} // namespace Kratos
