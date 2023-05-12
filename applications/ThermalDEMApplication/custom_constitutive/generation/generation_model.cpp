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
#include "generation_model.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  GenerationModel::GenerationModel() {}
  GenerationModel::~GenerationModel() {}

  //------------------------------------------------------------------------------------------------------------
  void GenerationModel::SetHeatGenerationMechanismInProperties(Properties::Pointer pProp, bool verbose) const {
    pProp->SetValue(GENERATION_MODEL_POINTER, this->CloneShared());
  }

  //------------------------------------------------------------------------------------------------------------
  double GenerationModel::ComputeHeatGeneration(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return 0.0;
  }

} // namespace Kratos
