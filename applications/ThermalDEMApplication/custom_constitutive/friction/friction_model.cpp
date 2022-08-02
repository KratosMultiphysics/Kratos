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
#include "friction_model.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  FrictionModel::FrictionModel() {}
  FrictionModel::~FrictionModel() {}

  //------------------------------------------------------------------------------------------------------------
  void FrictionModel::SetHeatGenerationMechanismInProperties(Properties::Pointer pProp, bool verbose) const {
    pProp->SetValue(FRICTION_MODEL_POINTER, this->CloneShared());
  }

  //------------------------------------------------------------------------------------------------------------
  double FrictionModel::ComputeHeatGeneration(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double FrictionModel::ComputePartitionCoeff(ThermalSphericParticle* particle) {
    return 0.0;
  }

} // namespace Kratos
