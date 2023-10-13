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
#include "direct_conduction_model.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  DirectConductionModel::DirectConductionModel() {}
  DirectConductionModel::~DirectConductionModel() {}

  //------------------------------------------------------------------------------------------------------------
  void DirectConductionModel::SetHeatExchangeMechanismInProperties(Properties::Pointer pProp, bool verbose) const {
    pProp->SetValue(DIRECT_CONDUCTION_MODEL_POINTER, this->CloneShared());
  }

  //------------------------------------------------------------------------------------------------------------
  double DirectConductionModel::GetSearchDistance(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double DirectConductionModel::ComputeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return 0.0;
  }

} // namespace Kratos
