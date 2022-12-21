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
#include "indirect_conduction_model.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  IndirectConductionModel::IndirectConductionModel() {}
  IndirectConductionModel::~IndirectConductionModel() {}

  //------------------------------------------------------------------------------------------------------------
  void IndirectConductionModel::SetHeatExchangeMechanismInProperties(Properties::Pointer pProp, bool verbose) const {
    pProp->SetValue(INDIRECT_CONDUCTION_MODEL_POINTER, this->CloneShared());
  }

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionModel::GetSearchDistance(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionModel::ComputeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return 0.0;
  }

} // namespace Kratos
