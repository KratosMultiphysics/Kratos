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
#include "heat_exchange_mechanism.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  HeatExchangeMechanism::HeatExchangeMechanism() {}
  HeatExchangeMechanism::~HeatExchangeMechanism() {}

  //------------------------------------------------------------------------------------------------------------
  void HeatExchangeMechanism::SetHeatExchangeMechanismInProperties(Properties::Pointer pProp, bool verbose) const {}

  //------------------------------------------------------------------------------------------------------------
  double HeatExchangeMechanism::GetSearchDistance(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double HeatExchangeMechanism::ComputeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double HeatExchangeMechanism::FinalizeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return 0.0;
  }

} // namespace Kratos
