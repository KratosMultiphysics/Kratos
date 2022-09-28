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
#include "radiation_model.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  RadiationModel::RadiationModel() {}
  RadiationModel::~RadiationModel() {}

  //------------------------------------------------------------------------------------------------------------
  void RadiationModel::SetHeatExchangeMechanismInProperties(Properties::Pointer pProp, bool verbose) const {
    pProp->SetValue(RADIATION_MODEL_POINTER, this->CloneShared());
  }

  //------------------------------------------------------------------------------------------------------------
  double RadiationModel::GetSearchDistance(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return particle->GetParticleRadius() * (r_process_info[MAX_RADIATION_DISTANCE]);
  }

  //------------------------------------------------------------------------------------------------------------
  double RadiationModel::ComputeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double RadiationModel::FinalizeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  void RadiationModel::AccumulateEnvironmentTemperature(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {}

} // namespace Kratos
