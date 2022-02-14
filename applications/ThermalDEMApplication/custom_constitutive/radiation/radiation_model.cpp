//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics ThermalDEM Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Rafael Rangel (rrangel@cimne.upc.edu)
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
    pProp->SetValue(DEM_RADIATION_MODEL_POINTER, this->CloneShared());
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
