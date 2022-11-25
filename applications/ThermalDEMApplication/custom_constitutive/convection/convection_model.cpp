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
#include "convection_model.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  ConvectionModel::ConvectionModel() {}
  ConvectionModel::~ConvectionModel() {}

  //------------------------------------------------------------------------------------------------------------
  void ConvectionModel::SetHeatExchangeMechanismInProperties(Properties::Pointer pProp, bool verbose) const {
    pProp->SetValue(CONVECTION_MODEL_POINTER, this->CloneShared());
  }

  //------------------------------------------------------------------------------------------------------------
  double ConvectionModel::ComputeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    const double surface_area       = particle->GetParticleSurfaceArea();
    const double char_length        = particle->GetParticleCharacteristicLength();
    const double fluid_conductivity = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    const double temp_grad          = r_process_info[FLUID_TEMPERATURE] - particle->GetParticleTemperature();

    // Compute Nusselt number
    const double Nu = ComputeNusseltNumber(r_process_info, particle);

    // Compute heat flux
    return (Nu * fluid_conductivity / char_length) * surface_area * temp_grad;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ConvectionModel::ComputeNusseltNumber(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return 0.0;
  }

} // namespace Kratos
