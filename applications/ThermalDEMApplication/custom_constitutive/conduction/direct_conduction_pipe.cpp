//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \
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
#include "direct_conduction_pipe.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  DirectConductionPipe::DirectConductionPipe() {}
  DirectConductionPipe::~DirectConductionPipe() {}

  //------------------------------------------------------------------------------------------------------------
  double DirectConductionPipe::ComputeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Check for contact
    if (!particle->mNeighborInContact)
      return 0.0;

    const double Rc        = particle->mContactRadiusAdjusted;
    const double d         = particle->mNeighborDistanceAdjusted;
    const double kavg      = particle->ComputeAverageConductivity();
    const double temp_grad = particle->GetNeighborTemperature() - particle->GetParticleTemperature();

    return kavg * (Globals::Pi * Rc * Rc) * temp_grad / d;

    KRATOS_CATCH("")
  }

} // namespace Kratos
