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
#include "nusselt_whitaker.h"

namespace Kratos {
  //------------------------------------------------------------------------------------------------------------
  NusseltWhitaker::NusseltWhitaker() {}
  NusseltWhitaker::~NusseltWhitaker() {}

  //------------------------------------------------------------------------------------------------------------
  double NusseltWhitaker::ComputeNusseltNumber(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    const double Pr = particle->ComputePrandtlNumber(r_process_info);
    const double Re = particle->ComputeReynoldNumber(r_process_info);

    // Assumption: temperature-dependent viscosity at particle surface is negleted
    return 2.0 + (0.4 * pow(Re,0.5) + 0.06 * pow(Re,2.0/3.0)) * pow(Pr,0.4);

    KRATOS_CATCH("")
  }

} // namespace Kratos
