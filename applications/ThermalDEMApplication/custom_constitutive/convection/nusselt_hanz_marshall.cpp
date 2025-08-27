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
#include "nusselt_hanz_marshall.h"

namespace Kratos {
  //------------------------------------------------------------------------------------------------------------
  NusseltHanzMarshall::NusseltHanzMarshall() {}
  NusseltHanzMarshall::~NusseltHanzMarshall() {}

  //------------------------------------------------------------------------------------------------------------
  double NusseltHanzMarshall::ComputeNusseltNumber(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    const double Pr = particle->ComputePrandtlNumber(r_process_info);
    const double Re = particle->ComputeReynoldNumber(r_process_info);

    return 2.0 + 0.6 * pow(Re,0.5) * pow(Pr,1.0/3.0);

    KRATOS_CATCH("")
  }

} // namespace Kratos
