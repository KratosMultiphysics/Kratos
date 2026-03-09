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
#include "nusselt_li_mason.h"

namespace Kratos {
  //------------------------------------------------------------------------------------------------------------
  NusseltLiMason::NusseltLiMason() {}
  NusseltLiMason::~NusseltLiMason() {}

  //------------------------------------------------------------------------------------------------------------
  double NusseltLiMason::ComputeNusseltNumber(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    const double Pr  = particle->ComputePrandtlNumber(r_process_info);
    const double Re  = particle->ComputeReynoldNumber(r_process_info);
    const double por = r_process_info[AVERAGE_POROSITY];

    // Assumption: exponent "m = 4.75" recommended for dense systems (3.50 is recommended for dilute systems)
    const double m = 4.75;

    if (Re < 200.0)
      return 2.0 + 0.6 * pow(por,m) * pow(Re,0.5) * pow(Pr,1.0/3.0);

    else if (Re < 1500.0)
      return 2.0 + pow(por,m) * (0.5 * pow(Re,0.5) + 0.02 * pow(Re,0.8)) * pow(Pr,1.0/3.0);

    else
      return 2.0 + 0.000045 * pow(por,m) * pow(Re,1.8);

    KRATOS_CATCH("")
  }

} // namespace Kratos
