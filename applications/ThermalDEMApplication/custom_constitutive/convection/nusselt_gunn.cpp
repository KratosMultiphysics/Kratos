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
#include "nusselt_gunn.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  NusseltGunn::NusseltGunn() {}
  NusseltGunn::~NusseltGunn() {}

  //------------------------------------------------------------------------------------------------------------
  double NusseltGunn::ComputeNusseltNumber(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    const double Pr  = particle->ComputePrandtlNumber(r_process_info);
    const double Re  = particle->ComputeReynoldNumber(r_process_info);
    const double por = r_process_info[AVERAGE_POROSITY];
    
    return (7.0 - 10.0 * por + 5.0 * por * por) * (1.0 + 0.7 * pow(Re,0.2) * pow(Pr,1.0/3.0)) + (1.33 - 2.4 * por + 1.2 * por * por) * pow(Re,0.7) * pow(Pr,1.0/3.0);

    KRATOS_CATCH("")
  }

} // namespace Kratos
