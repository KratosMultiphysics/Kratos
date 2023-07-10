//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

// System includes

// External includes
#include "DEM_application_variables.h"

// Project includes
#include "thermal_forward_euler_scheme.h"
#include "thermal_dem_application_variables.h"

namespace Kratos {

  //------------------------------------------------------------------------------------------------------------
  ThermalForwardEulerScheme::ThermalForwardEulerScheme() {}
  ThermalForwardEulerScheme::~ThermalForwardEulerScheme() {}

  //------------------------------------------------------------------------------------------------------------
  void ThermalForwardEulerScheme::SetThermalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
    pProp->SetValue(THERMAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalForwardEulerScheme::UpdateTemperature(Node& i, const double delta_t, const double c) {
    // Particle properties
    const double q = i.FastGetSolutionStepValue(HEATFLUX);
    const double m = i.FastGetSolutionStepValue(NODAL_MASS);
    const double t = i.FastGetSolutionStepValue(TEMPERATURE);

    // Compute new temperature
    const double temp_incr = delta_t * q / (m * c);
    const double temp_new  = t + temp_incr;

    // Set new temperature
    i.FastGetSolutionStepValue(TEMPERATURE) = temp_new;
  }

} // namespace Kratos
