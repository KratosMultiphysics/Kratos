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
#include "thermal_dem_integration_scheme.h"
#include "thermal_dem_application_variables.h"

namespace Kratos {

  //------------------------------------------------------------------------------------------------------------
  ThermalDEMIntegrationScheme::ThermalDEMIntegrationScheme() {}
  ThermalDEMIntegrationScheme::~ThermalDEMIntegrationScheme() {}

  //------------------------------------------------------------------------------------------------------------
  void ThermalDEMIntegrationScheme::SetThermalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
    pProp->SetValue(THERMAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalDEMIntegrationScheme::UpdateTemperature(Node& i, const double delta_t, const double c) {}

} // namespace Kratos
