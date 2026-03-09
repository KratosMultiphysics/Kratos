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
#include "numerical_integration_method.h"
#include "thermal_dem_application_variables.h"

namespace Kratos {
  //------------------------------------------------------------------------------------------------------------
  NumericalIntegrationMethod::NumericalIntegrationMethod() {
    CleanParameters();
  }

  NumericalIntegrationMethod::~NumericalIntegrationMethod() {}

  //------------------------------------------------------------------------------------------------------------
  void NumericalIntegrationMethod::SetNumericalIntegrationMethodInProperties(Properties::Pointer pProp, bool verbose) const {
    pProp->SetValue(NUMERICAL_INTEGRATION_METHOD_POINTER, this->CloneShared());
  }

  //------------------------------------------------------------------------------------------------------------
  double NumericalIntegrationMethod::SolveIntegral(void) {
    return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  void NumericalIntegrationMethod::CleanParameters(void) {
    mLimMin     = 0.0;
    mLimMax     = 0.0;
    mCoord      = 0.0;
    mTol        = 0.000001;
    mParams.p1  = 0.0;
    mParams.p2  = 0.0;
    mParams.p3  = 0.0;
    mParams.p4  = 0.0;
    mParams.p5  = 0.0;
    mParams.p6  = 0.0;
    mParams.p7  = 0.0;
    mParams.p8  = 0.0;
    mParams.p9  = 0.0;
    mParams.p10 = 0.0;
  }

} // namespace Kratos
