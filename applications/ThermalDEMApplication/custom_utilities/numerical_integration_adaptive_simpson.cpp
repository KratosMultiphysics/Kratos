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
#include "numerical_integration_adaptive_simpson.h"
#include "thermal_dem_application_variables.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  AdaptiveSimpsonQuadrature::AdaptiveSimpsonQuadrature():NumericalIntegrationMethod() {}
  AdaptiveSimpsonQuadrature::~AdaptiveSimpsonQuadrature() {}

  //------------------------------------------------------------------------------------------------------------
  void AdaptiveSimpsonQuadrature::SetNumericalIntegrationMethodInProperties(Properties::Pointer pProp, bool verbose) const {
    pProp->SetValue(NUMERICAL_INTEGRATION_METHOD_POINTER, this->CloneShared());
  }

  //------------------------------------------------------------------------------------------------------------
  double AdaptiveSimpsonQuadrature::SolveIntegral(void) {
    KRATOS_TRY

    // Initialization
    mCoord = mLimMin;
    double fa = mpEvalIntegrand(this);

    mCoord = mLimMax;
    double fb = mpEvalIntegrand(this);

    mCoord = (mLimMin + mLimMax) / 2.0;
    double fc = mpEvalIntegrand(this);

    // Get tolerance
    constexpr double eps = std::numeric_limits<double>::epsilon();
    if (mTol < 10.0 * eps)
      mTol = 10.0 * eps;

    // Solve integral recursively with adaptive Simpson quadrature
    return RecursiveIntegration(mLimMin, mLimMax, fa, fb, fc);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double AdaptiveSimpsonQuadrature::RecursiveIntegration(double a, double b, double fa, double fb, double fc) {
    KRATOS_TRY

    // TODO: in order to catch possible erros that can occur in singularities,
    //       add a min value for subdivision size (to contain machine representable point) and a max number of function evaluation (+- 10000).

    double c = (a + b) / 2.0;

    mCoord = (a + c) / 2.0;
    double fd = mpEvalIntegrand(this);

    mCoord = (c + b) / 2.0;
    double fe = mpEvalIntegrand(this);

    double I1 = (b - a) / 6.0  * (fa + 4.0 * fc + fb);
    double I2 = (b - a) / 12.0 * (fa + 4.0 * fd + 2.0 * fc + 4.0 * fe + fb);

    if (std::abs(I2 - I1) <= mTol) {
      return I2 + (I2 - I1) / 15.0;
    }
    else { // sub-divide interval recursively
      double Ia = RecursiveIntegration(a, c, fa, fc, fd);
      double Ib = RecursiveIntegration(c, b, fc, fb, fe);
      return Ia + Ib;
    }

    KRATOS_CATCH("")
  }

} // namespace Kratos
