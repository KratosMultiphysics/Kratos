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
#include "real_contact_model.h"
#include "thermal_dem_application_variables.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  RealContactModel::RealContactModel() {}
  RealContactModel::~RealContactModel() {}

  //------------------------------------------------------------------------------------------------------------
  void RealContactModel::SetRealContactModelInProperties(Properties::Pointer pProp, bool verbose) const {
    pProp->SetValue(REAL_CONTACT_MODEL_POINTER, this->CloneShared());
  }

  //------------------------------------------------------------------------------------------------------------
  void RealContactModel::AdjustContact(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {}

} // namespace Kratos
