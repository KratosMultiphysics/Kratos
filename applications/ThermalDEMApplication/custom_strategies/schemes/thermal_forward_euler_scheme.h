//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes
#include "includes/define.h"
#include "includes/model_part.h"

// Project includes
#include "thermal_dem_integration_scheme.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) ThermalForwardEulerScheme : public ThermalDEMIntegrationScheme
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(ThermalForwardEulerScheme);

      // Constructor / Destructor
      ThermalForwardEulerScheme();
      virtual ~ThermalForwardEulerScheme();

      // Public methods
      void SetThermalIntegrationSchemeInProperties (Properties::Pointer pProp, bool verbose = true) const override;
      void UpdateTemperature                       (Node& i, const double delta_t, const double c) override;

      // Clone
      ThermalDEMIntegrationScheme* CloneRaw() const override {
        ThermalDEMIntegrationScheme* cloned_scheme(new ThermalForwardEulerScheme(*this));
        return cloned_scheme;
      }

      ThermalDEMIntegrationScheme::Pointer CloneShared() const override {
        ThermalDEMIntegrationScheme::Pointer cloned_scheme(new ThermalForwardEulerScheme(*this));
        return cloned_scheme;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "ThermalForwardEulerScheme";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "ThermalForwardEulerScheme"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      ThermalForwardEulerScheme& operator=(ThermalForwardEulerScheme const& rOther) {return *this;}
      ThermalForwardEulerScheme(ThermalForwardEulerScheme const& rOther) {*this = rOther;}

  }; // Class ThermalForwardEulerScheme

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    ThermalForwardEulerScheme& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const ThermalForwardEulerScheme& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos
