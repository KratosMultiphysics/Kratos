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

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) ThermalDEMIntegrationScheme
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(ThermalDEMIntegrationScheme);

      // Constructor / Destructor
      ThermalDEMIntegrationScheme();
      virtual ~ThermalDEMIntegrationScheme();

      // Public methods
      virtual void SetThermalIntegrationSchemeInProperties (Properties::Pointer pProp, bool verbose = true) const;
      virtual void UpdateTemperature                       (Node& i, const double delta_t, const double c);

      // Clone
      virtual ThermalDEMIntegrationScheme* CloneRaw() const {
        ThermalDEMIntegrationScheme* cloned_scheme(new ThermalDEMIntegrationScheme(*this));
        return cloned_scheme;
      }

      virtual ThermalDEMIntegrationScheme::Pointer CloneShared() const {
        ThermalDEMIntegrationScheme::Pointer cloned_scheme(new ThermalDEMIntegrationScheme(*this));
        return cloned_scheme;
      }

      // Turn back information as a string
      virtual std::string Info() const {
        std::stringstream buffer;
        buffer << "ThermalDEMIntegrationScheme";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const { rOStream << "ThermalDEMIntegrationScheme"; }
      virtual void PrintData(std::ostream& rOStream) const {}

    private:

      // Assignment operator / Copy constructor
      ThermalDEMIntegrationScheme& operator=(ThermalDEMIntegrationScheme const& rOther) {return *this;}
      ThermalDEMIntegrationScheme(ThermalDEMIntegrationScheme const& rOther) {*this = rOther;}

      // Serializer
      friend class Serializer;

      virtual void save(Serializer& rSerializer) const {
        //rSerializer.save("MyMemberName",myMember);
      }

      virtual void load(Serializer& rSerializer) {
        //rSerializer.load("MyMemberName",myMember);
      }

  }; // Class ThermalDEMIntegrationScheme

  // input stream function
  inline std::istream& operator >> (std::istream& rIStream, ThermalDEMIntegrationScheme& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator << (std::ostream& rOStream, const ThermalDEMIntegrationScheme& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos
