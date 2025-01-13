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
  // Forward declarations
  class ThermalSphericParticle;

  class KRATOS_API(THERMAL_DEM_APPLICATION) HeatExchangeMechanism
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(HeatExchangeMechanism);

      // Constructor / Destructor
      HeatExchangeMechanism();
      virtual ~HeatExchangeMechanism();

      // Public methods
      virtual void   SetHeatExchangeMechanismInProperties (Properties::Pointer pProp, bool verbose = true) const;
      virtual double GetSearchDistance                    (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);
      virtual double ComputeHeatFlux                      (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);
      virtual double FinalizeHeatFlux                     (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);
      virtual double ComputeEffectiveThermalConductivity  (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);

      // Clone
      virtual HeatExchangeMechanism* CloneRaw() const {
        HeatExchangeMechanism* cloned_mechanism(new HeatExchangeMechanism(*this));
        return cloned_mechanism;
      }

      virtual HeatExchangeMechanism::Pointer CloneShared() const {
        HeatExchangeMechanism::Pointer cloned_mechanism(new HeatExchangeMechanism(*this));
        return cloned_mechanism;
      }

      // Turn back information as a string
      virtual std::string Info() const {
        std::stringstream buffer;
        buffer << "HeatExchangeMechanism";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const { rOStream << "HeatExchangeMechanism"; }
      virtual void PrintData(std::ostream& rOStream) const {}

    private:

      // Assignment operator / Copy constructor
      HeatExchangeMechanism& operator=(HeatExchangeMechanism const& rOther) {return *this;}
      HeatExchangeMechanism(HeatExchangeMechanism const& rOther) {*this = rOther;}

      // Serializer
      friend class Serializer;

      virtual void save(Serializer& rSerializer) const {
        //rSerializer.save("MyMemberName",myMember);
      }

      virtual void load(Serializer& rSerializer) {
        //rSerializer.load("MyMemberName",myMember);
      }

  }; // Class HeatExchangeMechanism

  // input stream function
  inline std::istream& operator >> (std::istream& rIStream, HeatExchangeMechanism& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator << (std::ostream& rOStream, const HeatExchangeMechanism& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos
