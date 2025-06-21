//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#pragma once

// System includes

// External includes

// Project includes
#include "custom_constitutive/heat_exchange_mechanism.h"
#include "custom_elements/thermal_spheric_particle.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) DirectConductionModel : public HeatExchangeMechanism
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(DirectConductionModel);

      // Constructor / Destructor
      DirectConductionModel();
      virtual ~DirectConductionModel();

      // Public methods
      void   SetHeatExchangeMechanismInProperties (Properties::Pointer pProp, bool verbose = true) const override;
      double GetSearchDistance                    (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;
      double ComputeHeatFlux                      (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;

      // Clone
      HeatExchangeMechanism* CloneRaw() const override {
        HeatExchangeMechanism* cloned_model(new DirectConductionModel(*this));
        return cloned_model;
      }

      HeatExchangeMechanism::Pointer CloneShared() const override {
        HeatExchangeMechanism::Pointer cloned_model(new DirectConductionModel(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "DirectConductionModel";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "DirectConductionModel"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      DirectConductionModel& operator=(DirectConductionModel const& rOther) {return *this;}
      DirectConductionModel(DirectConductionModel const& rOther) {*this = rOther;}

  }; // Class DirectConductionModel

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    DirectConductionModel& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const DirectConductionModel& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos
