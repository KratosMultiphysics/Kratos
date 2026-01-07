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
#include "direct_conduction_model.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) DirectConductionPipe : public DirectConductionModel
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(DirectConductionPipe);

      // Constructor / Destructor
      DirectConductionPipe();
      virtual ~DirectConductionPipe();

      // Public methods
      double ComputeHeatFlux                     (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;
      double ComputeEffectiveThermalConductivity (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;

      // Clone
      HeatExchangeMechanism* CloneRaw() const override {
        HeatExchangeMechanism* cloned_model(new DirectConductionPipe(*this));
        return cloned_model;
      }

      HeatExchangeMechanism::Pointer CloneShared() const override {
        HeatExchangeMechanism::Pointer cloned_model(new DirectConductionPipe(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "DirectConductionPipe";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "DirectConductionPipe"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      DirectConductionPipe& operator=(DirectConductionPipe const& rOther) {return *this;}
      DirectConductionPipe(DirectConductionPipe const& rOther) {*this = rOther;}

  }; // Class DirectConductionPipe

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    DirectConductionPipe& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const DirectConductionPipe& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos
