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
#include "indirect_conduction_model.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) IndirectConductionVargas : public IndirectConductionModel
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(IndirectConductionVargas);

      // Constructor / Destructor
      IndirectConductionVargas();
      virtual ~IndirectConductionVargas();

      // Public methods
      double ComputeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;

      // Clone
      HeatExchangeMechanism* CloneRaw() const override {
        HeatExchangeMechanism* cloned_model(new IndirectConductionVargas(*this));
        return cloned_model;
      }

      HeatExchangeMechanism::Pointer CloneShared() const override {
        HeatExchangeMechanism::Pointer cloned_model(new IndirectConductionVargas(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "IndirectConductionVargas";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "IndirectConductionVargas"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      IndirectConductionVargas& operator=(IndirectConductionVargas const& rOther) {return *this;}
      IndirectConductionVargas(IndirectConductionVargas const& rOther) {*this = rOther;}

  }; // Class IndirectConductionVargas

    // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    IndirectConductionVargas& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const IndirectConductionVargas& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos
