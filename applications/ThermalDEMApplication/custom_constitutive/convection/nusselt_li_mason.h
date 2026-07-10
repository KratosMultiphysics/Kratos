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
#include "convection_model.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) NusseltLiMason : public ConvectionModel
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(NusseltLiMason);

      // Constructor / Destructor
      NusseltLiMason();
      virtual ~NusseltLiMason();

      // Public methods
      double ComputeNusseltNumber(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;

      // Clone
      HeatExchangeMechanism* CloneRaw() const override {
        HeatExchangeMechanism* cloned_model(new NusseltLiMason(*this));
        return cloned_model;
      }

      HeatExchangeMechanism::Pointer CloneShared() const override {
        HeatExchangeMechanism::Pointer cloned_model(new NusseltLiMason(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "NusseltLiMason";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "NusseltLiMason"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      NusseltLiMason& operator=(NusseltLiMason const& rOther) {return *this;}
      NusseltLiMason(NusseltLiMason const& rOther) {*this = rOther;}

  }; // Class NusseltLiMason

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    NusseltLiMason& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const NusseltLiMason& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos
