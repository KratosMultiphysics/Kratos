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
  class KRATOS_API(THERMAL_DEM_APPLICATION) IndirectConductionVoronoiB : public IndirectConductionModel
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(IndirectConductionVoronoiB);

      // Constructor / Destructor
      IndirectConductionVoronoiB();
      virtual ~IndirectConductionVoronoiB();

      // Public methods
      double GetSearchDistance          (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;
      double ComputeHeatFlux            (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;
      double SphereWallCoeff            (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);
      double SphereSphereMonoSizeCoeff  (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);
      double SphereSphereMultiSizeCoeff (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);

      // Clone
      HeatExchangeMechanism* CloneRaw() const override {
        HeatExchangeMechanism* cloned_model(new IndirectConductionVoronoiB(*this));
        return cloned_model;
      }

      HeatExchangeMechanism::Pointer CloneShared() const override {
        HeatExchangeMechanism::Pointer cloned_model(new IndirectConductionVoronoiB(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "IndirectConductionVoronoiB";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream & rOStream) const override { rOStream << "IndirectConductionVoronoiB"; }
      virtual void PrintData(std::ostream & rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      IndirectConductionVoronoiB& operator=(IndirectConductionVoronoiB const& rOther) {return *this;}
      IndirectConductionVoronoiB(IndirectConductionVoronoiB const& rOther) {*this = rOther;}

  }; // Class IndirectConductionVoronoiB

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    IndirectConductionVoronoiB& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const IndirectConductionVoronoiB& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos
