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
#include "generation_model.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) GenerationDissipation : public GenerationModel
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(GenerationDissipation);

      // Constructor / Destructor
      GenerationDissipation();
      virtual ~GenerationDissipation();

      // Public methods
      double ComputeHeatGeneration (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;
      double ComputePartitionCoeff (ThermalSphericParticle* particle);
      void   FillHeatMap (const ProcessInfo& r_process_info,
                          ThermalSphericParticle* particle,
                          const double time,
                          const double heat_gen_damping_pp,
                          const double heat_gen_damping_pw,
                          const double heat_gen_sliding_pp,
                          const double heat_gen_sliding_pw,
                          const double heat_gen_rolling_pp,
                          const double heat_gen_rolling_pw);

      // Clone
      HeatGenerationMechanism* CloneRaw() const override {
        HeatGenerationMechanism* cloned_model(new GenerationDissipation(*this));
        return cloned_model;
      }

      HeatGenerationMechanism::Pointer CloneShared() const override {
        HeatGenerationMechanism::Pointer cloned_model(new GenerationDissipation(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "GenerationDissipation";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "GenerationDissipation"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      GenerationDissipation& operator=(GenerationDissipation const& rOther) {return *this;}
      GenerationDissipation(GenerationDissipation const& rOther) {*this = rOther;}

  }; // Class GenerationDissipation

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    GenerationDissipation& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const GenerationDissipation& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos
