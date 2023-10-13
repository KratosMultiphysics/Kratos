//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(GENERATION_MODEL_H_INCLUDED)
#define GENERATION_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/heat_generation_mechanism.h"
#include "custom_elements/thermal_spheric_particle.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) GenerationModel : public HeatGenerationMechanism
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(GenerationModel);

      // Constructor / Destructor
      GenerationModel();
      virtual ~GenerationModel();

      // Public methods
      void           SetHeatGenerationMechanismInProperties (Properties::Pointer pProp, bool verbose = true) const override;
      double         ComputeHeatGeneration                  (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;

      // Clone
      HeatGenerationMechanism* CloneRaw() const override {
        HeatGenerationMechanism* cloned_model(new GenerationModel(*this));
        return cloned_model;
      }

      HeatGenerationMechanism::Pointer CloneShared() const override {
        HeatGenerationMechanism::Pointer cloned_model(new GenerationModel(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "GenerationModel";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "GenerationModel"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      GenerationModel& operator=(GenerationModel const& rOther) {return *this;}
      GenerationModel(GenerationModel const& rOther) {*this = rOther;}

  }; // Class GenerationModel

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    GenerationModel& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const GenerationModel& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos

#endif // GENERATION_MODEL_H_INCLUDED
