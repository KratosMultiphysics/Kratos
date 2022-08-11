//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(FRICTION_MODEL_H_INCLUDED)
#define FRICTION_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/heat_generation_mechanism.h"
#include "custom_elements/thermal_spheric_particle.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) FrictionModel : public HeatGenerationMechanism
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(FrictionModel);

      // Constructor / Destructor
      FrictionModel();
      virtual ~FrictionModel();

      // Public methods
      void           SetHeatGenerationMechanismInProperties (Properties::Pointer pProp, bool verbose = true) const override;
      double         ComputeHeatGeneration                  (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;
      virtual double ComputePartitionCoeff                  (ThermalSphericParticle* particle);

      // Clone
      HeatGenerationMechanism* CloneRaw() const override {
        HeatGenerationMechanism* cloned_model(new FrictionModel(*this));
        return cloned_model;
      }

      HeatGenerationMechanism::Pointer CloneShared() const override {
        HeatGenerationMechanism::Pointer cloned_model(new FrictionModel(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "FrictionModel";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "FrictionModel"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      FrictionModel& operator=(FrictionModel const& rOther) {return *this;}
      FrictionModel(FrictionModel const& rOther) {*this = rOther;}

  }; // Class FrictionModel

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    FrictionModel& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const FrictionModel& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos

#endif // FRICTION_MODEL_H_INCLUDED
