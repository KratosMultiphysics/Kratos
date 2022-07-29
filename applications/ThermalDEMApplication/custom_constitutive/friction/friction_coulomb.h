//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(FRICTION_MODEL_COULOMB_H_INCLUDED)
#define FRICTION_MODEL_COULOMB_H_INCLUDED

// System includes

// External includes

// Project includes
#include "friction_model.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) FrictionCoulomb : public FrictionModel
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(FrictionCoulomb);

      // Constructor / Destructor
      FrictionCoulomb();
      virtual ~FrictionCoulomb();

      // Public methods
      double ComputeHeatGeneration (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;
      double ComputePartitionCoeff (ThermalSphericParticle* particle) override;

      // Clone
      HeatGenerationMechanism* CloneRaw() const override {
        HeatGenerationMechanism* cloned_model(new FrictionCoulomb(*this));
        return cloned_model;
      }

      HeatGenerationMechanism::Pointer CloneShared() const override {
        HeatGenerationMechanism::Pointer cloned_model(new FrictionCoulomb(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "FrictionCoulomb";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "FrictionCoulomb"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      FrictionCoulomb& operator=(FrictionCoulomb const& rOther) {return *this;}
      FrictionCoulomb(FrictionCoulomb const& rOther) {*this = rOther;}

  }; // Class FrictionCoulomb

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    FrictionCoulomb& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const FrictionCoulomb& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos

#endif // FRICTION_MODEL_COULOMB_H_INCLUDED
