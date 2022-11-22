//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(RADIATION_MODEL_CONTINUUM_KRAUSE_H_INCLUDED)
#define RADIATION_MODEL_CONTINUUM_KRAUSE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "radiation_model.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) RadiationContinuumKrause : public RadiationModel
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(RadiationContinuumKrause);

      // Constructor / Destructor
      RadiationContinuumKrause();
      virtual ~RadiationContinuumKrause();

      // Public methods
      double ComputeHeatFlux                  (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;
      double FinalizeHeatFlux                 (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;
      void   AccumulateEnvironmentTemperature (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;

      // Clone
      HeatExchangeMechanism* CloneRaw() const override {
        HeatExchangeMechanism* cloned_model(new RadiationContinuumKrause(*this));
        return cloned_model;
      }

      HeatExchangeMechanism::Pointer CloneShared() const override {
        HeatExchangeMechanism::Pointer cloned_model(new RadiationContinuumKrause(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "RadiationContinuumKrause";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "RadiationContinuumKrause"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      RadiationContinuumKrause& operator=(RadiationContinuumKrause const& rOther) {return *this;}
      RadiationContinuumKrause(RadiationContinuumKrause const& rOther) {*this = rOther;}

  }; // Class RadiationContinuumKrause

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    RadiationContinuumKrause& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const RadiationContinuumKrause& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos

#endif // RADIATION_MODEL_CONTINUUM_KRAUSE_H_INCLUDED
