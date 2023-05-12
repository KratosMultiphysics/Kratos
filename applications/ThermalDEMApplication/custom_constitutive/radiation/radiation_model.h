//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(RADIATION_MODEL_H_INCLUDED)
#define RADIATION_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/heat_exchange_mechanism.h"
#include "custom_elements/thermal_spheric_particle.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) RadiationModel : public HeatExchangeMechanism
  {
    public:

      // Definitions
      #define STEFAN_BOLTZMANN 5.670374419e-8

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(RadiationModel);

      // Constructor / Destructor
      RadiationModel();
      virtual ~RadiationModel();

      // Public methods
      void         SetHeatExchangeMechanismInProperties (Properties::Pointer pProp, bool verbose = true) const override;
      double       GetSearchDistance                    (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;
      double       ComputeHeatFlux                      (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;
      double       FinalizeHeatFlux                     (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;
      virtual void AccumulateEnvironmentTemperature     (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);

      // Clone
      HeatExchangeMechanism* CloneRaw() const override {
        HeatExchangeMechanism* cloned_model(new RadiationModel(*this));
        return cloned_model;
      }

      HeatExchangeMechanism::Pointer CloneShared() const override {
        HeatExchangeMechanism::Pointer cloned_model(new RadiationModel(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "RadiationModel";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "RadiationModel"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      RadiationModel& operator=(RadiationModel const& rOther) {return *this;}
      RadiationModel(RadiationModel const& rOther) {*this = rOther;}

  }; // Class RadiationModel

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    RadiationModel& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const RadiationModel& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos

#endif // RADIATION_MODEL_H_INCLUDED
