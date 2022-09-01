//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(DIRECT_CONDUCTION_MODEL_BOB_COMPLETE_H_INCLUDED)
#define DIRECT_CONDUCTION_MODEL_BOB_COMPLETE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "direct_conduction_model.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) DirectConductionBOBComplete : public DirectConductionModel
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(DirectConductionBOBComplete);

      // Definitions
      #define ETA_MIN  1.0
      #define ETA_MAX  100.0
      #define LAMB_MIN 0.01
      #define LAMB_MAX 100.0

      // Constructor / Destructor
      DirectConductionBOBComplete();
      virtual ~DirectConductionBOBComplete();

      // Public methods
      double GetSearchDistance        (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;
      double ComputeHeatFlux          (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;
      double ComputeHeatTransferCoeff (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);
      double ContactCoeff             (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);
      double SeparatedCoeff           (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);
      double ContactCoeffEtaMin       (const double Reff, const double kf, const double eta, const double kpf);
      double ContactCoeffEtaMax       (const double Reff, const double kf, const double eta, const double kpf);
      double SeparatedCoeffLambMin    (const double Reff, const double kf, const double kpf2);
      double SeparatedCoeffLambMax    (const double Reff, const double kf, const double Rcyl, const double Ds);

      // Clone
      HeatExchangeMechanism* CloneRaw() const override {
        HeatExchangeMechanism* cloned_model(new DirectConductionBOBComplete(*this));
        return cloned_model;
      }

      HeatExchangeMechanism::Pointer CloneShared() const override {
        HeatExchangeMechanism::Pointer cloned_model(new DirectConductionBOBComplete(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "DirectConductionBOBComplete";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "DirectConductionBOBComplete"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      DirectConductionBOBComplete& operator=(DirectConductionBOBComplete const& rOther) {return *this;}
      DirectConductionBOBComplete(DirectConductionBOBComplete const& rOther) {*this = rOther;}

  }; // Class DirectConductionBOBComplete

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    DirectConductionBOBComplete& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const DirectConductionBOBComplete& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos

#endif // DIRECT_CONDUCTION_MODEL_BOB_COMPLETE_H_INCLUDED
