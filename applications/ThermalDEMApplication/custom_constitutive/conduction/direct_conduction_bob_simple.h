//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(DIRECT_CONDUCTION_MODEL_BOB_SIMPLE_H_INCLUDED)
#define DIRECT_CONDUCTION_MODEL_BOB_SIMPLE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "direct_conduction_model.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) DirectConductionBOBSimple : public DirectConductionModel
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(DirectConductionBOBSimple);

      // Constructor / Destructor
      DirectConductionBOBSimple();
      virtual ~DirectConductionBOBSimple();

      // Public methods
      double ComputeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;

      // Clone
      HeatExchangeMechanism* CloneRaw() const override {
        HeatExchangeMechanism* cloned_model(new DirectConductionBOBSimple(*this));
        return cloned_model;
      }

      HeatExchangeMechanism::Pointer CloneShared() const override {
        HeatExchangeMechanism::Pointer cloned_model(new DirectConductionBOBSimple(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "DirectConductionBOBSimple";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "DirectConductionBOBSimple"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      DirectConductionBOBSimple& operator=(DirectConductionBOBSimple const& rOther) {return *this;}
      DirectConductionBOBSimple(DirectConductionBOBSimple const& rOther) {*this = rOther;}

  }; // Class DirectConductionBOBSimple

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    DirectConductionBOBSimple& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const DirectConductionBOBSimple& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos

#endif // DIRECT_CONDUCTION_MODEL_BOB_SIMPLE_H_INCLUDED
