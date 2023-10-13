//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(NUSSELT_NUMBER_HANZ_MARSHALL_H_INCLUDED)
#define NUSSELT_NUMBER_HANZ_MARSHALL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "convection_model.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) NusseltHanzMarshall : public ConvectionModel
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(NusseltHanzMarshall);

      // Constructor / Destructor
      NusseltHanzMarshall();
      virtual ~NusseltHanzMarshall();

      // Public methods
      double ComputeNusseltNumber(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;

      // Clone
      HeatExchangeMechanism* CloneRaw() const override {
        HeatExchangeMechanism* cloned_model(new NusseltHanzMarshall(*this));
        return cloned_model;
      }

      HeatExchangeMechanism::Pointer CloneShared() const override {
        HeatExchangeMechanism::Pointer cloned_model(new NusseltHanzMarshall(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "NusseltHanzMarshall";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "NusseltHanzMarshall"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      NusseltHanzMarshall& operator=(NusseltHanzMarshall const& rOther) {return *this;}
      NusseltHanzMarshall(NusseltHanzMarshall const& rOther) {*this = rOther;}

  }; // Class NusseltHanzMarshall

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    NusseltHanzMarshall& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const NusseltHanzMarshall& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos

#endif // NUSSELT_NUMBER_HANZ_MARSHALL_H_INCLUDED
