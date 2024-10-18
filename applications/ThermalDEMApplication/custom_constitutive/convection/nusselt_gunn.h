//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(NUSSELT_NUMBER_GUNN_H_INCLUDED)
#define NUSSELT_NUMBER_GUNN_H_INCLUDED

// System includes

// External includes

// Project includes
#include "convection_model.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) NusseltGunn : public ConvectionModel
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(NusseltGunn);

      // Constructor / Destructor
      NusseltGunn();
      virtual ~NusseltGunn();

      // Public methods
      double ComputeNusseltNumber(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;

      // Clone
      HeatExchangeMechanism* CloneRaw() const override {
        HeatExchangeMechanism* cloned_model(new NusseltGunn(*this));
        return cloned_model;
      }

      HeatExchangeMechanism::Pointer CloneShared() const override {
        HeatExchangeMechanism::Pointer cloned_model(new NusseltGunn(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "NusseltGunn";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "NusseltGunn"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      NusseltGunn& operator=(NusseltGunn const& rOther) {return *this;}
      NusseltGunn(NusseltGunn const& rOther) {*this = rOther;}

  }; // Class NusseltGunn

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    NusseltGunn& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const NusseltGunn& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos

#endif // NUSSELT_NUMBER_GUNN_H_INCLUDED
