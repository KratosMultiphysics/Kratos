//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(REAL_CONTACT_MODEL_MORRIS_AREA_H_INCLUDED)
#define REAL_CONTACT_MODEL_MORRIS_AREA_H_INCLUDED

// System includes

// External includes

// Project includes
#include "real_contact_model.h"
#include "custom_elements/thermal_spheric_particle.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) RealContactMorrisArea : public RealContactModel
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(RealContactMorrisArea);

      // Constructor / Destructor
      RealContactMorrisArea();
      virtual ~RealContactMorrisArea();

      // Public methods
      void AdjustContact(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;

      // Clone
      RealContactModel* CloneRaw() const override {
        RealContactModel* cloned_model(new RealContactMorrisArea(*this));
        return cloned_model;
      }

      RealContactModel::Pointer CloneShared() const override {
        RealContactModel::Pointer cloned_model(new RealContactMorrisArea(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "RealContactMorrisArea";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "RealContactMorrisArea"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      RealContactMorrisArea& operator=(RealContactMorrisArea const& rOther) {return *this;}
      RealContactMorrisArea(RealContactMorrisArea const& rOther) {*this = rOther;}

  }; // Class RealContactMorrisArea

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    RealContactMorrisArea& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const RealContactMorrisArea& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos

#endif // REAL_CONTACT_MODEL_MORRIS_AREA_H_INCLUDED
