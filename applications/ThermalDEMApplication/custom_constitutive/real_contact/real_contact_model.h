//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(REAL_CONTACT_MODEL_H_INCLUDED)
#define REAL_CONTACT_MODEL_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes
#include "includes/define.h"
#include "includes/model_part.h"

// Project includes

namespace Kratos
{
  // Forward declarations
  class ThermalSphericParticle;

  class KRATOS_API(THERMAL_DEM_APPLICATION) RealContactModel
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(RealContactModel);

      // Constructor / Destructor
      RealContactModel();
      virtual ~RealContactModel();

      // Public methods
      virtual void SetRealContactModelInProperties (Properties::Pointer pProp, bool verbose = true) const;
      virtual void AdjustContact                   (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);

      // Clone
      virtual RealContactModel* CloneRaw() const {
        RealContactModel* cloned_model(new RealContactModel(*this));
        return cloned_model;
      }

      virtual RealContactModel::Pointer CloneShared() const {
        RealContactModel::Pointer cloned_model(new RealContactModel(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const {
        std::stringstream buffer;
        buffer << "RealContactModel";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const { rOStream << "RealContactModel"; }
      virtual void PrintData(std::ostream& rOStream) const {}

    private:

      // Assignment operator / Copy constructor
      RealContactModel& operator=(RealContactModel const& rOther) {return *this;}
      RealContactModel(RealContactModel const& rOther) {*this = rOther;}

      // Serializer
      friend class Serializer;

      virtual void save(Serializer& rSerializer) const {
        //rSerializer.save("MyMemberName",myMember);
      }

      virtual void load(Serializer& rSerializer) {
        //rSerializer.load("MyMemberName",myMember);
      }

  }; // Class RealContactModel

  // input stream function
  inline std::istream& operator >> (std::istream& rIStream, RealContactModel& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator << (std::ostream& rOStream, const RealContactModel& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos

#endif // REAL_CONTACT_MODEL_H_INCLUDED
