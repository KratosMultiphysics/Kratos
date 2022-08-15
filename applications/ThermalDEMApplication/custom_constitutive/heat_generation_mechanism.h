//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(HEAT_GENERATION_MECHANISM_H_INCLUDED)
#define HEAT_GENERATION_MECHANISM_H_INCLUDED

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

  class KRATOS_API(THERMAL_DEM_APPLICATION) HeatGenerationMechanism
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(HeatGenerationMechanism);

      // Constructor / Destructor
      HeatGenerationMechanism();
      virtual ~HeatGenerationMechanism();

      // Public methods
      virtual void   SetHeatGenerationMechanismInProperties (Properties::Pointer pProp, bool verbose = true) const;
      virtual double ComputeHeatGeneration                  (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);
      virtual double FinalizeHeatGeneration                 (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);

      // Clone
      virtual HeatGenerationMechanism* CloneRaw() const {
        HeatGenerationMechanism* cloned_mechanism(new HeatGenerationMechanism(*this));
        return cloned_mechanism;
      }

      virtual HeatGenerationMechanism::Pointer CloneShared() const {
        HeatGenerationMechanism::Pointer cloned_mechanism(new HeatGenerationMechanism(*this));
        return cloned_mechanism;
      }

      // Turn back information as a string
      virtual std::string Info() const {
        std::stringstream buffer;
        buffer << "HeatGenerationMechanism";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const { rOStream << "HeatGenerationMechanism"; }
      virtual void PrintData(std::ostream& rOStream) const {}

    private:

      // Assignment operator / Copy constructor
      HeatGenerationMechanism& operator=(HeatGenerationMechanism const& rOther) {return *this;}
      HeatGenerationMechanism(HeatGenerationMechanism const& rOther) {*this = rOther;}

      // Serializer
      friend class Serializer;

      virtual void save(Serializer& rSerializer) const {
        //rSerializer.save("MyMemberName",myMember);
      }

      virtual void load(Serializer& rSerializer) {
        //rSerializer.load("MyMemberName",myMember);
      }

  }; // Class HeatGenerationMechanism

  // input stream function
  inline std::istream& operator >> (std::istream& rIStream, HeatGenerationMechanism& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator << (std::ostream& rOStream, const HeatGenerationMechanism& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos

#endif // HEAT_GENERATION_MECHANISM_H_INCLUDED
