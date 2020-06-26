//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_CYLINDER_PARTICLE_H_INCLUDED )
#define  KRATOS_CYLINDER_PARTICLE_H_INCLUDED



// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "discrete_element.h"
#include "spheric_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"


namespace Kratos
{

  class KRATOS_API(DEM_APPLICATION) CylinderParticle : public SphericParticle
    {
    public:

      KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(CylinderParticle);

      typedef GlobalPointersVector<Element> ParticleWeakVectorType;  //M: l'he afegit jo.. esta be aquesta?
      typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
      typedef GlobalPointersVector<Element >::iterator ParticleWeakIteratorType;


      CylinderParticle( IndexType NewId, GeometryType::Pointer pGeometry );
      CylinderParticle( IndexType NewId, NodesArrayType const& ThisNodes);
      CylinderParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

      /// Destructor.
      virtual ~CylinderParticle();

      double CalculateVolume() override;
      double CalculateMomentOfInertia() override;

      void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) override;
      void Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& r_process_info) override;
      void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info) override;
      void Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& r_process_info) override;



      /// Turn back information as a string.
      virtual std::string Info() const override
      {
          std::stringstream buffer;
          buffer << "CylinderParticle" ;
          return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "CylinderParticle";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const override {}


    protected:

      CylinderParticle();

    private:


      friend class Serializer;

      virtual void save(Serializer& rSerializer) const override
      {
          KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle );
      }

      virtual void load(Serializer& rSerializer) override
      {
          KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle );
      }



    }; // Class SphericParticle


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
                    CylinderParticle& rThis){ return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
                    const CylinderParticle& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SPHERIC_PARTICLE_H_INCLUDED  defined


