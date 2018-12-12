//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_CYLINDER_CONTINUUM_PARTICLE_H_INCLUDED )
#define  KRATOS_CYLINDER_CONTINUUM_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>
// Project includes
#include "includes/define.h"
#include "spheric_continuum_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"

namespace Kratos
{

  class KRATOS_API(DEM_APPLICATION) CylinderContinuumParticle: public SphericContinuumParticle
    {
    public:

      KRATOS_CLASS_POINTER_DEFINITION(CylinderContinuumParticle);

      typedef WeakPointerVector<Element> ParticleWeakVectorType;  //M: l'he afegit jo.. esta be aquesta?
      typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
      typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;

      CylinderContinuumParticle( IndexType NewId, GeometryType::Pointer pGeometry );
      CylinderContinuumParticle( IndexType NewId, NodesArrayType const& ThisNodes);
      CylinderContinuumParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

      /// Destructor.
      virtual ~CylinderContinuumParticle();


      virtual std::string Info() const override
      {
          std::stringstream buffer;
          buffer << "CylinderContinuumParticle" ;
          return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "CylinderContinuumParticle";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const override {}

      void ContactAreaWeighting() override;

    protected:

      CylinderContinuumParticle();

      double CalculateVolume() override;
      double CalculateMomentOfInertia() override;
      void AddContributionToRepresentativeVolume(const double distance, const double radius_sum, const double contact_area) override ;


    private:

      friend class Serializer;

      virtual void save(Serializer& rSerializer) const override
      {
          KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DiscreteElement );
      }

      virtual void load(Serializer& rSerializer) override
      {
          KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DiscreteElement );
      }

    };


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
                    CylinderContinuumParticle& rThis){ return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
                    const CylinderContinuumParticle& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

}  // namespace Kratos.

#endif // KRATOS_SPHERIC_PARTICLE_H_INCLUDED  defined


