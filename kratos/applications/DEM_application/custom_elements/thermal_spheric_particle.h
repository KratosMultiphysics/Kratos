//   
//   Project Name:        ThermalD       
//   Last Modified by:    $Author: Ferran Arrufat $
//   Date:                $Date: 2015-02-01  $
//   Revision:            $Revision: 1.0.0.0 $
//
//


#if !defined(KRATOS_THERMAL_SPHERIC_PARTICLE_H_INCLUDED )
#define  KRATOS_THERMAL_SPHERIC_PARTICLE_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "spheric_particle.h"
#include "spheric_continuum_particle.h"
#include "Particle_Contact_Element.h"
#include "containers/vector_component_adaptor.h"
#include "discrete_element.h"


namespace Kratos
{
  template< class TBaseElement >  
  class KRATOS_API(DEM_APPLICATION) ThermalSphericParticle : public TBaseElement
    {
    public:
      
      /// Pointer definition of ThermalSphericParticle
      KRATOS_CLASS_POINTER_DEFINITION(ThermalSphericParticle);
      
      typedef WeakPointerVector<Element> ParticleWeakVectorType;
      typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
      typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;
      
      typedef Node <3> NodeType;
      typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
      typedef std::size_t IndexType;
      typedef Geometry<Node < 3 > > GeometryType;
      typedef Properties PropertiesType;
      
      using TBaseElement::GetGeometry;
      using TBaseElement::GetProperties;
      using TBaseElement::mNeighbourElements;
      using TBaseElement::GetRadius;
      using TBaseElement::GetMass;

      /// Default constructor. 
      ThermalSphericParticle():TBaseElement(){};
      ThermalSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry):TBaseElement(NewId, pGeometry){};
      ThermalSphericParticle(IndexType NewId, NodesArrayType const& ThisNodes):TBaseElement(NewId, ThisNodes){};
      ThermalSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties):TBaseElement(NewId, pGeometry, pProperties){};

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {          
        return Element::Pointer(new ThermalSphericParticle<TBaseElement>(NewId, GetGeometry().Create(ThisNodes), pProperties));
      };
                   
      /// Destructor.
      virtual ~ThermalSphericParticle(){};
       
    
      double GetTemperature();     
      void ComputeConductiveHeatFlux(const ProcessInfo& r_process_info);   
      void ComputeConvectiveHeatFlux(const ProcessInfo& r_process_info);  
      void CalculateRightHandSide(ProcessInfo& r_current_process_info,
                                  double dt, 
                                  const array_1d<double,3>& gravity,
                                  int search_control);  
      void FinalizeSolutionStep(ProcessInfo& r_process_info); 
      void UpdateTemperature(const ProcessInfo& r_process_info); 
          
    
      /// Turn back information as a string.
      virtual std::string Info() const
      {
        std::stringstream buffer;
        buffer << "ThermalSphericParticle" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "ThermalSphericParticle";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const {}
      
    protected:
       
       virtual void CustomInitialize();	
              
       //thermal sphere neighbor information
        
        double mTemperature;
        double mConductiveHeatFlux;
        double mThermalConductivity;
        double mSpecificHeat;       
    

    private:

      friend class Serializer;

      virtual void save(Serializer& rSerializer) const
      {
          KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle );

      }

      virtual void load(Serializer& rSerializer)
      {
          KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle );

      }
            
    }; 
 
  /// input stream function
   /* template<TBaseElement>
    inline std::istream& operator >> (std::istream& rIStream, ThermalSphericParticle<TBaseElement>& rThis){ return rIStream;}
    
  /// output stream function
    template<TBaseElement>
    inline std::ostream& operator << (std::ostream& rOStream, const ThermalSphericParticle<TBaseElement>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
  
}  // namespace Kratos.

#endif // KRATOS_THERMAL_SPHERIC_PARTICLE_H_INCLUDED  defined 


