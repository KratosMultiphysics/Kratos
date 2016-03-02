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
    
  class ThermalSphericParticle : public SphericContinuumParticle
    {
    public:
      
      /// Pointer definition of ThermalSphericParticle
      KRATOS_CLASS_POINTER_DEFINITION(ThermalSphericParticle);

      typedef WeakPointerVector<Element> ParticleWeakVectorType;
      typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
      typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;

      /// Default constructor. 
      ThermalSphericParticle( IndexType NewId, GeometryType::Pointer pGeometry );
      ThermalSphericParticle( IndexType NewId, NodesArrayType const& ThisNodes);
      ThermalSphericParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );
      
      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;
         
      /// Destructor.
      virtual ~ThermalSphericParticle();
       
    
      double GetTemperature();     
      void ComputeConductiveHeatFlux(const ProcessInfo& r_process_info);   
      void ComputeConvectiveHeatFlux(const ProcessInfo& r_process_info);  
      void CalculateRightHandSide(VectorType& r_right_hand_side_vector, 
                                  ProcessInfo& r_current_process_info,
                                  double dt, 
                                  const array_1d<double,3>& gravity,
                                  int search_control);  
      void FinalizeSolutionStep(const ProcessInfo& r_process_info); 
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

      
//      
//      //member variables DEM_CONTINUUM
//      int mContinuumGroup;
//      std::vector<SphericContinuumParticle*> mContinuumIniNeighbourElements;
//      std::vector<ParticleContactElement*> mBondElements;
//      std::vector<int> mIniNeighbourIds;
//      
//      
//      //member variables DEM_THERMAL
//      int mThermalGroup;
//      std::vector<ThermalSphericParticle*> mThermalIniNeighbourElements;
//      std::vector<ParticleContactElement*> mBondElements;
//      std::vector<int> mIniNeighbourIds;
            
      
      
      
      ///@}
      
    protected:

       ThermalSphericParticle();
       
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
      
      
      
      /*
      /// Assignment operator.
      SphericContinuumParticle& operator=(SphericContinuumParticle const& rOther)
      {
    return *this;
      }

      /// Copy constructor.
      SphericContinuumParticle(SphericContinuumParticle const& rOther)
      {
    *this = rOther;
      }
      */
        
      ///@}        
    }; 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
                    ThermalSphericParticle& rThis){ return rIStream;}
  
  
  

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
                    const ThermalSphericParticle& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_THERMAL_SPHERIC_PARTICLE_H_INCLUDED  defined 


