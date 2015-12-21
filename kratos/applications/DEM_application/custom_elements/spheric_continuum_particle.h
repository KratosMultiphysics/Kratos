//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//


#if !defined(KRATOS_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED )
#define  KRATOS_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "spheric_particle.h"
#include "Particle_Contact_Element.h"
#include "../custom_constitutive/DEM_continuum_constitutive_law.h"
#include "containers/vector_component_adaptor.h"



namespace Kratos
{
  
  class SphericContinuumParticle : public SphericParticle
    {
    public:
      
      /// Pointer definition of SphericContinuumParticle
      KRATOS_CLASS_POINTER_DEFINITION(SphericContinuumParticle);

      typedef WeakPointerVector<Element> ParticleWeakVectorType;
      typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
      typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;

      /// Default constructor. 
      SphericContinuumParticle( IndexType NewId, GeometryType::Pointer pGeometry );
      SphericContinuumParticle( IndexType NewId, NodesArrayType const& ThisNodes);
      SphericContinuumParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );
      
      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;
         
      /// Destructor.
      virtual ~SphericContinuumParticle();
       
      void FullInitialize(const ProcessInfo& rCurrentProcessInfo);
      void SetInitialSphereContacts(ProcessInfo& rCurrentProcessInfo);
      void SetInitialFemContacts();
      void CreateContinuumConstitutiveLaws(ProcessInfo& rCurrentProcessInfo);
      void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);    
      virtual void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);     
      virtual void ContinuumSphereMemberDeclarationFirstStep(const ProcessInfo& rCurrentProcessInfo); 
      void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo);      
      void ComputeNewNeighboursHistoricalData(std::vector<unsigned int>&                  mTempNeighboursIds, 
                                              std::vector<array_1d<double, 3> >& mTempNeighbourElasticContactForces,
                                              std::vector<array_1d<double, 3> >& mTempNeighbourTotalContactForces,
                                              std::vector<SphericParticle*>&     mTempNeighbourElements,
                                              std::vector<double>&               mTempNeighboursDelta,
                                              std::vector<int>&                  mTempNeighboursFailureId,
                                              std::vector<int>&                  mTempNeighboursMapping,
                                              std::vector<int>&                  mTempContNeighboursMapping) ;
      virtual void ComputeNewRigidFaceNeighboursHistoricalData();      
      
      virtual void CalculateMeanContactArea(const bool has_mpi, const ProcessInfo& rCurrentProcessInfo, const bool first);      
      virtual void CalculateOnContactElements(unsigned int neighbour_iterator_id, size_t i_neighbour_count, int mapping, double LocalElasticContactForce[3], 
                                              double contact_sigma, double contact_tau, double failure_criterion_state, double acumulated_damage, int time_steps);
           
      virtual void AddPoissonContribution( const double equiv_poisson, double LocalCoordSystem[3][3], double& normal_force, double calculation_area);      
      virtual void ContactAreaWeighting();

      /// Turn back information as a string.
      virtual std::string Info() const
      {
        std::stringstream buffer;
        buffer << "SphericContinuumParticle" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "SphericContinuumParticle";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const {}

      
      
      //member variables DEM_CONTINUUM
      int mContinuumGroup;
      std::vector<SphericContinuumParticle*> mContinuumIniNeighbourElements;
      std::vector<Particle_Contact_Element*> mBondElements;
      std::vector<int> mIniNeighbourIds;
      std::vector<int> mNeighbourFailureId;
      std::vector<int> mIniNeighbourFailureId;
      std::vector<double> mNeighbourDelta;
      std::vector<int> mMappingNewCont;
      std::vector<int> mMappingNewIni;

      ///@}
      
    protected:

        SphericContinuumParticle();
        
        double AreaDebugging(const ProcessInfo& rCurrentProcessInfo); //MSIMSI DEBUG                
        void SymmetrizeTensor(const ProcessInfo& rCurrentProcessInfo );
        virtual void CustomInitialize();	
        virtual double GetInitialDeltaWithFEM(int index);      
        void ComputeAdditionalForces(array_1d<double, 3>& additionally_applied_force, array_1d<double, 3>& additionally_applied_moment, ProcessInfo& rCurrentProcessInfo, const array_1d<double,3>& gravity);
        virtual void ComputeBallToBallContactForce(array_1d<double, 3>& rElasticForce,
                                           array_1d<double, 3>& rContactForce, 
                                           array_1d<double, 3>& InitialRotaMoment, 
                                           ProcessInfo& rCurrentProcessInfo, 
                                           double dt,
                                           const bool multi_stage_RHS);         

        void ComputePressureForces(array_1d<double, 3>& externally_applied_force, ProcessInfo& rCurrentProcessInfo);
        void ApplyLocalMomentsDamping(const ProcessInfo& rCurrentProcessInfo );
        void CharacteristicParticleFailureId(const ProcessInfo& rCurrentProcessInfo );                
        void ComputeParticleBlockContactForce(const ProcessInfo& rCurrentProcessInfo);
        void ComputeParticleRotationSpring();
        void ComputeParticleSurfaceContactForce(ProcessInfo& rCurrentProcessInfo);
        void ComputeParticleRotationSpring_TRIAL(const ProcessInfo& rCurrentProcessInfo); //provisional                
                
        //double mFinalPressureTime;
        //double mFinalSimulationTime;
     
        int*  mSkinSphere; 

        //sphere neighbour information
        
        std::vector<double>         mContIniNeighArea;        
        std::vector<int>            mIniNeighbourToIniContinuum;        
        std::vector<double>         mIniNeighbourDelta;
                  
        //fem neighbour information
        std::vector<double>         mFemNeighbourDelta;                        
        std::vector<int>            mFemIniNeighbourIds;
        std::vector<double>         mFemIniNeighbourDelta;
        std::vector<int>            mFemMappingNewIni;
                        
        std::vector<Kratos::DEMContinuumConstitutiveLaw::Pointer> mContinuumConstitutiveLawArray;
        
                                    
        //FOR DEM_FEM APP        
        void ComputeParticleBlockContactForce_With_Rotation();
        void ComputeParticleBlockContactForce_Without_Rotation();
        void FindContactFaceOfBlockForParticle(ParticleWeakIteratorType rObj_2, int & RightFace, double LocalCoordSystem[3][3], double Coeff[4],double &DistPToB);       

        double distances_squared;
      
        
    private:


      friend class Serializer;

      virtual void save(Serializer& rSerializer) const
      {
          KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle );
          //rSerializer.save("mContinuumGroup",mContinuumGroup);
          //rSerializer.save("mIniNeighbourIds",mIniNeighbourIds);
          //rSerializer.save("mSymmStressTensor",mSymmStressTensor);
      }

      virtual void load(Serializer& rSerializer)
      {
          KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle );
          //rSerializer.load("mContinuumGroup",mContinuumGroup);
          //rSerializer.load("mIniNeighbourIds",mIniNeighbourIds);
          //rSerializer.load("mSymmStressTensor",mSymmStressTensor);
          mContinuumGroup        = this->GetGeometry()[0].FastGetSolutionStepValue(COHESIVE_GROUP);  
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
        
    }; // Class SphericContinuumParticle 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
                    SphericContinuumParticle& rThis){ return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
                    const SphericContinuumParticle& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED  defined 


