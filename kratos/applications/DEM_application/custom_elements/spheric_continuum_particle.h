//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
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
#include "containers/vector_component_adaptor.h"


namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
  ///@{

  ///@name Kratos Globals
  ///@{ 
  
  ///@} 
  ///@name Type Definitions
  ///@{ 
  
  ///@} 
  ///@name  Enum's
  ///@{
      
  ///@}
  ///@name  Functions 
  ///@{
      
  ///@}
  ///@name Kratos Classes
  ///@{
  
  /// Short class definition.
  /** Detail class definition.
  */
  class SphericContinuumParticle : public SphericParticle
    {
    public:


      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of SphericContinuumParticle
      KRATOS_CLASS_POINTER_DEFINITION(SphericContinuumParticle);

      typedef WeakPointerVector<Element> ParticleWeakVectorType;
      typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
      typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;
      
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor. 
      SphericContinuumParticle( IndexType NewId, GeometryType::Pointer pGeometry );
      SphericContinuumParticle( IndexType NewId, NodesArrayType const& ThisNodes);
      SphericContinuumParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );
      
      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;
         
      /// Destructor.
      virtual ~SphericContinuumParticle();

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
       


      void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);
    
      void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);
      
      void ContinuumSphereMemberDeclarationFirstStep(const ProcessInfo& rCurrentProcessInfo);
 
      void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo);
      
      void ComputeNewNeighboursHistoricalData(std::vector<int>& mTempNeighboursIds, std::vector<array_1d<double, 3> >& mTempNeighbourElasticContactForces,
                                                       std::vector<array_1d<double, 3> >& mTempNeighbourTotalContactForces);
      void ComputeNewRigidFaceNeighboursHistoricalData();
      
      virtual void NonlinearNormalForceCalculation(double LocalElasticContactForce[3], double kn1, double kn2, double distance, double max_dist, double initial_dist) ;


      virtual void EvaluateFailureCriteria(double LocalElasticContactForce[3], double ShearForceNow, double corrected_area, int i_neighbour_count, double& contact_sigma, double& contact_tau, double& failure_criterion_state, bool& sliding, int mapping);
      
      virtual void CalculateOnContactElements(unsigned int neighbour_iterator_id, size_t i_neighbour_count, int mapping, double LocalElasticContactForce[3], 
                                              double contact_sigma, double contact_tau, double failure_criterion_state, double acumulated_damage, int time_steps);

      virtual void ComputeStressStrain(   double mStressTensor[3][3],
                                          ProcessInfo& rCurrentProcessInfo,
                                          double& rRepresentative_Volume);
      
            
      virtual void StressTensorOperations(double mStressTensor[3][3],
                                          double GlobalElasticContactForce[3],
                                          array_1d<double,3> &other_to_me_vect,
                                          const double &distance,
                                          const double &radius_sum,
                                          const double &corrected_area,
                                          //ParticleWeakIteratorType neighbour_iterator,
                                          SphericParticle* neighbour_iterator,
                                          ProcessInfo& rCurrentProcessInfo,
                                          double &rRepresentative_Volume);
      
      virtual void AddPoissonContribution(double LocalCoordSystem[3][3],
                                          double GlobalContactForce[3],
                                          double GlobalElasticContactForce[3],
                                          //double ViscoDampingGlobalContactForce[3],
                                          array_1d<double, 3>& rContactForce,
                                          array_1d<double,3>& damp_forces);

      
      //virtual void InitializeContactElements(ParticleWeakIteratorType neighbour_iterator, double& corrected_area);

      ///***********************************************////////////// AIXO ES DECLARA AKI O LA INITIALITZACIÓ.

      //l'he de definir aqui aquest???
 
       //std::size_t& GetNumberOfNeighbours(){return(GetGeometry()(0)->FastGetSolutionStepValue(NUMBER_OF_NEIGHBOURS));};

       
       //double mInitialDelta;
       //vector<int> mVectorContactFailureId;
       //int mContactFailureId;

       //vector< double > mVectorContactInitialDelta; R: cal cridar-ho cada cop per no fer copia!!
       //double mContactInitialDelta;

       //vector<array_1d<double,3> > mVectorContactForces;
       //array_1d<double,3>& mContactForces;



       //auxiliar variables
/*
       double mOld_Displacement_X;
       double mOld_Displacement_Y;
       double mOld_Displacement_Z;
       double mDisplacement_X;
       double mDisplacement_Y;
       double mDisplacement_Z;
*/
      // std::vector<double>   mInitialDelta = GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_CONTACT_INITIAL_DELTA);
      // std::vector<int>  mContactFailureId = GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_CONTACT_FAILURE_ID);
       ///@}
      ///@name Access
      ///@{ 
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

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
      
            
      ///@}      
      ///@name Friends
      ///@{
      std::vector<SphericContinuumParticle*> mContinuumIniNeighbourElements;
      //std::vector<int> mContinuumIniNeighbourIds;
      std::vector<Particle_Contact_Element*> mBondElements;
      
      //member variables DEM_CONTINUUM
      int mContinuumGroup;
            
      ///@}
      
    protected:

       SphericContinuumParticle();

        void SetInitialSphereContacts();
        void CreateContinuumConstitutiveLaws();
        void SetInitialFemContacts();
        void NeighNeighMapping( ProcessInfo& rCurrentProcessInfo ); //MSIMSI DEBUG
        void CheckPairWiseBreaking(); //MSIMSI DEBUG
        double AreaDebugging(const ProcessInfo& rCurrentProcessInfo); //MSIMSI DEBUG
        
        virtual void ContactAreaWeighting2D();
        void ContactAreaWeighting3D();
        void SymmetrizeTensor(const ProcessInfo& rCurrentProcessInfo );
        
        virtual void CustomInitialize();
	
        virtual double GetInitialDelta(int index);
      
        void ComputeAdditionalForces(array_1d<double, 3>& additionally_applied_force, array_1d<double, 3>& additionally_applied_moment, ProcessInfo& rCurrentProcessInfo);
        void ComputeBallToBallContactForce(   //array_1d<double, 3>& rContactForce, 
                                                    //array_1d<double, 3>& rContactMoment, 
                                                    array_1d<double, 3>& rElasticForce, 
                                                    array_1d<double, 3>& InitialRotaMoment, 
                                                    ProcessInfo& rCurrentProcessInfo, 
                                                    const bool multi_stage_RHS);         
        //virtual void ComputeBallToSurfaceContactForce(array_1d<double, 3>& rContactForce, array_1d<double, 3>& rContactMoment, array_1d<double, 3>& InitialRotaMoment, array_1d<double, 3>& MaxRotaMoment, ProcessInfo& rCurrentProcessInfo);
        //MSIMSI 6 aixo hauria de cridar el del basic o cal ke sigui del continu?
        
        void ComputePressureForces(array_1d<double, 3>& externally_applied_force, ProcessInfo& rCurrentProcessInfo);
        void PlasticityAndDamage1D(double LocalElasticContactForce[3], double kn, double equiv_young, double indentation, double corrected_area, double radius_sum_i, double& failure_criterion_state, double& acumulated_damage, int i_neighbour_count, int mapping_new_cont, int mapping_new_ini, int time_steps);
        
        //void ApplyLocalForcesDamping(const ProcessInfo& rCurrentProcessInfo );
        void ApplyLocalMomentsDamping(const ProcessInfo& rCurrentProcessInfo );
        void CharacteristicParticleFailureId(const ProcessInfo& rCurrentProcessInfo );
        //void CalculateInitialLocalAxes(const ProcessInfo& rCurrentProcessInfo );
        //void CalculateLocalAxes(const ProcessInfo& rCurrentProcessInfo );
        
        
        
        void ComputeParticleBlockContactForce(const ProcessInfo& rCurrentProcessInfo);
        void ComputeParticleRotationSpring();
        void ComputeParticleSurfaceContactForce(ProcessInfo& rCurrentProcessInfo);
        void ComputeParticleRotationSpring_TRIAL(const ProcessInfo& rCurrentProcessInfo); //provisional                
                
        //DEMPACK
        int mDempack;
        double mDempack_damping;
        double mDempack_global_damping;
        vector< array_1d<double, 4> > mHistory;
        double mNcstr1_el;
        double mNcstr2_el;
        double mYoungPlastic;
        double mPlasticityLimit;
        double mDamageMaxDisplacementFactor;
     
        double mGamma1;
        double mGamma2;
        double mGamma3;
        double mMaxDef;
        
        
        double mStressTensor[3][3]; 
        double mSymmStressTensor[3][3]; 
        bool mContinuumSimulationOption;
        bool mContactMeshOption;
        int mTriaxialOption;
        int mStressStrainOption;

        double mFinalPressureTime;
        double mFinalSimulationTime;
     
        int mpCaseOption;
        int  mFailureId;
        int*  mSkinSphere;
        
        const double *mpCurrentTime;
   
        int mFailureCriterionOption;
        
        double mTension;
        double mCohesion;
        double mSectionalInertia;
        
        
        double mTensionLimit;
        double mTauZero;
        double mContactInternalFriccion;
        double mTanContactInternalFriccion;
        double mSinContactInternalFriccion;
        double mCosContactInternalFriccion;

        //sphere neighbour information
        
        Vector mcont_ini_neigh_area;
        std::vector<int> mIniNeighbourIds;
        Vector mIniNeighbourDelta;

        vector<int> mIniNeighbourFailureId;
        vector<int> mIniNeighbourToIniContinuum;
        
        std::vector<int> mMapping_New_Ini;
        std::vector<int> mMapping_New_Cont;
        std::vector<double> mNeighbourDelta;
        std::vector<int> mNeighbourFailureId;
        ParticleWeakVectorType mTempNeighbours;
        
        std::vector<double>               mTempNeighboursDelta;
        std::vector<int>                  mTempNeighboursFailureId;
        std::vector<int>                  mTempNeighboursMapping;
        std::vector<int>                  mTempContNeighboursMapping;
        
        Vector mHistDist;
  
        
        //fem neighbour information
        std::vector<double>               mFemNeighbourDelta;
        std::vector<double>               mFemTempNeighboursDelta;
        
        std::vector<DEMWall*>            mFemTempNeighbours;
        
        std::vector<int>                  mFemIniNeighbourIds;
        Vector                            mFemIniNeighbourDelta;
        std::vector<int>                  mFemMappingNewIni;
        std::vector<int>                  mFemTempNeighboursMapping;
        std::vector<array_1d<double, 3> > mFemTempNeighboursContactForces;
        std::vector<int> mFemTempNeighboursIds;
        
        
        //std::vector<DEMContinuumConstitutiveLaw::Pointer> mContinuumConstitutiveLawArray;
        //std::vector<DEMContinuumConstitutiveLaw::Pointer> mContinuumConstitutiveLawArray;

        
        //Non-linear
         double mN1;
         double mN2;
         double mN3;
         double mC1;
         double mC2;
         double mC3;
                            
        
        //FOR DEM_FEM APP
        
        void ComputeParticleBlockContactForce_With_Rotation();
        void ComputeParticleBlockContactForce_Without_Rotation();
        void FindContactFaceOfBlockForParticle(ParticleWeakIteratorType rObj_2, int & RightFace, double LocalCoordSystem[3][3], double Coeff[4],double &DistPToB);

        
        //double mDampType;
        //double mTimeStep;




        //std::vector<double> mForce;


        ///@name Protected static Member Variables
      ///@{ 
        
        
      ///@} 
      ///@name Protected member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Protected Operators
      ///@{ 
        
        
      ///@} 
      ///@name Protected Operations
      ///@{ 
        
        
      ///@} 
      ///@name Protected  Access 
      ///@{ 
        
        
      ///@}      
      ///@name Protected Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Protected LifeCycle 
      ///@{ 
      
            
      ///@}
      
    private:


      ///@name Static Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
     

       
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
        
        
      ///@} 
      ///@name Private  Access 
      ///@{ 
        
        
      ///@}    
      ///@name Private Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{


      ///@}
      ///@name Serialization
      ///@{

      friend class Serializer;

      virtual void save(Serializer& rSerializer) const
      {
          KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle );
          rSerializer.save("mContinuumGroup",mContinuumGroup);
          //rSerializer.save("mContinuumIniNeighbourElements",mContinuumIniNeighbourElements);
          //rSerializer.save("mBondElements",mBondElements);                              
      }

      virtual void load(Serializer& rSerializer)
      {
          KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle );
          rSerializer.load("mContinuumGroup",mContinuumGroup);
          //rSerializer.load("mContinuumIniNeighbourElements",mContinuumIniNeighbourElements);
          //rSerializer.load("mBondElements",mBondElements);
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


