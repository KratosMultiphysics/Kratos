/* 
 * File:   spheric_particle.h
 * Author: gcasas
 *
 * Created on 16 de septiembre de 2011, 11:35
 */

#ifndef _SPHERIC_PARTICLE_HERTZIAN_H
#define	_SPHERIC_PARTICLE_HERTZIAN_H

#include "utilities/timer.h"
// System includes
#include <math.h>
// External includes

// Project includes
#include "includes/model_part.h"

//Database includes

namespace Kratos{

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
class SphericHertzianParticle: public IndexedObject{

public:

    ///@name Type Definitions
    ///@{

    typedef std::vector<double> DistanceVectorType;
    typedef std::vector<double>::iterator DistanceIteratorType;
    /// Pointer definition of SphericHertzianParticle
    KRATOS_CLASS_POINTER_DEFINITION(SphericHertzianParticle);
    typedef WeakPointerVector<SphericHertzianParticle > ParticleWeakVectorType;
    typedef WeakPointerVector<SphericHertzianParticle >::iterator ParticleWeakIteratorType;

   //Cfeng:For Shear Contact Force,020312
    typedef std::vector< array_1d<double,3> > ContactForceVectorType;
    typedef std::vector< int > ContactFailureIdType;
    typedef std::vector< double > ContactInitialDeltaType;
    typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SphericHertzianParticle():IndexedObject(0){
        mpCenterNode = Node<3>::Pointer(new Node<3>(1,0.0,0.0,0.0));
        mMaterial = 1;
        mContinuumGroup = 0;
        mFailureId = 1; // detahced by default.
        mRadius = 1.0;
        mDensity = 1.0;
        mYoung = 1000.0;
        mPoisson = 0.5;
        mYoungStar = mYoung / (1.0 - mPoisson * mPoisson);
        mRestitutionCoef = 1.0;
        mZeta = 0.1;
        mProximity_Tol = 0.000000001;
        mMass = 1.33333333333333333333333 * M_PI * mDensity * mRadius * mRadius * mRadius;
        mTolerance = 0;
        mTension = 0.0;
        mFriction = 20.0;
        mCohesion = 0.0;
        mLocalDampRatio = 0.2;
	
       
        };

    SphericHertzianParticle(double tol, Node<3>::Pointer center):IndexedObject(center->Id()){
        mpCenterNode = center;
        mMaterial = mpCenterNode->GetSolutionStepValue(PARTICLE_MATERIAL);
        mContinuumGroup = mpCenterNode->GetSolutionStepValue(PARTICLE_CONTINUUM);
        mFailureId = !(mContinuumGroup); // if ContinuumGroup != 0 --> mFailureId = 0; mFailureId is 1 when mContinuumGroup=0;
        mRadius = mpCenterNode->GetSolutionStepValue(RADIUS);
        mDensity = mpCenterNode->GetSolutionStepValue(PARTICLE_DENSITY);
        mYoung = mpCenterNode->GetSolutionStepValue(YOUNG_MODULUS);
        mPoisson = mpCenterNode->GetSolutionStepValue(POISSON_RATIO);
        mYoungStar = mYoung / (1.0 - mPoisson * mPoisson);
        mRestitutionCoef = mpCenterNode->GetSolutionStepValue(PARTICLE_COEF_RESTITUTION);
        if(mRestitutionCoef < 0.0){mZeta = mpCenterNode->GetSolutionStepValue(PARTICLE_ZETA);}
        else{
            if(mRestitutionCoef > 0.9999){mZeta = 0.0;}
            else{
                mZeta = -sqrt(5.0) * log(mRestitutionCoef) / (sqrt(log(mRestitutionCoef) * log(mRestitutionCoef) + M_PI * M_PI));
                }
            }
        mMass = 1.33333333333333333333333 * M_PI * mDensity * mRadius * mRadius * mRadius;
        mProximity_Tol = tol * mRadius;
        mInitialPosition = mpCenterNode->GetInitialPosition().Coordinates();
        mOldPosition = mInitialPosition + mpCenterNode->FastGetSolutionStepValue(DISPLACEMENT, 1);
        mTolerance = 0; //M: no se com passarli...

        mTension = mpCenterNode->GetSolutionStepValue(PARTICLE_TENSION);
        mFriction = mpCenterNode->GetSolutionStepValue(PARTICLE_FRICTION);
        mCohesion = mpCenterNode->GetSolutionStepValue(PARTICLE_COHESION);
        mLocalDampRatio = mpCenterNode->GetSolutionStepValue(PARTICLE_LOCAL_DAMP_RATIO);
      
        };

  /// Destructor.
    virtual ~SphericHertzianParticle(){};

  //  Member variables

  ///@}
  ///@name Operators
  ///@{

    double operator()(SphericHertzianParticle const& p1, SphericHertzianParticle const& p2 ){
        double dist_squared = ((*(p1.mpCenterNode))[0] - (*(p2.mpCenterNode))[0]) * ((*(p1.mpCenterNode))[0] - (*(p2.mpCenterNode))[0])  + ((*(p1.mpCenterNode))[1] - (*(p2.mpCenterNode))[1]) * ((*(p1.mpCenterNode))[1] - (*(p2.mpCenterNode))[1]);
        double inter_particle_signed_distance = sqrt(dist_squared) - p1.mRadius - p2.mRadius;
        return (inter_particle_signed_distance);
        };

    double & operator[](std::size_t dimension){
        return ((*mpCenterNode)[dimension]);
        }

    double & operator[](std::size_t dimension) const{
        return ((*mpCenterNode)[dimension]);
        }

    ///@}
    ///@name Operations
    ///@{

    void SetInitialContacts();
    void AddContinuumContacts();
   /*
    void CaseSelector(int type_id, int damp_id, double dt, array_1d<double, 3 >& gravity);
    void NeighbourTypeSelector (SphericHertzianParticle::ParticleWeakVectorType& rGenericNeighbours, int type_id, int damp_id, int caseIdentifier, double dt, array_1d<double,3>& gravity);
    */
   
    void ComputeForcesGeneral ( int type_id, int damp_id, double dt_input, array_1d<double,3>& gravity);

    void AfterForceCalculation(int damp_id);

////Cfeng: Add some functions to calculate the direction of the contact plane
    void ComputeContactLocalCoordSystem(double NormalDirection[3], double LocalCoordSystem[3][3]);
    void norm(double Vector[3]);
    void VectorGlobal2Local(double LocalCoordSystem[3][3], double GlobalVector[3], double LocalVector[3]);
    void VectorLocal2Global(double LocalCoordSystem[3][3], double LocalVector[3], double GlobalVector[3]);
    double DotProduct(double Vector1[3], double Vector2[3]);
    void CrossProduct(double u[3], double v[3], double ReturnVector[3]);


	///@}
    ///@name Access
    ///@{
    double& GetRadius(){return (mRadius);};
    int& GetMaterial(){return (mMaterial);};
    int& GetContinuumGroup(){return (mContinuumGroup);};
     /////////////////////////////////////M: sta repetit com PARTICLE_FAILURE_ID
    double& GetMass(){return (mMass);};
    double& GetDensity(){return (mDensity);};
    double& GetYoung(){return (mYoung);};
    double& GetYoungStar(){return (mYoungStar);};
    double& GetPoisson(){return (mPoisson);};
    double& GetRestitutionCoef(){return (mRestitutionCoef);};
    double& GetZeta(){return (mZeta);};
    double& GetNumberOfNeighbours(){return(mpCenterNode->FastGetSolutionStepValue(NUMBER_OF_NEIGHBOURS));};
    array_1d<double,3>& GetDisplacement(){return (mpCenterNode->FastGetSolutionStepValue(DISPLACEMENT));};
    array_1d<double,3>& GetInitialPosition(){return (mInitialPosition);};
    array_1d<double,3>& GetPosition(){return (mpCenterNode->Coordinates());};
    array_1d<double,3>& GetVelocity(){return (mpCenterNode->FastGetSolutionStepValue(VELOCITY));};
    array_1d<double,3>& GetForce(){return (mpCenterNode->FastGetSolutionStepValue(FORCE));};
    Node<3>::Pointer& GetPointerToCenterNode(){return(mpCenterNode);};
    ParticleWeakVectorType& GetNeighbours(){return(mNeighbours);};
    ParticleWeakVectorType& GetInitialNeighbours(){return(mInitialNeighbours);};
    std::vector<double>& GetInitialDelta(){return(mInitialDelta);};
    DistanceVectorType& GetDistancesToNeighbours(){return(mDistancesToNeighbours);};
    DistanceVectorType& GetDistancesToContacts(){return(mDistancesToContacts);};
    void GetTangentialDisplacementOfNeighbours(DistanceVectorType& vector){};

    ContactForceVectorType&  GetContactForces(){return(mContactForces);};
    ContactFailureIdType& GetContactFailureId() {return(mContactFailureId);};
    ContactInitialDeltaType& GetContactInitialDelta(){return (mContactInitialDelta);};
        

    double& GetTension(){return (mpCenterNode->FastGetSolutionStepValue(PARTICLE_TENSION));};
    double& GetFriction(){return (mpCenterNode->FastGetSolutionStepValue(PARTICLE_FRICTION));};
    double& GetCohesion(){return (mpCenterNode->FastGetSolutionStepValue(PARTICLE_COHESION));};
    double& GetLocalDampRatio(){return (mpCenterNode->FastGetSolutionStepValue(PARTICLE_LOCAL_DAMP_RATIO));};
    double& GetFailureId(){return (mpCenterNode->FastGetSolutionStepValue(PARTICLE_FAILURE_ID));};/////////////////////////////////////M: sta repetit com mFaiuireId


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{
    /// Turn back information as a string.
    virtual std::string Info() const{
        return "SphericHertzianParticle";
        }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const{
        rOStream << "SphericHertzianParticle";
        }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const{
        //       rOStream << Perfix << "Center" << SearchUtils::PointerDistance(mPointsBegin, mPointsEnd) << "] : ";
        //       for(IteratorType i = mPointsBegin ; i != mPointsEnd ; i++)
        //       rOStream << **i << "    ";
        //       rOStream << std::endl;
        }

    SphericHertzianParticle(const SphericHertzianParticle& rOtherParticle){
        this->SetId(rOtherParticle.Id());
        this->mpCenterNode = rOtherParticle.mpCenterNode;
        this->mMaterial = rOtherParticle.mMaterial;
        this->mContinuumGroup =rOtherParticle.mContinuumGroup;
        this->mRadius = rOtherParticle.mRadius;
        this->mFailureId = rOtherParticle.mFailureId;
        this->mMass = rOtherParticle.mMass;
        this->mYoung = rOtherParticle.mYoung;
        this->mPoisson = rOtherParticle.mPoisson;
        this->mYoungStar = rOtherParticle.mYoungStar;
        this->mRestitutionCoef = rOtherParticle.mRestitutionCoef;
        this->mZeta = rOtherParticle.mZeta;
        this->mProximity_Tol = rOtherParticle.mProximity_Tol;

        this->mTension = rOtherParticle.mTension;
        this->mCohesion = rOtherParticle.mCohesion;
        this->mFriction = rOtherParticle.mFriction;
        this->mLocalDampRatio = rOtherParticle.mLocalDampRatio;
        this->mFailureId = rOtherParticle.mFailureId;
        }

    /// Print information about this object.
    //     virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    //     virtual void PrintData(std::ostream& rOStream) const;

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
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
    double mRadius;
    int mMaterial;
    int mContinuumGroup;
    double mMass;
    double mYoung;
    double mYoungStar;
    double mPoisson;
    double mRestitutionCoef;
    double mZeta;
    double mDensity;
    double mProximity_Tol;
    double mTolerance;
    Node<3>::Pointer mpCenterNode;
    array_1d<double, 3> mInitialPosition;
    array_1d<double, 3> mPosition;
    array_1d<double, 3> mOldPosition;
    array_1d<double, 3> mVelocity;
    array_1d<double, 3> mForce;
    ParticleWeakVectorType mNeighbours;
    ParticleWeakVectorType mInitialNeighbours;
    ParticleWeakVectorType mGenericNeighbours;
    std::vector<double> mInitialDelta;
   
    ContactForceVectorType mContactForces;
    ContactFailureIdType mContactFailureId;
    ContactInitialDeltaType mContactInitialDelta;

    
    double mCohesion;
    double mFriction;
    double mTension;
    double mLocalDampRatio;
    int mFailureId;      // the dominating type of failure of the contacts.
    
    DistanceVectorType mDistancesToNeighbours;
    DistanceVectorType mDistancesToContacts;

    //M: nomes per probes.
    double mTestVariable;
    double maxForce;
    double minForce;

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

    /// Assignment operator.
    SphericHertzianParticle& operator=(const SphericHertzianParticle& rOtherParticle);
    // Copy constructor.
    //      SphericParticleHertzian(const SphericParticleHertzian& rOtherParticle);


    ///@}

}; // Class SphericParticleHertzian

    ///@}

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, SphericHertzianParticle& rThis);

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const SphericHertzianParticle& rThis){
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);
        return rOStream;
        }
    ///@}
}
#endif	/* _SPHERIC_PARTICLE_HERTZIAN_H */

