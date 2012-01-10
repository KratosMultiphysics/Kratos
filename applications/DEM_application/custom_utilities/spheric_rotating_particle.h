/* 
 * File:   spheric_particle.h
 * Author: gcasas
 *
 * Created on 16 de septiembre de 2011, 11:35
 */

#ifndef _SPHERIC_ROTATING_PARTICLE_H
#define	_SPHERIC_ROTATING_PARTICLE_H

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
class SphericRotatingParticle: public IndexedObject{

public:

    ///@name Type Definitions
    ///@{

    typedef std::vector<double> DistanceVectorType;
    typedef std::vector<array_1d<double, 3 > > TangDisplacementsVectorType;
    typedef std::vector<double > ForceNormsVectorType;
    typedef std::vector<bool > BoolSlidingContactsVectorType;
    typedef std::vector<double>::iterator DistanceIteratorType;
    typedef std::vector<array_1d<double, 3 > >::iterator TangDisplacementsIteratorType;
    typedef std::vector<double >::iterator ForceNormsIteratorType;
    typedef std::vector<bool >::iterator BoolSlidingContactsIteratorType;
    /// Pointer definition of SphericRotatingParticle
    KRATOS_CLASS_POINTER_DEFINITION(SphericRotatingParticle);
    typedef WeakPointerVector<SphericRotatingParticle > ParticleWeakVectorType;
    typedef WeakPointerVector<SphericRotatingParticle >::iterator ParticleWeakIteratorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SphericRotatingParticle():IndexedObject(0){
        mpCenterNode = Node<3>::Pointer(new Node<3>(1,0.0,0.0,0.0));
        mRadius = 1.0;
        mDensity = 1.0;
        mYoung = 1000.0;
        mPoisson = 0.5;
        mNormalStiffness = 1000.0;
        mTangentialStiffness = mNormalStiffness * (1.0 - mYoung) / (1 - 0.5 * mYoung);
        mStaticFriction = 0.5;
        mDynamicFriction = 0.3;
        mRestitutionCoef = 1.0;
        mZeta = 0.0; //mZeta = 1 for critical damping
        mProximity_Tol = 0.00000001;
        mMass = 1.33333333333333333333333 * M_PI * mDensity * mRadius * mRadius * mRadius;
        mInertia = 0.4 * mMass * mRadius * mRadius;
        mOldPosition = ZeroVector(0);
        };

    SphericRotatingParticle(double tol, Node<3>::Pointer center):IndexedObject(center->Id()){
        mpCenterNode = center;
        mRadius = mpCenterNode->GetSolutionStepValue(RADIUS);
        mDensity = mpCenterNode->GetSolutionStepValue(PARTICLE_DENSITY);
        mYoung = mpCenterNode->GetSolutionStepValue(YOUNG_MODULUS);
        mPoisson = mpCenterNode->GetSolutionStepValue(POISSON_RATIO);
        mNormalStiffness = mpCenterNode->GetSolutionStepValue(PARTICLE_STIFFNESS);
        mTangentialStiffness = mNormalStiffness * (1.0 - mYoung) / (1 - 0.5 * mYoung); //Elastic solid mechanics analysis of Mindlin (1949)
        mStaticFriction = mpCenterNode->GetSolutionStepValue(PARTICLE_STATIC_FRICTION_COEF);
        mDynamicFriction = mpCenterNode->GetSolutionStepValue(PARTICLE_DYNAMIC_FRICTION_COEF);
        mRestitutionCoef = mpCenterNode->GetSolutionStepValue(PARTICLE_COEF_RESTITUTION);
        if(mRestitutionCoef < 0.0){mZeta = mpCenterNode->GetSolutionStepValue(PARTICLE_ZETA);}
        else{
            double beta = M_PI / log(mRestitutionCoef);
            if(mRestitutionCoef > 0.9999){mZeta = 0.0;}
            else{
                mZeta = sqrt(1 / (1 + beta * beta)); //From hertzian theory
                }
            }
        mMass = 1.33333333333333333333333 * M_PI * mDensity * mRadius * mRadius * mRadius;
        mInertia = 0.4 * mMass * mRadius * mRadius;
        mProximity_Tol = tol * mRadius;
        mInitialPosition = mpCenterNode->GetInitialPosition().Coordinates();
        mOldPosition = mInitialPosition + mpCenterNode->FastGetSolutionStepValue(DISPLACEMENT, 1);
        };

  /// Destructor.
    virtual ~SphericRotatingParticle(){};

  //  Member variables

  ///@}
  ///@name Operators
  ///@{

    double operator()(SphericRotatingParticle const& p1, SphericRotatingParticle const& p2 ){
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

    void ComputeForcesOnCenterNode(double dt, array_1d<double, 3 >& gravity);
    void UpdateContactsList();
    void ComputeNewTangentialDisplacements();
    ///@}
    ///@name Access
    ///@{
    double& GetNumberOfNeighbours(){return(mpCenterNode->FastGetSolutionStepValue(NUMBER_OF_NEIGHBOURS));};
    int GetNumberOfContactingNeighbours(){return(mContactingNeighbours.size());};
    double& GetRadius(){return (mRadius);};
    double& GetMass(){return (mMass);};
    double& GetInertia(){return (mInertia);};
    double& GetDensity(){return (mDensity);};
    double& GetDynamicFriction(){return (mDynamicFriction);};
    double& GetStaticFriction(){return (mStaticFriction);};
    double& GetNormalStiffness(){return (mNormalStiffness);};
    double& GetTangentialStiffness(){return (mTangentialStiffness);};
    double& GetYoung(){return (mYoung);};
    double& GetPoisson(){return (mPoisson);};
    double& GetRestitutionCoef(){return (mRestitutionCoef);};
    double& GetZeta(){return (mZeta);};
    double& GetProxTol(){return (mProximity_Tol);};
    array_1d<double,3>& GetDisplacement(){return (mpCenterNode->FastGetSolutionStepValue(DISPLACEMENT));};
    array_1d<double,3>& GetInitialPosition(){return (mInitialPosition);};
    array_1d<double,3>& GetPosition(){return (mpCenterNode->Coordinates());};
    array_1d<double,3>& GetOldPosition(){return (mOldPosition);};
    array_1d<double,3>& GetRotation(){return (mpCenterNode->FastGetSolutionStepValue(ROTATION, 0));};
    array_1d<double,3>& GetOldRotation(){return (mpCenterNode->FastGetSolutionStepValue(ROTATION, 1));};
    array_1d<double,3>& GetVelocity(){return (mpCenterNode->FastGetSolutionStepValue(VELOCITY));};
    array_1d<double,3>& GetOldVelocity(){return (mpCenterNode->FastGetSolutionStepValue(VELOCITY, 1));};
    array_1d<double,3>& GetAngularVelocity(){return (mpCenterNode->FastGetSolutionStepValue(ANGULAR_VELOCITY));};
    array_1d<double,3>& GetForce(){return (mpCenterNode->FastGetSolutionStepValue(FORCE));};
    array_1d<double,3>& GetMoment(){return (mpCenterNode->FastGetSolutionStepValue(MOMENT));};
    Node<3>::Pointer& GetPointerToCenterNode(){return(mpCenterNode);};
    ParticleWeakVectorType& GetNeighbours(){return(mNeighbours);};
    ParticleWeakVectorType& GetContactingNeighbours(){return(mContactingNeighbours);};
    DistanceVectorType& GetDistancesToNeighbours(){return(mDistancesToNeighbours);};
    DistanceVectorType& GetDistancesToContacts(){return(mDistancesToContacts);};
    void ClearTangentialDisplacementOfNeighbours(){mTangentialDisplacementsOfNeighbours.clear();};
    void GetTangentialDisplacementOfNeighbours(TangDisplacementsVectorType& vector){vector = mTangentialDisplacementsOfNeighbours;};
    void GetMaxTangForcesSquared(ForceNormsVectorType& vector){vector = mMaxTangForces;};
    void GetSlidingFlagsOfNeighbours(BoolSlidingContactsVectorType& vector){vector = mSlidingFlagsOfNeighbours;};

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{
    /// Turn back information as a string.
    virtual std::string Info() const{
        return "SphericRotatingParticle";
        }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const{
        rOStream << "SphericRotatingParticle";
        }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const{
        //       rOStream << Perfix << "Center" << SearchUtils::PointerDistance(mPointsBegin, mPointsEnd) << "] : ";
        //       for(IteratorType i = mPointsBegin ; i != mPointsEnd ; i++)
        //       rOStream << **i << "    ";
        //       rOStream << std::endl;
        }

    SphericRotatingParticle(const SphericRotatingParticle& rOtherParticle){
        this->SetId(rOtherParticle.Id());
        this->mpCenterNode = rOtherParticle.mpCenterNode;
        this->mRadius = rOtherParticle.mRadius;
        this->mMass = rOtherParticle.mMass;
        this->mStaticFriction = rOtherParticle.mStaticFriction;
        this->mDynamicFriction = rOtherParticle.mDynamicFriction;
        this->mInertia = rOtherParticle.mInertia;
        this->mNormalStiffness = rOtherParticle.mNormalStiffness;
        this->mTangentialStiffness = rOtherParticle.mTangentialStiffness;
        this->mYoung = rOtherParticle.mYoung;
        this->mPoisson = rOtherParticle.mPoisson;
        this->mRestitutionCoef = rOtherParticle.mRestitutionCoef;
        this->mZeta = rOtherParticle.mZeta;
        this->mProximity_Tol = rOtherParticle.mProximity_Tol;
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
    double mMass;
    double mInertia;
    double mNormalStiffness;
    double mTangentialStiffness;
    double mYoung;
    double mPoisson;
    double mRestitutionCoef;
    double mZeta;
    double mDensity;
    double mStaticFriction;
    double mDynamicFriction;
    double mProximity_Tol;
    Node<3>::Pointer mpCenterNode;
    array_1d<double, 3> mInitialPosition;
    array_1d<double, 3> mPosition;
    array_1d<double, 3> mOldPosition;
    array_1d<double, 3> mVelocity;
    array_1d<double, 3> mForce;
    ParticleWeakVectorType mNeighbours;
    ParticleWeakVectorType mContactingNeighbours;
    DistanceVectorType mDistancesToNeighbours;
    DistanceVectorType mDistancesToContacts;
    TangDisplacementsVectorType mTangentialDisplacementsOfNeighbours;
    ForceNormsVectorType mMaxTangForces;
    BoolSlidingContactsVectorType mSlidingFlagsOfNeighbours;

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
    SphericRotatingParticle& operator=(const SphericRotatingParticle& rOtherParticle);
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
    inline std::istream& operator >> (std::istream& rIStream, SphericRotatingParticle& rThis);

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const SphericRotatingParticle& rThis){
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);
        return rOStream;
        }
    ///@}
}
#endif	/* _SPHERIC_ROTATING_PARTICLE_H */

