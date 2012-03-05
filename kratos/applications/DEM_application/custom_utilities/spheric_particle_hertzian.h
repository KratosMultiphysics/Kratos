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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SphericHertzianParticle():IndexedObject(0){
        mpCenterNode = Node<3>::Pointer(new Node<3>(1,0.0,0.0,0.0));
        mMaterial = 1;
        mRadius = 1.0;
        mDensity = 1.0;
        mYoung = 1000.0;
        mPoisson = 0.5;
        mYoungStar = mYoung / (1.0 - mPoisson * mPoisson);
        mRestitutionCoef = 1.0;
        mZeta = 0.1;
        mProximity_Tol = 0.000000001;
        mMass = 1.33333333333333333333333 * M_PI * mDensity * mRadius * mRadius * mRadius;
        };

    SphericHertzianParticle(double tol, Node<3>::Pointer center):IndexedObject(center->Id()){
        mpCenterNode = center;
        mMaterial = mpCenterNode->GetSolutionStepValue(PARTICLE_MATERIAL);
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

    void ComputeForcesOnCenterNode(double dt, array_1d<double, 3 >& gravity);
    ///@}
    ///@name Access
    ///@{
    double& GetRadius(){return (mRadius);};
    int& GetMaterial(){return (mMaterial);};
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
    DistanceVectorType& GetDistancesToNeighbours(){return(mDistancesToNeighbours);};
    DistanceVectorType& GetDistancesToContacts(){return(mDistancesToContacts);};
    void GetTangentialDisplacementOfNeighbours(DistanceVectorType& vector){};

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
        this->mRadius = rOtherParticle.mRadius;
        this->mMass = rOtherParticle.mMass;
        this->mYoung = rOtherParticle.mYoung;
        this->mPoisson = rOtherParticle.mPoisson;
        this->mYoungStar = rOtherParticle.mYoungStar;
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
    int mMaterial;
    double mMass;
    double mYoung;
    double mYoungStar;
    double mPoisson;
    double mRestitutionCoef;
    double mZeta;
    double mDensity;
    double mProximity_Tol;
    Node<3>::Pointer mpCenterNode;
    array_1d<double, 3> mInitialPosition;
    array_1d<double, 3> mPosition;
    array_1d<double, 3> mOldPosition;
    array_1d<double, 3> mVelocity;
    array_1d<double, 3> mForce;
    ParticleWeakVectorType mNeighbours;
    DistanceVectorType mDistancesToNeighbours;
    DistanceVectorType mDistancesToContacts;

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

