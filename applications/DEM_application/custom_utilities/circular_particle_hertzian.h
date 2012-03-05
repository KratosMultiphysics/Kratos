//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: G.Casas $
//   Date:                $Date: 2011-6-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//

#ifndef _CIRCULAR_PARTICLE_HERTZIAN_H
#define	_CIRCULAR_PARTICLE_HERTZIAN_H

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
class CircularHertzianParticle: public IndexedObject{

public:

    ///@name Type Definitions
    ///@{

    typedef std::vector<double> DistanceVectorType;
    typedef std::vector<double>::iterator DistanceIteratorType;
    /// Pointer definition of CircularHertzianParticle
    KRATOS_CLASS_POINTER_DEFINITION(CircularHertzianParticle);
    typedef WeakPointerVector<CircularHertzianParticle > ParticleWeakVectorType;
    typedef WeakPointerVector<CircularHertzianParticle >::iterator ParticleWeakIteratorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CircularHertzianParticle():IndexedObject(0){
        mpCenterNode = Node<3>::Pointer(new Node<3>(1,0.0,0.0,0.0));
        mMaterial = 1;
        mRadius = 1.0;
        mDensity = 1.0;
        mYoung = 1000.0;
        mPoisson = 0.5;
        mYoungStar = mYoung / (1.0 - mPoisson * mPoisson);
        mZeta = 0.1;
        mProximity_Tol = 0.00001;
        mMass = M_PI * mDensity * mRadius * mRadius;
        };

    CircularHertzianParticle(double tol, Node<3>::Pointer center):IndexedObject(center->Id()){
        mpCenterNode = center;
        mMaterial = mpCenterNode->GetSolutionStepValue(PARTICLE_MATERIAL);
        mRadius = mpCenterNode->GetSolutionStepValue(RADIUS);
        mDensity = mpCenterNode->GetSolutionStepValue(PARTICLE_DENSITY);
        mYoung = mpCenterNode->GetSolutionStepValue(YOUNG_MODULUS);
        mPoisson = mpCenterNode->GetSolutionStepValue(POISSON_RATIO);
        mYoungStar = mYoung / (1.0 - mPoisson * mPoisson);
        mZeta = mpCenterNode->GetSolutionStepValue(PARTICLE_ZETA);
        mMass = M_PI * mDensity * mRadius * mRadius;
        mProximity_Tol = tol * mRadius;
        mInitialPosition = mpCenterNode->GetInitialPosition().Coordinates();
        mOldPosition = mInitialPosition + mpCenterNode->FastGetSolutionStepValue(DISPLACEMENT, 1);
        };

  /// Destructor.
    virtual ~CircularHertzianParticle(){};

  //  Member variables

  ///@}
  ///@name Operators
  ///@{

    double operator()( CircularHertzianParticle const& p1, CircularHertzianParticle const& p2 ){
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
    int& GetNumberOfNeighbours(){return(mNumberOfNeighbours);};
    int& GetMaterial(){return (mMaterial);};
    double& GetRadius(){return (mRadius);};
    double& GetMass(){return (mMass);};
    double& GetDensity(){return (mDensity);};
    double& GetYoung(){return (mYoung);};
    double& GetYoungStar(){return (mYoungStar);};
    double& GetPoisson(){return (mPoisson);};
    double& GetZeta(){return (mZeta);};
    array_1d<double,3>& GetDisplacement(){return (mpCenterNode->FastGetSolutionStepValue(DISPLACEMENT));};
    array_1d<double,3>& GetInitialPosition(){return (mInitialPosition);};
    array_1d<double,3>& GetPosition(){return (mpCenterNode->Coordinates());};
    array_1d<double,3>& GetVelocity(){return (mpCenterNode->FastGetSolutionStepValue(VELOCITY));};
    array_1d<double,3>& GetForce(){return (mpCenterNode->FastGetSolutionStepValue(FORCE));};
    Node<3>::Pointer& GetPointerToCenterNode(){return(mpCenterNode);};
    ParticleWeakVectorType& GetNeighbours(){return(mNeighbours);};
    DistanceVectorType& GetDistancesToNeighbours(){return(mDistancesToNeighbours);};
    DistanceVectorType& GetDistancesToContacts(){return(mDistancesToContacts);};
    void ClearTangentialDisplacementOfNeighbours(){};
    void GetTangentialDisplacementOfNeighbours(DistanceVectorType& vector){};

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{
    /// Turn back information as a string.
    virtual std::string Info() const{
        return "CircularHertzianParticle";
        }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const{
        rOStream << "CircularHertzianParticle";
        }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const{
        //       rOStream << Perfix << "Center" << SearchUtils::PointerDistance(mPointsBegin, mPointsEnd) << "] : ";
        //       for(IteratorType i = mPointsBegin ; i != mPointsEnd ; i++)
        //       rOStream << **i << "    ";
        //       rOStream << std::endl;
        }

    CircularHertzianParticle(const CircularHertzianParticle& rOtherParticle){
        this->SetId(rOtherParticle.Id());
        this->mpCenterNode = rOtherParticle.mpCenterNode;
        this->mMaterial = rOtherParticle.mMaterial;
        this->mRadius = rOtherParticle.mRadius;
        this->mMass = rOtherParticle.mMass;
        this->mYoung = rOtherParticle.mYoung;
        this->mYoungStar = rOtherParticle.mYoungStar;
        this->mPoisson = rOtherParticle.mPoisson;
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
    int mNumberOfNeighbours;
    int mMaterial;
    double mRadius;
    double mMass;
    double mYoung;
    double mYoungStar;
    double mPoisson;
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
    CircularHertzianParticle& operator=(const CircularHertzianParticle& rOtherParticle);
    // Copy constructor.
    //      CircularHertzianParticle(const CircularHertzianParticle& rOtherParticle);


    ///@}

}; // Class CircularHertzianParticle

    ///@}

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, CircularHertzianParticle& rThis);

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const CircularHertzianParticle& rThis){
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);
        return rOStream;
        }
    ///@}

}  // namespace Kratos.
#endif // _CIRCULAR_PARTICLE_HERTZIAN_H  defined


