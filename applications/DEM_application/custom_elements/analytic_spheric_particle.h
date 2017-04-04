//
// Author: Guillermo Casas gcasas@cimne.upc.edu
//


#if !defined(KRATOS_ANALYTIC_SPHERIC_PARTICLE_H_INCLUDED)
#define  KRATOS_ANALYTIC_SPHERIC_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include "spheric_particle.h"


namespace Kratos
{

class KRATOS_API(DEM_APPLICATION) AnalyticSphericParticle : public SphericParticle
{
public:

/// Pointer definition of AnalyticSphericParticle
KRATOS_CLASS_POINTER_DEFINITION(AnalyticSphericParticle);

typedef WeakPointerVector<Condition> ConditionWeakVectorType;
typedef WeakPointerVector<Condition >::iterator ConditionWeakIteratorType;

typedef WeakPointerVector<Element> ParticleWeakVectorType;
typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;

/// Default constructor.
AnalyticSphericParticle( IndexType NewId, GeometryType::Pointer pGeometry );
AnalyticSphericParticle( IndexType NewId, NodesArrayType const& ThisNodes);
AnalyticSphericParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

/// Destructor.
virtual ~AnalyticSphericParticle(){};

/// Turn back information as a string.
std::string Info() const override
{
std::stringstream buffer;
buffer << "AnalyticSphericParticle" ;
return buffer.str();
}

/// Print information about this object.
void PrintInfo(std::ostream& rOStream) const override {rOStream << "AnalyticSphericParticle";}

/// Print object's data.
void PrintData(std::ostream& rOStream) const override {}

protected:

AnalyticSphericParticle();

private:

friend class Serializer;
std::vector<bool> NeighboursContactStatus;
unsigned int mNumberOfCollidingSpheres;
/*
4 is taken as the maximum number of particles simultaneously coming into contact
with this sphere. Whenever more than 4 particles happen to come into contact at
the same time step, which should be extremely rare in collisional regime, the extra
collisions will not be recorded in the present step, but will in forecomming steps
(with the corresponding slightly increased error in measurement) unless an
extraordinarily short contact time (of only 1 time step) takes place.
*/
array_1d<int, 4> mCollidingIds;
array_1d<double, 4> mCollidingRadii;
array_1d<double, 4> mCollidingNormalVelocities;
array_1d<double, 4> mCollidingTangentialVelocities;

void save(Serializer& rSerializer) const override
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DiscreteElement );
    rSerializer.save("mRadius",mRadius);
    rSerializer.save("mSearchRadius", mSearchRadius);
    rSerializer.save("mSearchRadiusWithFem", mSearchRadiusWithFem);
    rSerializer.save("mRealMass",mRealMass);
    rSerializer.save("mClusterId",mClusterId);
    rSerializer.save("mBoundDeltaDispSq",mBoundDeltaDispSq);
    rSerializer.save("HasStressTensor", (int)this->Is(DEMFlags::HAS_STRESS_TENSOR));
    if (this->Is(DEMFlags::HAS_STRESS_TENSOR)){
        rSerializer.save("mSymmStressTensor", mSymmStressTensor);
    }
}

void load(Serializer& rSerializer) override
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DiscreteElement );
    rSerializer.load("mRadius",mRadius);
    rSerializer.load("mSearchRadius", mSearchRadius);
    rSerializer.load("mSearchRadiusWithFem", mSearchRadiusWithFem);
    rSerializer.load("mRealMass",mRealMass);
    rSerializer.load("mClusterId",mClusterId);
    rSerializer.load("mBoundDeltaDispSq",mBoundDeltaDispSq); 
    int aux_int=0;
    rSerializer.load("HasStressTensor", aux_int);
    if(aux_int) this->Set(DEMFlags::HAS_STRESS_TENSOR, true);
    if (this->Is(DEMFlags::HAS_STRESS_TENSOR)){
        mStressTensor  = new Matrix(3,3);
        *mStressTensor = ZeroMatrix(3,3);
        mSymmStressTensor  = new Matrix(3,3);
        *mSymmStressTensor = ZeroMatrix(3,3);
        rSerializer.load("mSymmStressTensor", mSymmStressTensor);
    }
}

}; // Class AnalyticSphericParticle

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
            AnalyticSphericParticle& rThis){ return rIStream;}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
            const AnalyticSphericParticle& rThis)
{
rThis.PrintInfo(rOStream);
rOStream << std::endl;
rThis.PrintData(rOStream);

return rOStream;
}

}  // namespace Kratos.

#endif // KRATOS_ANALYTIC_SPHERIC_PARTICLE_H_INCLUDED  defined
