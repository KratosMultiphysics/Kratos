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
typedef SphericParticle BaseType;
typedef BaseType::ParticleDataBuffer BaseBufferType;
typedef std::unique_ptr<BaseType::ParticleDataBuffer> BaseBufferPointerType;

/// Default constructor.
AnalyticSphericParticle();
AnalyticSphericParticle( IndexType NewId, GeometryType::Pointer pGeometry);
AnalyticSphericParticle( IndexType NewId, NodesArrayType const& ThisNodes);
AnalyticSphericParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);
AnalyticSphericParticle(Element::Pointer p_spheric_particle);

Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

/// Destructor.
virtual ~AnalyticSphericParticle(){}

AnalyticSphericParticle& operator=(const AnalyticSphericParticle& rOther);

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

int GetNumberOfCollisions();
int GetNumberOfCollisionsWithFaces();
int GetNumberOfCollisionsWithEdges();

static const unsigned int mMaxCollidingSpheres = 4;
static const unsigned int mMaxCollidingFaceSpheres = 4;

array_1d<int, 4> &GetCollidingIds();
array_1d<double, 4> &GetCollidingNormalRelativeVelocity();
array_1d<double, 4> &GetCollidingTangentialRelativeVelocity();
array_1d<int, 4> &GetCollidingFaceIds();
array_1d<double, 4> &GetCollidingFaceNormalRelativeVelocity();
array_1d<double, 4> &GetCollidingFaceTangentialRelativeVelocity();
array_1d<double, 4> &GetCollidingLinearImpulse();

/*
void GetCollidingNormalRelativeVelocity(array_1d<double, 4>& colliding_normal_vel);
void GetCollidingTangentialRelativeVelocity(array_1d<double, 4>& colliding_tangential_vel);
void GetCollidingFaceIds(array_1d<int, 4>& colliding_ids);
void GetCollidingFaceNormalRelativeVelocity(array_1d<double, 4>& colliding_normal_vel);
void GetCollidingFaceTangentialRelativeVelocity(array_1d<double, 4>& colliding_tangential_vel);
void GetCollidingLinearImpulse(array_1d<double, 4>& colliding_linear_impulse);
*/

protected:

class ParticleDataBuffer: public SphericParticle::ParticleDataBuffer
{
public:

ParticleDataBuffer(SphericParticle* p_this_particle): SphericParticle::ParticleDataBuffer(p_this_particle){}

virtual ~ParticleDataBuffer(){}

std::vector<int> mCurrentContactingNeighbourIds;
std::vector<int> mCurrentContactingFaceNeighbourIds;
};

std::unique_ptr<SphericParticle::ParticleDataBuffer> CreateParticleDataBuffer(SphericParticle* p_this_particle) override
{
    ClearImpactMemberVariables();
    return std::unique_ptr<SphericParticle::ParticleDataBuffer>(new ParticleDataBuffer(p_this_particle));
}

void PushBackIdToContactingNeighbours(BaseBufferType &data_buffer, int id);

void PushBackIdToContactingFaceNeighbours(BaseBufferType & data_buffer, int p_wall_id);


void ClearNeighbours(BaseBufferType & data_buffer);

private:

friend class Serializer;
std::vector<bool> NeighboursContactStatus;
unsigned int mNumberOfCollidingSpheres;
unsigned int mNumberOfCollidingSpheresWithFaces;
unsigned int mNumberOfCollidingSpheresWithEdges;
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
array_1d<double, 4> mCollidingLinearImpulse;
std::vector<int> mContactingNeighbourIds;
//std::unique_ptr<ParticleDataBuffer> mpDataBuffer;

array_1d<int, 4> mCollidingFaceIds;
array_1d<double, 4> mCollidingFaceNormalVelocities;
array_1d<double, 4> mCollidingFaceTangentialVelocities;
array_1d<double, 4> mCollidingFaceSecondTangentialVelocities;
array_1d<char, 4> mCollidingFaceCollisionTypes;
std::vector<int> mContactingFaceNeighbourIds;


ParticleDataBuffer* GetPointerToDerivedDataBuffer(BaseBufferType& data_buffer)
{
  BaseBufferType *p_raw_data_buffer = &data_buffer;
  return static_cast<ParticleDataBuffer*>(p_raw_data_buffer);
}

void ClearImpactMemberVariables();

void FinalizeForceComputation(BaseBufferType & data_buffer) override;

void EvaluateBallToBallForcesForPositiveIndentiations(SphericParticle::ParticleDataBuffer & data_buffer,
                                                       const ProcessInfo& r_process_info,
                                                       double LocalElasticContactForce[3],
                                                       double DeltDisp[3],
                                                       double LocalDeltDisp[3],
                                                       double RelVel[3],
                                                       const double indentation,
                                                       double ViscoDampingLocalContactForce[3],
                                                       double& cohesive_force,
                                                       SphericParticle* p_neighbour_element,
                                                       bool& sliding,
                                                       double LocalCoordSystem[3][3],
                                                       double OldLocalCoordSystem[3][3],
                                                       array_1d<double, 3>& neighbour_elastic_contact_force) override;

void ComputeBallToRigidFaceContactForce(SphericParticle::ParticleDataBuffer & data_buffer,
                                                         array_1d<double, 3>& r_elastic_force,
                                                         array_1d<double, 3>& r_contact_force,
                                                         double& RollingResistance,
                                                         array_1d<double, 3>& rigid_element_force,
                                                         ProcessInfo& r_process_info,                                                        
                                                         int search_control) override;


bool IsNewNeighbour(const int neighbour_id);

bool IsNewFaceNeighbour(const int neighbour_id);

void RecordNewImpact(BaseType::ParticleDataBuffer & data_buffer);

void RecordNewFaceImpact(BaseType::ParticleDataBuffer & data_buffer);

void save(Serializer& rSerializer) const override
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle);
}

void load(Serializer& rSerializer) override
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle);
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
