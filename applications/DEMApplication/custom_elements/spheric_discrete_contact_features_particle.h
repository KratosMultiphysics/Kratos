//
// Author: Joaquín Irazábal jirazabal@cimne.upc.edu
//


#if !defined(KRATOS_SPHERIC_DISCRETE_CONTACT_FEATURES_PARTICLE_H_INCLUDED)
#define  KRATOS_SPHERIC_DISCRETE_CONTACT_FEATURES_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include "spheric_particle.h"

namespace Kratos
{

class KRATOS_API(DEM_APPLICATION) SphericDiscreteContactFeaturesParticle : public SphericParticle
{
public:

/// Pointer definition of SphericDiscreteContactFeaturesParticle
KRATOS_CLASS_POINTER_DEFINITION(SphericDiscreteContactFeaturesParticle);

/// Default constructor.
SphericDiscreteContactFeaturesParticle();
SphericDiscreteContactFeaturesParticle(IndexType NewId, GeometryType::Pointer pGeometry);
SphericDiscreteContactFeaturesParticle(IndexType NewId, NodesArrayType const& ThisNodes);
SphericDiscreteContactFeaturesParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);
SphericDiscreteContactFeaturesParticle(Element::Pointer p_spheric_particle);

Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

/// Destructor.
virtual ~SphericDiscreteContactFeaturesParticle(){}

SphericDiscreteContactFeaturesParticle& operator=(const SphericDiscreteContactFeaturesParticle& rOther);

/// Turn back information as a string.
std::string Info() const override
{
std::stringstream buffer;
buffer << "SphericDiscreteContactFeaturesParticle" ;
return buffer.str();
}

/// Print information about this object.
void PrintInfo(std::ostream& rOStream) const override {rOStream << "SphericDiscreteContactFeaturesParticle";}

/// Print object's data.
void PrintData(std::ostream& rOStream) const override {}

std::vector<double> mNeighbourContactRadius;
std::vector<double> mNeighbourRigidContactRadius;
std::vector<double> mNeighbourIndentation;
std::vector<double> mNeighbourRigidIndentation;
std::vector<double> mNeighbourContactStress;
std::vector<double> mNeighbourRigidContactStress;

private:

friend class Serializer;

void ComputeNewNeighboursHistoricalData(DenseVector<int>& temp_neighbours_ids,
                                        std::vector<array_1d<double, 3> >& temp_neighbour_elastic_contact_forces) override;

void ComputeNewRigidFaceNeighboursHistoricalData() override;

void save(Serializer& rSerializer) const override
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle);
}

void load(Serializer& rSerializer) override
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle);
}

}; // Class SphericDiscreteContactFeaturesParticle

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
            SphericDiscreteContactFeaturesParticle& rThis){ return rIStream;}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
            const SphericDiscreteContactFeaturesParticle& rThis)
{
rThis.PrintInfo(rOStream);
rOStream << std::endl;
rThis.PrintData(rOStream);

return rOStream;
}

}  // namespace Kratos.

#endif // KRATOS_SPHERIC_DISCRETE_CONTACT_FEATURES_PARTICLE_H_INCLUDED  defined
