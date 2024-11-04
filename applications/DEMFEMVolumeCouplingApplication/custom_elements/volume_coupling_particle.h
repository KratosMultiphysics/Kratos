

#if !defined(KRATOS_VOLUME_COUPLING_PARTICLE_H_INCLUDED)
#define  KRATOS_VOLUME_COUPLING_PARTICLE_H_INCLUDED

// System includes
// #include <string>
// #include <iostream>

// Project includes
// #include "includes/define.h"
// #include "discrete_element.h"
// #include "custom_utilities/AuxiliaryFunctions.h"
// #include "custom_constitutive/DEM_discontinuum_constitutive_law.h"
// #include "custom_constitutive/DEM_rolling_friction_model.h"
// #include "custom_conditions/RigidFace.h"
// #include "custom_conditions/dem_wall.h"
// #include "custom_strategies/schemes/dem_integration_scheme.h"
// #include "includes/kratos_export_api.h"
// #include "custom_utilities/properties_proxies.h"
// #include "includes/kratos_flags.h"
#include "custom_elements/spheric_particle.h"


namespace Kratos
{

class KRATOS_API(DEMFEM_VOLUME_COUPLING_APPLICATION) VolumeCouplingParticle : public SphericParticle
{


public:
   
VolumeCouplingParticle(IndexType NewId, GeometryType::Pointer pGeometry);
VolumeCouplingParticle();// Default constructor needed for serialization
VolumeCouplingParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);
Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

//virtual void ComputeAdditionalForces(array_1d<double, 3>& externally_applied_force, array_1d<double, 3>& externally_applied_moment, const ProcessInfo& r_process_info, const array_1d<double,3>& gravity) override;

protected:

virtual void EvaluateBallToRigidFaceForcesForPositiveIndentations(SphericParticle::ParticleDataBuffer &data_buffer,
                                                                   const int rigid_neighbour_index,
                                                                   const double DeltVel[3],
                                                                   const ProcessInfo& r_process_info,
                                                                   double OldLocalElasticContactForce[3],
                                                                   double LocalElasticContactForce[3],
                                                                   double LocalDeltDisp[3],
                                                                   const double indentation,
                                                                   const double  previous_indentation,
                                                                   double ViscoDampingLocalContactForce[3],
                                                                   double& cohesive_force,
                                                                   Condition* const wall,
                                                                   bool& sliding) override;

double GetMass () override;  

virtual void ComputeBallToRigidFaceContactForceAndMoment(
    SphericParticle::ParticleDataBuffer & data_buffer,
    array_1d<double, 3>& r_elastic_force,
    array_1d<double, 3>& r_contact_force,
    array_1d<double, 3>& rigid_element_force,
    const ProcessInfo& r_process_info) override;

void Initialize(const ProcessInfo& r_process_info) override;


virtual void EvaluateBallToBallForcesForPositiveIndentiations(SphericParticle::ParticleDataBuffer & data_buffer,
                                                            const ProcessInfo& r_process_info,
                                                            double LocalElasticContactForce[3],
                                                            double DeltDisp[3],
                                                            double LocalDeltDisp[3],
                                                            double RelVel[3],
                                                            double indentation,
                                                            double ViscoDampingLocalContactForce[3],
                                                            double& cohesive_force,
                                                            SphericParticle* element2,
                                                            bool& sliding,
                                                            double LocalCoordSystem[3][3],
                                                            double OldLocalCoordSystem[3][3],
                                                            array_1d<double, 3>& neighbour_elastic_contact_force) override;


   
};

}  // namespace Kratos.
#endif // KRATOS_VOLUME_COUPLING_PARTICLE_H_INCLUDED  defined

