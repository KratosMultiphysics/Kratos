//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

// Project includes
#include "custom_conditions/dem_wall.h"
#include "../custom_elements/spheric_particle.h"
#include "custom_utilities/GeometryFunctions.h"

namespace Kratos
{
	using namespace GeometryFunctions;

//***********************************************************************************
//***********************************************************************************


// Constructor

DEMWall::DEMWall()
{
}

// Constructor

DEMWall::DEMWall(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
}

// Constructor

DEMWall::DEMWall(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
    //setting up the nodal degrees of freedom
}

//***********************************************************************************
//***********************************************************************************

Condition::Pointer DEMWall::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new DEMWall(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//***********************************************************************************
//***********************************************************************************
// Destructor

DEMWall::~DEMWall()
{
}

//***********************************************************************************
//***********************************************************************************

void DEMWall::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR << "This function (DEMWall::Initialize) shouldn't be accessed, use derived class instead"<<std::endl;
}

//***********************************************************************************
//***********************************************************************************

void DEMWall::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& r_process_info) {

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * dim;

    if (rRightHandSideVector.size() != MatSize) {
        rRightHandSideVector.resize(MatSize, false);
    }
    rRightHandSideVector = ZeroVector(MatSize);

    std::vector<SphericParticle*>& vector_of_glued_particles = GetVectorOfGluedParticles();
    for (unsigned int i=0; i<vector_of_glued_particles.size(); i++) {
        SphericParticle* p_particle = vector_of_glued_particles[i];
        DEMIntegrationScheme& dem_scheme = p_particle->GetTranslationalIntegrationScheme();
        GluedToWallScheme* p_glued_scheme = dynamic_cast<GluedToWallScheme*>(&dem_scheme);
        #ifdef KRATOS_DEBUG
        Condition* p_condition = p_glued_scheme->pGetCondition();
        if(p_condition != this) {
            KRATOS_ERROR << "Inconsistency in the pointers to faces of the glued spheres!! Condition with id: " <<this->Id()<<" used a sphere with Id: "<<p_particle->Id()<<" glued to it, but the sphere was actually glued to a Condition with Id: "<<p_condition->Id()<<std::endl;
        }
        #endif
        array_1d<double, 3> force = ZeroVector(3);
        std::vector<double> weights_vector(number_of_nodes, 0.0);
        noalias(force) = p_particle->GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);
        Vector& r_shape_functions_values = p_glued_scheme->GetShapeFunctionsValues();
        for(size_t j=0; j<r_shape_functions_values.size(); j++) {
            weights_vector[j] = r_shape_functions_values[j];
        }
        for (unsigned int k=0; k< number_of_nodes; k++) {
            unsigned int w =  k * dim;
            for(size_t l=0; l<dim; l++) {
                rRightHandSideVector[w + l] += force[l] * weights_vector[k];
            }
        }

        //AddForcesDueToTorque(rRightHandSideVector, r_shape_functions_values, weights_vector, force, p_particle);
    }
    std::vector<SphericParticle*>& rNeighbours = this->mNeighbourSphericParticles;

    for (unsigned int i=0; i<rNeighbours.size(); i++) {
        if(rNeighbours[i]->Is(BLOCKED)) continue; //Inlet Generator Spheres are ignored when integrating forces.
        array_1d<double, 3> force = ZeroVector(3);
        std::vector<double> weights_vector(number_of_nodes, 0.0);
        ComputeForceAndWeightsOfSphereOnThisFace(rNeighbours[i], force, weights_vector);

        for (unsigned int k=0; k< number_of_nodes; k++) {
            unsigned int w =  k * dim;
            for(size_t l=0; l<dim; l++) {
                rRightHandSideVector[w + l] += -force[l] * weights_vector[k];
            }
        }
    }
}

void DEMWall::CalculateElasticForces(VectorType& rElasticForces, const ProcessInfo& r_process_info)
{

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * dim;

    if (rElasticForces.size() != MatSize)
    {
        rElasticForces.resize(MatSize, false);
    }
    rElasticForces = ZeroVector(MatSize);

    std::vector<SphericParticle*>& rNeighbours = this->mNeighbourSphericParticles;

    for (unsigned int i=0; i<rNeighbours.size(); i++)
    {

        if(rNeighbours[i]->Is(BLOCKED)) continue; //Inlet Generator Spheres are ignored when integrating forces.

        std::vector<DEMWall*>& rRFnei = rNeighbours[i]->mNeighbourRigidFaces;

        for (unsigned int i_nei = 0; i_nei < rRFnei.size(); i_nei++)
        {
            const int& contact_type = rNeighbours[i]->mContactConditionContactTypes[i_nei];

            if ( (rRFnei[i_nei] != nullptr) && ( rRFnei[i_nei]->Id() == this->Id() ) && (contact_type > 0 ) ) {
                const array_1d<double, 4>& weights_vector = rNeighbours[i]->mContactConditionWeights[i_nei];
                double ContactElasticForce[3] = {0.0};

                const array_1d<double, 3>& neighbour_rigid_faces_elastic_contact_force = rNeighbours[i]->mNeighbourRigidFacesElasticContactForce[i_nei];
                ContactElasticForce[0] = neighbour_rigid_faces_elastic_contact_force[0];
                ContactElasticForce[1] = neighbour_rigid_faces_elastic_contact_force[1];
                ContactElasticForce[2] = neighbour_rigid_faces_elastic_contact_force[2];

                for (unsigned int k=0; k< number_of_nodes; k++) {
                    unsigned int w =  k * dim;
                    for(size_t l=0; l<dim; l++) {
                        rElasticForces[w + l] += -ContactElasticForce[l] * weights_vector[k];
                    }
                }
            }//if the condition neighbour of my sphere neighbour is myself.
        }//Loop spheres neighbours (condition)
    }//Loop condition neighbours (spheres)
}

void DEMWall::ComputeForceAndWeightsOfSphereOnThisFace(SphericParticle* p_particle, array_1d<double, 3>& force, std::vector<double>& weights_vector) {
    if(p_particle->Is(DEMFlags::STICKY)) return;

    std::vector<DEMWall*>& rRFnei = p_particle->mNeighbourRigidFaces;

    for (unsigned int i_nei = 0; i_nei < rRFnei.size(); i_nei++) {
        int Contact_Type = p_particle->mContactConditionContactTypes[i_nei];

        if ( (rRFnei[i_nei] == this) && (Contact_Type > 0 ) ) {
            for(size_t i=0; i<weights_vector.size(); i++) weights_vector[i] = p_particle->mContactConditionWeights[i_nei][i];
            const array_1d<double, 3>& neighbour_rigid_faces_contact_force = p_particle->mNeighbourRigidFacesTotalContactForce[i_nei];
            noalias(force) = neighbour_rigid_faces_contact_force;
        }//if the condition neighbour of my sphere neighbour is myself.
    }

}

void DEMWall::GetDeltaDisplacement( array_1d<double, 3> & delta_displacement, int inode)
{
  delta_displacement = this->GetGeometry()[inode].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
}

void DEMWall::InitializeSolutionStep(const ProcessInfo& r_process_info){
}


void DEMWall::CalculateNormal(array_1d<double, 3>& rnormal){

  KRATOS_ERROR << "This function (DEMWall::CalculateNormal) shouldn't be accessed, use derived class instead"<<std::endl;
}

 void DEMWall::AddExplicitContribution(const VectorType& rRHS,
                         const Variable<VectorType>& rRHSVariable,
                         const Variable<array_1d<double,3> >& rDestinationVariable,
                         const ProcessInfo& r_process_info)
{
    KRATOS_TRY


    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

   if( rDestinationVariable == EXTERNAL_FORCE )
      {

    for(unsigned int i=0; i< number_of_nodes; i++)
      {
        int index = dimension * i;

        GetGeometry()[i].SetLock();

        array_1d<double, 3 > &ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
        for(unsigned int j=0; j<dimension; j++) {
            ExternalForce[j] += rRHS[index + j];
        }

        GetGeometry()[i].UnSetLock();
      }
      }

    if( rDestinationVariable == FORCE_RESIDUAL )
      {

    for(unsigned int i=0; i< number_of_nodes; i++)
      {
        int index = dimension * i;

        GetGeometry()[i].SetLock();

        array_1d<double, 3 > &ForceResidual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
        for(unsigned int j=0; j<dimension; j++)
          {
        ForceResidual[j] += rRHS[index + j];
          }

        GetGeometry()[i].UnSetLock();
      }
      }

    KRATOS_CATCH( "" )
}

double DEMWall::GetYoung() const                    { return GetProperties()[YOUNG_MODULUS];    }
double DEMWall::GetPoisson() const                  { return GetProperties()[POISSON_RATIO];    }

void DEMWall::FinalizeSolutionStep(const ProcessInfo& r_process_info)
{

}

void DEMWall::GetRightHadSideVector(std::vector<array_1d <double, 3> >& rRightHandSideVector)
{

  for(unsigned int a = 0; a< mRightHandSideVector.size(); a++)
  {
    for(unsigned int b = 0; b< 3; b++)
    {
      rRightHandSideVector[a][b] = mRightHandSideVector[a][b];
    }
  }

}

void DEMWall::SetRightHadSideVector(const std::vector<array_1d <double, 3> >& rRightHandSideVector)
{
  for(unsigned int a = 0; a< mRightHandSideVector.size(); a++)
  {
    for(unsigned int b = 0; b< 3; b++)
    {
      mRightHandSideVector[a][b] = rRightHandSideVector[a][b];
    }
  }

}

void DEMWall::AddToRightHadSideVector(const std::vector<array_1d <double, 3> >& rRightHandSideVector)
{
  for(unsigned int a = 0; a< mRightHandSideVector.size(); a++)
  {
    for(unsigned int b = 0; b< 3; b++)
    {
      mRightHandSideVector[a][b] += rRightHandSideVector[a][b];
    }
  }

}

//***********************************************************************************
//***********************************************************************************

} // Namespace Kratos.
