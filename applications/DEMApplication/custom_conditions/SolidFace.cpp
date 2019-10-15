//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

// Project includes
#include "custom_conditions/SolidFace.h"
#include "custom_elements/spheric_particle.h"
#include "custom_conditions/dem_wall.h"
#include "custom_utilities/GeometryFunctions.h"

namespace Kratos {

    using namespace GeometryFunctions;

//***********************************************************************************
//***********************************************************************************

// Constructor

SolidFace3D::SolidFace3D() {}

// Constructor

SolidFace3D::SolidFace3D(IndexType NewId, GeometryType::Pointer pGeometry) : DEMWall(NewId, pGeometry) {}

// Constructor

SolidFace3D::SolidFace3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
           : DEMWall(NewId, pGeometry, pProperties) {
    //setting up the nodal degrees of freedom
    }

//***********************************************************************************
//***********************************************************************************

Condition::Pointer SolidFace3D::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new SolidFace3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//***********************************************************************************
//***********************************************************************************
// Destructor

SolidFace3D::~SolidFace3D() {}


//***********************************************************************************
//***********************************************************************************

void SolidFace3D::Initialize(const ProcessInfo& rCurrentProcessInfo) {
    if (! rCurrentProcessInfo[IS_RESTARTED]){
        const unsigned int number_of_nodes = GetGeometry().size();

        for (unsigned int i=0; i< number_of_nodes; i++)
        {
            this->GetGeometry()[i].FastGetSolutionStepValue(NON_DIMENSIONAL_VOLUME_WEAR) = 0.0;
            this->GetGeometry()[i].FastGetSolutionStepValue(IMPACT_WEAR) = 0.0;
        }
    }
}

//***********************************************************************************
//***********************************************************************************

void SolidFace3D::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                         ProcessInfo& r_process_info)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * 3;

    if (rRightHandSideVector.size() != MatSize)
    {
        rRightHandSideVector.resize(MatSize, false);
    }
    rRightHandSideVector = ZeroVector(MatSize);

    std::vector<SphericParticle*>& rNeighbours = this->mNeighbourSphericParticles;

    for (unsigned int i=0; i<rNeighbours.size(); i++)
    {
        if(rNeighbours[i]->Is(BLOCKED)) continue; //Inlet Generator Spheres are ignored when integrating forces.

        std::vector<DEMWall*>& rRFnei = rNeighbours[i]->mNeighbourRigidFaces;

        for (unsigned int i_nei = 0; i_nei < rRFnei.size(); i_nei++)
        {
            int Contact_Type = rNeighbours[i]->mContactConditionContactTypes[i_nei];

            if ( ( rRFnei[i_nei]->Id() == this->Id() ) && (Contact_Type > 0 ) )
            {

                const array_1d<double, 4>& weights_vector = rNeighbours[i]->mContactConditionWeights[i_nei];
                double weight = 0.0;

                double ContactForce[3] = {0.0};

                const array_1d<double, 3>& neighbour_rigid_faces_contact_force = rNeighbours[i]->mNeighbourRigidFacesTotalContactForce[i_nei];

                ContactForce[0] = neighbour_rigid_faces_contact_force[0];
                ContactForce[1] = neighbour_rigid_faces_contact_force[1];
                ContactForce[2] = neighbour_rigid_faces_contact_force[2];

                for (unsigned int k=0; k< number_of_nodes; k++)
                {
                    weight = weights_vector[k];

                    unsigned int w =  k * 3;

                    rRightHandSideVector[w + 0] += -ContactForce[0] * weight;
                    rRightHandSideVector[w + 1] += -ContactForce[1] * weight;
                    rRightHandSideVector[w + 2] += -ContactForce[2] * weight;
                }

            }//if the condition neighbour of my sphere neighbour is myself.
        }//Loop spheres neighbours (condition)
    }//Loop condition neighbours (spheres)
}//CalculateRightHandSide

void SolidFace3D::GetDeltaDisplacement( array_1d<double, 3> & delta_displacement, int inode){

  array_1d<double, 3> displacement_now = this->GetGeometry()[inode].FastGetSolutionStepValue(DISPLACEMENT);
  array_1d<double, 3> displacement_old = this->GetGeometry()[inode].FastGetSolutionStepValue(DISPLACEMENT,1);
  delta_displacement = displacement_now - displacement_old;

}

void SolidFace3D::CalculateNormal(array_1d<double, 3>& rnormal){

    array_1d<double, 3> v1, v2;

    v1[0] = GetGeometry()[1].X() - GetGeometry()[0].X();
    v1[1] = GetGeometry()[1].Y() - GetGeometry()[0].Y();
    v1[2] = GetGeometry()[1].Z() - GetGeometry()[0].Z();

    v2[0] = GetGeometry()[2].X() - GetGeometry()[0].X();
    v2[1] = GetGeometry()[2].Y() - GetGeometry()[0].Y();
    v2[2] = GetGeometry()[2].Z() - GetGeometry()[0].Z();

    MathUtils<double>::CrossProduct(rnormal, v1, v2);

    rnormal /= MathUtils<double>::Norm3(rnormal);
}


void SolidFace3D::FinalizeSolutionStep(ProcessInfo& r_process_info)
{


}

//***********************************************************************************
//***********************************************************************************

} // Namespace Kratos.
