// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "DEM_application.h"
#include "DEM_D_linear_viscous_Coulomb_CL.h"

namespace Kratos
{


void DEM_D_linear_viscous_Coulomb::Initialize(const ProcessInfo& rCurrentProcessInfo){}

DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_linear_viscous_Coulomb::Clone() const
{
    DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_linear_viscous_Coulomb(*this));
    return p_clone;
}

void DEM_D_linear_viscous_Coulomb::SetConstitutiveLawInProperties(Properties::Pointer pProp) const
{
    std::cout << " Assigning DEM_D_linear_viscous_Coulomb to properties " << pProp->Id() << std::endl;
    pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER,this->Clone());
}


void DEM_D_linear_viscous_Coulomb::CalculateContactForces(double LocalElasticContactForce[3],
                                                          double indentation,
                                                          double kn_el,
                                                          double LocalDeltDisp[3],
                                                          double kt_el,
                                                          int& neighbour_failure_id,
                                                          double equiv_tg_of_fri_ang)
{

    /* Translational Forces */

    if  (indentation > 0.0 || (neighbour_failure_id == 0) )
    {
        NormalForceCalculation(LocalElasticContactForce,kn_el,indentation);

    } //if compression or cohesive contact

    //    //Tangential. With degradation:

    //    double degradation = 1.0;
    //    LocalElasticContactForce[0] += - degradation*kt_el * LocalDeltDisp[0];  // 0: first tangential
    //    LocalElasticContactForce[1] += - degradation*kt_el * LocalDeltDisp[1];  // 1: second tangential

    double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
                                +   LocalElasticContactForce[1] * LocalElasticContactForce[1]);


    /* Tangential Friction for broken bonds */  //dempack and kdem do the same.

    if ( neighbour_failure_id != 0 )
        //degut als canvis de DEMPACK hi ha hagut una modificació, ara despres de trencar es fa akest maping de maxima tangencial que és correcte!
    {
        double Frictional_ShearForceMax = equiv_tg_of_fri_ang * LocalElasticContactForce[2];

        if (Frictional_ShearForceMax < 0.0)
        {
            Frictional_ShearForceMax = 0.0;

        }


        //        failure_criterion_state = 1.0;

        if( (ShearForceNow >  Frictional_ShearForceMax) && (ShearForceNow != 0.0) )
        {
            LocalElasticContactForce[0] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[0];
            LocalElasticContactForce[1] = (Frictional_ShearForceMax / ShearForceNow )* LocalElasticContactForce[1];
            //            sliding = true;

        }

    }

} //DEM_D_linear_viscous_Coulomb  CalculateContactForces






void DEM_D_linear_viscous_Coulomb::NormalForceCalculation(double LocalElasticContactForce[3],
                                                          double kn,
                                                          double indentation){
    LocalElasticContactForce[2] = kn * indentation;

} //NormalForceCalculation





} /* namespace Kratos.*/
