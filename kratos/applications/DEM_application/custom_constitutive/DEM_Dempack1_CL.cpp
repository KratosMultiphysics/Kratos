// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "DEM_application.h"
#include "DEM_Dempack1_CL.h"

namespace Kratos
{


void DEM_Dempack1::Initialize(const ProcessInfo& rCurrentProcessInfo)
{

    mN1 = rCurrentProcessInfo[SLOPE_FRACTION_N1];
    mN2 = rCurrentProcessInfo[SLOPE_FRACTION_N2];
    mN3 = rCurrentProcessInfo[SLOPE_FRACTION_N3];
    mC1 = rCurrentProcessInfo[SLOPE_LIMIT_COEFF_C1]*1e6;
    mC2 = rCurrentProcessInfo[SLOPE_LIMIT_COEFF_C2]*1e6;
    mC3 = rCurrentProcessInfo[SLOPE_LIMIT_COEFF_C3]*1e6;
    mYoungPlastic = rCurrentProcessInfo[YOUNG_MODULUS_PLASTIC];
    mPlasticityLimit = rCurrentProcessInfo[PLASTIC_YIELD_STRESS]*1e6;
    mDamageMaxDisplacementFactor = rCurrentProcessInfo[DAMAGE_FACTOR];
    mTauZero                       = rCurrentProcessInfo[CONTACT_TAU_ZERO]*1e6;
    mTanContactInternalFriccion    = rCurrentProcessInfo[CONTACT_INTERNAL_FRICC];
    double atanInternalFriccion    = atan(mTanContactInternalFriccion);
    mSinContactInternalFriccion    = sin(atanInternalFriccion);
    mCosContactInternalFriccion    = cos(atanInternalFriccion);
    mFailureCriterionOption        = rCurrentProcessInfo[FAILURE_CRITERION_OPTION];
    mTensionLimit                  = rCurrentProcessInfo[CONTACT_SIGMA_MIN]*1e6; //N/m2

}

DEMContinuumConstitutiveLaw::Pointer DEM_Dempack1::Clone() const
{
    DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_Dempack1(*this));
    return p_clone;
}



void DEM_Dempack1::SetConstitutiveLawInProperties(Properties::Pointer pProp) const
{
    std::cout << " Assigning DEM_Dempack1 to properties " << pProp->Id() << std::endl;
    pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER,this->Clone());
}


void DEM_Dempack1::CalculateContactForces(double mRadius,
                                          double mSqrtOfRealMass,
                                          double other_radius,
                                          double otherSqrtMass,
                                          double distance,
                                          double initial_delta,
                                          int& neighbour_failure_id,
                                          ProcessInfo& rCurrentProcessInfo,
                                          PropertiesProxy *myProperties,
                                          PropertiesProxy *neighbourProperties,
                                          int mapping_new_ini,
                                          int mapping_new_cont,
                                          unsigned int i_neighbour_count,
                                          double LocalElasticContactForce[3],
                                          double ViscoDampingLocalContactForce[3],
                                          double LocalDeltDisp[3],
                                          Vector mcont_ini_neigh_area,
                                          array_1d<double, 4 > &mHistory_mapping_new_cont,
                                          double mDempack_damping,
                                          int mDampType,
                                          int mIniNeighbourFailureId_mapping_new_ini,
                                          double LocalCoordSystem[3][3],
                                          double RelVel[3])

{

    KRATOS_TRY

    const double &other_sqrt_of_mass      = otherSqrtMass;
    double radius_sum                     = mRadius + other_radius;
    double radius_sum_i                   = 1.0 / radius_sum;
    // double equiv_radius                   = 2.0 * mRadius * other_radius * radius_sum_i;
    double initial_dist                   = (radius_sum - initial_delta);
    double initial_dist_i                 = 1.0 / initial_dist;
    double indentation                    = initial_dist - distance;

    double rmin = mRadius;
    if(other_radius<mRadius) rmin = other_radius;

    // double calculation_area               = 0.25*M_PI * equiv_radius * equiv_radius;
    double equiv_mass                     = mSqrtOfRealMass * other_sqrt_of_mass;
    double myYoung                        = myProperties->GetYoung();
    double myPoisson                      = myProperties->GetPoisson();

    double myTgOfFrictionAngle            = myProperties->GetTgOfFrictionAngle();

    double equiv_young;
    double equiv_poisson;
    double equiv_visco_damp_coeff_normal;
    double equiv_visco_damp_coeff_tangential;

    double other_tg_of_fri_angle     = neighbourProperties->GetTgOfFrictionAngle();
    double equiv_tg_of_fri_ang;
    equiv_tg_of_fri_ang               = 0.5 * (myTgOfFrictionAngle + other_tg_of_fri_angle);


    bool sliding = false;

    double contact_tau = 0.0;
    double contact_sigma = 0.0;
    double failure_criterion_state = 0.0;
    double acumulated_damage = 0.0;

    // Getting neighbour properties
    double other_young               = neighbourProperties->GetYoung();
    double other_poisson             = neighbourProperties->GetPoisson();
    //double other_ln_of_restit_coeff  = neighbourProperties->GetLnOfRestitCoeff();

    equiv_young                       = 2.0 * myYoung * other_young / (myYoung + other_young);
    
    if((myPoisson + other_poisson)!= 0.0) {
        equiv_poisson                 = 2.0 * myPoisson * other_poisson / (myPoisson + other_poisson);
    } else {
        equiv_poisson = 0.0;
    }
    //equiv_ln_of_restit_coeff          = 0.5 * (myLnOfRestitCoeff + other_ln_of_restit_coeff);

    double calculation_area     = KRATOS_M_PI*rmin*rmin;
    double equiv_shear          = equiv_young/(2.0*(1+equiv_poisson));
    double kn_el                = equiv_young*calculation_area*initial_dist_i;
    double kt_el                = equiv_shear*calculation_area*initial_dist_i;
    double aux_norm_to_tang     = sqrt(kt_el / kn_el);


    // Translational Forces

    if  (indentation > 0.0 || (neighbour_failure_id == 0) )
    {

        PlasticityAndDamage(LocalElasticContactForce,
                            kn_el,
                            equiv_young,
                            indentation,
                            calculation_area,
                            radius_sum_i,
                            failure_criterion_state,
                            acumulated_damage,
                            neighbour_failure_id,
                            mapping_new_cont,
                            mapping_new_ini,
                            rCurrentProcessInfo [TIME_STEPS],
                            mHistory_mapping_new_cont,
                            mIniNeighbourFailureId_mapping_new_ini);


    } // if compression cohesive contact


    //Tangential. With degradation:

    double degradation = 1.0;

    if(indentation >= 0.0 ) //COMPRESSION
    {

        degradation = mHistory_mapping_new_cont[3];

    }
    else
    {

        degradation = (1.0 -  mHistory_mapping_new_cont[2]);

    }


    LocalElasticContactForce[0] += - degradation*kt_el * LocalDeltDisp[0];  // 0: first tangential
    LocalElasticContactForce[1] += - degradation*kt_el * LocalDeltDisp[1];  // 1: second tangential


    double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
                                +   LocalElasticContactForce[1] * LocalElasticContactForce[1]);

    //evaluating Failure for the continuum contacts/

    if(neighbour_failure_id == 0)
    {

        //                  mNeighbourFailureId[i_neighbour_count] = 2;  //shear in compression
        //                  mNeighbourFailureId[i_neighbour_count] = 3;  //shear failure tension
        //                  mNeighbourFailureId[i_neighbour_count] = 4;  //tension failure
        //                  mNeighbourFailureId[i_neighbour_count] = 12; //both shear and tension

        EvaluateFailureCriteria(LocalElasticContactForce,
                                ShearForceNow,
                                calculation_area,
                                i_neighbour_count,
                                contact_sigma,
                                contact_tau,
                                failure_criterion_state,
                                sliding,
                                mapping_new_ini,
                                mFailureCriterionOption,
                                mTauZero,
                                mTanContactInternalFriccion,
                                mSinContactInternalFriccion,
                                mCosContactInternalFriccion,
                                mIniNeighbourFailureId_mapping_new_ini,
                                neighbour_failure_id);
    }


    // Tangential Friction for broken bonds // dempack and kdem do the same.

    if ( neighbour_failure_id != 0 )
    {
        double Frictional_ShearForceMax = equiv_tg_of_fri_ang * LocalElasticContactForce[2];

        if (Frictional_ShearForceMax < 0.0)
        {
            Frictional_ShearForceMax = 0.0;

        }

        failure_criterion_state = 1.0;

        if( (ShearForceNow >  Frictional_ShearForceMax) && (ShearForceNow != 0.0) )
        {
            LocalElasticContactForce[0] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[0];
            LocalElasticContactForce[1] = (Frictional_ShearForceMax / ShearForceNow )* LocalElasticContactForce[1];
            sliding = true;

        }

    }


    //DAMPING:
    
    //double ViscoDampingLocalContactForce[3]    = {0.0};

    if  (indentation > 0.0 || (neighbour_failure_id == 0) )
    {

        double LocalRelVel[3] = {0.0};

        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, RelVel, LocalRelVel);

        equiv_visco_damp_coeff_normal = mDempack_damping*2.0*sqrt(kn_el/(mSqrtOfRealMass * mSqrtOfRealMass + other_sqrt_of_mass * other_sqrt_of_mass))*equiv_mass;   // := 2d0* sqrt ( kn_el*(m1*m2)/(m1+m2) )

        equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang;

        CalculateViscoDamping(LocalRelVel,
                              ViscoDampingLocalContactForce,
                              indentation,
                              equiv_visco_damp_coeff_normal,
                              equiv_visco_damp_coeff_tangential,
                              sliding,
                              mDampType);
    }


    KRATOS_CATCH("")

}


void DEM_Dempack1::PlasticityAndDamage(double LocalElasticContactForce[3],
                                       double kn_el,
                                       double equiv_young,
                                       double indentation,
                                       double calculation_area,
                                       double radius_sum_i,
                                       double& failure_criterion_state,
                                       double& acumulated_damage,
                                       int& neighbour_failure_id,
                                       int mapping_new_cont,
                                       int mapping_new_ini,
                                       int time_steps,
                                       array_1d<double, 4 > &mHistory_mapping_new_cont,
                                       int mIniNeighbourFailureId_mapping_new_ini){
    //VARIABLES

    // a guardar:

    double kn_b = kn_el / mN1;
    double kn_c = kn_el / mN2;
    double kn_d = kn_el / mN3;
    double kp_el = mYoungPlastic/equiv_young * kn_el;
    double Yields_el = mPlasticityLimit * calculation_area;

    double Ncstr1_el = mC1 * calculation_area;
    double Ncstr2_el = mC2 * calculation_area;
    double Ncstr3_el = mC3 * calculation_area;
    double Ntstr_el  = mTensionLimit * calculation_area;
    double u_max = mHistory_mapping_new_cont[0];

    double& fn = LocalElasticContactForce[2]; //[2] means 'normal' contact force


    if( indentation >= 0.0 ) //COMPRESSION
    {

        fn = kn_el * indentation;

        double u_ela1 = Ncstr1_el/kn_el;
        double u_ela2 = u_ela1 + (Ncstr2_el-Ncstr1_el)/(kn_b);
        double u_ela3 = u_ela2 + (Ncstr3_el-Ncstr2_el)/(kn_c);

        if ( ( indentation > u_max ) || ( time_steps <= 1) )  // maximum historical intentation OR first step  MSIMSI 0

        {

            mHistory_mapping_new_cont[0]  = indentation;    // Guarda el treshold del màxim desplaçament


            if (indentation > u_ela3) //4rt tram
            {

                fn = Ncstr3_el + ( indentation - u_ela3 )*kn_d;
                mHistory_mapping_new_cont[3] = kn_d/kn_el;

            }
            else if (indentation > u_ela2) //3r tram
            {

                fn = Ncstr2_el + ( indentation - u_ela2 )*kn_c;
                mHistory_mapping_new_cont[3] = kn_c/kn_el;

            }
            else
            {
                if( indentation > u_ela1) //2n tram
                {
                    fn = Ncstr1_el + (indentation - u_ela1)*kn_b;
                    mHistory_mapping_new_cont[3] = kn_b/kn_el;

                }

            }

            mHistory_mapping_new_cont[1] = fn; //actualitzar la força màxima a compressió.

        }

        else //Per sota del màxim.
        {

            if(mHistory_mapping_new_cont[1] > 0.0)  //Màxim en compressió.
            {

                double u_plas;        //MSIMSI 2 akesta operació de saber quant val la u_plastica es fa cada pas de temps i en realitat es fixe sempre.

                if(Yields_el <= Ncstr1_el) //si el punt de plastificació està en la primera rama elastica.
                {
                    u_plas = Yields_el/kn_el;

                }
                else
                {
                    if(Yields_el <= Ncstr2_el) //si està en la segona...
                    {
                        u_plas = u_ela1 + (Yields_el-Ncstr1_el)/(kn_b);

                    }
                    else if(Yields_el <= Ncstr3_el) //si està en la tercera...
                    {
                        u_plas = u_ela2 + (Yields_el-Ncstr2_el)/(kn_c);

                    }

                    else //en la quarta
                    {

                        u_plas = u_ela3 + (Yields_el-Ncstr3_el)/(kn_d);
                    }

                }


                if ( u_plas < u_max ) //si nosaltres estem per sota del maxim pero ja estem plastificant
                {
                    fn = mHistory_mapping_new_cont[1] - kp_el*(u_max - indentation); // Esta en zona de descarga plastica (pot estar en carga/descarga)
                    mHistory_mapping_new_cont[3] = kp_el/kn_el;


                }
                else                                   // Esta en zona descarga elastica, ens despreocupem de la plasticitat
                {

                    if ( indentation > u_ela3)  //en la 4a ramma
                    {
                        fn = Ncstr3_el + (indentation - u_ela3)*kn_d;

                    }

                    else if ( indentation > u_ela2)  //en la 3a ramma
                    {
                        fn = Ncstr2_el + (indentation - u_ela2)*kn_c;

                    }

                    else
                    {
                        if(indentation > u_ela1)  //en la 2a rama
                        {
                            fn = Ncstr1_el + (indentation-u_ela1)*kn_b;
                        }

                    }

                }

            } //si tenim precàrrega en compressió.

        }//Per sota del màxim.


    } //Compression


    else //tension
    {
        fn = kn_el * indentation;

        double u1 = Ntstr_el / kn_el;

        double u2 = u1*(1+ mDamageMaxDisplacementFactor);

        if(fabs(indentation) > u2)                  // FULL DAMAGE
        {
            neighbour_failure_id = 4; //tension failure
            mIniNeighbourFailureId_mapping_new_ini = 4;
            acumulated_damage = 1.0;
            fn = 0.0;
        }
        else
        {
            if (fabs(indentation) > u1)
            {
                double u_frac = (fabs(indentation) - u1)/(u2 - u1);
                //failure_criterion_state = fabs(indentation)/u2;
                acumulated_damage = u_frac;


                if (u_frac > mHistory_mapping_new_cont[2])
                {
                    mHistory_mapping_new_cont[2] = u_frac;
                }

            }

            fn = indentation * kn_el*(1.0 -  mHistory_mapping_new_cont[2]);  // normal adhesive force (gap +)

        }

    }//Tension

}


void DEM_Dempack1::EvaluateFailureCriteria(double LocalElasticContactForce[3],
                                           double ShearForceNow,
                                           double calculation_area,
                                           int i_neighbour_count,
                                           double& contact_sigma,
                                           double& contact_tau,
                                           double& failure_criterion_state,
                                           bool& sliding,
                                           int mapping_new_ini,
                                           int mFailureCriterionOption,
                                           double mTauZero,
                                           double mTanContactInternalFriccion,
                                           double mSinContactInternalFriccion,
                                           double mCosContactInternalFriccion,
                                           int mIniNeighbourFailureId_mapping_new_ini,
                                           int& neighbour_failure_id){

    //(1) MOHR-COULOMB FAILURE: (we don't consider rotational spring!!!!! here) need to be thought.

    if (mFailureCriterionOption==1)  //MOHR-COULOMB
    {
        contact_tau = ShearForceNow/(calculation_area);
        contact_sigma = LocalElasticContactForce[2]/(calculation_area);

        double sigma_max, sigma_min;

        if (LocalElasticContactForce[2]>=0)
        {
            sigma_max = contact_sigma;
            sigma_min = 0;
        }
        else
        {
            sigma_max = 0;
            sigma_min = contact_sigma;
        }

        //change into principal stresses

        double centre = 0.5*(sigma_max + sigma_min);
        double radius = sqrt( (sigma_max - centre)*(sigma_max - centre) + contact_tau*contact_tau   ) ;

        double sigma_I = centre + radius;
        double sigma_II = centre - radius;

        // Check:

        double distance_to_failure = ( mTauZero/(mTanContactInternalFriccion) + centre )*mSinContactInternalFriccion;

        failure_criterion_state = radius/distance_to_failure;


        if ( sigma_I - sigma_II >= 2*mTauZero*mCosContactInternalFriccion + (sigma_I + sigma_II)*mSinContactInternalFriccion )
        {
            //breaks

            neighbour_failure_id= 5; //mohr coulomb
            mIniNeighbourFailureId_mapping_new_ini= 5;
            failure_criterion_state = 1.0;
            sliding = true ;

        }


    } //MOHR-COULOMB

    ///(2) UNCOUPLED FRACTURE

    if (mFailureCriterionOption==2)//UNCOUPLED FRACTURE and DEMPACK
    {

        contact_tau = ShearForceNow/(calculation_area);

        contact_sigma = LocalElasticContactForce[2]/(calculation_area);

        //double mTauZero = 0.5*sqrt(mCompressionLimit*mTensionLimit);

        if (LocalElasticContactForce[2]>=0)
        {
            double tau_strength = mTauZero+mTanContactInternalFriccion*contact_sigma;

            failure_criterion_state = contact_tau/tau_strength;

            if(contact_tau>tau_strength)
            {
                neighbour_failure_id= 2; //shear in compression
                mIniNeighbourFailureId_mapping_new_ini= 2;
                failure_criterion_state = 1.0;
                sliding = true;
            }
        } //positive sigmas

        else //negative sigmas
        {

            double tau_strength = mTauZero;

            failure_criterion_state = GeometryFunctions::max(contact_tau/tau_strength, -contact_sigma/mTensionLimit) ;

            if(contact_tau > tau_strength)
            {
                neighbour_failure_id= 3;  //shear failure tension
                mIniNeighbourFailureId_mapping_new_ini= 3;
                sliding = true;
                failure_criterion_state = 1.0;

                //Amb Dempack la fractura tracció és el limit del dany i es mira al calcul de forces...


                /*

                            if(contact_sigma<-mTensionLimit && mElasticityType<2)
                            {
                                neighbour_failure_id= 12; //both shear and tension
                                mIniNeighbourFailureId[ mapping_new_ini ] = 12;
                                failure_criterion_state = 1.0;
                            } //both shear and tension
                          */


            }

            /*
                        else if (contact_sigma<-mTensionLimit && mElasticityType<2)
                        {
                            neighbour_failure_id= 4; //tension failure
                            mIniNeighbourFailureId[ mapping_new_ini ] = 4;
                            sliding = true;
                            failure_criterion_state = 1.0;

                        }*/


        } //negative values of sigma

    } //UNCOUPLED FRACTURE

}



} /* namespace Kratos.*/
