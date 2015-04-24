// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "DEM_application.h"
#include "DEM_ExponentialHC_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    
    
    
    
    void DEM_ExponentialHC::Initialize(const ProcessInfo& rCurrentProcessInfo) {
        
    mHistoryMaxInd              = 0.0; //maximum indentation achieved
    mHistoryMaxForce            = 0.0; //maximum force achieved
    mHistoryDamage              = 0.0; //cumulated_damage
    mHistoryDegradation         = 1.0; //degradation factor for G reducing in Dempack;
    mHistoryDisp                = 0.0; //displacement;
    mHistoryShearFlag           = 0.0; //superado el limite de cortante;       
    }

    DEMContinuumConstitutiveLaw::Pointer DEM_ExponentialHC::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_ExponentialHC(*this));
        return p_clone;
    }

    void DEM_ExponentialHC::SetConstitutiveLawInProperties(Properties::Pointer pProp) const {
        std::cout << " Assigning DEM_ExponentialHC to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }
    
    void DEM_ExponentialHC::CalculateContactArea(double mRadius, double other_radius, double &calculation_area){
        double rmin = mRadius;
        if (other_radius < mRadius) rmin = other_radius;    
        calculation_area = KRATOS_M_PI * rmin*rmin;
    }
    
    void DEM_ExponentialHC::CalculateElasticConstants(double &kn_el, 
                                                double &kt_el, 
                                                double initial_dist, 
                                                double equiv_young, 
                                                double equiv_poisson, 
                                                double calculation_area ){
        double equiv_shear = equiv_young / (2.0 * (1 + equiv_poisson));
        kn_el = equiv_young * calculation_area/initial_dist;
        kt_el = equiv_shear * calculation_area/initial_dist;
    }
    
    void DEM_ExponentialHC::CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
                                                double &equiv_visco_damp_coeff_tangential,
                                                SphericContinuumParticle* element1,
                                                SphericContinuumParticle* element2,
                                                double kn_el) {

        double aux_norm_to_tang = 0.0;
        const double mRealMass = element1->GetMass();
        const double &other_real_mass = element2->GetMass();
        const double mDempack_local_damping = element1->GetProperties()[DEMPACK_LOCAL_DAMPING]; 

        equiv_visco_damp_coeff_normal = mDempack_local_damping * 2.0 * sqrt(kn_el / (mRealMass + other_real_mass)) * (sqrt(mRealMass * other_real_mass)); // := 2d0* sqrt ( kn_el*(m1*m2)/(m1+m2) )
        equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang; // Dempack no ho fa servir...
    }
    
    
    void DEM_ExponentialHC::CalculateNormalForces(double LocalElasticContactForce[3],
                                            const double kn_el,
                                            double equiv_young,
                                            double indentation,
                                            double calculation_area,
                                            double& acumulated_damage,
                                            SphericContinuumParticle* element1,
                                            SphericContinuumParticle* element2,
                                            int &mNeighbourFailureId_count,
                                            int &mIniNeighbourFailureId_mapping,
                                            int time_steps) {
                       
        const double mN1 = element1->GetProperties()[SLOPE_FRACTION_N1];
        const double mN2 = element1->GetProperties()[SLOPE_FRACTION_N2];
        const double mN3 = element1->GetProperties()[SLOPE_FRACTION_N3];
        const double mC1 = element1->GetProperties()[SLOPE_LIMIT_COEFF_C1]*1e6;
        const double mC2 = element1->GetProperties()[SLOPE_LIMIT_COEFF_C2]*1e6;
        const double mC3 = element1->GetProperties()[SLOPE_LIMIT_COEFF_C3]*1e6;
        const double mYoungPlastic = element1->GetProperties()[YOUNG_MODULUS_PLASTIC];
        const double mPlasticityLimit = element1->GetProperties()[PLASTIC_YIELD_STRESS]*1e6;
        const double mDamageMaxDisplacementFactor = element1->GetProperties()[DAMAGE_FACTOR];
        const double mTensionLimit = element1->GetProperties()[CONTACT_SIGMA_MIN]*1e6; //N/m2
        
        double kn_b = kn_el / mN1;
        double kn_c = kn_el / mN2;
        double kn_d = kn_el / mN3;
        double kp_el = mYoungPlastic / equiv_young * kn_el;
        double Yields_el = mPlasticityLimit * calculation_area;

        double Ncstr1_el = mC1 * calculation_area;
        double Ncstr2_el = mC2 * calculation_area;
        double Ncstr3_el = mC3 * calculation_area;
        double Ntstr_el = mTensionLimit * calculation_area;
        double u_max = mHistoryMaxInd;

        double& fn = LocalElasticContactForce[2]; //[2] means 'normal' contact force                

        if (indentation >= 0.0) { //COMPRESSION

            fn = kn_el * indentation;
            double u_ela1 = Ncstr1_el / kn_el;
            ;
            double u_ela2 = u_ela1 + (Ncstr2_el - Ncstr1_el) / (kn_b);
            double u_ela3 = u_ela2 + (Ncstr3_el - Ncstr2_el) / (kn_c);

            if ((indentation > u_max) || (time_steps <= 1)) { //maximum historical intentation OR first step  MSIMSI 0

                mHistoryMaxInd = indentation; // Guarda el threshold del màxim desplaçament

                if (indentation > u_ela3) { //4rt tram

                    fn = Ncstr3_el + (indentation - u_ela3) * kn_d;
                    mHistoryDegradation = kn_d / kn_el;
                } else if (indentation > u_ela2) {
                    //3r tram
                    fn = Ncstr2_el + (indentation - u_ela2) * kn_c;
                    mHistoryDegradation = kn_c / kn_el;
                } else {

                    if (indentation > u_ela1) { //2n tram

                        fn = Ncstr1_el + (indentation - u_ela1) * kn_b;
                        mHistoryDegradation = kn_b / kn_el;
                    }
                }
                mHistoryMaxForce = fn;                                 //actualitzar la força màxima a compressió.
            } else {                                            //Per sota del màxim.

                if (mHistoryMaxForce > 0.0) {                          //Màxim en compressió. 

                    double u_plas; //MSIMSI 2 akesta operació de saber quant val la u_plastica es fa cada pas de temps i en realitat es fixe sempre.
                    if (Yields_el <= Ncstr1_el) { //si el punt de plastificació està en la primera rama elastica.

                        u_plas = Yields_el / kn_el;
                    } else {

                        if (Yields_el <= Ncstr2_el) //si està en la segona...
                        {
                            u_plas = u_ela1 + (Yields_el - Ncstr1_el) / (kn_b);
                        } else if (Yields_el <= Ncstr3_el) { //si està en la tercera...

                            u_plas = u_ela2 + (Yields_el - Ncstr2_el) / (kn_c);
                        } else { //en la quarta                   
                            u_plas = u_ela3 + (Yields_el - Ncstr3_el) / (kn_d);
                        }
                    }
                    if (u_plas < u_max) { //si nosaltres estem per sota del maxim pero ja estem plastificant 

                        fn = mHistoryMaxForce - kp_el * (u_max - indentation); // Esta en zona de descarga plastica (pot estar en carga/descarga)
                        mHistoryDegradation = kp_el / kn_el;
                    } else { // Esta en zona descarga elastica, ens despreocupem de la plasticitat                 
                        if (indentation > u_ela3) { //en la 4a ramma                   
                            fn = Ncstr3_el + (indentation - u_ela3) * kn_d;
                        } else if (indentation > u_ela2) { //en la 3a ramma                    
                            fn = Ncstr2_el + (indentation - u_ela2) * kn_c;
                        } else {
                            if (indentation > u_ela1) { //en la 2a rama                     
                                fn = Ncstr1_el + (indentation - u_ela1) * kn_b;
                            }
                        }
                    }
                } //si tenim precàrrega en compressió.              
            } //Per sota del màxim.
        }//Compression
        else { //tension      
            fn = kn_el * indentation;

            double u1 = Ntstr_el / kn_el;
            double u2 = u1 * (1 + mDamageMaxDisplacementFactor);

            if (fabs(indentation) > u2) { // FULL DAMAGE 

                mNeighbourFailureId_count = 4; //tension failure
                mIniNeighbourFailureId_mapping = 4;
                acumulated_damage = 1.0;
                fn = 0.0;
            } else {
                if (fabs(indentation) > u1) {
                    double u_frac = (fabs(indentation) - u1) / (u2 - u1);
                    //failure_criterion_state = fabs(indentation)/u2;
                    acumulated_damage = u_frac;

                    if (u_frac > mHistoryDamage) {
                        mHistoryDamage = u_frac;
                    }
                    double kn_damage = u1 / (fabs(indentation)) * kn_el * (1.0 - mHistoryDamage); // normal adhesive force (gap +)
                    fn = kn_damage * indentation;
                    //fn = indentation * kn_el*(1.0 -  mHistory[mapping_new_cont][2]);  // normal adhesive force (gap +)
                }

            }
        } //Tension    
    }

    void DEM_ExponentialHC::CalculateTangentialForces(double LocalElasticContactForce[3],
                                                double LocalDeltDisp[3],
                                                double kt_el,
                                                double indentation,
                                                double calculation_area,
                                                double& failure_criterion_state,
                                                SphericContinuumParticle* element1,
                                                SphericContinuumParticle* element2,
                                                int &mNeighbourFailureId_count,
                                                int &mIniNeighbourFailureId_mapping,
                                                bool& sliding,
                                                int search_control,
                                                vector<int>& search_control_vector,
                                                double mapping_new_cont) {
                
        const double other_tg_of_fri_angle = element2->GetTgOfFrictionAngle();        
        const double myTgOfFrictionAngle = element1->GetTgOfFrictionAngle();
        double contact_tau = 0.0;
        double contact_sigma = 0.0;
        const double mTensionLimit = element1-> GetProperties()[CONTACT_SIGMA_MIN]*1e6; //N/m2    revisar header
        const double mTauZero = element1-> GetProperties()[CONTACT_TAU_ZERO]*1e6;
        const double mInternalFriccion = element1-> GetProperties()[CONTACT_INTERNAL_FRICC];
        const double mShearEnergyCoef = element1-> GetProperties()[SHEAR_ENERGY_COEF]; // *************

        double degradation = 1.0; //Tangential. With degradation:

        if (mapping_new_cont != -1) {
            if (indentation >= 0.0) { //COMPRESSION              
                degradation = mHistoryDegradation;
            } else {
                degradation = (1.0 - mHistoryDamage);
            }
        }

//        degradation = 1.0; // revisar degradation = 1.0 hardcoded   ////////////////

        if (mNeighbourFailureId_count == 0) {
            if (mHistoryShearFlag == 0.0) {

                LocalElasticContactForce[0] += -degradation * kt_el * LocalDeltDisp[0]; // 0: first tangential
                LocalElasticContactForce[1] += -degradation * kt_el * LocalDeltDisp[1]; // 1: second tangential  
            }
            double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
                    + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

            contact_tau = ShearForceNow / calculation_area;
            contact_sigma = LocalElasticContactForce[2] / calculation_area;

            double tau_strength = mTauZero;

            if (contact_sigma >= 0) {
                tau_strength = mTauZero + mInternalFriccion * contact_sigma;
            }

            if (contact_tau > tau_strength) {
                mHistoryShearFlag = 1.0;
            }

            if (mHistoryShearFlag != 0.0) {

                double increment_disp = sqrt(LocalDeltDisp[0] * LocalDeltDisp[0] + LocalDeltDisp[1] * LocalDeltDisp[1]);
                mHistoryDisp += increment_disp;

                double u1_tau = tau_strength * calculation_area / kt_el;
                double damage_tau = 1.0;

                if (mShearEnergyCoef != 0.0) {
                    damage_tau = mHistoryDisp / (u1_tau * (mShearEnergyCoef));
                }

                double aux = (1.0 - damage_tau) * (tau_strength / contact_tau);

                LocalElasticContactForce[0] = aux * LocalElasticContactForce[0];
                LocalElasticContactForce[1] = aux * LocalElasticContactForce[1];

                failure_criterion_state = 1.0;

                //failure_criterion_state = (1+mShearEnergyCoef*damage_tau) /(1+mShearEnergyCoef);		    
                //if(contact_sigma<0){ failure_criterion_state = GeometryFunctions::max( failure_criterion_state, -contact_sigma / mTensionLimit ); }		    

                if (damage_tau >= 1.0) {
                    mNeighbourFailureId_count = 2; // shear
                    mIniNeighbourFailureId_mapping = 2;
                    //failure_criterion_state = 1.0;
                    sliding = true;
                }
            } else {

                if (contact_sigma < 0) {
                    failure_criterion_state = GeometryFunctions::max(contact_tau / tau_strength, -contact_sigma / mTensionLimit);
                } else failure_criterion_state = contact_tau / tau_strength;
                if (failure_criterion_state > 1.0) failure_criterion_state = 1.0;
            }

            /*mContinuumConstitutiveLawArray[mapping_new_cont]-> EvaluateFailureCriteria(contact_sigma, contact_tau, failure_criterion_state, sliding,
                                          mFailureCriterionOption, mTauZero, mInternalFriccion, mSinContactInternalFriccion, mCosContactInternalFriccion,
                                          mNeighbourFailureId_count, mIniNeighbourFailureId_mapping, mTensionLimit);      
             */
        }
        if (search_control == 0) {
            if (mNeighbourFailureId_count != 0) {
                search_control_vector[OpenMPUtils::ThisThread()] = 1;
            }
        }

        /* Tangential Friction for broken bonds */ //dempack and kdem do the same.

        if (mNeighbourFailureId_count != 0) //*   //degut als canvis de DEMPACK hi ha hagut una modificació, ara despres de trencar es fa akest maping de maxima tangencial que és correcte!
        {

            LocalElasticContactForce[0] += -degradation * kt_el * LocalDeltDisp[0]; // 0: first tangential
            LocalElasticContactForce[1] += -degradation * kt_el * LocalDeltDisp[1]; // 1: second tangential  

            double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
                    + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

            double Frictional_ShearForceMax = (0.5 * (myTgOfFrictionAngle + other_tg_of_fri_angle)) * LocalElasticContactForce[2];
            if (Frictional_ShearForceMax < 0.0) {
                Frictional_ShearForceMax = 0.0;
            }

            failure_criterion_state = 1.0;
            if ((ShearForceNow > Frictional_ShearForceMax) && (ShearForceNow != 0.0)) {
                LocalElasticContactForce[0] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[0];
                LocalElasticContactForce[1] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[1];
                sliding = true;
            }
        }
    }        

} /* namespace Kratos.*/


//    void SphericContinuumParticle::NonlinearNormalForceCalculation(double LocalElasticContactForce[3], double kn1, double kn2, double distance, double max_dist, double initial_dist) {
//
//        LocalElasticContactForce[2] = kn1 * (initial_dist - max_dist) + kn2 * (max_dist - distance);
//    }


//
//        mGamma1                         = rCurrentProcessInfo[DONZE_G1];
//        mGamma2                         = rCurrentProcessInfo[DONZE_G2];
//        mGamma3                         = rCurrentProcessInfo[DONZE_G3];
//        mMaxDef                         = rCurrentProcessInfo[DONZE_MAX_DEF];




//sabemos que double indentation = initial_dist - distance; //#1
           
//double current_def = (initial_dist - distance) / initial_dist;
//double kn1 = kn_el;
//double kn2 = kn1 * mGamma3 + kn1 * mGamma1 * exp(mGamma2 * (current_def - mMaxDef));
//if (kn2 > kn1) {kn2 = kn1;}
//double max_dist = initial_dist * (1 - mMaxDef);
//
//initial_dist ???
//double kn_plas = kn1; // modificable en el futuro con un input

double &fn = LocalElasticContactForce[2] q inicialmente entra como 0.0

if indentation >= 0 {
    fn = kn1 * (initial_dist - distance); // = kn1 * indentation; fuerza en parte lineal
    
    if (distance < mHistDist) //// indentation >= u_max   vemos si esta en carga, comparando la distancia actual entre particulas con la maxima historica
    {
        mHistDist = distance;    
        
        if (distance < max_dist) //// indentation >= C1*Area/kn_el   se supera el limite para el cambio de pendiente.
        {
           fn = kn1 * (initial_dist - max_dist) + kn2 * (max_dist - distance);
        }
    }
}
double mHistoryMaxForce = fn; //guardamos el maximo historico de fuerza fn

else // esta en descarga, la distancia entre particulas aumenta respecto la historica, current_dist > mHistDist
{
    if (hist_fn > 0); //  fuerza normal esta en el rango de compresion de la curva
    {
        plast_dist = max_dist; // initial_dist*(1-plast_def) distancia associada al valor de plast_def impuesto, por ahora coincide con el cambio de pendiente
        if (plast_dist > mHistDist) // mientras se este por encima de la maxima historica, estamos en plasticidad.
        {
            fn = hist_fn + kn_el_plas * (mHistDist - distance); // en descarga: 500 - kn_plas(10 - 12). en carga 500 - kn_el(10 - 11) pero con distance > mHistDistance
        } 
        else // esta en descarga pero no en la zona plastica, descarga por la linea elastica
        {
            if (distance < max_dist) // se supera el limite para le cambio de pendiente, mientras plast_dist=max_dist nunca pasara, nunca descargara por kn2
            {
                NonlinearNormalForceCalculation(LocalElasticContactForce, kn1, kn2, distance, max_dist, initial_dist);
            } 
            else // descarga por la primera rama elastica
            {
                LocalElasticContactForce[2] = kn1 * (initial_dist - distance); // fuerza en parte lineal
            }
        }
    }
}

            

