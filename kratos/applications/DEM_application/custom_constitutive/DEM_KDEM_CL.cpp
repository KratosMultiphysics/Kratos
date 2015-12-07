// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "DEM_application.h"
#include "DEM_KDEM_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    void DEM_KDEM::Initialize(const ProcessInfo& rCurrentProcessInfo) {

        KRATOS_TRY
        mHistoryMaxInd = 0.0; //maximum indentation achieved
        mHistoryMaxForce = 0.0; //maximum force achieved
        mHistoryDamage = 0.0; //cumulated_damage
        mHistoryDegradation = 1.0; //degradation factor for G reducing in Dempack;
        mHistoryDisp = 0.0; //displacement;
        mHistoryShearFlag = 0.0; //superado el limite de cortante;  
        KRATOS_CATCH("")  
    }

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM(*this));
        return p_clone;
    }

    void DEM_KDEM::SetConstitutiveLawInProperties(Properties::Pointer pProp) const {
        std::cout << "Assigning DEM_KDEM to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    void DEM_KDEM::CalculateContactArea(double radius, double other_radius, double &calculation_area) {

        KRATOS_TRY
        double radius_sum = radius + other_radius;
        double equiv_radius = 2.0 * radius * other_radius / radius_sum;
        calculation_area = 0.25 * KRATOS_M_PI * equiv_radius * equiv_radius;
        KRATOS_CATCH("")  
    }

    void DEM_KDEM::CalculateElasticConstants(double& kn_el, double& kt_el, double initial_dist, double equiv_young,
                                             double equiv_poisson, double calculation_area) {
        
        KRATOS_TRY
        double equiv_shear = equiv_young / (2.0 * (1 + equiv_poisson));
        kn_el = equiv_young * calculation_area / initial_dist;
        kt_el = equiv_shear * calculation_area / initial_dist;
        KRATOS_CATCH("")  
    }

    void DEM_KDEM::CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
            double &equiv_visco_damp_coeff_tangential,
            SphericContinuumParticle* element1,
            SphericContinuumParticle* element2,
            double kn_el,
            double kt_el) {

        KRATOS_TRY 
        const double my_mass    = element1->GetMass();
        const double other_mass = element2->GetMass();
        const double equiv_mass = 1.0 / (1.0/my_mass + 1.0/other_mass);        
        
        const double my_coeff_rest    = element1->GetProperties()[COEFFICIENT_OF_RESTITUTION];
        const double other_coeff_rest = element2->GetProperties()[COEFFICIENT_OF_RESTITUTION];
        const double equiv_coeff_rest = 0.5 * (my_coeff_rest + other_coeff_rest);
        
        double epsilon = log(equiv_coeff_rest) / sqrt(KRATOS_M_PI * KRATOS_M_PI + log(equiv_coeff_rest) * log(equiv_coeff_rest));
              
        equiv_visco_damp_coeff_normal     = -2.0 * epsilon * sqrt(equiv_mass * kn_el);
        equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal;

    KRATOS_CATCH("")      
    }

    void DEM_KDEM::CalculateForces(ProcessInfo& rCurrentProcessInfo,
                                   double LocalElasticContactForce[3],
            double LocalDeltDisp[3],
            const double kn_el,
            double kt_el,
            double& failure_criterion_state,
            double equiv_young,
            double indentation,
            double calculation_area,
            double& acumulated_damage,
            SphericContinuumParticle* element1,
            SphericContinuumParticle* element2,
            int i_neighbour_count,
            int time_steps,
            bool& sliding,
            int search_control,
            vector<int>& search_control_vector) {
        
        KRATOS_TRY
        CalculateNormalForces(LocalElasticContactForce,
                kn_el,
                equiv_young,
                indentation,
                calculation_area,
                acumulated_damage,
                element1,
                element2,
                i_neighbour_count,
                time_steps);

        CalculateTangentialForces(LocalElasticContactForce,
                LocalDeltDisp,
                kt_el,
                indentation,
                calculation_area,
                failure_criterion_state,
                element1,
                element2,
                i_neighbour_count,
                sliding,
                search_control,
                search_control_vector);

    KRATOS_CATCH("")      
    }

    void DEM_KDEM::CalculateNormalForces(double LocalElasticContactForce[3],
            const double kn_el,
            double equiv_young,
            double indentation,
            double calculation_area,
            double& acumulated_damage,
            SphericContinuumParticle* element1,
            SphericContinuumParticle* element2,
            int i_neighbour_count,
            int time_steps) {

        KRATOS_TRY

        int &mapping_new_ini = element1->mMapping_New_Ini[i_neighbour_count];
        int &mNeighbourFailureId_count = element1->mNeighbourFailureId[i_neighbour_count];
        int &mIniNeighbourFailureId_mapping = element1->mIniNeighbourFailureId[mapping_new_ini];

        Properties& element1_props = element1->GetProperties();
        Properties& element2_props = element2->GetProperties();
        double mN1;
        double mN2;
        double mN3;
        double mC1;
        double mC2;
        double mC3;
        double mYoungPlastic;
        double mPlasticityLimit;
        double mDamageMaxDisplacementFactor;
        double mTensionLimit;

        if (&element1_props == &element2_props) {

             mN1 = element1_props[SLOPE_FRACTION_N1];
             mN2 = element1_props[SLOPE_FRACTION_N2];
             mN3 = element1_props[SLOPE_FRACTION_N3];
             mC1 = element1_props[SLOPE_LIMIT_COEFF_C1]*1e6;
             mC2 = element1_props[SLOPE_LIMIT_COEFF_C2]*1e6;
             mC3 = element1_props[SLOPE_LIMIT_COEFF_C3]*1e6;
             mYoungPlastic = element1_props[YOUNG_MODULUS_PLASTIC];
             mPlasticityLimit = element1_props[PLASTIC_YIELD_STRESS]*1e6;
             mDamageMaxDisplacementFactor = element1_props[DAMAGE_FACTOR];
             mTensionLimit = element1_props[CONTACT_SIGMA_MIN]*1e6; //N/m2
        }

        else {

            mN1 = 0.5*(element1_props[SLOPE_FRACTION_N1] + element2_props[SLOPE_FRACTION_N1] );
            mN2 = 0.5*(element1_props[SLOPE_FRACTION_N2] + element2_props[SLOPE_FRACTION_N2] );
            mN3 = 0.5*(element1_props[SLOPE_FRACTION_N3] + element2_props[SLOPE_FRACTION_N3] );
            mC1 = 0.5*1e6*(element1_props[SLOPE_LIMIT_COEFF_C1] + element2_props[SLOPE_LIMIT_COEFF_C1]);
            mC2 = 0.5*1e6*(element1_props[SLOPE_LIMIT_COEFF_C2] + element2_props[SLOPE_LIMIT_COEFF_C2]);
            mC3 = 0.5*1e6*(element1_props[SLOPE_LIMIT_COEFF_C3] + element2_props[SLOPE_LIMIT_COEFF_C3]);
            mYoungPlastic = 0.5*(element1_props[YOUNG_MODULUS_PLASTIC] + element2_props[YOUNG_MODULUS_PLASTIC]);
            mPlasticityLimit = 0.5*1e6*(element1_props[PLASTIC_YIELD_STRESS] + element2_props[PLASTIC_YIELD_STRESS]);
            mDamageMaxDisplacementFactor = 0.5*(element1_props[DAMAGE_FACTOR] + element2_props[DAMAGE_FACTOR]);
            mTensionLimit = 0.5*1e6*(element1_props[CONTACT_SIGMA_MIN] + element2_props[CONTACT_SIGMA_MIN]); //N/m2
        }

        const double kn_b = kn_el / mN1;
        const double kn_c = kn_el / mN2;
        const double kn_d = kn_el / mN3;
        const double kp_el = mYoungPlastic / equiv_young * kn_el;
        const double Yields_el = mPlasticityLimit * calculation_area;

        const double Ncstr1_el = mC1 * calculation_area;
        const double Ncstr2_el = mC2 * calculation_area;
        const double Ncstr3_el = mC3 * calculation_area;
        const double Ntstr_el = mTensionLimit * calculation_area;
        double u_max = mHistoryMaxInd;

        double& fn = LocalElasticContactForce[2]; //[2] means 'normal' contact force                

        if (indentation >= 0.0) { //COMPRESSION

            fn = kn_el * indentation;
            double u_ela1 = Ncstr1_el / kn_el;
            double u_ela2 = u_ela1 + (Ncstr2_el - Ncstr1_el) / (kn_b);
            double u_ela3 = u_ela2 + (Ncstr3_el - Ncstr2_el) / (kn_c);

            if ((indentation > u_max) || (time_steps <= 1)) { //maximum historical indentation OR first step  MSIMSI 0

                mHistoryMaxInd = indentation; // Guarda el threshold del màxim desplaçament

                if (indentation > u_ela3) { //4rt tram
                    fn = Ncstr3_el + (indentation - u_ela3) * kn_d;
                    mHistoryDegradation = kn_d / kn_el;
                } else if (indentation > u_ela2) {//3r tram
                    fn = Ncstr2_el + (indentation - u_ela2) * kn_c;
                    mHistoryDegradation = kn_c / kn_el;
                } else {
                    if (indentation > u_ela1) { //2n tram
                        fn = Ncstr1_el + (indentation - u_ela1) * kn_b;
                        mHistoryDegradation = kn_b / kn_el;
                    }
                }
                mHistoryMaxForce = fn; //actualitzar la força màxima a compressió.
            } else { //Per sota del màxim.
                if (mHistoryMaxForce > 0.0) { //Màxim en compressió. 

                    double u_plas; //MSIMSI 2 akesta operació de saber quant val la u_plastica es fa cada pas de temps i en realitat es fixe sempre.
                    if (Yields_el <= Ncstr1_el) { //si el punt de plastificació està en la primera rama elastica.

                        u_plas = Yields_el / kn_el;
                    } else {
                        if (Yields_el <= Ncstr2_el) { //si està en la segona...
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
    KRATOS_CATCH("")      
    }

    void DEM_KDEM::CalculateTangentialForces(double LocalElasticContactForce[3],
            double LocalDeltDisp[3],
            double kt_el,
            double indentation,
            double calculation_area,
            double& failure_criterion_state,
            SphericContinuumParticle* element1,
            SphericContinuumParticle* element2,
            int i_neighbour_count,
            bool& sliding,
            int search_control,
            vector<int>& search_control_vector) {

        KRATOS_TRY

        int &mapping_new_ini = element1->mMapping_New_Ini[i_neighbour_count];
        int &mapping_new_cont = element1->mMapping_New_Cont[i_neighbour_count];

        int &mNeighbourFailureId_count = element1->mNeighbourFailureId[i_neighbour_count];
        int &mIniNeighbourFailureId_mapping = element1->mIniNeighbourFailureId[mapping_new_ini];

        Properties& element1_props = element1->GetProperties();
        Properties& element2_props = element2->GetProperties();

        double contact_tau = 0.0;
        double contact_sigma = 0.0;
        double mTensionLimit;
        double mTauZero;
        double mInternalFriccion;
        double mShearEnergyCoef;


        if (&element1_props == &element2_props ) {

            mTensionLimit = 1e6 * element1_props[CONTACT_SIGMA_MIN];
            mTauZero = 1e6 * element1_props[CONTACT_TAU_ZERO];
            mInternalFriccion = element1_props[CONTACT_INTERNAL_FRICC];
            mShearEnergyCoef = element1_props[SHEAR_ENERGY_COEF];

        }
        else {

            mTensionLimit = 0.5*1e6*(element1_props[CONTACT_SIGMA_MIN] + element2_props[CONTACT_SIGMA_MIN]);
            mTauZero = 0.5*1e6*(element1_props[CONTACT_TAU_ZERO] + element2_props[CONTACT_TAU_ZERO]);
            mInternalFriccion = 0.5*(element1_props[CONTACT_INTERNAL_FRICC] + element2_props[CONTACT_INTERNAL_FRICC]);
            mShearEnergyCoef = 0.5*(element1_props[SHEAR_ENERGY_COEF] + element2_props[SHEAR_ENERGY_COEF]);

        }

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
        }
        
        if (search_control == 0) {
            if (mNeighbourFailureId_count != 0) {
                search_control_vector[OpenMPUtils::ThisThread()] = 1;
            }
        }

        /* Tangential Friction for broken bonds //dempack and kdem do the same.
	
        if (mNeighbourFailureId_count != 0)   //degut als canvis de DEMPACK hi ha hagut una modificació, ara despres de trencar es fa akest maping de maxima tangencial que és correcte!
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
	*/
    KRATOS_CATCH("")      
    }

} /* namespace Kratos.*/
