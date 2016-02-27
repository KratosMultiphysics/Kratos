// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "DEM_application.h"
#include "DEM_ExponentialHC_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    void DEM_ExponentialHC::Initialize(const ProcessInfo& r_process_info) {
        

        KRATOS_TRY
        mHistoryMaxInd = 0.0; //maximum indentation achieved
        mHistoryMaxForce = 0.0; //maximum force achieved
        mHistoryDamage = 0.0; //cumulated_damage
        mHistoryDegradation = 1.0; //degradation factor for G reducing in Dempack;
        mHistoryDisp = 0.0; //displacement;
        mHistoryShearFlag = 0.0; //superado el limite de cortante;   

        mGamma1 = 0.0;
        mGamma2 = 0.0;
        mGamma3 = 0.0;
        mMaxDef = 0.0;
        
        KRATOS_CATCH("")  

    }

    DEMContinuumConstitutiveLaw::Pointer DEM_ExponentialHC::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_ExponentialHC(*this));
        return p_clone;
    }

    void DEM_ExponentialHC::SetConstitutiveLawInProperties(Properties::Pointer pProp) const {
        std::cout << "Assigning DEM_ExponentialHC to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    void DEM_ExponentialHC::CalculateContactArea(double radius, double other_radius, double &calculation_area) {
        
        KRATOS_TRY 
        double rmin = radius;
        if (other_radius < radius) rmin = other_radius;
        calculation_area = KRATOS_M_PI * rmin*rmin;
        KRATOS_CATCH("")  
    }

    void DEM_ExponentialHC::CalculateElasticConstants(double &kn_el,
            double &kt_el,
            double initial_dist,
            double equiv_young,
            double equiv_poisson,
            double calculation_area) {
        
        KRATOS_TRY 
        double equiv_shear = equiv_young / (2.0 * (1 + equiv_poisson));
        kn_el = equiv_young * calculation_area / initial_dist;
        kt_el = equiv_shear * calculation_area / initial_dist;
        KRATOS_CATCH("")  
    }

    void DEM_ExponentialHC::CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
            double &equiv_visco_damp_coeff_tangential,
            SphericContinuumParticle* element1,
            SphericContinuumParticle* element2,
            double kn_el,
            double kt_el) {

        KRATOS_TRY 
        double aux_norm_to_tang = 0.0;
        const double mRealMass = element1->GetMass();
        const double &other_real_mass = element2->GetMass();
        const double mCoefficientOfRestitution = element1->GetProperties()[COEFFICIENT_OF_RESTITUTION];

        equiv_visco_damp_coeff_normal = (1-mCoefficientOfRestitution) * 2.0 * sqrt(kn_el / (mRealMass + other_real_mass)) * (sqrt(mRealMass * other_real_mass)); // := 2d0* sqrt ( kn_el*(m1*m2)/(m1+m2) )
        equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang; // Dempack no ho fa servir...
        
        KRATOS_CATCH("")  
    }

    void DEM_ExponentialHC::CalculateForces(ProcessInfo& r_process_info,
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

    void DEM_ExponentialHC::CalculateTangentialForces(double LocalElasticContactForce[3],
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


        //int &mapping_new_ini = element1->mMappingNewIni[i_neighbour_count];
        //int &mapping_new_cont = element1->mMappingNewCont[i_neighbour_count];

        int &mNeighbourFailureId_count = element1->mIniNeighbourFailureId[i_neighbour_count];
        //int &mIniNeighbourFailureId_mapping = element1->mIniNeighbourFailureId[mapping_new_ini];


        const double other_tg_of_fri_angle = element2->GetTgOfFrictionAngle();
        const double myTgOfFrictionAngle = element1->GetTgOfFrictionAngle();
        double contact_tau = 0.0;
        double contact_sigma = 0.0;
        const double mTensionLimit = element1-> GetProperties()[CONTACT_SIGMA_MIN]*1e6; //N/m2    revisar header
        const double mTauZero = element1-> GetProperties()[CONTACT_TAU_ZERO]*1e6;
        const double mInternalFriccion = element1-> GetProperties()[CONTACT_INTERNAL_FRICC];
        const double mShearEnergyCoef = element1-> GetProperties()[SHEAR_ENERGY_COEF]; // *************

        double degradation = 1.0; //Tangential. With degradation:

        //if (mapping_new_cont != -1) {
        if (i_neighbour_count < int(element1->mContinuumInitialNeighborsSize)) {
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
                    //mIniNeighbourFailureId_mapping = 2;
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
    KRATOS_CATCH("")      
    }

    void DEM_ExponentialHC::CalculateNormalForces(double LocalElasticContactForce[3],
            const double kn_el,
            double equiv_young,
            double indentation,
            double calculation_area,
            double& acumulated_damage,
            SphericContinuumParticle* element1,
            SphericContinuumParticle* element2,
            int i_neighbour_count,
            int time_steps) {

        //        mGamma1 = r_process_info[DONZE_G1];
        //        mGamma2 = r_process_info[DONZE_G2];
        //        mGamma3 = r_process_info[DONZE_G3];
        //        mMaxDef = r_process_info[DONZE_MAX_DEF];
        
        
        KRATOS_TRY

        mGamma1 = 0.2;
        mGamma2 = 16;
        mGamma3 = 0.275;
        mMaxDef = 0.002;

        //int &mapping_new_ini = element1->mMappingNewIni[i_neighbour_count];
        int &mNeighbourFailureId_count = element1->mIniNeighbourFailureId[i_neighbour_count];
        //int &mIniNeighbourFailureId_mapping = element1->mIniNeighbourFailureId[mapping_new_ini];
        double &mNeighbourDelta_count = element1->mIniNeighbourDelta[i_neighbour_count];

        const double mDamageMaxDisplacementFactor = element1->GetProperties()[DAMAGE_FACTOR];
        const double mTensionLimit = element1->GetProperties()[CONTACT_SIGMA_MIN]*1e6; //N/m2
        const double &other_radius = element2->GetRadius();
        const double my_radius = element1->GetRadius();
        const double initial_delta = mNeighbourDelta_count; //*
        const double initial_dist = (other_radius + my_radius - initial_delta);
        double current_def = indentation / initial_dist;
        double u_max = mHistoryMaxInd;
        double kn_plas = kn_el; // modificable en el futuro con un input         double kn_plas = mYoungPlastic / equiv_young * kn_el;
        double Ntstr_el = mTensionLimit * calculation_area;
        double &fn = LocalElasticContactForce[2];
        double kn_exp = kn_el * mGamma3 + kn_el * mGamma1 * exp(mGamma2 * (current_def - mMaxDef));
        if (kn_exp > kn_el) {
            kn_exp = kn_el;
        }



        if (indentation >= 0.0) { // COMPRESSION STAGE
            fn = kn_el * indentation; //  fuerza en parte lineal

            if ((indentation > u_max) || (time_steps <= 1)) { //maximum historical vemos si esta en carga, comparando la indentation actual con la maxima historica
                {
                    mHistoryMaxInd = indentation; // Guarda el threshold de la maxima indentation

                    if (indentation > initial_dist * mMaxDef) //// dempack_indentation >= C1*Area/kn_el   se supera el limite para el cambio de pendiente.
                    {
                        fn = kn_el * initial_dist * mMaxDef + kn_exp * (indentation - initial_dist * mMaxDef);
                    }
                }
                mHistoryMaxForce = fn; //actualitzar la força màxima a compressió.
            } else { ////Per sota del màxim. esta en descarga, la distancia entre particulas aumenta respecto la historica, current_dist > mHistDist          
                if (mHistoryMaxForce > 0.0) { //Màxim en compressió. 

                    double u_plas; //  por ahora coincide con el cambio de pendiente
                    if (indentation <= initial_dist * mMaxDef) { // //si el punt de plastificació està en la primera rama elastica..

                        u_plas = indentation;
                    } else {u_plas = initial_dist * mMaxDef + fn/kn_exp;
                    }
                    if (u_plas < u_max) { //si nosaltres estem per sota del maxim pero ja estem plastificant 
                        fn = mHistoryMaxForce - kn_plas * (u_max - indentation); // Esta en zona de descarga plastica (pot estar en carga/descarga)
                        mHistoryDegradation = kn_plas / kn_el; // used in tangential degradation
                    } else {
                        if (indentation > initial_dist * mMaxDef) {
                            //                            por encima de la rama elastica
                            fn = kn_el * initial_dist * mMaxDef + kn_exp * (indentation - initial_dist * mMaxDef);
                        }
                    }
                }
            }
        } else { //tension   same as dempack atm
            fn = kn_el * indentation;
            double u1 = Ntstr_el / kn_el;
            double u2 = u1 * (1 + mDamageMaxDisplacementFactor);

            if (fabs(indentation) > u2) { // FULL DAMAGE 
                mNeighbourFailureId_count = 4; //tension failure
                //mIniNeighbourFailureId_mapping = 4;
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

} /* namespace Kratos.*/









//////////////////////////////////////////////////////////////////////////////////////////
//void SphericContinuumParticle::NonlinearNormalForceCalculation(double LocalElasticContactForce[3],
//        double kn_el,
//        double distance,
//        double max_dist,
//        double initial_dist) {
//
//    //mGamma1 = r_process_info[DONZE_G1];
//    //mGamma2 = r_process_info[DONZE_G2];
//    //mGamma3 = r_process_info[DONZE_G3];
//    //mMaxDef = r_process_info[DONZE_MAX_DEF];
//
//    //sabemos que double indentation = initial_dist - distance;
//
//    double current_def = (initial_dist - distance) / initial_dist;
//    double kn1 = kn_el;
//    double kn2 = kn1 * mGamma3 + kn1 * mGamma1 * exp(mGamma2 * (current_def - mMaxDef));
//    if (kn2 > kn1) {
//        kn2 = kn1;
//    }
//    double max_dist = initial_dist * (1 - mMaxDef);
//
//    //initial_dist ???
//    double kn_plas = kn1; // modificable en el futuro con un input
//
//    double &fn = LocalElasticContactForce[2] q inicialmente entra como 0.0
//
//    if indentation >= 0 {
//        fn = kn1 * (initial_dist - distance); // = kn1 * indentation; fuerza en parte lineal
//
//        if (distance < mHistDist) //// indentation >= u_max   vemos si esta en carga, comparando la distancia actual entre particulas con la maxima historica
//        {
//            mHistDist = distance;
//
//            if (distance < max_dist) //// indentation >= C1*Area/kn_el   se supera el limite para el cambio de pendiente.
//            {
//                fn = kn1 * (initial_dist - max_dist) + kn2 * (max_dist - distance);
//            }
//        }
//    }
//    double mHistoryMaxForce = fn; //guardamos el maximo historico de fuerza fn
//
//    else // esta en descarga, la distancia entre particulas aumenta respecto la historica, current_dist > mHistDist
//    {
//        if (hist_fn > 0); //  fuerza normal esta en el rango de compresion de la curva
//        {
//            plast_dist = max_dist; // initial_dist*(1-plast_def) distancia associada al valor de plast_def impuesto, por ahora coincide con el cambio de pendiente
//            if (plast_dist > mHistDist) // mientras se este por encima de la maxima historica, estamos en plasticidad.
//            {
//                fn = hist_fn + kn_el_plas * (mHistDist - distance); // en descarga: 500 - kn_plas(10 - 12). en carga 500 - kn_el(10 - 11) pero con distance > mHistDistance
//            } else // esta en descarga pero no en la zona plastica, descarga por la linea elastica
//            {
//                if (distance < max_dist) // se supera el limite para le cambio de pendiente, mientras plast_dist=max_dist nunca pasara, nunca descargara por kn2
//                {
//                    NonlinearNormalForceCalculation(LocalElasticContactForce, kn1, kn2, distance, max_dist, initial_dist);
//                } else // descarga por la primera rama elastica
//                {
//                    LocalElasticContactForce[2] = kn1 * (initial_dist - distance); // fuerza en parte lineal
//                }
//            }
//        }
//    }
//}












