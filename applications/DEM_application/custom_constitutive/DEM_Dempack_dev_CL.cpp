// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "DEM_Dempack_dev_CL.h"
//#include "custom_elements/spheric_particle.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_Dempack_dev::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_Dempack_dev(*this));
        return p_clone;
    }

    void DEM_Dempack_dev::SetConstitutiveLawInProperties(Properties::Pointer pProp) const {
        std::cout << "Assigning DEM_Dempack_dev to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }


    void DEM_Dempack_dev::GetContactArea(const double radius,
                                     const double other_radius,
                                     const Vector& vector_of_initial_areas,
                                     const int neighbour_position,
                                     double& calculation_area)
    {
        if (vector_of_initial_areas.size()) calculation_area = vector_of_initial_areas[neighbour_position];
        else CalculateContactArea(radius, other_radius, calculation_area);
        //CalculateContactArea(radius, other_radius, calculation_area);

    }

    void DEM_Dempack_dev::CalculateElasticConstants(double &kn_el,
                                                double &kt_el,
                                                double initial_dist,
                                                double equiv_young,
                                                double equiv_poisson,
                                                double calculation_area,
                                                SphericContinuumParticle* element1,
                                                SphericContinuumParticle* element2) {
        
        KRATOS_TRY 
        double equiv_shear = equiv_young / (2.0 * (1 + equiv_poisson));
        kn_el = equiv_young * calculation_area / initial_dist;
        kt_el = equiv_shear * calculation_area / initial_dist;
        std::ofstream outputfile("knkt.txt", std::ios_base::out | std::ios_base::app);
        outputfile << kn_el << " " << kt_el << "\n";
        outputfile.close();

        KRATOS_CATCH("")  
    }



    void DEM_Dempack_dev::CalculateNormalForces(double LocalElasticContactForce[3],
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

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];

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



        double rmin = element1->GetRadius();
        const double other_radius = element2->GetRadius();
        if (other_radius < rmin) rmin = other_radius;
        double effective_calculation_area = KRATOS_M_PI * rmin*rmin;

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
        } else {

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
        const double Yields_el = mPlasticityLimit * effective_calculation_area;

        const double Ncstr1_el = mC1 * effective_calculation_area;
        const double Ncstr2_el = mC2 * effective_calculation_area;
        const double Ncstr3_el = mC3 * effective_calculation_area;
        const double Ntstr_el = mTensionLimit * effective_calculation_area;
        double u_max = mHistoryMaxInd;

        double& fn = LocalElasticContactForce[2]; //[2] means 'normal' contact force

        if (indentation >= 0.0) { //COMPRESSION

            fn = kn_el * indentation;
            double u_ela1 = Ncstr1_el / kn_el;
            double u_ela2 = u_ela1 + (Ncstr2_el - Ncstr1_el) / (kn_b);
            double u_ela3 = u_ela2 + (Ncstr3_el - Ncstr2_el) / (kn_c);

            if ((indentation > u_max) || (time_steps <= 1)) { //maximum historical intentation OR first step  MSIMSI 0

                mHistoryMaxInd = indentation;               // Guarda el threshold del màxim desplaçament

                if (indentation > u_ela3) {                 //4rt tram
                    fn = Ncstr3_el + (indentation - u_ela3) * kn_d;
                    mHistoryDegradation = kn_d / kn_el;
                } else if (indentation > u_ela2) {          //3r tram
                    fn = Ncstr2_el + (indentation - u_ela2) * kn_c;
                    mHistoryDegradation = kn_c / kn_el;
                } else {
                    if (indentation > u_ela1) {             //2n tram
                        fn = Ncstr1_el + (indentation - u_ela1) * kn_b;
                        mHistoryDegradation = kn_b / kn_el;
                    }
                }
                mHistoryMaxForce = fn;                      //actualitzar la força màxima a compressió.
            } else {                                        //Per sota del màxim.
                if (mHistoryMaxForce > 0.0) {               //Màxim en compressió.

                    double u_plas;                          //MSIMSI 2 akesta operació de saber quant val la u_plastica es fa cada pas de temps i en realitat es fixe sempre.
                    if (Yields_el <= Ncstr1_el) {           //si el punt de plastificació està en la primera rama elastica.

                        u_plas = Yields_el / kn_el;
                    } else {
                        if (Yields_el <= Ncstr2_el) {       //si està en la segona...
                            u_plas = u_ela1 + (Yields_el - Ncstr1_el) / (kn_b);
                        } else if (Yields_el <= Ncstr3_el) { //si està en la tercera...
                            u_plas = u_ela2 + (Yields_el - Ncstr2_el) / (kn_c);
                        } else {                            //en la quarta
                            u_plas = u_ela3 + (Yields_el - Ncstr3_el) / (kn_d);
                        }
                    }
                    if (u_plas < u_max) {                   //si nosaltres estem per sota del maxim pero ja estem plastificant
                        fn = mHistoryMaxForce - kp_el * (u_max - indentation); // Esta en zona de descarga plastica (pot estar en carga/descarga)
                        mHistoryDegradation = kp_el / kn_el;
                    } else {                                // Esta en zona descarga elastica, ens despreocupem de la plasticitat
                        if (indentation > u_ela3) {         // en la 4a ramma
                            fn = Ncstr3_el + (indentation - u_ela3) * kn_d;
                        } else if (indentation > u_ela2) {  // en la 3a ramma
                            fn = Ncstr2_el + (indentation - u_ela2) * kn_c;
                        } else {
                            if (indentation > u_ela1) {     // en la 2a rama
                                fn = Ncstr1_el + (indentation - u_ela1) * kn_b;
                            }
                        }
                    }
                }       //si tenim precàrrega en compressió.
            }           //Per sota del màxim.
        }               //Compression
        else {          //tension
            fn = kn_el * indentation;
            double u1 = Ntstr_el / kn_el;
            double u2 = u1 * (1 + mDamageMaxDisplacementFactor);

            if (failure_type == 0) {

                if (fabs(indentation) > u2) {   // FULL DAMAGE
                    failure_type = 4;           //tension failure
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
            }
            else {fn = 0.0;}
        } //Tension

    KRATOS_CATCH("")
    }


      /////// check old tangential forces calculations with unstable rotations
/*
    void DEM_Dempack_dev::CalculateTangentialForces(double OldLocalElasticContactForce[3],
            double LocalElasticContactForce[3],
            double LocalElasticExtraContactForce[3],
            double LocalCoordSystem[3][3],
            double LocalDeltDisp[3],
            const double kt_el,
            const double equiv_shear,
            double& contact_sigma,
            double& contact_tau,
            double indentation,
            double calculation_area,
            double& failure_criterion_state,
            SphericContinuumParticle* element1,
            SphericContinuumParticle* element2,
            int i_neighbour_count,
            bool& sliding,
            int search_control,
            vector<int>& search_control_vector,
            const ProcessInfo& r_process_info) {

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];


        double rmin = element1->GetRadius();
        const double other_radius = element2->GetRadius();
        if (other_radius < rmin) rmin = other_radius;
        double effective_calculation_area = KRATOS_M_PI * rmin*rmin;

        Properties& element1_props = element1->GetProperties();
        Properties& element2_props = element2->GetProperties();

        if (r_process_info[SHEAR_STRAIN_PARALLEL_TO_BOND_OPTION]) {
            AddContributionOfShearStrainParallelToBond(OldLocalElasticContactForce, LocalElasticExtraContactForce, LocalCoordSystem, kt_el, calculation_area,  element1, element2);
        }

        double mTensionLimit;
        double mTauZero;
        double mInternalFriccion;
        double mShearEnergyCoef;

        if (&element1_props == &element2_props) {

            mTensionLimit = element1_props[CONTACT_SIGMA_MIN]*1e6;
            mTauZero = element1_props[CONTACT_TAU_ZERO]*1e6;
            mInternalFriccion = element1_props[CONTACT_INTERNAL_FRICC];
            mShearEnergyCoef = element1_props[SHEAR_ENERGY_COEF];
        }
        else{

            mTensionLimit = 0.5*1e6*(element1_props[CONTACT_SIGMA_MIN] + element2_props[CONTACT_SIGMA_MIN]);
            mTauZero = 0.5*1e6*(element1_props[CONTACT_TAU_ZERO] + element2_props[CONTACT_TAU_ZERO]);
            mInternalFriccion = 0.5*(element1_props[CONTACT_INTERNAL_FRICC] + element2_props[CONTACT_INTERNAL_FRICC]);
            mShearEnergyCoef = 0.5*(element1_props[SHEAR_ENERGY_COEF] + element2_props[SHEAR_ENERGY_COEF]);
        }

        double degradation = 1.0; //Tangential. With degradation:

        if (i_neighbour_count < int(element1->mContinuumInitialNeighborsSize)) {
            if (indentation >= 0.0) { //COMPRESSION
                degradation = mHistoryDegradation;
            } else {
                degradation = (1.0 - mHistoryDamage);
            }
        }

        if (failure_type == 0) { //This means it has not broken

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
                    failure_type = 2; // shear
                    //failure_criterion_state = 1.0;    //
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
            if (failure_type != 0) {
                search_control_vector[OpenMPUtils::ThisThread()] = 1;
            }
        }
    }

*/


    void DEM_Dempack_dev::CalculateTangentialForces(double OldLocalElasticContactForce[3],
            double LocalElasticContactForce[3],
            double LocalElasticExtraContactForce[3],
            double LocalCoordSystem[3][3],
            double LocalDeltDisp[3],
            const double kt_el,
            const double equiv_shear,
            double& contact_sigma,
            double& contact_tau,
            double indentation,
            double calculation_area,
            double& failure_criterion_state,
            SphericContinuumParticle* element1,
            SphericContinuumParticle* element2,
            int i_neighbour_count,
            bool& sliding,
            int search_control,
            vector<int>& search_control_vector,
            const ProcessInfo& r_process_info) {

        KRATOS_TRY
        const int time_steps = r_process_info[TIME_STEPS];
        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
        LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - kt_el * LocalDeltDisp[0]; // 0: first tangential
        LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - kt_el * LocalDeltDisp[1]; // 1: second tangential

        if (r_process_info[SHEAR_STRAIN_PARALLEL_TO_BOND_OPTION]) {
            AddContributionOfShearStrainParallelToBond(OldLocalElasticContactForce, LocalElasticExtraContactForce, element1->mNeighbourElasticExtraContactForces[i_neighbour_count], LocalCoordSystem, kt_el, calculation_area,  element1, element2);
        }

        double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
                                  + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

        double rmin = element1->GetRadius();
        const double other_radius = element2->GetRadius();
        if (other_radius < rmin) rmin = other_radius;
        double effective_calculation_area = KRATOS_M_PI * rmin*rmin;

        if (failure_type == 0) { // This means it has not broken
            Properties& element1_props = element1->GetProperties();
            Properties& element2_props = element2->GetProperties();
            const double mTauZero = 0.5 * 1e6 * (element1_props[CONTACT_TAU_ZERO] + element2_props[CONTACT_TAU_ZERO]);
            const double mInternalFriction = 0.5 * (element1_props[CONTACT_INTERNAL_FRICC] + element2_props[CONTACT_INTERNAL_FRICC]);

            contact_tau = ShearForceNow / calculation_area;
            contact_sigma = LocalElasticContactForce[2] / calculation_area;

            double tau_strength = mTauZero;

            if (contact_sigma >= 0) {
                tau_strength = mTauZero + mInternalFriction * contact_sigma;
            }

            if (contact_tau > tau_strength) {
                failure_type = 2; // shear
            }
        }
        else {
            const double equiv_tg_of_fri_ang = 0.5 * (element1->GetTgOfFrictionAngle() + element2->GetTgOfFrictionAngle());
            double Frictional_ShearForceMax = equiv_tg_of_fri_ang * LocalElasticContactForce[2];
            if (Frictional_ShearForceMax < 0.0) {
                Frictional_ShearForceMax = 0.0;
            }

            if ((ShearForceNow > Frictional_ShearForceMax) && (ShearForceNow != 0.0)) {
                LocalElasticContactForce[0] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[0];
                LocalElasticContactForce[1] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[1];
                sliding = true;
            }
        }

        KRATOS_CATCH("")
    }

    
    void DEM_Dempack_dev::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                    SphericContinuumParticle* neighbor,
                                                    double equiv_young,
                                                    double distance,
                                                    double calculation_area,
                                                    double LocalCoordSystem[3][3],
                                                    double ElasticLocalRotationalMoment[3],
                                                    double ViscoLocalRotationalMoment[3],
                                                    double equiv_poisson,
                                                    double indentation) {}//ComputeParticleRotationalMoments




    void DEM_Dempack_dev::AddContributionOfShearStrainParallelToBond(double OldLocalElasticContactForce[3],
                                                              double LocalElasticExtraContactForce[3],
                                                              array_1d<double, 3>& OldElasticExtraContactForce,
                                                              double LocalCoordSystem[3][3],
                                                              const double kt_el,
                                                              const double calculation_area,
                                                              SphericContinuumParticle* element1,
                                                              SphericContinuumParticle* element2) {



        if (element1->mSymmStressTensor == NULL /*|| element1->mOldSymmStressTensor == NULL*/) return;
        if(element1->IsSkin() || element2->IsSkin()) return;
        //unsigned int element_id = element1->Id();
        //unsigned int element_id2 = element2->Id();

        double average_stress_tensor[3][3];

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                average_stress_tensor[i][j]     = 0.5 * ((*(element1->mSymmStressTensor   ))(i,j) + (*(element2->mSymmStressTensor   ))(i,j));
            }
        }

        double current_sigma_local[3][3];
        GeometryFunctions::TensorGlobal2Local(LocalCoordSystem, average_stress_tensor, current_sigma_local);

        const double current_tangential_stress1 = current_sigma_local[0][2];
        const double current_tangential_stress2 = current_sigma_local[1][2];

        LocalElasticExtraContactForce[0] -= calculation_area * current_tangential_stress1;
        LocalElasticExtraContactForce[1] -= calculation_area * current_tangential_stress2;

        double tangential_delta_acumulated_0 = OldLocalElasticContactForce[0] / kt_el;
        double tangential_delta_acumulated_1 = OldLocalElasticContactForce[1] / kt_el;

        LocalElasticExtraContactForce[0] += kt_el * tangential_delta_acumulated_0;
        LocalElasticExtraContactForce[1] += kt_el * tangential_delta_acumulated_1;

    }

} /* namespace Kratos.*/
