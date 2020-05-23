//Authors: M.A. Celigueta and G. Casas (CIMNE)
//Date: January 2016

#include "DEM_D_Bentonite_Colloid_CL.h"
#include "custom_elements/spheric_particle.h"

double ToThePower(double base, unsigned int n)
{
    if (n == 0) return 1.0;
    double x = base;

    while (n > 1){
        x *= base;
        n -= 1;
    }

    return  x;
}

double GetDebyeLength(double cation_concentration)
{
    return 7.34e7 * std::sqrt(cation_concentration);
}

namespace Kratos {
    // Hard-coded values for the moment; some should probably be nodal
    DEM_D_Bentonite_Colloid::DEM_D_Bentonite_Colloid(){
        mA_H = 10e-19;
        mD_p = 2.0e-7; // particle diameter; it whould be equal for both particles or the third law of Newton will be violated
        mA_p = 0.25 * Globals::Pi * mD_p * mD_p;
        mThickness = 1.0e-9;
        mDDLCoefficient = 1.5e5;
        mEquivRadius = mD_p / Globals::Pi; // this is the "coin" equivalent radius
    }

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Bentonite_Colloid::Clone() const
    {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Bentonite_Colloid(*this));
        return p_clone;
    }

    void DEM_D_Bentonite_Colloid::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        KRATOS_INFO("DEM") << "Assigning DEM_D_Bentonite_Colloid to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    std::string DEM_D_Bentonite_Colloid::GetTypeOfLaw() {
        std::string type_of_law = "Linear";
        return type_of_law;
    }

    /////////////////////////
    // DEM-DEM INTERACTION //
    /////////////////////////

    void DEM_D_Bentonite_Colloid::CalculateForces(const ProcessInfo& r_process_info,
                                                  const double OldLocalContactForce[3],
                                                  double LocalElasticContactForce[3],
                                                  double LocalDeltDisp[3],
                                                  double LocalRelVel[3],
                                                  double indentation,
                                                  double previous_indentation,
                                                  double ViscoDampingLocalContactForce[3],
                                                  double& cohesive_force,
                                                  SphericParticle* element1,
                                                  SphericParticle* element2,
                                                  bool& sliding,
                                                  double LocalCoordSystem[3][3]){

        //InitializeContact(element1, element2, indentation);
//G
        //LocalElasticContactForce[2]  = CalculateNormalForce(indentation);
        if ((element2->Is(BLOCKED) && element1->Is(NEW_ENTITY)) || (element2->Is(NEW_ENTITY) && element1->Is(BLOCKED))){ // you are contacting an injector
//            const array_1d<double, 3>& global_force = element2->GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE);
//            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, global_force, LocalElasticContactForce);
//            if (element2->Is(ACTIVE) && element1->IsNot(NEW_ENTITY)){ // it has a particle inside, which should already be doing the pushing, or else it is an injected particle, which already has had its force imposed
//                LocalElasticContactForce[2] = 0.0;
//
//            }

//            else if (element2->Is(ACTIVE) && element1->Is(NEW_ENTITY)){
//                const array_1d<double, 3>& global_force = element2->GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE);
//                //GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, global_force, LocalElasticContactForce);
//            }

//            else { // it must push to keep neighbours away and allow a new particle to be inserted inside
//                const double cation_concentration = element1->GetGeometry()[0].FastGetSolutionStepValue(CATION_CONCENTRATION);
//                LocalElasticContactForce[2] = CalculateNormalForce(5e-8, cation_concentration);
//            }

        }

        else { // you are contacting a regular ball, do normal ball-to-ball force evaluation
            const double distance = element1->GetInteractionRadius() + element2->GetInteractionRadius() - indentation;
            const double cation_concentration = element1->GetGeometry()[0].FastGetSolutionStepValue(CATION_CONCENTRATION);
            LocalElasticContactForce[0] = 0.0;
            LocalElasticContactForce[1] = 0.0;
            LocalElasticContactForce[2] = CalculateNormalForce(distance, cation_concentration);
        }
//Z
        cohesive_force              = CalculateCohesiveNormalForce(element1, element2, indentation);
        CalculateViscoDampingForce(LocalRelVel, ViscoDampingLocalContactForce, element1, element2);

        double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];

        /*if (normal_contact_force < 0.0) {
            normal_contact_force = 0.0;
            ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
        }*/

        CalculateTangentialForceWithNeighbour(normal_contact_force, OldLocalContactForce, LocalElasticContactForce, ViscoDampingLocalContactForce, LocalDeltDisp,
                                        sliding, element1, element2, indentation, previous_indentation);

    }

    void DEM_D_Bentonite_Colloid::CalculateViscoDampingForce(double LocalRelVel[3],
                                                             double ViscoDampingLocalContactForce[3],
                                                             SphericParticle* const element1,
                                                             SphericParticle* const element2) {

//        const double my_mass    = element1->GetMass();
//        const double other_mass = element2->GetMass();

//        const double equiv_mass = 1.0 / (1.0/my_mass + 1.0/other_mass);

//        const double my_gamma    = element1->GetProperties()[DAMPING_GAMMA];
//        const double other_gamma = element2->GetProperties()[DAMPING_GAMMA];
//        const double equiv_gamma = 0.5 * (my_gamma + other_gamma);
//        const double equiv_visco_damp_coeff_normal     = 2.0 * equiv_gamma * sqrt(equiv_mass * mKn);
//        const double equiv_visco_damp_coeff_tangential = 2.0 * equiv_gamma * sqrt(equiv_mass * mKt);

//        ViscoDampingLocalContactForce[0] = - equiv_visco_damp_coeff_tangential * LocalRelVel[0];
//        ViscoDampingLocalContactForce[1] = - equiv_visco_damp_coeff_tangential * LocalRelVel[1];
//        ViscoDampingLocalContactForce[2] = - equiv_visco_damp_coeff_normal     * LocalRelVel[2];

        ViscoDampingLocalContactForce[0] = 0.0;
        ViscoDampingLocalContactForce[1] = 0.0;
        ViscoDampingLocalContactForce[2] = 0.0;
    }

    /////////////////////////
    // DEM-FEM INTERACTION //
    /////////////////////////

    void DEM_D_Bentonite_Colloid::CalculateForcesWithFEM(ProcessInfo& r_process_info,
                                                         const double OldLocalContactForce[3],
                                                         double LocalElasticContactForce[3],
                                                         double LocalDeltDisp[3],
                                                         double LocalRelVel[3],
                                                         double indentation,
                                                         double previous_indentation,
                                                         double ViscoDampingLocalContactForce[3],
                                                         double& cohesive_force,
                                                         SphericParticle* const element,
                                                         Condition* const wall,
                                                         bool& sliding) {

        //InitializeContactWithFEM(element, wall, indentation);
        const double distance = element->GetInteractionRadius() - indentation;
        const double cation_concentration = element->GetGeometry()[0].FastGetSolutionStepValue(CATION_CONCENTRATION);

        const double smoother = 1.0;//std::max(1.0, 9.0 * indentation / (element1->GetInteractionRadius() + element2->GetInteractionRadius()));
        LocalElasticContactForce[2] = smoother * CalculateNormalForce(distance, cation_concentration);
        //LocalElasticContactForce[2]  = smoother * (CalculateVanDerWaalsForce(distance) - 0.0000001 * CalculateVanDerWaalsForce(pow(distance, 1.2))) ;

        //CalculateViscoDampingForceWithFEM(LocalRelVel, ViscoDampingLocalContactForce, element, wall);
    }

    template<class NeighbourClassType>
    void DEM_D_Bentonite_Colloid::CalculateTangentialForceWithNeighbour(const double normal_contact_force,
                                                                          const double OldLocalContactForce[3],
                                                                          double LocalElasticContactForce[3],
                                                                          double ViscoDampingLocalContactForce[3],
                                                                          const double LocalDeltDisp[3],
                                                                          bool& sliding,
                                                                          SphericParticle* const element,
                                                                          NeighbourClassType* const neighbour,
                                                                          double indentation,
                                                                          double previous_indentation) {

//        LocalElasticContactForce[0] = 0.0;
//        LocalElasticContactForce[1] = 0.0;
//        ViscoDampingLocalContactForce[0] = 0.0;
//        ViscoDampingLocalContactForce[1] = 0.0;
    }

    void DEM_D_Bentonite_Colloid::CalculateViscoDampingForceWithFEM(double LocalRelVel[3],
                                                                    double ViscoDampingLocalContactForce[3],
                                                                    SphericParticle* const element,
                                                                    Condition* const wall) {

        const double my_mass    = element->GetMass();
        const double gamma = element->GetProperties()[DAMPING_GAMMA];
        const double normal_damping_coefficient     = 2.0 * gamma * sqrt(my_mass * mKn);
        const double tangential_damping_coefficient = 2.0 * gamma * sqrt(my_mass * mKt);

        ViscoDampingLocalContactForce[0] = - tangential_damping_coefficient * LocalRelVel[0];
        ViscoDampingLocalContactForce[1] = - tangential_damping_coefficient * LocalRelVel[1];
        ViscoDampingLocalContactForce[2] = - normal_damping_coefficient     * LocalRelVel[2];

    }

    double DEM_D_Bentonite_Colloid::CalculateNormalForce(const double distance, const double cation_concentration){

        const double F_vdW = CalculateVanDerWaalsForce(distance);
        const double F_DDL = CalculateDiffuseDoubleLayerForce(distance, cation_concentration);
        return F_vdW + F_DDL;
    }

    double DEM_D_Bentonite_Colloid::CalculateCohesiveNormalForce(SphericParticle* const element1, SphericParticle* const element2, const double indentation){
        return 0.0;
    }

    double DEM_D_Bentonite_Colloid::CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, Condition* const wall, const double indentation){
        return 0.0;
    }

    double DEM_D_Bentonite_Colloid::CalculateVanDerWaalsForce(const double distance)
    {
        return - mA_p * mA_H / (6.0 * Globals::Pi) * (1.0 / ToThePower(distance, 3) - 2.0 / ToThePower(distance + mThickness, 3) + 1.0 / ToThePower(distance + 2 * mThickness, 3));
    }

    double DEM_D_Bentonite_Colloid::CalculateDiffuseDoubleLayerForce(const double distance, const double cation_concentration)
    {
        return mA_p * mDDLCoefficient * cation_concentration * exp(- GetDebyeLength(cation_concentration) * distance);
    }

} // namespace Kratos
