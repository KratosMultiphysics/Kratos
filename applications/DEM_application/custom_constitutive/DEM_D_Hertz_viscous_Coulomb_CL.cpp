// Authors: M.A. Celigueta and S. Latorre (CIMNE)
// Date: October 2015

#include "DEM_D_Hertz_viscous_Coulomb_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    void DEM_D_Hertz_viscous_Coulomb::Initialize(const ProcessInfo& r_process_info) {}

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Hertz_viscous_Coulomb::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Hertz_viscous_Coulomb(*this));
        return p_clone;
    }

    void DEM_D_Hertz_viscous_Coulomb::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) const {
        if (verbose) KRATOS_INFO("DEM") << "Assigning DEM_D_Hertz_viscous_Coulomb to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    std::string DEM_D_Hertz_viscous_Coulomb::GetTypeOfLaw() {
        std::string type_of_law = "Hertz";
        return type_of_law;
    }

    /////////////////////////
    // DEM-DEM INTERACTION //
    /////////////////////////

    void DEM_D_Hertz_viscous_Coulomb::InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {
        //Get equivalent Radius
        const double my_radius       = element1->GetRadius();
        const double other_radius    = element2->GetRadius();
        const double radius_sum      = my_radius + other_radius;
        const double radius_sum_inv  = 1.0 / radius_sum;
        const double equiv_radius    = my_radius * other_radius * radius_sum_inv;

        //Get equivalent Young's Modulus
        const double my_young        = element1->GetYoung();
        const double other_young     = element2->GetYoung();
        const double my_poisson      = element1->GetPoisson();
        const double other_poisson   = element2->GetPoisson();
        const double equiv_young     = my_young * other_young / (other_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - other_poisson * other_poisson));

        //Get equivalent Shear Modulus
        const double my_shear_modulus = 0.5 * my_young / (1.0 + my_poisson);
        const double other_shear_modulus = 0.5 * other_young / (1.0 + other_poisson);
        const double equiv_shear = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - other_poisson)/other_shear_modulus);

        //Normal and Tangent elastic constants
        const double sqrt_equiv_radius_and_indentation = sqrt(equiv_radius * indentation);
        mKn = 2.0 * equiv_young * sqrt_equiv_radius_and_indentation;
        mKt = 4.0 * equiv_shear * mKn / equiv_young;
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateForces(const ProcessInfo& r_process_info,
                                                      const double OldLocalElasticContactForce[3],
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
                                                      double LocalCoordSystem[3][3]) {

        InitializeContact(element1, element2, indentation);

        LocalElasticContactForce[2]  = CalculateNormalForce(element1, element2, indentation, LocalCoordSystem);
        cohesive_force               = CalculateCohesiveNormalForce(element1, element2, indentation);

        CalculateViscoDampingForce(LocalRelVel, ViscoDampingLocalContactForce, element1, element2);

        double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];

        if (normal_contact_force < 0.0) {
            normal_contact_force = 0.0;
            ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
        }
        /* //COATING
        array_1d<double, 3> other_vel = ZeroVector(3);
        if (element2->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)) {
            other_vel = element2->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
            element2->GetGeometry()[0].Set(TO_ERASE); // This line applies only for coating simulations
        }
        double other_vel_module = sqrt(other_vel[0] * other_vel[0] + other_vel[1] * other_vel[1] + other_vel[2] * other_vel[2]);
        element1->GetGeometry()[0].FastGetSolutionStepValue(SPRAYED_MATERIAL) += 0.00001 * other_vel_module;
        */
        double AuxElasticShearForce;
        double MaximumAdmisibleShearForce;

        CalculateTangentialForceWithNeighbour(normal_contact_force, OldLocalElasticContactForce, LocalElasticContactForce, ViscoDampingLocalContactForce, LocalDeltDisp,
                                              sliding, element1, element2, indentation, previous_indentation, AuxElasticShearForce, MaximumAdmisibleShearForce);

        double& elastic_energy = element1->GetElasticEnergy();
        CalculateElasticEnergyDEM(elastic_energy, indentation, LocalElasticContactForce);

        if(sliding && MaximumAdmisibleShearForce != 0.0){
            double& inelastic_frictional_energy = element1->GetInelasticFrictionalEnergy();
            CalculateInelasticFrictionalEnergyDEM(inelastic_frictional_energy, AuxElasticShearForce, LocalElasticContactForce);
        }

        double& inelastic_viscodamping_energy = element1->GetInelasticViscodampingEnergy();
        CalculateInelasticViscodampingEnergyDEM(inelastic_viscodamping_energy, ViscoDampingLocalContactForce, LocalDeltDisp);

    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateViscoDampingForce(double LocalRelVel[3],
                                                                 double ViscoDampingLocalContactForce[3],
                                                                 SphericParticle* const element1,
                                                                 SphericParticle* const element2) {

        const double my_mass    = element1->GetMass();
        const double other_mass = element2->GetMass();

        const double equiv_mass = 1.0 / (1.0/my_mass + 1.0/other_mass);

        const double my_gamma    = element1->GetProperties()[DAMPING_GAMMA];
        const double other_gamma = element2->GetProperties()[DAMPING_GAMMA];
        const double equiv_gamma = 0.5 * (my_gamma + other_gamma);

        const double equiv_visco_damp_coeff_normal     = 2.0 * equiv_gamma * sqrt(equiv_mass * mKn);
        const double equiv_visco_damp_coeff_tangential = 2.0 * equiv_gamma * sqrt(equiv_mass * mKt);

        ViscoDampingLocalContactForce[0] = - equiv_visco_damp_coeff_tangential * LocalRelVel[0];
        ViscoDampingLocalContactForce[1] = - equiv_visco_damp_coeff_tangential * LocalRelVel[1];
        ViscoDampingLocalContactForce[2] = - equiv_visco_damp_coeff_normal     * LocalRelVel[2];
    }

    /////////////////////////
    // DEM-FEM INTERACTION //
    /////////////////////////
    
    void DEM_D_Hertz_viscous_Coulomb::InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta) {
        //Get effective Radius
        const double my_radius           = element->GetRadius(); //Get equivalent Radius
        const double effective_radius    = my_radius - ini_delta;

        //Get equivalent Young's Modulus
        const double my_young            = element->GetYoung();
        const double walls_young         = wall->GetProperties()[YOUNG_MODULUS];
        const double my_poisson          = element->GetPoisson();
        const double walls_poisson       = wall->GetProperties()[POISSON_RATIO];
        const double equiv_young         = my_young * walls_young / (walls_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - walls_poisson * walls_poisson));

        //Get equivalent Shear Modulus
        const double my_shear_modulus    = 0.5 * my_young / (1.0 + my_poisson);
        const double walls_shear_modulus = 0.5 * walls_young / (1.0 + walls_poisson);
        const double equiv_shear         = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - walls_poisson)/walls_shear_modulus);

        //Normal and Tangent elastic constants
        const double sqrt_equiv_radius_and_indentation = sqrt(effective_radius * indentation);
        mKn = 2.0 * equiv_young * sqrt_equiv_radius_and_indentation;
        mKt = 4.0 * equiv_shear * mKn / equiv_young;
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateForcesWithFEM(ProcessInfo& r_process_info,
                                                             const double OldLocalElasticContactForce[3],
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

        InitializeContactWithFEM(element, wall, indentation);

        LocalElasticContactForce[2] = CalculateNormalForce(element, wall, indentation);
        cohesive_force              = CalculateCohesiveNormalForceWithFEM(element, wall, indentation);

        CalculateViscoDampingForceWithFEM(LocalRelVel, ViscoDampingLocalContactForce, element, wall);

        double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];

        if (normal_contact_force < 0.0) {
            normal_contact_force = 0.0;
            ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
        }

        double AuxElasticShearForce;
        double MaximumAdmisibleShearForce;

        CalculateTangentialForceWithNeighbour(normal_contact_force, OldLocalElasticContactForce, LocalElasticContactForce, ViscoDampingLocalContactForce, LocalDeltDisp,
                                              sliding, element, wall, indentation, previous_indentation, AuxElasticShearForce, MaximumAdmisibleShearForce);

        double& elastic_energy = element->GetElasticEnergy();
        CalculateElasticEnergyFEM(elastic_energy, indentation, LocalElasticContactForce);//MSIMSI

        if(sliding && MaximumAdmisibleShearForce != 0.0){
            double& inelastic_frictional_energy = element->GetInelasticFrictionalEnergy();
            CalculateInelasticFrictionalEnergyFEM(inelastic_frictional_energy, AuxElasticShearForce, LocalElasticContactForce);
        }

        double& inelastic_viscodamping_energy = element->GetInelasticViscodampingEnergy();
        CalculateInelasticViscodampingEnergyFEM(inelastic_viscodamping_energy, ViscoDampingLocalContactForce, LocalDeltDisp);
    }

    template<class NeighbourClassType>

    void DEM_D_Hertz_viscous_Coulomb::CalculateTangentialForceWithNeighbour(const double normal_contact_force,
                                                                            const double OldLocalElasticContactForce[3],
                                                                            double LocalElasticContactForce[3],
                                                                            double ViscoDampingLocalContactForce[3],
                                                                            const double LocalDeltDisp[3],
                                                                            bool& sliding,
                                                                            SphericParticle* const element,
                                                                            NeighbourClassType* const neighbour,
                                                                            double indentation,
                                                                            double previous_indentation,
                                                                            double& AuxElasticShearForce,
                                                                            double& MaximumAdmisibleShearForce) {

        LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - mKt * LocalDeltDisp[0];
        LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - mKt * LocalDeltDisp[1];

        if (previous_indentation > indentation) {
            const double minoring_factor = sqrt (indentation / previous_indentation);
            LocalElasticContactForce[0] = OldLocalElasticContactForce[0] * minoring_factor - mKt * LocalDeltDisp[0];
            LocalElasticContactForce[1] = OldLocalElasticContactForce[1] * minoring_factor - mKt * LocalDeltDisp[1];
        }

        AuxElasticShearForce = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

        const double my_tg_of_friction_angle        = element->GetTgOfFrictionAngle();
        const double neighbour_tg_of_friction_angle = neighbour->GetProperties()[FRICTION];
        const double equiv_tg_of_fri_ang            = 0.5 * (my_tg_of_friction_angle + neighbour_tg_of_friction_angle);    
        
        MaximumAdmisibleShearForce = normal_contact_force * equiv_tg_of_fri_ang;

        const double tangential_contact_force_0 = LocalElasticContactForce[0] + ViscoDampingLocalContactForce[0];
        const double tangential_contact_force_1 = LocalElasticContactForce[1] + ViscoDampingLocalContactForce[1];

        const double ActualTotalShearForce = sqrt(tangential_contact_force_0 * tangential_contact_force_0 + tangential_contact_force_1 * tangential_contact_force_1);

        if (ActualTotalShearForce > MaximumAdmisibleShearForce) {

            const double ActualElasticShearForce = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

            const double dot_product = LocalElasticContactForce[0] * ViscoDampingLocalContactForce[0] + LocalElasticContactForce[1] * ViscoDampingLocalContactForce[1];
            const double ViscoDampingLocalContactForceModule = sqrt(ViscoDampingLocalContactForce[0] * ViscoDampingLocalContactForce[0] +\
                                                                    ViscoDampingLocalContactForce[1] * ViscoDampingLocalContactForce[1]);

            if (dot_product >= 0.0) {

                if (ActualElasticShearForce > MaximumAdmisibleShearForce) {
                    const double fraction = MaximumAdmisibleShearForce / ActualElasticShearForce;
                    LocalElasticContactForce[0]      = LocalElasticContactForce[0] * fraction;
                    LocalElasticContactForce[1]      = LocalElasticContactForce[1] * fraction;
                    ViscoDampingLocalContactForce[0] = 0.0;
                    ViscoDampingLocalContactForce[1] = 0.0;
                }
                else {
                    const double ActualViscousShearForce = MaximumAdmisibleShearForce - ActualElasticShearForce;
                    const double fraction = ActualViscousShearForce / ViscoDampingLocalContactForceModule;
                    ViscoDampingLocalContactForce[0] *= fraction;
                    ViscoDampingLocalContactForce[1] *= fraction;
                }
            }
            else {
                if (ViscoDampingLocalContactForceModule >= ActualElasticShearForce) {
                    const double fraction = (MaximumAdmisibleShearForce + ActualElasticShearForce) / ViscoDampingLocalContactForceModule;
                    ViscoDampingLocalContactForce[0] *= fraction;
                    ViscoDampingLocalContactForce[1] *= fraction;
                }
                else {
                    const double fraction = MaximumAdmisibleShearForce / ActualElasticShearForce;
                    LocalElasticContactForce[0]      = LocalElasticContactForce[0] * fraction;
                    LocalElasticContactForce[1]      = LocalElasticContactForce[1] * fraction;
                    ViscoDampingLocalContactForce[0] = 0.0;
                    ViscoDampingLocalContactForce[1] = 0.0;
                }
            }
            sliding = true;
        }
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateViscoDampingForceWithFEM(double LocalRelVel[3],
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

    double DEM_D_Hertz_viscous_Coulomb::CalculateNormalForce(const double indentation) {
        return 0.666666666666666666667 * mKn * indentation;
    }

    double DEM_D_Hertz_viscous_Coulomb::CalculateCohesiveNormalForce(SphericParticle* const element1, SphericParticle* const element2, const double indentation){
        return 0.0;
    }
    
    double DEM_D_Hertz_viscous_Coulomb::CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, Condition* const wall, const double indentation){
        return 0.0;
    }

    //MSIMSI
    void DEM_D_Hertz_viscous_Coulomb::CalculateElasticEnergyDEM(double& elastic_energy, double indentation, double LocalElasticContactForce[3])
    {
        double normal_elastic         = 0.20*LocalElasticContactForce[2]*indentation; //each ball in a contact with another ball receives half the contact energy
        double tangential_elastic     = 0.25*(LocalElasticContactForce[0]*LocalElasticContactForce[0]+LocalElasticContactForce[1]*LocalElasticContactForce[1])/mKt; //each ball in a contact with another ball receives half the contact energy
        elastic_energy += normal_elastic;
        elastic_energy += tangential_elastic;
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateInelasticFrictionalEnergyDEM(double& inelastic_frictional_energy, double& AuxElasticShearForce, double LocalElasticContactForce[3])
    {
        double frictional_energy = 0.25*((AuxElasticShearForce*AuxElasticShearForce)-(LocalElasticContactForce[0]*LocalElasticContactForce[0]+LocalElasticContactForce[1]*LocalElasticContactForce[1]))/mKt;
        inelastic_frictional_energy += frictional_energy;
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateInelasticViscodampingEnergyDEM(double& inelastic_viscodamping_energy, double ViscoDampingLocalContactForce[3], double LocalDeltDisp[3])
    {
        double viscodamping_energy = 0.50*sqrt(ViscoDampingLocalContactForce[0]*ViscoDampingLocalContactForce[0]*LocalDeltDisp[0]*LocalDeltDisp[0]+ViscoDampingLocalContactForce[1]*ViscoDampingLocalContactForce[1]*LocalDeltDisp[1]*LocalDeltDisp[1]+ViscoDampingLocalContactForce[2]*ViscoDampingLocalContactForce[2]*LocalDeltDisp[2]*LocalDeltDisp[2]);
        inelastic_viscodamping_energy += viscodamping_energy;
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateElasticEnergyFEM(double& elastic_energy, double indentation, double LocalElasticContactForce[3])
    {
        double normal_elastic     = 0.40*LocalElasticContactForce[2]*indentation; //each ball in a contact with a wall receives all the contact energy
        double tangential_elastic = 0.50*(LocalElasticContactForce[0]*LocalElasticContactForce[0]+LocalElasticContactForce[1]*LocalElasticContactForce[1])/mKt;
        elastic_energy += normal_elastic;
        elastic_energy += tangential_elastic;
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateInelasticFrictionalEnergyFEM(double& inelastic_frictional_energy, double& AuxElasticShearForce, double LocalElasticContactForce[3])
    {
        double frictional_energy = 0.50*((AuxElasticShearForce*AuxElasticShearForce)-(LocalElasticContactForce[0]*LocalElasticContactForce[0]+LocalElasticContactForce[1]*LocalElasticContactForce[1]))/mKt;
        inelastic_frictional_energy += frictional_energy;
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateInelasticViscodampingEnergyFEM(double& inelastic_viscodamping_energy, double ViscoDampingLocalContactForce[3], double LocalDeltDisp[3])
    {
        double viscodamping_energy = sqrt(ViscoDampingLocalContactForce[0]*ViscoDampingLocalContactForce[0]*LocalDeltDisp[0]*LocalDeltDisp[0]+ViscoDampingLocalContactForce[1]*ViscoDampingLocalContactForce[1]*LocalDeltDisp[1]*LocalDeltDisp[1]+ViscoDampingLocalContactForce[2]*ViscoDampingLocalContactForce[2]*LocalDeltDisp[2]*LocalDeltDisp[2]);
        inelastic_viscodamping_energy += viscodamping_energy;
    }

} // namespace Kratos
