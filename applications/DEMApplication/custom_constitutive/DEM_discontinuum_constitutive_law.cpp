// Last modified by: S. Latorre (CIMNE)
// Date: October 2015

#include "DEM_discontinuum_constitutive_law.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    DEMDiscontinuumConstitutiveLaw::DEMDiscontinuumConstitutiveLaw() {
        //KRATOS_INFO("DEM") << " DEMDiscontinuumConstitutiveLaw constructor..." << std::endl;

    } // Class DEMDiscontinuumConstitutiveLaw

    DEMDiscontinuumConstitutiveLaw::DEMDiscontinuumConstitutiveLaw(const DEMDiscontinuumConstitutiveLaw &rReferenceDiscontinuumConstitutiveLaw) {
        //KRATOS_INFO("DEM") << " DEMDiscontinuumConstitutiveLaw copy constructor..." << std::endl;
    }

    //    DEMDiscontinuumConstitutiveLaw::DEMDiscontinuumConstitutiveLaw( const DEMDiscontinuumConstitutiveLaw &rReferenceDiscontinuumConstitutiveLaw){
    //        //KRATOS_INFO("DEM") << " DEMDiscontinuumConstitutiveLaw copy constructor..." << std::endl;
    //    }

    void DEMDiscontinuumConstitutiveLaw::Initialize(const ProcessInfo& r_process_info) {
    }

    void DEMDiscontinuumConstitutiveLaw::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        //if (verbose) KRATOS_INFO("DEM") << "Assigning DEMDiscontinuumConstitutiveLaw to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
    }

    void DEMDiscontinuumConstitutiveLaw::Check(Properties::Pointer pProp) const {
        if(!pProp->Has(FRICTION)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable FRICTION should be present in the properties when using DEMDiscontinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(FRICTION) = 0.0;
        }
        if(!pProp->Has(YOUNG_MODULUS)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable YOUNG_MODULUS should be present in the properties when using DEMDiscontinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(YOUNG_MODULUS) = 0.0;
        }
        if(!pProp->Has(POISSON_RATIO)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable POISSON_RATIO should be present in the properties when using DEMDiscontinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(POISSON_RATIO) = 0.0;
        }
        if(!pProp->Has(COEFFICIENT_OF_RESTITUTION)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable COEFFICIENT_OF_RESTITUTION should be present in the properties when using DEMDiscontinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(COEFFICIENT_OF_RESTITUTION) = 0.0;
        }
    }

    std::string DEMDiscontinuumConstitutiveLaw::GetTypeOfLaw() {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::GetTypeOfLaw) should not be called.","")
        std::string type_of_law = "";
        return type_of_law;
    }

    DEMDiscontinuumConstitutiveLaw::Pointer DEMDiscontinuumConstitutiveLaw::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEMDiscontinuumConstitutiveLaw(*this));
        return p_clone;
    }

    DEMDiscontinuumConstitutiveLaw::~DEMDiscontinuumConstitutiveLaw() {
        //KRATOS_INFO("DEM") << "Law destructor..." ;
    }

    void DEMDiscontinuumConstitutiveLaw::CalculateContactArea(double radius, double other_radius, double &calculation_area) {

        KRATOS_TRY
        double rmin = radius;
        if (other_radius < radius) rmin = other_radius;
        calculation_area = Globals::Pi * rmin*rmin;
        KRATOS_CATCH("")
    }

    void DEMDiscontinuumConstitutiveLaw::CalculateElasticConstants(double &kn_el,
            double &kt_el,
            double initial_dist,
            double equiv_young,
            double equiv_poisson,
            double calculation_area,
            SphericParticle* element1,
            SphericParticle* element2) {

        KRATOS_TRY
        double equiv_shear = equiv_young / (2.0 * (1.0 + equiv_poisson));
        kn_el = equiv_young * calculation_area / initial_dist;
        kt_el = equiv_shear * calculation_area / initial_dist;
        KRATOS_CATCH("")
    }


    double DEMDiscontinuumConstitutiveLaw::LocalPeriod(const int i, SphericParticle* element1,
                                               SphericParticle* element2) {

        double myYoung = element1->GetYoung();
        double other_young = element2->GetYoung();
        double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);

        const double my_radius = element1->GetRadius();
        const double other_radius = element2->GetRadius();
        double calculation_area = 0;
        CalculateContactArea(my_radius, other_radius, calculation_area);

        double radius_sum = my_radius + other_radius;
        const double radius_sum_inv  = 1.0 / radius_sum;
        const double equiv_radius    = my_radius * other_radius * radius_sum_inv;

        const double modified_radius = equiv_radius * 0.31225; // sqrt(alpha * (2.0 - alpha)) = 0.31225
        double kn = equiv_young * Globals::Pi * modified_radius;        // 2.0 * equiv_young * sqrt_equiv_radius;

        const double mRealMass = element1->GetMass();
        const double other_real_mass = element2->GetMass();
        double equiv_mass = (mRealMass*other_real_mass)/(mRealMass+other_real_mass);

        // calculation of damping gamma
        const double my_gamma    = element1->GetProperties()[DAMPING_GAMMA];
        const double other_gamma = element2->GetProperties()[DAMPING_GAMMA];
        const double friction_coeff = element1->GetProperties()[FRICTION];
        const double equiv_gamma = 0.5 * (my_gamma + other_gamma);
        const double viscous_damping_coeff     = 2.0 * equiv_gamma * sqrt(equiv_mass * kn);
        double rescaled_damping = viscous_damping_coeff/(2*equiv_mass);
        double sqr_period = sqrt(1+friction_coeff*friction_coeff) * kn / equiv_mass - rescaled_damping*rescaled_damping;
        return sqr_period;
    }

    void DEMDiscontinuumConstitutiveLaw::InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double ini_delta) {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::InitializeContact) should not be called.","")
    }

    void DEMDiscontinuumConstitutiveLaw::InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta) {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::InitializeContactWithFEM) should not be called.","")
    }

    void DEMDiscontinuumConstitutiveLaw::GetContactStiffness(SphericParticle* const element1, SphericParticle* const element2, const double ini_delta, double& kn,double& kt){

      InitializeContact(element1, element2, ini_delta);
      kn = mKn;
      kt = mKt;

    }

    void DEMDiscontinuumConstitutiveLaw::CalculateForces(const ProcessInfo& r_process_info,
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
                                                         bool& sliding, double LocalCoordSystem[3][3]) {

        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::CalculateForces) should not be called.","")
    }

    void DEMDiscontinuumConstitutiveLaw::CalculateElasticEnergy(double& normal_elastic_energy,
                                                                double indentation,
                                                                double& cohesive_force,
                                                                SphericParticle* element1,
                                                                SphericParticle* element2) {


    }

    void DEMDiscontinuumConstitutiveLaw::CalculateForcesWithFEM(ProcessInfo& r_process_info,
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

        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::CalculateForcesWithFEM) should not be called.","")

    }

    double DEMDiscontinuumConstitutiveLaw::CalculateNormalForce(const double indentation) {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::CalculateNormalForce) should not be called.","")
    }


    double DEMDiscontinuumConstitutiveLaw::CalculateNormalForce(SphericParticle* const element1, SphericParticle* const element2, const double indentation,
        double LocalCoordSystem[3][3]) {
        return CalculateNormalForce(indentation);
    }

    double DEMDiscontinuumConstitutiveLaw::CalculateNormalForce(SphericParticle* const element, Condition* const wall, const double indentation){
        return CalculateNormalForce(indentation);
    }

    double DEMDiscontinuumConstitutiveLaw::CalculateCohesiveNormalForce(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::CalculateCohesiveNormalForce) should not be called.","")
        return 0.0;
    }

    double DEMDiscontinuumConstitutiveLaw::CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, Condition* const wall, const double indentation){
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::CalculateCohesiveNormalForceWithFEM) should not be called.","")
        return 0.0;
    }

    void DEMDiscontinuumConstitutiveLaw::CalculateViscoDamping(double LocalRelVel[3],
            double ViscoDampingLocalContactForce[3],
            double indentation,
            double equiv_visco_damp_coeff_normal,
            double equiv_visco_damp_coeff_tangential,
            bool sliding) {

        //*** component-wise since localContactForce and RelVel have in principle no relationship.
        // The visco force can be higher than the contact force only if they go to the same direction. (in my opinion)
        // But in opposite direction the visco damping can't overpass the force...

        KRATOS_TRY

        ViscoDampingLocalContactForce[2] = -equiv_visco_damp_coeff_normal * LocalRelVel[2];

        if (sliding == false) { //only applied when no sliding to help to the regularized friction law or the spring convergence
            ViscoDampingLocalContactForce[0] = -equiv_visco_damp_coeff_tangential * LocalRelVel[0];
            ViscoDampingLocalContactForce[1] = -equiv_visco_damp_coeff_tangential * LocalRelVel[1];
        }

        KRATOS_CATCH("")
    }

    void DEMDiscontinuumConstitutiveLaw::CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
            double &equiv_visco_damp_coeff_tangential,
            SphericParticle* element1,
            SphericParticle* element2,
            double kn_el,
            double kt_el) {

        KRATOS_TRY
        double aux_norm_to_tang = 0.0;
        const double my_mass = element1->GetMass();
        const double &other_real_mass = element2->GetMass();
        const double mCoefficientOfRestitution = element1->GetProperties()[COEFFICIENT_OF_RESTITUTION];

        equiv_visco_damp_coeff_normal = (1-mCoefficientOfRestitution) * 2.0 * sqrt(kn_el / (my_mass + other_real_mass)) * (sqrt(my_mass * other_real_mass)); // := 2d0* sqrt ( kn_el*(m1*m2)/(m1+m2) )
        equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang; // Dempack no ho fa servir...
        KRATOS_CATCH("")
    }

} // KRATOS
