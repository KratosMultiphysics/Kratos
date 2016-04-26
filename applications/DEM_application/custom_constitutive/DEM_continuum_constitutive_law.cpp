#include "DEM_application.h"
#include "DEM_continuum_constitutive_law.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::DEMContinuumConstitutiveLaw() {
        //std::cout << " DEMContinuumConstitutiveLaw constructor.." << std::endl;

    } // Class DEMContinuumConstitutiveLaw

    DEMContinuumConstitutiveLaw::DEMContinuumConstitutiveLaw(const DEMContinuumConstitutiveLaw &rReferenceContinuumConstitutiveLaw) {
        //std::cout << " DEMContinuumConstitutiveLaw copy constructor.." << std::endl;
    }
    
    std::string DEMContinuumConstitutiveLaw::GetTypeOfLaw() {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMContinuumConstitutiveLaw::GetTypeOfLaw) should not be called.","")
        std::string type_of_law = "";
        return type_of_law;
    }

    void DEMContinuumConstitutiveLaw::Initialize() {
    }

    void DEMContinuumConstitutiveLaw::SetConstitutiveLawInProperties(Properties::Pointer pProp) const {
        std::cout << "Assigning DEMContinuumConstitutiveLaw to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    DEMContinuumConstitutiveLaw::Pointer DEMContinuumConstitutiveLaw::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEMContinuumConstitutiveLaw(*this));
        return p_clone;
    }

    DEMContinuumConstitutiveLaw::~DEMContinuumConstitutiveLaw() { //std::cout << "Law destructor..." ;
    }

    void DEMContinuumConstitutiveLaw::CalculateViscoDamping(double LocalRelVel[3],
                                                            double ViscoDampingLocalContactForce[3],
                                                            double indentation,
                                                            double equiv_visco_damp_coeff_normal,
                                                            double equiv_visco_damp_coeff_tangential,
                                                            bool& sliding,
                                                            int failure_id) {

        //*** component-wise since localContactForce and RelVel have in principle no relationship.
        // The visco force can be higher than the contact force only if they go to the same direction. (in my opinion)
        // But in opposite direction the visco damping can't overpass the force...

        KRATOS_TRY  

        ViscoDampingLocalContactForce[2] = -equiv_visco_damp_coeff_normal * LocalRelVel[2];

        if (sliding == false) { //only applied when no sliding to help the regularized friction law or the spring convergence
            ViscoDampingLocalContactForce[0] = -equiv_visco_damp_coeff_tangential * LocalRelVel[0];
            ViscoDampingLocalContactForce[1] = -equiv_visco_damp_coeff_tangential * LocalRelVel[1];
        }

        KRATOS_CATCH("")
    }
    
    void DEMContinuumConstitutiveLaw::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                                       SphericContinuumParticle* neighbor,
                                                                       double equiv_young,
                                                                       double distance,
                                                                       double calculation_area,
                                                                       double LocalCoordSystem[3][3],
                                                                       double ElasticLocalRotationalMoment[3],
                                                                       double ViscoLocalRotationalMoment[3]) {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMContinuumConstitutiveLaw::ComputeParticleRotationalMoments) should not be called.","")
    }
    
    void DEMContinuumConstitutiveLaw::AddPoissonContribution(const double equiv_poisson, double LocalCoordSystem[3][3], double& normal_force, double calculation_area, Matrix* mSymmStressTensor, 
                                                             SphericContinuumParticle* element1, SphericContinuumParticle* element2){    
    }


    double DEMContinuumConstitutiveLaw::LocalMaxSearchDistance(const int i,
                                          SphericContinuumParticle* element1,
                                          SphericContinuumParticle* element2) {

        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMContinuumConstitutiveLaw::LocalMaxSearchDistance) should not be called.","")

    }

    double DEMContinuumConstitutiveLaw::LocalPeriod(const int i, SphericContinuumParticle* element1,
                                                                 SphericContinuumParticle* element2) {

        // calculation of equivalent young modulus
        double myYoung = element1->GetYoung();
        double other_young = element2->GetYoung();
        double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);

        const double my_radius = element1->GetRadius();
        const double other_radius = element2->GetRadius();
        double calculation_area = 0;
        CalculateContactArea(my_radius, other_radius, calculation_area);

        double radius_sum = my_radius + other_radius;
        double initial_delta = element1->GetInitialDelta(i);
        double initial_dist = radius_sum - initial_delta;

        // calculation of elastic constants
        double kn_el = equiv_young * calculation_area / initial_dist;

        //mDensity = element1_props[PARTICLE_DENSITY];
        //other_density = element2_props[PARTICLE_DENSITY];
        //double m1 = 4/3 * KRATOS_M_PI * my_radius * my_radius * my_radius * mDensity;
        //double m2 = 4/3 * KRATOS_M_PI * other_radius * other_radius * other_radius * other_density;

        const double mRealMass = element1->GetMass();  // { mRealMass = real_mass;  GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = real_mass;}
        const double &other_real_mass = element2->GetMass();
        const double mCoefficientOfRestitution = element1->GetProperties()[COEFFICIENT_OF_RESTITUTION];
        const double mOtherCoefficientOfRestitution = element2->GetProperties()[COEFFICIENT_OF_RESTITUTION];
        const double equiv_coefficientOfRestitution = 0.5 * (mCoefficientOfRestitution + mOtherCoefficientOfRestitution);
        // calculation of damping gamma

        double equiv_mass = (mRealMass*other_real_mass)/(mRealMass+other_real_mass);
        double viscous_damping_coeff = (1-equiv_coefficientOfRestitution) * 2.0 * sqrt(kn_el * equiv_mass);

        double rescaled_damping = viscous_damping_coeff/(2*equiv_mass);

        double sqr_period = kn_el / equiv_mass - rescaled_damping*rescaled_damping;
        return sqr_period;

    }



} //kratos
