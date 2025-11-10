//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//
// References using this model:
// R.L. Rangel et al. (2024). Multiscale data-driven modeling of the thermomechanical behavior of granular media with thermal expansion effects. Computers and Geotechnics, 176:106789.
// R.L. Rangel et al. (2024). A continuum-discrete multiscale methodology using machine learning for thermal analysis of granular media. Computers and Geotechnics, 168:106118.
// Guo & Zhao (2014). coupled FEM/DEM approach for hierarchical multiscale modelling of granular media. IJNME, 99(11):789-818.
//

#include "DEM_D_Linear_Simple_Coulomb_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

  DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Linear_Simple_Coulomb::Clone() const {
    DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Linear_Simple_Coulomb(*this));
    return p_clone;
  }

  std::unique_ptr<DEMDiscontinuumConstitutiveLaw> DEM_D_Linear_Simple_Coulomb::CloneUnique() {
    return Kratos::make_unique<DEM_D_Linear_Simple_Coulomb>();
  }

  std::string DEM_D_Linear_Simple_Coulomb::GetTypeOfLaw() {
    std::string type_of_law = "Linear";
    return type_of_law;
  }

  void DEM_D_Linear_Simple_Coulomb::Check(Properties::Pointer pProp) const {
    if (!pProp->Has(STATIC_FRICTION)) {
      if (!pProp->Has(FRICTION)) { //deprecated since April 6th, 2020
        KRATOS_WARNING("DEM") << std::endl;
        KRATOS_WARNING("DEM") << "WARNING: Variable STATIC_FRICTION or FRICTION should be present in the properties when using DEMDiscontinuumConstitutiveLaw. 0.0 value assigned by default." << std::endl;
        KRATOS_WARNING("DEM") << std::endl;
        pProp->GetValue(STATIC_FRICTION) = 0.0;
      }
      else {
        pProp->GetValue(STATIC_FRICTION) = pProp->GetValue(FRICTION);
      }
    }
  }

  void DEM_D_Linear_Simple_Coulomb::CalculateForces(const ProcessInfo& r_process_info,
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
    // Get properties of the two particles
    const double r1 = element1->GetRadius();
    const double r2 = element2->GetRadius();
    const double reff = r1 * r2 / (r1 + r2);
    const double m1 = element1->GetMass();
    const double m2 = element2->GetMass();
    const double meff = m1 * m2 / (m1 + m2);
    const double v1 = element1->GetPoisson();
    const double v2 = element2->GetPoisson();
    const double E1 = element1->GetYoung();
    const double E2 = element2->GetYoung();
    const double Eeff = 1.0 / ((1.0 - v1*v1) / E1 + (1.0 - v2*v2) / E2);
    const double G1 = 0.5 * E1 / (1.0 + v1);
    const double G2 = 0.5 * E2 / (1.0 + v2);
    const double Geff = 1.0 / ((2.0 - v1) / G1 + (2.0 - v2) / G2);
    Properties& contact_props = element1->GetProperties().GetSubProperties(element2->GetProperties().Id());
    const double phi = contact_props[DAMPING_GAMMA];
    const double friction_angle_tg = std::tan(contact_props[STATIC_FRICTION]);
    
    // Compute normal force (elastic)
    const double Kn = 2.0 * Eeff * sqrt(reff * indentation);
    const double Fne = (2.0/3.0) * Kn * indentation;

    // Compute normal force (viscous)
    double Fnv = -(2.0 * phi * sqrt(meff * Kn)) * LocalRelVel[2];

    // Check for artificial cohesion
    double Fn = Fne + Fnv;
    if (Fn < 0.0) {
        Fn = 0.0;
        Fnv = -Fne;
    }

    // Compute tangential force
    const double Kt = 4.0 * Geff * Kn / Eeff;
    double Ft_y = OldLocalElasticContactForce[0] - Kt * LocalDeltDisp[0];
    double Ft_z = OldLocalElasticContactForce[1] - Kt * LocalDeltDisp[1];
    double Ft = sqrt(Ft_y * Ft_y + Ft_z * Ft_z);

    // Check for Coulomb condition
    const double Ftmax = Fn * friction_angle_tg;
    if (Ft > Ftmax) {
        const double fraction = Ftmax / Ft;
        Ft_y *= fraction;
        Ft_z *= fraction;
        Ft = sqrt(Ft_y * Ft_y + Ft_z * Ft_z);
    }

    // Store forces in their respective arrays ([0]: tangent in local y, [1]: tangent in local z, [2]: normal)
    LocalElasticContactForce[0] = Ft_y;
    LocalElasticContactForce[1] = Ft_z;
    LocalElasticContactForce[2] = Fne;
    ViscoDampingLocalContactForce[0] = 0.0; // No tangential viscous force
    ViscoDampingLocalContactForce[1] = 0.0; // No tangential viscous force
    ViscoDampingLocalContactForce[2] = Fnv;

    // Calculate elastic energy (each particle in a contact with another particle receives half the contact energy)
    double& elastic_energy = element1->GetElasticEnergy();
    elastic_energy += 0.20 * Fne * indentation; // normal component
    elastic_energy += 0.25 * Ft * Ft / Kt; // tangential component
  }

  void DEM_D_Linear_Simple_Coulomb::CalculateForcesWithFEM(const ProcessInfo& r_process_info,
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
    // Get properties of particle and wall
    const double reff = element->GetRadius();
    const double meff = element->GetMass();
    const double v1 = element->GetPoisson();
    const double v2 = wall->GetProperties()[POISSON_RATIO];
    const double E1 = element->GetYoung();
    const double E2 = wall->GetProperties()[YOUNG_MODULUS];
    const double Eeff = 1.0 / ((1.0 - v1*v1) / E1 + (1.0 - v2*v2) / E2);
    const double G1 = 0.5 * E1 / (1.0 + v1);
    const double G2 = 0.5 * E2 / (1.0 + v2);
    const double Geff = 1.0 / ((2.0 - v1) / G1 + (2.0 - v2) / G2);
    Properties& properties_of_this_contact = element->GetProperties().GetSubProperties(wall->GetProperties().Id());
    const double phi = properties_of_this_contact[DAMPING_GAMMA];
    const double friction_angle_tg = std::tan(properties_of_this_contact[STATIC_FRICTION]);

    // Compute normal force (elastic)
    const double Kn = 2.0 * Eeff * sqrt(reff * indentation);
    const double Fne = (2.0/3.0) * Kn * indentation;

    // Compute normal force (viscous)
    double Fnv = -(2.0 * phi * sqrt(meff * Kn)) * LocalRelVel[2];

    // Check for artificial cohesion
    double Fn = Fne + Fnv;
    if (Fn < 0.0) {
        Fn = 0.0;
        Fnv = -Fne;
    }

    // Compute tangential force
    const double Kt = 4.0 * Geff * Kn / Eeff;
    double Ft_y = OldLocalElasticContactForce[0] - Kt * LocalDeltDisp[0];
    double Ft_z = OldLocalElasticContactForce[1] - Kt * LocalDeltDisp[1];
    double Ft = sqrt(Ft_y * Ft_y + Ft_z * Ft_z);

    // Check for Coulomb condition
    const double Ftmax = Fn * friction_angle_tg;
    if (Ft > Ftmax) {
        const double fraction = Ftmax / Ft;
        Ft_y *= fraction;
        Ft_z *= fraction;
        Ft = sqrt(Ft_y * Ft_y + Ft_z * Ft_z);
    }

    // Store forces in their respective arrays ([0]: tangent in local y, [1]: tangent in local z, [2]: normal)
    LocalElasticContactForce[0] = Ft_y;
    LocalElasticContactForce[1] = Ft_z;
    LocalElasticContactForce[2] = Fne;
    ViscoDampingLocalContactForce[0] = 0.0; // No tangential viscous force
    ViscoDampingLocalContactForce[1] = 0.0; // No tangential viscous force
    ViscoDampingLocalContactForce[2] = Fnv;

    // Calculate elastic energy (each particle in a contact with a wall receives all the contact energy)
    double& elastic_energy = element->GetElasticEnergy();
    elastic_energy += 0.40 * Fne * indentation; // normal component
    elastic_energy += 0.50 * Ft * Ft / Kt; // tangential component
    
  }
} // namespace Kratos
