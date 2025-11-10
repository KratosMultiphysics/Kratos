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

#pragma once

#include <string>
#include "../custom_elements/spheric_particle.h"
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos {

  class KRATOS_API(DEM_APPLICATION) DEM_D_Linear_Simple_Coulomb : public DEMDiscontinuumConstitutiveLaw {

  public:

    using DEMDiscontinuumConstitutiveLaw::CalculateNormalForce;

    KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Linear_Simple_Coulomb);

    DEM_D_Linear_Simple_Coulomb() {}
    ~DEM_D_Linear_Simple_Coulomb() {}

    DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;
    std::unique_ptr<DEMDiscontinuumConstitutiveLaw> CloneUnique() override;

    std::string GetTypeOfLaw() override;
    virtual void Check(Properties::Pointer pProp) const override;

    void CalculateForces(const ProcessInfo& r_process_info,
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
                         double LocalCoordSystem[3][3]) override;

    void CalculateForcesWithFEM(const ProcessInfo& r_process_info,
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
                                bool& sliding) override;

  private:

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
    }

    virtual void load(Serializer& rSerializer) override {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
    }

  };
}
