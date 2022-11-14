//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(DEM_D_LINEAR_ZHAO_CL_H_INCLUDED)
#define DEM_D_LINEAR_ZHAO_CL_H_INCLUDED

#include <string>
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos {

  class SphericParticle;
  class KRATOS_API(DEM_APPLICATION) DEM_D_Linear_Zhao : public DEMDiscontinuumConstitutiveLaw {
  
  public:
    
    using DEMDiscontinuumConstitutiveLaw::CalculateNormalForce;
    
    KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Linear_Zhao);
    
    DEM_D_Linear_Zhao() {}
    ~DEM_D_Linear_Zhao() {}
    
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

    template <class NeighbourClassType>
    void CalculateTangentialForceWithNeighbour(const double normal_contact_force,
                                               const double OldLocalElasticContactForce[3],
                                               double LocalElasticContactForce[3],
                                               const double LocalDeltDisp[3],
                                               bool& sliding,
                                               SphericParticle* const element,
                                               NeighbourClassType* const neighbour);

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
#endif
