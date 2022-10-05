/////////////////////////////////////////////////
// Main author: Rafael Rangel (CIMNE)
// Email: rrangel@cimne.upc.edu
// Date: Oct 2022
/////////////////////////////////////////////////

#if !defined(DEM_ROLLING_FRICTION_MODEL_VISCOUS_TORQUE_H_INCLUDED)
#define DEM_ROLLING_FRICTION_MODEL_VISCOUS_TORQUE_H_INCLUDED

/* Project includes */
#include "../custom_elements/spheric_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "DEM_rolling_friction_model.h"

namespace Kratos
{
  class KRATOS_API(DEM_APPLICATION) DEMRollingFrictionModelViscousTorque : public DEMRollingFrictionModel
  {
  public:
    // Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(DEMRollingFrictionModelViscousTorque);

    // Clone
    virtual DEMRollingFrictionModel::Pointer Clone() const override;
    virtual std::unique_ptr<DEMRollingFrictionModel> CloneUnique() override;

    // Constructor / Destructor
    DEMRollingFrictionModelViscousTorque() {}
    ~DEMRollingFrictionModelViscousTorque() {}

    // Public methods
    void Check(Properties::Pointer pProp) const override;
    virtual void ComputeRollingFriction(SphericParticle* p_element, SphericParticle* p_neighbor, double LocalCoordSystem_2[3], double LocalContactForce[3], array_1d<double, 3>& mContactMoment)   override;
    virtual void ComputeRollingFrictionWithWall(SphericParticle* p_element, Condition* const wall, double LocalCoordSystem_2[3], double LocalContactForce[3], array_1d<double, 3>& mContactMoment) override;

  private:
    friend class Serializer;
    virtual void save(Serializer & rSerializer) const override {}
    virtual void load(Serializer & rSerializer) override {}
  };
} /* namespace Kratos */

#endif /* DEM_ROLLING_FRICTION_MODEL_VISCOUS_TORQUE_H_INCLUDED defined */