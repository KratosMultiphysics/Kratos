/////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: chengshun.shang1996@gmail.com
// Date: Aug 2022
/////////////////////////////////////////////////

#if !defined(DEM_ROLLING_FRICTION_MODEL_CONSTANT_TORQUE_H_INCLUDED)
#define DEM_ROLLING_FRICTION_MODEL_CONSTANT_TORQUE_H_INCLUDED

/* Project includes */
#include "../custom_elements/spheric_particle.h"
#include "DEM_rolling_friction_model.h"


namespace Kratos{

    class KRATOS_API(DEM_APPLICATION) DEMRollingFrictionModelConstantTorque : public DEMRollingFrictionModel{
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEMRollingFrictionModelConstantTorque);

        DEMRollingFrictionModelConstantTorque() {}

        void Check(Properties::Pointer pProp) const override;

        ~DEMRollingFrictionModelConstantTorque() {}

        virtual DEMRollingFrictionModel::Pointer Clone() const override;

        virtual std::unique_ptr<DEMRollingFrictionModel> CloneUnique() override;

        virtual void ComputeRollingFriction(SphericParticle* p_element, SphericParticle* p_neighbor, double LocalContactForce[3], array_1d<double, 3>& mContactMoment, double indentation) override;
        
        virtual void ComputeRollingFrictionWithWall(double LocalContactForce[3], SphericParticle* p_element, Condition* const wall, double indentation, array_1d<double, 3>& mContactMoment) override;
    
    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override {
                    //rSerializer.load("MyMemberName",myMember);
        }

    };

    
} /* namespace Kratos.*/

#endif /*DEM_ROLLING_FRICTION_MODEL_CONSTANT_TORQUE_H_INCLUDED defined*/