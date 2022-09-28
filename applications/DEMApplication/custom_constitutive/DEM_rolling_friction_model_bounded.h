/////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
//              Joaqu¨ªn Iraz¨¢bal
// Email: chengshun.shang1996@gmail.com
// Date: Sep 2022
/////////////////////////////////////////////////

#if !defined(DEM_ROLLING_FRICTION_MODEL_BOUNDED_H_INCLUDED)
#define DEM_ROLLING_FRICTION_MODEL_BOUNDED_H_INCLUDED

/* Project includes */
#include "../custom_elements/spheric_particle.h"
#include "DEM_rolling_friction_model.h"
#include "includes/kratos_flags.h"


namespace Kratos{

    class KRATOS_API(DEM_APPLICATION) DEMRollingFrictionModelBounded : public DEMRollingFrictionModel{
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEMRollingFrictionModelBounded);

        DEMRollingFrictionModelBounded() {}

        double mRollingResistance = 0.0;

        void Check(Properties::Pointer pProp) const override;

        bool CheckIfThisModelRequiresRecloningForEachNeighbour() override;

        ~DEMRollingFrictionModelBounded() {}

        virtual DEMRollingFrictionModel::Pointer Clone() const override;

        virtual std::unique_ptr<DEMRollingFrictionModel> CloneUnique() override;

        virtual void InitializeSolutionStep() override;

        virtual void ComputeRollingResistance(const double& NormalLocalContactForce, const double& equiv_rolling_friction_coeff, const unsigned int i) override;

        virtual void DoFinalOperations(SphericParticle* p_element, double dt, array_1d<double, 3>& mContactMoment) override;
    
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

#endif /*DEM_ROLLING_FRICTION_MODEL_BOUNDED_H_INCLUDED defined*/