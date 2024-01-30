//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Chengshun Shang (cshang@cimne.upc.edu)
//                 Joaquin Irazabal (jirazabal@cimne.upc.edu)
//

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

        DEMRollingFrictionModel::Pointer Clone() const override;

        std::unique_ptr<DEMRollingFrictionModel> CloneUnique() override;

        void InitializeSolutionStep() override;

        void ComputeRollingResistance(SphericParticle* p_element, SphericParticle* p_neighbor, double LocalContactForce[3]) override;

        void ComputeRollingResistanceWithWall(SphericParticle* p_element, Condition* const wall, double LocalContactForce[3]) override;

        void DoFinalOperations(SphericParticle* p_element, double dt, array_1d<double, 3>& mContactMoment) override;
    
    private:

        friend class Serializer;

        void save(Serializer& rSerializer) const override {
                    //rSerializer.save("MyMemberName",myMember);
        }

        void load(Serializer& rSerializer) override {
                    //rSerializer.load("MyMemberName",myMember);
        }

    };

    
} /* namespace Kratos.*/

#endif /*DEM_ROLLING_FRICTION_MODEL_BOUNDED_H_INCLUDED defined*/