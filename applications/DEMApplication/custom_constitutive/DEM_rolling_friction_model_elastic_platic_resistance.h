//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Chengshun Shang (cshang@cimne.upc.edu)
//

#if !defined(DEM_ROLLING_FRICTION_MODEL_ELASTIC_PLASTIC_RESISTANCE_H_INCLUDED)
#define DEM_ROLLING_FRICTION_MODEL_ELASTIC_PLASTIC_RESISTANCE_H_INCLUDED

/* Project includes */
#include "../custom_elements/spheric_particle.h"
#include "DEM_rolling_friction_model.h"


namespace Kratos{

    class KRATOS_API(DEM_APPLICATION) DEMRollingFrictionModelElasticPlasticResistance : public DEMRollingFrictionModel{
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEMRollingFrictionModelElasticPlasticResistance);

        DEMRollingFrictionModelElasticPlasticResistance() {}

        void Check(Properties::Pointer pProp) const override;

        ~DEMRollingFrictionModelElasticPlasticResistance() {}

        DEMRollingFrictionModel::Pointer Clone() const override;

        std::unique_ptr<DEMRollingFrictionModel> CloneUnique() override;

        void InitializeContact(SphericParticle* const p_element, SphericParticle* const p_neighbor, const double indentation) override;

        virtual void InitializeContactWall(SphericParticle* const element, Condition* const wall, const double indentation);

        void ComputeRollingFriction(SphericParticle* p_element, SphericParticle* p_neighbor, const ProcessInfo& r_process_info, double LocalContactForce[3], double indentation, array_1d<double, 3>& mContactMoment) override;
        
        void ComputeRollingFrictionWithWall(SphericParticle* p_element, Condition* const wall, const ProcessInfo& r_process_info, double LocalContactForce[3], double indentation, array_1d<double, 3>& mContactMoment, double LocalCoordSystem2[3]) override;

        void CalculateInelasticRollingResistanceEnergy(double& inelastic_rollingresistance_energy, const array_1d<double, 3>& rolling_friction_moment, const array_1d<double, 3>& relative_angular_velocity, double dt) override;

        void CalculateInelasticRollingResistanceEnergyWithWall(double& inelastic_rollingresistance_energy, const array_1d<double, 3>& rolling_friction_moment, const array_1d<double, 3>& relative_angular_velocity, double dt) override;
    
        double m_rolling_friction_moment[3] = {0.0};
        double mKt = 0.0;
    
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

#endif /*DEM_ROLLING_FRICTION_MODEL_ELASTIC_PLASTIC_RESISTANCE_H_INCLUDED defined*/