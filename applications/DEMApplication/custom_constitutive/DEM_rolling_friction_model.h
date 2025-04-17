//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Chengshun Shang (cshang@cimne.upc.edu)
//

#if !defined(DEM_ROLLING_FRICTION_MODEL_H_INCLUDED)
#define DEM_ROLLING_FRICTION_MODEL_H_INCLUDED

// System includes
#include <string>
#include <iostream>

/* Project includes */
#include "includes/define.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "includes/serializer.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_elements/discrete_element.h"
#include "custom_elements/Particle_Contact_Element.h"
#include "containers/array_1d.h"
#include "custom_conditions/RigidFace.h"
#include "custom_conditions/dem_wall.h"

namespace Kratos{

    class Properties;
    class SphericParticle; // forward declaration of spheric cont particle

    class KRATOS_API(DEM_APPLICATION) DEMRollingFrictionModel {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEMRollingFrictionModel);

        DEMRollingFrictionModel();

        virtual bool CheckIfThisModelRequiresRecloningForEachNeighbour() {
            return true;
        }

        virtual void SetAPrototypeOfThisInProperties(Properties::Pointer pProp, bool verbose = true);

        virtual void Check(Properties::Pointer pProp) const;

        virtual ~DEMRollingFrictionModel();

        virtual DEMRollingFrictionModel::Pointer Clone() const;

        virtual std::unique_ptr<DEMRollingFrictionModel> CloneUnique();

        virtual void InitializeSolutionStep() {}

        virtual void InitializeContact(SphericParticle* const p_element, SphericParticle* const p_neighbor, const double indentation) {}

        virtual void ComputeRollingFriction(SphericParticle* p_element, SphericParticle* p_neighbor, const ProcessInfo& r_process_info, double LocalContactForce[3], double indentation, array_1d<double, 3>& mContactMoment, double LocalCoordSystem2[3]) {}
        
        virtual void ComputeRollingFrictionWithWall(SphericParticle* p_element, Condition* const wall, const ProcessInfo& r_process_info, double LocalContactForce[3], double indentation, array_1d<double, 3>& mContactMoment, double LocalCoordSystem2[3]) {}

        virtual void ComputeRollingResistance(SphericParticle* p_element, SphericParticle* p_neighbor, double LocalContactForce[3]) {}

        virtual void ComputeRollingResistanceWithWall(SphericParticle* p_element, Condition* const wall, double LocalContactForce[3]) {}

        virtual void DoFinalOperations(SphericParticle* p_element, double dt, array_1d<double, 3>& mContactMoment) {}

        virtual void CalculateInelasticRollingResistanceEnergy(double& inelastic_rollingresistance_energy, const array_1d<double, 3>& rolling_friction_moment, const array_1d<double, 3>& relative_angular_velocity, double dt) {}

        virtual void CalculateInelasticRollingResistanceEnergyWithWall(double& inelastic_rollingresistance_energy, const array_1d<double, 3>& rolling_friction_moment, const array_1d<double, 3>& relative_angular_velocity, double dt) {}
    
    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const {           
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) {
                    //rSerializer.load("MyMemberName",myMember);
        }
    };

    //This definition is done here to avoid recursive inclusion of header files
    KRATOS_DEFINE_APPLICATION_VARIABLE(DEM_APPLICATION, DEMRollingFrictionModel::Pointer, DEM_ROLLING_FRICTION_MODEL_POINTER)

} /* namespace Kratos.*/

#endif /*DEM_ROLLING_FRICTION_MODEL_H_INCLUDED defined*/