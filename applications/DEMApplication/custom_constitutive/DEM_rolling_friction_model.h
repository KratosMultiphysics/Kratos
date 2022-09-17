/////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: chengshun.shang1996@gmail.com
// Date: Aug 2022
/////////////////////////////////////////////////

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

        virtual void ComputeRollingFriction(SphericParticle* p_element, SphericParticle* p_neighbor, double LocalContactForce[3], array_1d<double, 3>& mContactMoment, double indentation);
        
        virtual void ComputeRollingFrictionWithWall(double LocalContactForce[3], SphericParticle* p_element, Condition* const wall, double indentation, array_1d<double, 3>& mContactMoment);

        virtual void InitializeSolutionStep() {}

        virtual void ComputeRollingResistance(const double& NormalLocalContactForce, const double& equiv_rolling_friction_coeff, const unsigned int i) {}

        virtual void DoFinalOperations(SphericParticle* p_element, double dt, array_1d<double, 3>& mContactMoment) {}
    
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