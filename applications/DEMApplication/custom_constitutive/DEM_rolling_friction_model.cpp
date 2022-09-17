/////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: chengshun.shang1996@gmail.com
// Date: Aug 2022
/////////////////////////////////////////////////

#include "custom_elements/spheric_particle.h"
#include "DEM_rolling_friction_model.h"

namespace Kratos{

    DEMRollingFrictionModel::DEMRollingFrictionModel() {}

    void DEMRollingFrictionModel::SetAPrototypeOfThisInProperties(Properties::Pointer pProp, bool verbose) {
        if (verbose) KRATOS_INFO("DEM")  << "Assigning " << pProp->GetValue(DEM_ROLLING_FRICTION_MODEL_NAME) << " to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_ROLLING_FRICTION_MODEL_POINTER, this->Clone());
        this->Check(pProp);
    }

    std::unique_ptr<DEMRollingFrictionModel> DEMRollingFrictionModel::CloneUnique() {
        KRATOS_ERROR << "This function (DEMRollingFrictionModel::CloneUnique) shouldn't be accessed, use derived class instead"<<std::endl;
    }

    DEMRollingFrictionModel::Pointer DEMRollingFrictionModel::Clone() const {
        DEMRollingFrictionModel::Pointer p_clone(new DEMRollingFrictionModel(*this));
        return p_clone;
    }

    DEMRollingFrictionModel::~DEMRollingFrictionModel() {}

    void DEMRollingFrictionModel::Check(Properties::Pointer pProp) const {
        KRATOS_ERROR << "This function (DEMRollingFrictionModel::Check) shouldn't be accessed, use derived class instead"<<std::endl;
    }

    void DEMRollingFrictionModel::ComputeRollingFriction(SphericParticle* p_element, SphericParticle* p_neighbor, double LocalContactForce[3], array_1d<double, 3>& mContactMoment, double indentation){
        //KRATOS_ERROR << "This function (DEMRollingFrictionModel::ComputeRollingFriction) shouldn't be accessed, use derived class instead"<<std::endl;
    }
        
    void DEMRollingFrictionModel::ComputeRollingFrictionWithWall(double LocalContactForce[3], SphericParticle* p_element, Condition* const wall, double indentation, array_1d<double, 3>& mContactMoment){
        //KRATOS_ERROR << "This function (DEMRollingFrictionModel::ComputeRollingFrictionWithWall) shouldn't be accessed, use derived class instead"<<std::endl;
    }

}/*namespace Kratos*/