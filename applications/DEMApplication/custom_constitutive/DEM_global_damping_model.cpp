//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

// Includes
#include "DEM_global_damping_model.h"

namespace Kratos
{
    void DEMGlobalDampingModel::SetGlobalDampingModelInProperties(Properties::Pointer pProp, bool verbose) {
        if (verbose) KRATOS_INFO("DEM")  << "Assigning " << pProp->GetValue(DEM_GLOBAL_DAMPING_MODEL_NAME) << " to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_GLOBAL_DAMPING_MODEL_POINTER, this->Clone());
    }

    DEMGlobalDampingModel::Pointer DEMGlobalDampingModel::Clone() const {
        DEMGlobalDampingModel::Pointer p_clone(new DEMGlobalDampingModel(*this));
        return p_clone;
    }

    std::unique_ptr<DEMGlobalDampingModel> DEMGlobalDampingModel::CloneUnique() {
        KRATOS_ERROR << "This function (DEMGlobalDampingModel::CloneUnique) shouldn't be accessed, use derived class instead"<<std::endl;
    }
}
