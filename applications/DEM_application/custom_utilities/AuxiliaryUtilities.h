#pragma once
// Author: Miguel Angel Celigueta maceli@cimne.upc.edu

#include <pybind11/pybind11.h>
#include "includes/define.h"
#include "includes/model_part.h"


namespace Kratos {
    
    class AuxiliaryUtilities {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(AuxiliaryUtilities);

        AuxiliaryUtilities() {};
        virtual ~AuxiliaryUtilities() {};

        pybind11::object GetIthSubModelPartIsForceIntegrationGroup(ModelPart& rParentModelPart, const int& required_i);
        pybind11::object GetIthSubModelPartName(ModelPart& rParentModelPart, const int& required_i);
        pybind11::object GetIthSubModelPartIdentifier(ModelPart& rParentModelPart, const int& required_i);        
        template<class T> pybind11::object GetIthSubModelPartData(ModelPart& rParentModelPart, const int& required_i, const Variable<T>& rVariable);
        ModelPart::NodesContainerType::Pointer GetIthSubModelPartNodes(ModelPart& rParentModelPart, const int& required_i);
        
    };
}
