// Author: Miguel Angel Celigueta maceli@cimne.upc.edu

#include "includes/define.h"
#include <boost/python.hpp>
#include "includes/model_part.h"

namespace Kratos {
    
    class AuxiliaryUtilities {

    public:
        AuxiliaryUtilities() {};
        virtual ~AuxiliaryUtilities() {};

        boost::python::object GetIthSubModelPartIsForceIntegrationGroup(ModelPart& rParentModelPart, const int& required_i);
        boost::python::object GetIthSubModelPartName(ModelPart& rParentModelPart, const int& required_i);
        boost::python::object GetIthSubModelPartIdentifier(ModelPart& rParentModelPart, const int& required_i);        
        template<class T> boost::python::object GetIthSubModelPartData(ModelPart& rParentModelPart, const int& required_i, const Variable<T>& rVariable);
        ModelPart::NodesContainerType::Pointer GetIthSubModelPartNodes(ModelPart& rParentModelPart, const int& required_i);
        
    };
}
