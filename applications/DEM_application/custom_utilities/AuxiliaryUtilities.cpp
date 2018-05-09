//
// Author: 
// Miguel Angel Celigueta maceli@cimne.upc.edu
//

#include "AuxiliaryUtilities.h"
#include "DEM_application_variables.h"

// System includes
#include <string>
#include <iostream>

namespace Kratos {
    pybind11::object AuxiliaryUtilities::GetIthSubModelPartIsForceIntegrationGroup(ModelPart& rParentModelPart, const int& required_i){
        KRATOS_TRY;                             
        int current_i = 0;
        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = rParentModelPart.SubModelPartsBegin(); sub_model_part != rParentModelPart.SubModelPartsEnd(); ++sub_model_part) {

            if(current_i == required_i) {
                const int is_force_integration_group = (*sub_model_part)[FORCE_INTEGRATION_GROUP];
                return pybind11::cast( (bool) (is_force_integration_group) );
            }
            current_i++;
        }

        KRATOS_THROW_ERROR(std::runtime_error, "The function GetIthSubModelPartIsForceIntegrationGroup required a position greater than the number of SubModelParts! The number of present submodelparts is ", rParentModelPart.NumberOfSubModelParts());
        
        KRATOS_CATCH("");        
    };
    
    pybind11::object AuxiliaryUtilities::GetIthSubModelPartName(ModelPart& rParentModelPart, const int& required_i){
        KRATOS_TRY;                             
        int current_i = 0;
        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = rParentModelPart.SubModelPartsBegin(); sub_model_part != rParentModelPart.SubModelPartsEnd(); ++sub_model_part) {

            if(current_i == required_i) {
                const std::string& name = (*sub_model_part).Name();
                return pybind11::cast(name);
            }
            current_i++;
        }

        KRATOS_THROW_ERROR(std::runtime_error, "The function GetIthSubModelPartName required a position greater than the number of SubModelParts! The number of present submodelparts is ", rParentModelPart.NumberOfSubModelParts());
        
        KRATOS_CATCH("");  
    };
    
    pybind11::object AuxiliaryUtilities::GetIthSubModelPartIdentifier(ModelPart& rParentModelPart, const int& required_i){
        KRATOS_TRY; 
        
        int current_i = 0;
        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = rParentModelPart.SubModelPartsBegin(); sub_model_part != rParentModelPart.SubModelPartsEnd(); ++sub_model_part) {

            if(current_i == required_i) {
                const std::string& identifier = (*sub_model_part)[IDENTIFIER];
                return pybind11::cast(identifier);
            }
            
            current_i++;
        }

        KRATOS_THROW_ERROR(std::runtime_error, "The function GetIthSubModelPartIdentifier required a position greater than the number of SubModelParts! The number of present submodelparts is ", rParentModelPart.NumberOfSubModelParts());
        
        KRATOS_CATCH("");  
    };
    
    template<class T>
    pybind11::object AuxiliaryUtilities::GetIthSubModelPartData(ModelPart& rParentModelPart, const int& required_i, const Variable<T>& rVariable){
        KRATOS_TRY; 
        
        int current_i = 0;
        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = rParentModelPart.SubModelPartsBegin(); sub_model_part != rParentModelPart.SubModelPartsEnd(); ++sub_model_part) {

            if(current_i == required_i) {
                const T& value = (*sub_model_part)[rVariable];
                return pybind11::cast(value);
            }
            
            current_i++;
        }

        KRATOS_THROW_ERROR(std::runtime_error, "The function GetIthSubModelPartData required a position greater than the number of SubModelParts! The number of present submodelparts is ", rParentModelPart.NumberOfSubModelParts());
        
        KRATOS_CATCH("");  
    }      
    
    // Explicit instantiations of the previous template function (otherwise it does not link):
    template pybind11::object AuxiliaryUtilities::GetIthSubModelPartData<int>(ModelPart& rParentModelPart, const int& required_i, const Variable<int>& rVariable);
    template pybind11::object AuxiliaryUtilities::GetIthSubModelPartData<double>(ModelPart& rParentModelPart, const int& required_i, const Variable<double>& rVariable);
    template pybind11::object AuxiliaryUtilities::GetIthSubModelPartData<array_1d<double,3> >(ModelPart& rParentModelPart, const int& required_i, const Variable<array_1d<double,3> >& rVariable);
    template pybind11::object AuxiliaryUtilities::GetIthSubModelPartData<std::string>(ModelPart& rParentModelPart, const int& required_i, const Variable<std::string>& rVariable);
       
    
    ModelPart::NodesContainerType::Pointer AuxiliaryUtilities::GetIthSubModelPartNodes(ModelPart& rParentModelPart, const int& required_i){
        KRATOS_TRY; 
        
        int current_i = 0;
        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = rParentModelPart.SubModelPartsBegin(); sub_model_part != rParentModelPart.SubModelPartsEnd(); ++sub_model_part) {

            if(current_i == required_i) {
                return (*sub_model_part).pNodes();
            }
            
            current_i++;
        }

        KRATOS_THROW_ERROR(std::runtime_error, "The function GetIthSubModelPartNodes required a position greater than the number of SubModelParts! The number of present submodelparts is ", rParentModelPart.NumberOfSubModelParts());
        
        KRATOS_CATCH("");  
        
    }       
    
} // Namespace Kratos
