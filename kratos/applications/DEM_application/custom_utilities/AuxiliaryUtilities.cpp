//
// Author: 
// Miguel Angel Celigueta maceli@cimne.upc.edu
//

// System includes
#include <string>
#include <iostream>

#include "AuxiliaryUtilities.h"
#include "DEM_application_variables.h"

namespace Kratos {
    boost::python::object AuxiliaryUtilities::GetIthSubModelPartIsForceIntegrationGroup(ModelPart& rParentModelPart, const int& required_i){
        KRATOS_TRY;                             
        int current_i = 0;
        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = rParentModelPart.SubModelPartsBegin(); sub_model_part != rParentModelPart.SubModelPartsEnd(); ++sub_model_part) {

            if(current_i == required_i) {
                const int is_force_integration_group = (*sub_model_part)[FORCE_INTEGRATION_GROUP];
                return boost::python::object( (bool) (is_force_integration_group) );
            }
            current_i++;
        }

        KRATOS_THROW_ERROR(std::runtime_error, "The function GetIthSubModelPartIsForceIntegrationGroup required a position greater than the number of SubModelParts! The number of present submodelparts is ", rParentModelPart.NumberOfSubModelParts());
        
        KRATOS_CATCH("");        
    };
    
    boost::python::object AuxiliaryUtilities::GetIthSubModelPartName(ModelPart& rParentModelPart, const int& required_i){
        KRATOS_TRY;                             
        int current_i = 0;
        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = rParentModelPart.SubModelPartsBegin(); sub_model_part != rParentModelPart.SubModelPartsEnd(); ++sub_model_part) {

            if(current_i == required_i) {
                const std::string& name = (*sub_model_part).Name();
                return boost::python::object(name);
            }
            current_i++;
        }

        KRATOS_THROW_ERROR(std::runtime_error, "The function GetIthSubModelPartName required a position greater than the number of SubModelParts! The number of present submodelparts is ", rParentModelPart.NumberOfSubModelParts());
        
        KRATOS_CATCH("");  
    };
    
    boost::python::object AuxiliaryUtilities::GetIthSubModelPartIdentifier(ModelPart& rParentModelPart, const int& required_i){
        KRATOS_TRY; 
        
        int current_i = 0;
        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = rParentModelPart.SubModelPartsBegin(); sub_model_part != rParentModelPart.SubModelPartsEnd(); ++sub_model_part) {

            if(current_i == required_i) {
                const std::string& identifier = (*sub_model_part)[IDENTIFIER];
                return boost::python::object(identifier);
            }
            
            current_i++;
        }

        KRATOS_THROW_ERROR(std::runtime_error, "The function GetIthSubModelPartIdentifier required a position greater than the number of SubModelParts! The number of present submodelparts is ", rParentModelPart.NumberOfSubModelParts());
        
        KRATOS_CATCH("");  
    };
    
    template<class T>
    boost::python::object AuxiliaryUtilities::GetIthSubModelPartData(ModelPart& rParentModelPart, const int& required_i, const Variable<T>& rVariable){
        KRATOS_TRY; 
        
        int current_i = 0;
        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = rParentModelPart.SubModelPartsBegin(); sub_model_part != rParentModelPart.SubModelPartsEnd(); ++sub_model_part) {

            if(current_i == required_i) {
                const T& value = (*sub_model_part)[rVariable];
                return boost::python::object(value);
            }
            
            current_i++;
        }

        KRATOS_THROW_ERROR(std::runtime_error, "The function GetIthSubModelPartData required a position greater than the number of SubModelParts! The number of present submodelparts is ", rParentModelPart.NumberOfSubModelParts());
        
        KRATOS_CATCH("");  
    }      
    
    // Explicit instantiations of the previous template function (otherwise it does not link):
    template boost::python::object AuxiliaryUtilities::GetIthSubModelPartData<int>(ModelPart& rParentModelPart, const int& required_i, const Variable<int>& rVariable);
    template boost::python::object AuxiliaryUtilities::GetIthSubModelPartData<double>(ModelPart& rParentModelPart, const int& required_i, const Variable<double>& rVariable);
    template boost::python::object AuxiliaryUtilities::GetIthSubModelPartData<array_1d<double,3> >(ModelPart& rParentModelPart, const int& required_i, const Variable<array_1d<double,3> >& rVariable);
    template boost::python::object AuxiliaryUtilities::GetIthSubModelPartData<std::string>(ModelPart& rParentModelPart, const int& required_i, const Variable<std::string>& rVariable);
       
    
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
