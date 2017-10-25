//  KratosShapeOptimizationApplication
//
//  License:		 BSD License
//					 license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//                   


#if !defined(KRATOS_REPLACE_ELEMENTS_AND_CONDITIONS_FOR_ADJOINT_PROBLEM_PROCESS_H_INCLUDED )
#define  KRATOS_REPLACE_ELEMENTS_AND_CONDITIONS_FOR_ADJOINT_PROBLEM_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for all processes in Kratos.
/** This function applies a constant value (and fixity) to all of the nodes in a given mesh
*/
class ReplaceElementsAndConditionsForAdjointProblemProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_DEFINE_LOCAL_FLAG(VARIABLE_IS_FIXED);


    /// Pointer definition of ReplaceElementsAndConditionsForAdjointProblemProcess
    KRATOS_CLASS_POINTER_DEFINITION(ReplaceElementsAndConditionsForAdjointProblemProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    ReplaceElementsAndConditionsForAdjointProblemProcess(ModelPart& model_part, 
                              Parameters Settings
                                   ) : Process(Flags()) , mr_model_part(model_part), mSettings( Settings)
    {
        KRATOS_TRY

//only include validation with c++11 since raw_literals do not exist in c++03
#if __cplusplus >= 201103L

        Parameters default_parameters( R"(
            {
                "Add_string": "NAME_OF_ADD_STRING",
                "Add_before_in_element_name": "ADD_STRING_BEFORE",
                "Add_before_in_condition_name": "ADD_STRING_BEFORE"
            }  )" );

        //now validate agains defaults -- this also ensures no type mismatch*/
        Settings.ValidateAndAssignDefaults(default_parameters);
#endif

        
        KRATOS_CATCH("")
    }
    

    
    /// Destructor.
    ~ReplaceElementsAndConditionsForAdjointProblemProcess() override {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    /// Execute method is used to execute the ReplaceElementsAndConditionsForAdjointProblemProcess algorithms.
    void Execute()  override
    {
        ModelPart& r_root_model_part = ObtainRootModelPart( mr_model_part );

        std::string sub_name_element = mSettings["Add_before_in_element_name"].GetString();
        std::string sub_name_condition = mSettings["Add_before_in_condition_name"].GetString();
        std::string adding_string = mSettings["Add_string"].GetString();
        
    #pragma omp parallel for                              //--> TODO: Check if this really works in parallel
        for(int i=0; i< (int)r_root_model_part.Elements().size(); i++)
        {
            ModelPart::ElementsContainerType::iterator it = r_root_model_part.ElementsBegin() + i;

            std::string element_name = it->Info();
            std::string::size_type position_ele = 0;
            std::string::size_type found_ele;
            found_ele = element_name.find(sub_name_element, position_ele);
            element_name.insert(found_ele, adding_string);

            if( !KratosComponents< Element >::Has( element_name ) )
                KRATOS_THROW_ERROR(std::invalid_argument, "Element name not found in KratosComponents< Element > -- name is ", element_name);
            const Element& rReferenceElement = KratosComponents<Element>::Get(element_name); //is this a problem to initialize the element in the loop?
            Element::Pointer p_element = rReferenceElement.Create(it->Id(), it->pGetGeometry(), it->pGetProperties());
 
            //deep copy elemental data
            p_element->Data() = it->Data();
            
            (*it.base()) = p_element;
        }
        
    #pragma omp parallel for                              //--> TODO: Check if this really works in parallel
        for(int i=0; i< (int)r_root_model_part.Conditions().size(); i++)
        {
            ModelPart::ConditionsContainerType::iterator it = r_root_model_part.ConditionsBegin() + i;

            std::string condition_name = it->Info();
            std::string::size_type position_cond = 0;
            std::string::size_type found_cond;
            found_cond = condition_name.find(sub_name_condition, position_cond);
            condition_name.insert(found_cond, adding_string);
      
            if( !KratosComponents< Condition >::Has( condition_name ) )
                KRATOS_THROW_ERROR(std::invalid_argument, "Condition name not found in KratosComponents< Condition > -- name is ", condition_name);
            const Condition& rReferenceCondition = KratosComponents<Condition>::Get(condition_name); //is this a problem to initialize the element in the loop?
      
            Condition::Pointer p_condition = rReferenceCondition.Create(it->Id(), it->pGetGeometry(), it->pGetProperties());
     
            //deep copy elemental data
            p_condition->Data() = it->Data();
            
            (*it.base()) = p_condition;
        }      
        
        //change the sons
        for (ModelPart::SubModelPartIterator i_sub_model_part = r_root_model_part.SubModelPartsBegin(); i_sub_model_part != r_root_model_part.SubModelPartsEnd(); i_sub_model_part++)
            UpdateSubModelPart( *i_sub_model_part, r_root_model_part );


    }




    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ReplaceElementsAndConditionsForAdjointProblemProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ReplaceElementsAndConditionsForAdjointProblemProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}
protected:
    
    ModelPart& mr_model_part;
    Parameters mSettings;

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{    
    ModelPart& ObtainRootModelPart( ModelPart& r_model_part )
    {
        if (r_model_part.IsSubModelPart())
            return ObtainRootModelPart(*r_model_part.GetParentModelPart());
        else
            return r_model_part;
    }
    
    void UpdateSubModelPart(ModelPart& r_model_part, ModelPart& r_root_model_part)
    {
        //change the model part itself
        #pragma omp parallel for
        for(int i=0; i< (int)r_model_part.Elements().size(); i++)
        {
            ModelPart::ElementsContainerType::iterator it = r_model_part.ElementsBegin() + i;
            
            (*it.base()) = r_root_model_part.Elements()(it->Id());
        }
        
        #pragma omp parallel for
        for(int i=0; i< (int)r_model_part.Conditions().size(); i++)
        {
            ModelPart::ConditionsContainerType::iterator it = r_model_part.ConditionsBegin() + i;
            
            (*it.base()) = r_root_model_part.Conditions()(it->Id());
        }
            
        //change the sons
        for (ModelPart::SubModelPartIterator i_sub_model_part = r_model_part.SubModelPartsBegin(); i_sub_model_part != r_model_part.SubModelPartsEnd(); i_sub_model_part++)
            UpdateSubModelPart( *i_sub_model_part, r_root_model_part );

    }

    /// Assignment operator.
    ReplaceElementsAndConditionsForAdjointProblemProcess& operator=(ReplaceElementsAndConditionsForAdjointProblemProcess const& rOther);

    /// Copy constructor.
    //ReplaceElementsAndConditionsForAdjointProblemProcess(ReplaceElementsAndConditionsForAdjointProblemProcess const& rOther);


    ///@}

}; // Class ReplaceElementsAndConditionsForAdjointProblemProcess

KRATOS_CREATE_LOCAL_FLAG(ReplaceElementsAndConditionsForAdjointProblemProcess,VARIABLE_IS_FIXED, 0);

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ReplaceElementsAndConditionsForAdjointProblemProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ReplaceElementsAndConditionsForAdjointProblemProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_REPLACE_ELEMENTS_AND_CONDITIONS_FOR_ADJOINT_PROBLEM_PROCESS_H_INCLUDED  defined 


