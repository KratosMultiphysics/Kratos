//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                    
//

#if !defined(KRATOS_REPLACE_ELEMENTS_AND_CONDITIONS_PROCESS_H_INCLUDED )
#define  KRATOS_REPLACE_ELEMENTS_AND_CONDITIONS_PROCESS_H_INCLUDED



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
class ReplaceElementsAndConditionsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_DEFINE_LOCAL_FLAG(VARIABLE_IS_FIXED);


    /// Pointer definition of ReplaceElementsAndConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ReplaceElementsAndConditionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    ReplaceElementsAndConditionsProcess(ModelPart& model_part, 
                              Parameters Settings
                                   ) : Process(Flags()) , mr_model_part(model_part), mSettings( Settings)
    {
        KRATOS_TRY

//only include validation with c++11 since raw_literals do not exist in c++03
#if __cplusplus >= 201103L

        Parameters default_parameters( R"(
            {
                "element_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "condition_name": "PLEASE_PRESCRIBE_VARIABLE_NAME"
            }  )" );

        //some vvalues need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        //so that an error is thrown if they don't exist
        if( !KratosComponents< Element >::Has( Settings["element_name"].GetString() ) )
                KRATOS_THROW_ERROR(std::invalid_argument, "Element name not found in KratosComponents< Element > -- name is ", Settings["element_name"].GetString());
        if( !KratosComponents< Condition >::Has( Settings["condition_name"].GetString() ) )
                KRATOS_THROW_ERROR(std::invalid_argument, "Condition name not found in KratosComponents< Condition > -- name is ", Settings["condition_name"].GetString());        
        
        //now validate agains defaults -- this also ensures no type mismatch
        
        Settings.ValidateAndAssignDefaults(default_parameters);
#endif

        
        KRATOS_CATCH("")
    }
    

    
    /// Destructor.
    virtual ~ReplaceElementsAndConditionsProcess() {}


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


    /// Execute method is used to execute the ReplaceElementsAndConditionsProcess algorithms.
    virtual void Execute() 
    {
        ModelPart& r_root_model_part = ObtainRootModelPart( mr_model_part );
        
	KRATOS_WATCH(mSettings.PrettyPrintJsonString());
        const Element& rReferenceElement = KratosComponents<Element>::Get(mSettings["element_name"].GetString());
        const Condition& rReferenceCondition = KratosComponents<Condition>::Get(mSettings["condition_name"].GetString());
        
        #pragma omp parallel for
        for(int i=0; i< (int)r_root_model_part.Elements().size(); i++)
        {
            ModelPart::ElementsContainerType::iterator it = r_root_model_part.ElementsBegin() + i;
            
            Element::Pointer p_element = rReferenceElement.Create(it->Id(), it->pGetGeometry(), it->pGetProperties());
            
            (*it.base()) = p_element;

        }
        
        #pragma omp parallel for
        for(int i=0; i< (int)r_root_model_part.Conditions().size(); i++)
        {
            ModelPart::ConditionsContainerType::iterator it = r_root_model_part.ConditionsBegin() + i;
            
            Condition::Pointer p_condition = rReferenceCondition.Create(it->Id(), it->pGetGeometry(), it->pGetProperties());
            
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
    virtual std::string Info() const
    {
        return "ReplaceElementsAndConditionsProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ReplaceElementsAndConditionsProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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
    ReplaceElementsAndConditionsProcess& operator=(ReplaceElementsAndConditionsProcess const& rOther);

    /// Copy constructor.
    //ReplaceElementsAndConditionsProcess(ReplaceElementsAndConditionsProcess const& rOther);


    ///@}

}; // Class ReplaceElementsAndConditionsProcess

KRATOS_CREATE_LOCAL_FLAG(ReplaceElementsAndConditionsProcess,VARIABLE_IS_FIXED, 0);

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ReplaceElementsAndConditionsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ReplaceElementsAndConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_REPLACE_ELEMENTS_AND_CONDITIONS_PROCESS_H_INCLUDED  defined 


