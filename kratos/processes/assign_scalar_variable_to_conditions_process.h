//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Josep Maria Carbonell
//
#if !defined(KRATOS_ASSIGN_SCALAR_VARIABLE_TO_CONDITIONS_PROCESS_H_INCLUDED )
#define  KRATOS_ASSIGN_SCALAR_VARIABLE_TO_CONDITIONS_PROCESS_H_INCLUDED



// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for assigning a value to scalar variables or array_1d components processes in Kratos.
/** This function assigns a value to a variable belonging to all of the nodes in a given mesh
*/
class AssignScalarVariableToConditionsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssignScalarVariableToConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignScalarVariableToConditionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    AssignScalarVariableToConditionsProcess(ModelPart& model_part,
						 Parameters rParameters
						 ) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY
			 
        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "mesh_id": 0,
                "variable_name": "VARIABLE_NAME",
                "value" : 1.0
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mmesh_id       = rParameters["mesh_id"].GetInt();
        mvariable_name = rParameters["variable_name"].GetString();

	if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
        {
	  mdouble_value = rParameters["value"].GetDouble();
        }
        else if( KratosComponents< Variable<int> >::Has( mvariable_name ) ) //case of int variable
        {
	  mint_value = rParameters["value"].GetInt();
        }
        else if( KratosComponents< Variable<bool> >::Has( mvariable_name ) ) //case of bool variable
        {
	  mbool_value = rParameters["value"].GetBool();
        }
	else{
	  KRATOS_THROW_ERROR(std::runtime_error,"trying to set a variable that is not in the model_part - variable name is ",mvariable_name); 
	}
	
        KRATOS_CATCH("");
    }



    /// Destructor.
    virtual ~AssignScalarVariableToConditionsProcess() {}


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


    /// Execute method is used to execute the AssignScalarVariableToConditionsProcess algorithms.
    virtual void Execute() override
    {

        KRATOS_TRY;
 
	if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
        {
            InternalAssignValue<>(KratosComponents< Variable<double> >::Get(mvariable_name), mdouble_value);
        }
        else if( KratosComponents< Variable<int> >::Has( mvariable_name ) ) //case of int variable
        {
            InternalAssignValue<>(KratosComponents< Variable<int> >::Get(mvariable_name) , mint_value);
        }
        else if( KratosComponents< Variable<bool> >::Has( mvariable_name ) ) //case of bool variable
        {
            InternalAssignValue<>(KratosComponents< Variable<bool> >::Get(mvariable_name), mbool_value);
        }
        else
        {
            KRATOS_THROW_ERROR(std::logic_error, "Not able to set the variable. Attempting to set variable:",mvariable_name);
        }

        KRATOS_CATCH("");

    }
    
    virtual void ExecuteInitializeSolutionStep() override
    {
        Execute();
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
        return "AssignScalarVariableToConditionsProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "AssignScalarVariableToConditionsProcess";
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

    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{
    ///@}
    ///@name Protected Operators
    ///@{

    /// Copy constructor.
    AssignScalarVariableToConditionsProcess(AssignScalarVariableToConditionsProcess const& rOther);

    ///@}
    ///@name Protected Operations
    ///@{
    ///@}
    ///@name Protected  Access
    ///@{
    ///@}
    ///@name Protected Inquiry
    ///@{
    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:

    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mr_model_part;
    std::string mvariable_name;
    double mdouble_value;
    int mint_value;
    bool mbool_value;
    std::size_t mmesh_id;

    ///@}
    ///@name Private Operators
    ///@{

    template< class TVarType, class TDataType >
    void InternalAssignValue(TVarType& rVar, const TDataType value)
    {
        const int nconditions = mr_model_part.GetMesh(mmesh_id).Conditions().size();

        if(nconditions != 0)
        {
            ModelPart::ConditionsContainerType::iterator it_begin = mr_model_part.GetMesh(mmesh_id).ConditionsBegin();

             #pragma omp parallel for
            for(int i = 0; i<nconditions; i++)
            {
                ModelPart::ConditionsContainerType::iterator it = it_begin + i;

                it->SetValue(rVar, value);
            }
        }
    }

    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    AssignScalarVariableToConditionsProcess& operator=(AssignScalarVariableToConditionsProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AssignScalarVariableToConditionsProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignScalarVariableToConditionsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignScalarVariableToConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ASSIGN_SCALAR_VARIABLE_TO_CONDITIONS_PROCESS_H_INCLUDED  defined
