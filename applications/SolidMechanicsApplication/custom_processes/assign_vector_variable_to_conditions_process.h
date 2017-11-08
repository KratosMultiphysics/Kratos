//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_ASSIGN_VECTOR_VARIABLE_TO_CONDITIONS_PROCESS_H_INCLUDED )
#define  KRATOS_ASSIGN_VECTOR_VARIABLE_TO_CONDITIONS_PROCESS_H_INCLUDED



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
class AssignVectorVariableToConditionsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssignVectorVariableToConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignVectorVariableToConditionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    AssignVectorVariableToConditionsProcess(ModelPart& model_part,
					    Parameters rParameters) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY
			 
        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "variable_name": "VARIABLE_NAME",
                "value" : [0.0, 0.0, 0.0]
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mvariable_name = rParameters["variable_name"].GetString();

        if(KratosComponents< Variable<array_1d<double,3> > >::Has(mvariable_name) == false)
        {
            KRATOS_THROW_ERROR(std::runtime_error,"trying to set a variable that is not in the model_part - variable name is ",mvariable_name);
        }

        mvector_value[0] = rParameters["value"][0].GetDouble();
        mvector_value[1] = rParameters["value"][1].GetDouble();
        mvector_value[2] = rParameters["value"][2].GetDouble();

        KRATOS_CATCH("");
    }

    
    AssignVectorVariableToConditionsProcess(ModelPart& model_part,
					    const Variable<array_1d<double,3> >& rVariable,
					    const array_1d<double,3>& rvector_value) : Process() , mr_model_part(model_part), mvector_value(rvector_value)
    {
        KRATOS_TRY;

        mvariable_name = rVariable.Name();

	if( KratosComponents< Variable<array_1d<double,3> > >::Has( mvariable_name ) == false ) //case of array_1d variable
	  KRATOS_THROW_ERROR(std::runtime_error,"trying to set a variable that is not in the model_part - variable name is ",mvariable_name);

        KRATOS_CATCH("")
    }


    /// Destructor.
    virtual ~AssignVectorVariableToConditionsProcess() {}


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


    /// Execute method is used to execute the AssignVectorVariableToConditionsProcess algorithms.
    virtual void Execute() 
    {

        KRATOS_TRY
 
	InternalAssignValue(KratosComponents< Variable<array_1d<double,3> > >::Get(mvariable_name), mvector_value);

        KRATOS_CATCH("")

    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    virtual void ExecuteInitialize()
    {
    }

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    virtual void ExecuteBeforeSolutionLoop()
    {
    }


    /// this function will be executed at every time step BEFORE performing the solve phase
    virtual void ExecuteInitializeSolutionStep()
    {
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    virtual void ExecuteFinalizeSolutionStep()
    {
    }


    /// this function will be executed at every time step BEFORE  writing the output
    virtual void ExecuteBeforeOutputStep()
    {
    }


    /// this function will be executed at every time step AFTER writing the output
    virtual void ExecuteAfterOutputStep()
    {
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    virtual void ExecuteFinalize()
    {
        KRATOS_TRY

	array_1d<double,3> vector_value;
	vector_value.clear();
	InternalAssignValue(KratosComponents< Variable<array_1d<double,3> > >::Get(mvariable_name), vector_value);

        KRATOS_CATCH("")
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
        return "AssignVectorVariableToConditionsProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "AssignVectorVariableToConditionsProcess";
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
   
    /// Copy constructor.
    AssignVectorVariableToConditionsProcess(AssignVectorVariableToConditionsProcess const& rOther);

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
    array_1d<double,3> mvector_value;
  
    ///@}
    ///@name Private Operators
    ///@{
    
    void InternalAssignValue(const Variable<array_1d<double,3> >& rVariable,
			     const array_1d<double,3>& rvector_value)
    {
        const int nconditions = mr_model_part.GetMesh().Conditions().size();

        if(nconditions != 0)
        {
            ModelPart::ConditionsContainerType::iterator it_begin = mr_model_part.GetMesh().ConditionsBegin();

            #pragma omp parallel for
            for(int i = 0; i<nconditions; i++)
            {
                ModelPart::ConditionsContainerType::iterator it = it_begin + i;

                it->SetValue(rVariable, rvector_value);
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
    AssignVectorVariableToConditionsProcess& operator=(AssignVectorVariableToConditionsProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AssignVectorVariableToConditionsProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignVectorVariableToConditionsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignVectorVariableToConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ASSIGN_VECTOR_VARIABLE_TO_CONDITIONS_PROCESS_H_INCLUDED  defined
