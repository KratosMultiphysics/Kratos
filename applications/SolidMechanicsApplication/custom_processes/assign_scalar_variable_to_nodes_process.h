//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_ASSIGN_SCALAR_VARIABLE_TO_NODES_PROCESS_H_INCLUDED )
#define  KRATOS_ASSIGN_SCALAR_VARIABLE_TO_NODES_PROCESS_H_INCLUDED



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
class AssignScalarVariableToNodesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssignScalarVariableToNodesProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignScalarVariableToNodesProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    AssignScalarVariableToNodesProcess(ModelPart& model_part,
				       Parameters rParameters) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY
			 
        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "variable_name": "VARIABLE_NAME",
                "value" : 1.0
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mvariable_name = rParameters["variable_name"].GetString();

        if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(mvariable_name) ) //case of component variable
        {
            typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
            component_type var_component = KratosComponents< component_type >::Get(mvariable_name);

            if( model_part.GetNodalSolutionStepVariablesList().Has( var_component.GetSourceVariable() ) == false )
            {
                KRATOS_THROW_ERROR(std::runtime_error,"trying to set a variable that is not in the model_part - variable name is ",mvariable_name);
            }

            mdouble_value = rParameters["value"].GetDouble();
        }
        else if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
        {
            if( model_part.GetNodalSolutionStepVariablesList().Has( KratosComponents< Variable<double> >::Get( mvariable_name ) ) == false )
            {
                KRATOS_THROW_ERROR(std::runtime_error,"trying to set a variable that is not in the model_part - variable name is ",mvariable_name);
            }

            mdouble_value = rParameters["value"].GetDouble();
        }
        else if( KratosComponents< Variable<int> >::Has( mvariable_name ) ) //case of int variable
        {
            if( model_part.GetNodalSolutionStepVariablesList().Has( KratosComponents< Variable<int> >::Get( mvariable_name ) ) == false )
            {
                KRATOS_THROW_ERROR(std::runtime_error,"trying to set a variable that is not in the model_part - variable name is ",mvariable_name);
            }

            mint_value = rParameters["value"].GetInt();

        }
        else if( KratosComponents< Variable<bool> >::Has( mvariable_name ) ) //case of bool variable
        {
	    if( model_part.GetNodalSolutionStepVariablesList().Has( KratosComponents< Variable<bool> >::Get( mvariable_name ) ) == false )
            {
                KRATOS_THROW_ERROR(std::runtime_error,"trying to set a variable that is not in the model_part - variable name is ",mvariable_name);
            }

	    mbool_value = rParameters["value"].GetBool();
        }

        KRATOS_CATCH("");
    }

    AssignScalarVariableToNodesProcess(ModelPart& model_part,
				       const Variable<double>& rVariable,
				       const double double_value) : Process() , mr_model_part(model_part),mdouble_value(double_value), mint_value(0), mbool_value(false)
    {
        KRATOS_TRY;


        if( model_part.GetNodalSolutionStepVariablesList().Has( rVariable ) == false )
        {
                KRATOS_THROW_ERROR(std::runtime_error,"trying to set a variable that is not in the model_part - variable name is ",rVariable);
        }

        mvariable_name = rVariable.Name();

        KRATOS_CATCH("");
    }

    AssignScalarVariableToNodesProcess(ModelPart& model_part,
				       const VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >& rVariable,
				       const double double_value) : Process() , mr_model_part(model_part),mdouble_value(double_value), mint_value(0), mbool_value(false)
    {
        KRATOS_TRY;

        if( model_part.GetNodalSolutionStepVariablesList().Has( rVariable.GetSourceVariable() ) == false )
        {
                KRATOS_THROW_ERROR(std::runtime_error,"trying to set a variable that is not in the model_part - variable name is ",rVariable);
        }

        mvariable_name = rVariable.Name();

        KRATOS_CATCH("");
    }

    AssignScalarVariableToNodesProcess(ModelPart& model_part,
				       const Variable< int >& rVariable,
				       const int int_value) : Process() , mr_model_part(model_part),mdouble_value(0.0), mint_value(int_value), mbool_value(false)
    {
        KRATOS_TRY;

        if( model_part.GetNodalSolutionStepVariablesList().Has( rVariable ) == false )
        {
                KRATOS_THROW_ERROR(std::runtime_error,"Trying to set a variable that is not in the model_part - variable name is ",rVariable);
        }

        mvariable_name = rVariable.Name();

        KRATOS_CATCH("");
    }

    AssignScalarVariableToNodesProcess(ModelPart& model_part,
				       const Variable< bool >& rVariable,
				       const bool bool_value) : Process() , mr_model_part(model_part),mdouble_value(0.0), mint_value(0), mbool_value(bool_value)
    {
        KRATOS_TRY;


        if( model_part.GetNodalSolutionStepVariablesList().Has( rVariable ) == false )
        {
                KRATOS_THROW_ERROR(std::runtime_error,"Trying to set a variable that is not in the model_part - variable name is ",rVariable);
        }

        mvariable_name = rVariable.Name();

        KRATOS_CATCH("");
    }


    /// Destructor.
    virtual ~AssignScalarVariableToNodesProcess() {}


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


    /// Execute method is used to execute the AssignScalarVariableToNodesProcess algorithms.
    virtual void Execute() 
    {

        KRATOS_TRY;
 
        if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(mvariable_name) ) //case of component variable
        {
            typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
            component_type var_component = KratosComponents< component_type >::Get(mvariable_name);
            InternalAssignValue< component_type, double>(var_component , mdouble_value);
        }
        else if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
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
        return "AssignScalarVariableToNodesProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "AssignScalarVariableToNodesProcess";
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
    AssignScalarVariableToNodesProcess(AssignScalarVariableToNodesProcess const& rOther);

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

    ///@}
    ///@name Private Operators
    ///@{

    template< class TVarType, class TDataType >
    void InternalAssignValue(TVarType& rVar, const TDataType value)
    {
        const int nnodes = mr_model_part.Nodes().size();

        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh().NodesBegin();

             #pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                it->FastGetSolutionStepValue(rVar) = value;
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
    AssignScalarVariableToNodesProcess& operator=(AssignScalarVariableToNodesProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AssignScalarVariableToNodesProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignScalarVariableToNodesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignScalarVariableToNodesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ASSIGN_SCALAR_VARIABLE_TO_NODES_PROCESS_H_INCLUDED  defined
