//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_FREE_SCALAR_DOF_PROCESS_H_INCLUDED )
#define  KRATOS_FREE_SCALAR_DOF_PROCESS_H_INCLUDED


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

/// The base class for freeing scalar variable Dof or array_1d component Dof processes in Kratos.
/** This function free the variable dof belonging to all of the nodes in a given mesh
*/
class FreeScalarDofProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FreeScalarDofProcess
    KRATOS_CLASS_POINTER_DEFINITION(FreeScalarDofProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    FreeScalarDofProcess(ModelPart& model_part,
			Parameters rParameters
			) : Process() , mr_model_part(model_part)
    {
        KRATOS_TRY
			 
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id": 0,
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME"
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mmesh_id       = rParameters["mesh_id"].GetInt();
        mvariable_name = rParameters["variable_name"].GetString();

        if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(mvariable_name) ) //case of component variable
        {
            typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
            component_type var_component = KratosComponents< component_type >::Get(mvariable_name);

            if( model_part.GetNodalSolutionStepVariablesList().Has( var_component.GetSourceVariable() ) == false )
            {
                KRATOS_THROW_ERROR(std::runtime_error,"trying to set a variable that is not in the model_part - variable name is ",mvariable_name);
            }

        }
        else if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
        {
            if( model_part.GetNodalSolutionStepVariablesList().Has( KratosComponents< Variable<double> >::Get( mvariable_name ) ) == false )
            {
                KRATOS_THROW_ERROR(std::runtime_error,"trying to set a variable that is not in the model_part - variable name is ",mvariable_name);
            }

        }
        else if( KratosComponents< Variable<int> >::Has( mvariable_name ) ) //case of int variable
        {
            if( model_part.GetNodalSolutionStepVariablesList().Has( KratosComponents< Variable<int> >::Get( mvariable_name ) ) == false )
            {
                KRATOS_THROW_ERROR(std::runtime_error,"trying to set a variable that is not in the model_part - variable name is ",mvariable_name);
            }

        }
        else if( KratosComponents< Variable<bool> >::Has( mvariable_name ) ) //case of bool variable
        {
	  if( model_part.GetNodalSolutionStepVariablesList().Has(KratosComponents< Variable<bool> >::Get( mvariable_name ) ) == false )
            {
                KRATOS_THROW_ERROR(std::runtime_error,"trying to set a variable that is not in the model_part - variable name is ",mvariable_name);
            }
        }

        KRATOS_CATCH("");
    }

    FreeScalarDofProcess(ModelPart& model_part,
			const Variable<double>& rVariable,
			std::size_t mesh_id
			) : Process(), mr_model_part(model_part), mmesh_id(mesh_id)
    {
        KRATOS_TRY;


        if( model_part.GetNodalSolutionStepVariablesList().Has( rVariable ) == false )
        {
                KRATOS_THROW_ERROR(std::runtime_error,"trying to set a variable that is not in the model_part - variable name is ",rVariable);
        }

        mvariable_name = rVariable.Name();

        KRATOS_CATCH("");
    }

    FreeScalarDofProcess(ModelPart& model_part,
			const VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >& rVariable,
			std::size_t mesh_id
			) : Process(), mr_model_part(model_part), mmesh_id(mesh_id)
    {
        KRATOS_TRY;

        if( model_part.GetNodalSolutionStepVariablesList().Has( rVariable.GetSourceVariable() ) == false )
        {
                KRATOS_THROW_ERROR(std::runtime_error,"trying to set a variable that is not in the model_part - variable name is ",rVariable);
        }

	mvariable_name = rVariable.Name();

        KRATOS_CATCH("");
    }

    FreeScalarDofProcess(ModelPart& model_part,
			const Variable< int >& rVariable,
			std::size_t mesh_id
			) : Process(), mr_model_part(model_part), mmesh_id(mesh_id)
    {
        KRATOS_TRY;

        if( model_part.GetNodalSolutionStepVariablesList().Has( rVariable ) == false )
        {
                KRATOS_THROW_ERROR(std::runtime_error,"Trying to set a variable that is not in the model_part - variable name is ",rVariable);
        }

	mvariable_name = rVariable.Name();

        KRATOS_CATCH("");
    }

    FreeScalarDofProcess(ModelPart& model_part,
			const Variable< bool >& rVariable,
			std::size_t mesh_id
			) : Process(), mr_model_part(model_part), mmesh_id(mesh_id)
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
    virtual ~FreeScalarDofProcess() {}


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


    /// Execute method is used to execute the FreeScalarDofProcess algorithms.
    virtual void Execute() 
    {

        KRATOS_TRY;
 
        if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(mvariable_name) ) //case of component variable
        {
            typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
            component_type var_component = KratosComponents< component_type >::Get(mvariable_name);
            InternalFreeDof< component_type >(var_component);
        }
        else if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
        {
            InternalFreeDof<>(KratosComponents< Variable<double> >::Get(mvariable_name));
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
        return "FreeScalarDofProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "FreeScalarDofProcess";
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
    FreeScalarDofProcess(FreeScalarDofProcess const& rOther);

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
    std::size_t mmesh_id;

    ///@}
    ///@name Private Operators
    ///@{

    template< class TVarType >
    void InternalFreeDof(TVarType& rVar)
    {
        const int nnodes = mr_model_part.GetMesh(mmesh_id).Nodes().size();

        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh(mmesh_id).NodesBegin();

             #pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

		it->pAddDof(rVar)->FreeDof();
		//it->Free(rVar);
                
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
    FreeScalarDofProcess& operator=(FreeScalarDofProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class FreeScalarDofProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  FreeScalarDofProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FreeScalarDofProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FREE_SCALAR_DOF_PROCESS_H_INCLUDED  defined
