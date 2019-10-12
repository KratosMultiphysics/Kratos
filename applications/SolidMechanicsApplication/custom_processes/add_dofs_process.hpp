//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                June 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_ADD_DOFS_PROCESS_H_INCLUDED)
#define  KRATOS_ADD_DOFS_PROCESS_H_INCLUDED



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

/// The base class for fixing scalar variable Dof or array_1d component Dof processes in Kratos.
/** This function fix the variable dof belonging to all of the nodes in a given mesh
*/
class AddDofsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Variable<array_1d<double, 3> >                                    VectorVariableType;
    typedef Variable<double>                                                  ScalarVariableType;
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > >      ComponentType;

    /// Pointer definition of AddDofsProcess
    KRATOS_CLASS_POINTER_DEFINITION(AddDofsProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    AddDofsProcess(ModelPart& model_part,
		   Parameters rParameters
		   ) : Process() , mrModelPart(model_part)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variables_list": [],
                "reactions_list": []

            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

	// Check variables vs reactions consistency
	if( rParameters["variables_list"].size() != rParameters["reactions_list"].size() )
	  KRATOS_ERROR << "variables_list and reactions_list has not the same number of components "<<std::endl;


	for(unsigned int i=0; i<rParameters["variables_list"].size(); i++)
	  {
	    if( !rParameters["variables_list"][i].IsString() )
	      KRATOS_ERROR << "variables_list contains a non-string variable name "<<std::endl;

	    std::string variable_name = rParameters["variables_list"][i].GetString();

	    bool supplied_reaction = true;
	    if(rParameters["reactions_list"][i].IsNull())
	      supplied_reaction = false;

	    if( KratosComponents< VectorVariableType >::Has( variable_name ) ){ //case of array_1d (vector with components) variable

	      const VectorVariableType& VectorVariable = KratosComponents< VectorVariableType >::Get(variable_name);
	      if( model_part.GetNodalSolutionStepVariablesList().Has( VectorVariable ) == false ){
		KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is "<<variable_name<<std::endl;
	      }
	      else{
		for(unsigned int j=0; j<3; j++)
		  {
		    std::string component_name = variable_name;
		    component_name += ms_components[j];
		    const ComponentType& ComponentVariable = KratosComponents< ComponentType >::Get(component_name);

		    if(supplied_reaction){
		      std::string reaction_component_name = rParameters["reactions_list"][i].GetString();
		      reaction_component_name += ms_components[j];
		      const ComponentType& ReactionComponentVariable = KratosComponents< ComponentType >::Get(reaction_component_name);
		      m_component_variables_list.push_back(&ComponentVariable);
		      m_component_reactions_list.push_back(&ReactionComponentVariable);
		    }
		    else{
		      m_component_variables_no_reaction_list.push_back(&ComponentVariable);
		    }

		  }
	      }
	    }
	    else if( KratosComponents< ComponentType >::Has(variable_name) ){ //case of component variable

	      const ComponentType& ComponentVariable = KratosComponents< ComponentType >::Get(variable_name);

	      if( model_part.GetNodalSolutionStepVariablesList().Has( ComponentVariable.GetSourceVariable() ) == false ){

		KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is "<<variable_name<<std::endl;
	      }
	      else{

		if(supplied_reaction){
		  std::string reaction_name = rParameters["reactions_list"][i].GetString();
		  const ComponentType& ReactionComponentVariable = KratosComponents< ComponentType >::Get(reaction_name);
		  m_component_variables_list.push_back(&ComponentVariable);
		  m_component_reactions_list.push_back(&ReactionComponentVariable);
		}
		else{
		  m_component_variables_no_reaction_list.push_back(&ComponentVariable);
		}

	      }


	    }
	    else if( KratosComponents< ScalarVariableType >::Has( variable_name ) ){ //case of double variable

	      const ScalarVariableType& ScalarVariable = KratosComponents< ScalarVariableType >::Get( variable_name );
	      if( model_part.GetNodalSolutionStepVariablesList().Has( ScalarVariable ) ==  false ){
		KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is "<<variable_name<<std::endl;
	      }
	      else{

		if(supplied_reaction){
		  std::string reaction_name = rParameters["reactions_list"][i].GetString();
		  const ScalarVariableType& ReactionVariable = KratosComponents< ScalarVariableType >::Get(reaction_name);
		  m_scalar_variables_list.push_back(&ScalarVariable);
		  m_scalar_reactions_list.push_back(&ReactionVariable);
		}
		else{
		  m_scalar_variables_no_reaction_list.push_back(&ScalarVariable);
		}

	      }

	    }
	    else{
	      KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is "<<variable_name<<std::endl;
	    }
	  }


        KRATOS_CATCH("")
    }


    AddDofsProcess(ModelPart& model_part,
		   const pybind11::list& rVariablesList,
		   const pybind11::list& rReactionsList
		   ) : Process(), mrModelPart(model_part)
    {
        KRATOS_TRY

	unsigned int number_variables = len(rVariablesList);
	unsigned int number_reactions = len(rReactionsList);

	// Check variables vs reactions consistency
	if( number_variables != number_reactions )
	  KRATOS_ERROR << "variables_list and reactions_list has not the same number of components "<<std::endl;

	for(unsigned int i=0; i<number_variables; i++)
	  {

	    //std::string variable_name = boost::python::extract<std::string>(rVariablesList[i]);
	    //std::string reaction_name = boost::python::extract<std::string>(rReactionsList[i]);
	    std::string variable_name = pybind11::cast<std::string>(rVariablesList[i]);
	    std::string reaction_name = pybind11::cast<std::string>(rReactionsList[i]);

	    bool supplied_reaction = true;
	    if(reaction_name == "NOT_DEFINED")
	      supplied_reaction = false;

	    if( KratosComponents< VectorVariableType >::Has( variable_name ) ){ //case of array_1d (vector with components) variable

	      const VectorVariableType& VectorVariable = KratosComponents< VectorVariableType >::Get(variable_name);
	      if( model_part.GetNodalSolutionStepVariablesList().Has( VectorVariable ) == false ){
		KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is "<<variable_name<<std::endl;
	      }
	      else{
		for(unsigned int j=0; j<3; j++)
		  {
		    std::string component_name = variable_name;
		    component_name += ms_components[j];
		    const ComponentType& ComponentVariable = KratosComponents< ComponentType >::Get(component_name);

		    if(supplied_reaction){
		      std::string reaction_component_name = reaction_name;
		      reaction_component_name += ms_components[j];
		      const ComponentType& ReactionComponentVariable = KratosComponents< ComponentType >::Get(reaction_component_name);
		      m_component_variables_list.push_back(&ComponentVariable);
		      m_component_reactions_list.push_back(&ReactionComponentVariable);
		    }
		    else{
		      m_component_variables_no_reaction_list.push_back(&ComponentVariable);
		    }

		  }
	      }
	    }
	    else if( KratosComponents< ComponentType >::Has(variable_name) ){ //case of component variable

	      const ComponentType& ComponentVariable = KratosComponents< ComponentType >::Get(variable_name);

	      if( model_part.GetNodalSolutionStepVariablesList().Has( ComponentVariable.GetSourceVariable() ) == false ){

		KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is "<<variable_name<<std::endl;
	      }
	      else{

		if(supplied_reaction){
		  const ComponentType& ReactionComponentVariable = KratosComponents< ComponentType >::Get(reaction_name);
		  m_component_variables_list.push_back(&ComponentVariable);
		  m_component_reactions_list.push_back(&ReactionComponentVariable);
		}
		else{
		  m_component_variables_list.push_back(&ComponentVariable);
		}

	      }

	    }
	    else if( KratosComponents< ScalarVariableType >::Has( variable_name ) ){ //case of double variable

	      const ScalarVariableType& ScalarVariable = KratosComponents< ScalarVariableType >::Get( variable_name );
	      if( model_part.GetNodalSolutionStepVariablesList().Has( ScalarVariable ) ==  false ){
		KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is "<<variable_name<<std::endl;
	      }
	      else{

		if(supplied_reaction){
		  const ScalarVariableType& ReactionVariable = KratosComponents< ScalarVariableType >::Get(reaction_name);
		  m_scalar_variables_list.push_back(&ScalarVariable);
		  m_scalar_reactions_list.push_back(&ReactionVariable);
		}
		else{
		  m_scalar_variables_no_reaction_list.push_back(&ScalarVariable);
		}

	      }

	    }
	    else{
	      KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is "<<variable_name<<std::endl;
	    }
	  }

        KRATOS_CATCH("")
    }


    /// Destructor.
    ~AddDofsProcess() override {}


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


    /// Execute method is used to execute the AddDofsProcess algorithms.
    void Execute() override
    {

        KRATOS_TRY;

	int number_of_nodes = mrModelPart.NumberOfNodes();
	ModelPart::NodeConstantIterator nodes_begin = mrModelPart.NodesBegin();

	/*
	//1nd way: (fastest) generating the dofs for the initial node and add to others (still fails if a variable or a dof is set when mdpa is read)
	AddNodalDofs(nodes_begin);
	ModelPart::NodeType::DofsContainerType& reference_dofs = nodes_begin->GetDofs();
        #pragma omp parallel for
	for (int k=0; k<number_of_nodes; k++)
	  {
	    ModelPart::NodeConstantIterator it = nodes_begin + k;

	    for(ModelPart::NodeType::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
	      {
		it->pAddDof( *iii );
	      }
	  }
	*/

	//2nd way:  (faster)
        // #pragma omp parallel for
	for (int k=0; k<number_of_nodes; k++)
	  {
	    ModelPart::NodeConstantIterator it = nodes_begin + k;
	    AddNodalDofs(it);
	  }


	/*
	//3rt way: add dofs in the standard way one by one to all nodes  (slower)
	AddNodalDofs();
	*/

	//CheckNodalData(nodes_begin);

        KRATOS_CATCH("");

    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
    }

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    void ExecuteBeforeSolutionLoop() override
    {
    }


    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {
    }


    /// this function will be executed at every time step BEFORE  writing the output
    void ExecuteBeforeOutputStep() override
    {
    }


    /// this function will be executed at every time step AFTER writing the output
    void ExecuteAfterOutputStep() override
    {
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    void ExecuteFinalize() override
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
    std::string Info() const override
    {
        return "AddDofsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AddDofsProcess";
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

    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{
    ///@}
    ///@name Protected Operators
    ///@{

    /// Copy constructor.
    AddDofsProcess(AddDofsProcess const& rOther);

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

    ModelPart& mrModelPart;

    const std::vector<std::string> ms_components {"_X", "_Y", "_Z"};

    std::vector<ComponentType const *> m_component_variables_list;
    std::vector<ComponentType const *> m_component_reactions_list;
    std::vector<ComponentType const *> m_component_variables_no_reaction_list;

    std::vector<ScalarVariableType const *> m_scalar_variables_list;
    std::vector<ScalarVariableType const *> m_scalar_reactions_list;
    std::vector<ScalarVariableType const *> m_scalar_variables_no_reaction_list;


    ///@}
    ///@name Private Operators
    ///@{


    void AddNodalDofs()
    {
      KRATOS_TRY

	int number_of_nodes = mrModelPart.NumberOfNodes();
	ModelPart::NodeConstantIterator nodes_begin = mrModelPart.NodesBegin();

	for( unsigned int i=0; i < m_component_variables_list.size(); i++ )
	  {
            #pragma omp parallel for
	    for (int k=0; k<number_of_nodes; k++)
	      {
		ModelPart::NodeConstantIterator it = nodes_begin + k;
		it->pAddDof(*m_component_variables_list[i],*m_component_reactions_list[i]);
	      }
	  }

	for( unsigned int j=0; j < m_component_variables_no_reaction_list.size(); j++ )
	  {
           #pragma omp parallel for
	    for (int k=0; k<number_of_nodes; k++)
	      {
		ModelPart::NodeConstantIterator it = nodes_begin + k;
		it->pAddDof(*m_component_variables_no_reaction_list[j]);
	      }
	  }

	for( unsigned int l=0; l < m_scalar_variables_list.size(); l++ )
	  {
           #pragma omp parallel for
	    for (int k=0; k<number_of_nodes; k++)
	      {
		ModelPart::NodeConstantIterator it = nodes_begin + k;
		it->pAddDof(*m_scalar_variables_list[l],*m_scalar_reactions_list[l]);
	      }
	  }

	for( unsigned int m=0; m < m_scalar_variables_no_reaction_list.size(); m++ )
	  {
           #pragma omp parallel for
	    for (int k=0; k<number_of_nodes; k++)
	      {
		ModelPart::NodeConstantIterator it = nodes_begin + k;
		it->pAddDof(*m_scalar_variables_no_reaction_list[m]);
	      }
	  }

      KRATOS_CATCH(" ")
    }


    void AddNodalDofs( ModelPart::NodeConstantIterator& node_it )
    {
      KRATOS_TRY

      for( unsigned int i=0; i < m_component_variables_list.size(); i++ )
	{
	  node_it->pAddDof(*m_component_variables_list[i],*m_component_reactions_list[i]);
	}

      for( unsigned int j=0; j < m_component_variables_no_reaction_list.size(); j++ )
	{
	  node_it->pAddDof(*m_component_variables_no_reaction_list[j]);
	}

      for( unsigned int l=0; l < m_scalar_variables_list.size(); l++ )
	{
	  node_it->pAddDof(*m_scalar_variables_list[l],*m_scalar_reactions_list[l]);
	}

      for( unsigned int m=0; m < m_scalar_variables_no_reaction_list.size(); m++ )
	{
	  node_it->pAddDof(*m_scalar_variables_no_reaction_list[m]);
	}

      KRATOS_CATCH(" ")
    }


    void CheckNodalData( ModelPart::NodeConstantIterator& node_it )
    {
      KRATOS_TRY

      std::cout<<" CHECK VARIABLES LIST KEYS "<<std::endl;

      VariablesListDataValueContainer VariablesList = (node_it)->SolutionStepData();

      std::cout<<" list size "<<VariablesList.pGetVariablesList()->size()<<std::endl;
      std::cout<<" Variable: "<<(*VariablesList.pGetVariablesList())[0]<<std::endl;
      std::cout<<" end "<<std::endl;

      KRATOS_CATCH(" ")
    }

    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    AddDofsProcess& operator=(AddDofsProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AddDofsProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AddDofsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AddDofsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ADD_DOFS_PROCESS_H_INCLUDED  defined
