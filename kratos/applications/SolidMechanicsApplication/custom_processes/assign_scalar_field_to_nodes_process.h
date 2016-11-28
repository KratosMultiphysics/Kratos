//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_ASSIGN_SCALAR_FIELD_TO_NODES_PROCESS_H_INCLUDED )
#define  KRATOS_ASSIGN_SCALAR_FIELD_TO_NODES_PROCESS_H_INCLUDED



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
class AssignScalarFieldToNodesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssignScalarFieldToNodesProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignScalarFieldToNodesProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    AssignScalarFieldToNodesProcess(ModelPart& model_part,
					  PyObject* pPyObject,
					  const char* pPyMethodName,
					  const bool SpatialFieldFunction,
					  Parameters rParameters
					  ) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY
			 
        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "mesh_id": 0,
                "variable_name": "VARIABLE_NAME",
                "local_axes" : {}
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mmesh_id       = rParameters["mesh_id"].GetInt();
        mvariable_name = rParameters["variable_name"].GetString();

	mpPyObject      =  pPyObject;	
	mpPyMethodName  =  pPyMethodName;

	mIsSpatialField = SpatialFieldFunction;

	// Admissible values for local axes, are "empty" or 
        //"local_axes" :{
        //    "origin" : [0.0, 0.0, 0.0]
        //    "axes"   : [ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ] 
        //    }

	mHasLocalOrigin = false;
	if( rParameters["local_axes"].Has("origin") ){
	  mHasLocalOrigin = true;
	  mLocalOrigin.resize(3,false);
	  noalias(mLocalOrigin) = ZeroVector(3);
	  for( unsigned int i=0; i<3; i++)
	    mLocalOrigin[i] = rParameters["local_axes"]["origin"][i].GetDouble();
	}

	mHasLocalAxes = false;
	if( rParameters["local_axes"].Has("axes") ){
	  mHasLocalAxes = true;
	  mTransformationMatrix.resize(3,3,false);
	  noalias(mTransformationMatrix) = ZeroMatrix(3,3);	  
	  for( unsigned int i=0; i<3; i++)
	    for( unsigned int j=0; j<3; j++)
	      mTransformationMatrix(i,j) = rParameters["local_axes"]["axes"][i][j].GetDouble();
	}
	
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
        else
        {
            KRATOS_THROW_ERROR(std::logic_error, "Not able to set the variable type/name. Attempting to set variable:",mvariable_name);
        }

        KRATOS_CATCH("")
    }


    /// Destructor.
    virtual ~AssignScalarFieldToNodesProcess() {}


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


    /// Execute method is used to execute the AssignScalarFieldToNodesProcess algorithms.
    virtual void Execute() 
    {

        KRATOS_TRY;

        ProcessInfo& rCurrentProcessInfo = mr_model_part.GetProcessInfo();

	const double& rCurrentTime = rCurrentProcessInfo[TIME];
	
        if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(mvariable_name) ) //case of component variable
        {
            typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
            component_type var_component = KratosComponents< component_type >::Get(mvariable_name);
            InternalAssignValue< component_type>(var_component, rCurrentTime);
        }
        else if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
        {
	  InternalAssignValue<>(KratosComponents< Variable<double> >::Get(mvariable_name), rCurrentTime);
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
        return "AssignScalarFieldToNodesProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "AssignScalarFieldToNodesProcess";
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
    AssignScalarFieldToNodesProcess(AssignScalarFieldToNodesProcess const& rOther);

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

    PyObject* mpPyObject;  
    const char* mpPyMethodName;

    Vector mLocalOrigin;
    Matrix mTransformationMatrix;

    bool mIsSpatialField;
    
    bool mHasLocalOrigin;
    bool mHasLocalAxes;
    
    std::size_t mmesh_id;

    ///@}
    ///@name Private Operators
    ///@{

    void LocalAxesTransform(const double& rX_global, const double& rY_global, const double& rZ_global,
			    double& rx_local, double& ry_local, double& rz_local)
    {
      rx_local = rX_global;
      ry_local = rY_global;
      rz_local = rZ_global;

      if( mHasLocalOrigin  || mHasLocalAxes ){
      
	//implement global to local axes transformation
	Vector GlobalPosition(3);
	GlobalPosition[0] = rX_global;
	GlobalPosition[1] = rY_global;
	GlobalPosition[2] = rZ_global;

	if( mHasLocalOrigin )
	  GlobalPosition -= mLocalOrigin;

	Vector LocalPosition(3);
	noalias(LocalPosition) = ZeroVector(3);

	if( mHasLocalAxes )
	  noalias(LocalPosition) = prod(mTransformationMatrix,GlobalPosition);

	rx_local = LocalPosition[0];
	ry_local = LocalPosition[1];
	rz_local = LocalPosition[2];

      }
    }

    void CallFunction(const Node<3>::Pointer& pNode, const double& time, double& rValue)
    {
      
      if( mIsSpatialField ){

	double x = 0.0, y = 0.0, z = 0.0;
	
	LocalAxesTransform(pNode->X(), pNode->Y(), pNode->Z(), x, y, z);
    
	rValue = boost::python::call_method<double>(mpPyObject, mpPyMethodName, x, y, z, time);
	
      }
      else{
	
	rValue = boost::python::call_method<double>(mpPyObject, mpPyMethodName, 0.0, 0.0, 0.0, time);
	
      }
      
    }
    
    
    template< class TVarType>
    void InternalAssignValue(TVarType& rVariable, const double& rTime)
    {
        const int nnodes = mr_model_part.GetMesh(mmesh_id).Nodes().size();

	double Value = 0;
	
        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh(mmesh_id).NodesBegin();

            //#pragma omp parallel for  //it does not work in parallel
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

		this->CallFunction(*(it.base()), rTime, Value);
				
		it->FastGetSolutionStepValue(rVariable) = Value;
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
    AssignScalarFieldToNodesProcess& operator=(AssignScalarFieldToNodesProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AssignScalarFieldToNodesProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignScalarFieldToNodesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignScalarFieldToNodesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ASSIGN_SCALAR_FIELD_TO_NODES_PROCESS_H_INCLUDED  defined
