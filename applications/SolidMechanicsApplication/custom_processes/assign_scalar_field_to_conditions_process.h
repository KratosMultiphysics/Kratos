//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_ASSIGN_SCALAR_FIELD_TO_CONDITIONS_PROCESS_H_INCLUDED )
#define  KRATOS_ASSIGN_SCALAR_FIELD_TO_CONDITIONS_PROCESS_H_INCLUDED



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
class AssignScalarFieldToConditionsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssignScalarFieldToConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignScalarFieldToConditionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    AssignScalarFieldToConditionsProcess(ModelPart& model_part,
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
                "variable_name": "VARIABLE_NAME",
                "local_axes" : {}
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

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

	if( KratosComponents< Variable<Vector> >::Has( mvariable_name ) == false ) //case of double variable
        {
	  KRATOS_THROW_ERROR(std::runtime_error,"trying to set a variable that is not in the model_part - variable name is ",mvariable_name); 
	}
	
        KRATOS_CATCH("")
    }



    /// Destructor.
    virtual ~AssignScalarFieldToConditionsProcess() {}


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


    /// Execute method is used to execute the AssignScalarFieldToConditionsProcess algorithms.
    virtual void Execute() 
    {

        KRATOS_TRY

	ProcessInfo& rCurrentProcessInfo = mr_model_part.GetProcessInfo();

	const double& rCurrentTime = rCurrentProcessInfo[TIME];

 
	if( KratosComponents< Variable<Vector> >::Has( mvariable_name ) ) //case of double variable
        {
	  InternalAssignValue<>(KratosComponents< Variable<Vector> >::Get(mvariable_name), rCurrentTime);
        }
        else
        {
            KRATOS_THROW_ERROR(std::logic_error, "Not able to set the variable. Attempting to set variable:",mvariable_name);
        }

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
 
	if( KratosComponents< Variable<Vector> >::Has( mvariable_name ) ) //case of double variable
        {
	  Vector Value(3);
	  noalias(Value) = ZeroVector(3);
	  InternalAssignValue<>(KratosComponents< Variable<Vector> >::Get(mvariable_name), Value);
        }
        else
        {
            KRATOS_THROW_ERROR(std::logic_error, "Not able to set the variable. Attempting to set variable:",mvariable_name);
        }

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
        return "AssignScalarFieldToConditionsProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "AssignScalarFieldToConditionsProcess";
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
    AssignScalarFieldToConditionsProcess(AssignScalarFieldToConditionsProcess const& rOther);

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

    void CallFunction(const Condition::Pointer& pCondition, const double& time, Vector& rValue)
    {
      
      Condition::GeometryType& rConditionGeometry = pCondition->GetGeometry();
      unsigned int size = rConditionGeometry.size();
      
      rValue.resize(size,false);
      
      
      if( mIsSpatialField ){
	
	double x = 0, y = 0, z = 0;
	  
	for(unsigned int i=0; i<size; i++)
	  {
	    LocalAxesTransform(rConditionGeometry[i].X(), rConditionGeometry[i].Y(), rConditionGeometry[i].Z(), x, y, z);    
	    rValue[i] = boost::python::call_method<double>(mpPyObject, mpPyMethodName, x, y, z, time);
	  }
	
      }
      else{

	double value = boost::python::call_method<double>(mpPyObject, mpPyMethodName, 0.0, 0.0, 0.0, time);
	for(unsigned int i=0; i<size; i++)
	  {
	    rValue[i] = value;
	  }
	
      }
      
    }

    
    template< class TVarType >
    void InternalAssignValue(TVarType& rVar, const double& rTime)
    {
        const int nconditions = mr_model_part.GetMesh().Conditions().size();

	Vector Value;
	
        if(nconditions != 0)
        {
            ModelPart::ConditionsContainerType::iterator it_begin = mr_model_part.GetMesh().ConditionsBegin();

            // #pragma omp parallel for //it does not work in parallel
            for(int i = 0; i<nconditions; i++)
            {
                ModelPart::ConditionsContainerType::iterator it = it_begin + i;

		this->CallFunction(*(it.base()), rTime, Value);
		
                it->SetValue(rVar, Value);
            }
        }
    }

    template< class TVarType, class TDataType >
    void InternalAssignValue(TVarType& rVar, const TDataType value)
    {
      const int nconditions = mr_model_part.GetMesh().Conditions().size();

        if(nconditions != 0)
        {
            ModelPart::ConditionsContainerType::iterator it_begin = mr_model_part.GetMesh().ConditionsBegin();

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
    AssignScalarFieldToConditionsProcess& operator=(AssignScalarFieldToConditionsProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AssignScalarFieldToConditionsProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignScalarFieldToConditionsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignScalarFieldToConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ASSIGN_SCALAR_FIELD_TO_CONDITIONS_PROCESS_H_INCLUDED  defined
