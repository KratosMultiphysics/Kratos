//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_ASSIGN_VECTOR_FIELD_TO_CONDITIONS_PROCESS_H_INCLUDED)
#define  KRATOS_ASSIGN_VECTOR_FIELD_TO_CONDITIONS_PROCESS_H_INCLUDED



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
class KRATOS_API(SOLID_MECHANICS_APPLICATION) AssignVectorFieldToConditionsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssignVectorFieldToConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignVectorFieldToConditionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    AssignVectorFieldToConditionsProcess(ModelPart& model_part,
					 pybind11::object& rPyObject,
					 const std::string& rPyMethodName,
					 const bool SpatialFieldFunction,
					 Parameters rParameters
					 ) : Process(Flags()), mrModelPart(model_part), mPyObject(rPyObject), mPyMethodName(rPyMethodName), mIsSpatialField(SpatialFieldFunction)
    {
        KRATOS_TRY
			 
        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "variable_name": "VARIABLE_NAME",
                "value" : [0.0, 0.0, 0.0],
                "local_axes" : {}
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mvariable_name = rParameters["variable_name"].GetString();

	// Admissible values for local axes, are "empty" or 
        //"local_axes" :{
        //    "origin" : [0.0, 0.0, 0.0]
        //    "axes"   : [ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ] 
        //    }

	mHasLocalOrigin = false;
	if( rParameters["local_axes"].Has("origin") ){
	  mHasLocalOrigin = true;
	  mLocalOrigin.resize(3,false);
	  for( unsigned int i=0; i<3; i++)
	    mLocalOrigin[i] = rParameters["local_axes"]["origin"][i].GetDouble();
	}

	mHasLocalAxes = false;
	if( rParameters["local_axes"].Has("axes") ){
	  mHasLocalAxes = true;
	  mTransformationMatrix.resize(3,3,false);
	  for( unsigned int i=0; i<3; i++)
	    for( unsigned int j=0; j<3; j++)
	      mTransformationMatrix(i,j) = rParameters["local_axes"]["axes"][i][j].GetDouble();
	}
	
        if(KratosComponents< Variable<Vector> >::Has(mvariable_name) == false)
        {
            KRATOS_THROW_ERROR(std::runtime_error,"trying to set a variable that is not in the model_part - variable name is ",mvariable_name);
        }

        mvector_value[0] = rParameters["value"][0].GetDouble();
        mvector_value[1] = rParameters["value"][1].GetDouble();
        mvector_value[2] = rParameters["value"][2].GetDouble();

        KRATOS_CATCH("");
    }

    /// Destructor.
    virtual ~AssignVectorFieldToConditionsProcess() {}


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


    /// Execute method is used to execute the AssignVectorFieldToConditionsProcess algorithms.
    virtual void Execute() 
    {

        KRATOS_TRY

	ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();

	const double& rCurrentTime = rCurrentProcessInfo[TIME];

 	if( KratosComponents< Variable<Vector> >::Has( mvariable_name ) ) //case of vector variable
        {

	  Vector Value;
	  //KratosComponents< Variable<Vector> >, Vector
	  InternalAssignValue<>(KratosComponents< Variable<Vector> >::Get(mvariable_name), Value, rCurrentTime);

        }
	else if( KratosComponents< Variable<array_1d<double,3> > >::Has( mvariable_name ) ) //case of array_1d variable
        {

	  array_1d<double,3> Value;
	  //KratosComponents< Variable<array_1d<double,3> > >,array_1d<double,3> 
	  InternalAssignValue<>(KratosComponents< Variable<array_1d<double,3> > >::Get(mvariable_name), Value, rCurrentTime);

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

        KRATOS_TRY

	if( KratosComponents< Variable<Vector> >::Has( mvariable_name ) ) //case of vector variable
        {

	  Vector Value;
	  noalias(Value) = ZeroVector(3);
	  InternalAssignValue<>(KratosComponents< Variable<Vector> >::Get(mvariable_name), Value);

        }
	else if( KratosComponents< Variable<array_1d<double,3> > >::Has( mvariable_name ) ) //case of array_1d variable
        {

	  array_1d<double,3> Value;
	  Value.clear();
	  InternalAssignValue<>(KratosComponents< Variable<array_1d<double,3> > >::Get(mvariable_name), Value);

        }
	else
	{
	  KRATOS_THROW_ERROR(std::logic_error, "Not able to set the variable. Attempting to set variable:",mvariable_name);
        }
	

        KRATOS_CATCH("");
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
        return "AssignVectorFieldToConditionsProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "AssignVectorFieldToConditionsProcess";
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
    ///@name Protected Operations
    ///@{

    /// Copy constructor.
    AssignVectorFieldToConditionsProcess(AssignVectorFieldToConditionsProcess const& rOther);

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
    std::string mvariable_name;

    array_1d<double,3> mvector_value;

    pybind11::object mPyObject;
    std::string mPyMethodName;

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

    template<class TDataType >
    void CallFunction(const Condition::Pointer& pCondition, const double& time, TDataType& rValue)
    {

      Condition::GeometryType& rConditionGeometry = pCondition->GetGeometry();
      unsigned int size = rConditionGeometry.size();
      double value = 0;
      unsigned int counter = 0;
      
      rValue.resize(size*3,false);
      
      if( mIsSpatialField ){

	double x = 0.0, y = 0.0, z = 0.0;

	for(unsigned int i=0; i<size; i++)
	  {
	    LocalAxesTransform(rConditionGeometry[i].X(), rConditionGeometry[i].Y(), rConditionGeometry[i].Z(), x, y, z);

            value = mPyObject.attr(mPyMethodName.c_str())(x,y,z,time).cast<double>();
            
	    for(unsigned int j=0; j<3; j++)
	      {
		rValue[counter] = value * mvector_value[j];
		counter++;
	      }
	  }
	
      }
      else{

        value = mPyObject.attr(mPyMethodName.c_str())(0.0,0.0,0.0,time).cast<double>();
	for(unsigned int i=0; i<size; i++)
	  {
	    for(unsigned int j=0; j<3; j++)
	      {
		rValue[counter] = value * mvector_value[j];
		counter++;
	      }
	  }

	
      }
      
    }


    template< class TVarType, class TDataType >
    void InternalAssignValue(TVarType& rVariable, TDataType& Value, const double& rTime )
    {
        const int nconditions = mrModelPart.GetMesh().Conditions().size();
	
        if(nconditions != 0)
        {
            ModelPart::ConditionsContainerType::iterator it_begin = mrModelPart.GetMesh().ConditionsBegin();

            //#pragma omp parallel for  //it does not work in parallel
            for(int i = 0; i<nconditions; i++)
            {
                ModelPart::ConditionsContainerType::iterator it = it_begin + i;

		this->CallFunction<TDataType>(*(it.base()), rTime, Value);
		
                it->SetValue(rVariable, Value);
            }
        }

    }


    template< class TVarType, class TDataType >
    void InternalAssignValue(TVarType& rVar, const TDataType value)
    {
      const int nconditions = mrModelPart.GetMesh().Conditions().size();

        if(nconditions != 0)
        {
            ModelPart::ConditionsContainerType::iterator it_begin = mrModelPart.GetMesh().ConditionsBegin();

             #pragma omp parallel for
            for(int i = 0; i<nconditions; i++)
            {
                ModelPart::ConditionsContainerType::iterator it = it_begin + i;

                it->SetValue(rVar, value);
            }
        }
    }

    double CallPythonMethod(PyObject* pPyObject, const char* pPyMethodName,
                           double rX, double rY, double rZ, const double& rTime)
    {
        KRATOS_TRY

        if( PyObject_IsTrue(pPyObject) && PyCallable_Check(pPyObject) ){
	    
          KRATOS_INFO("pPyObject call") << " pPyObject exists and is callable " << std::endl;

          PyObject* pArgs = PyTuple_Pack(4,PyFloat_FromDouble(rX),PyFloat_FromDouble(rY),PyFloat_FromDouble(rZ),
                                           PyFloat_FromDouble(rTime));

          PyObject* pResult = PyObject_CallObject(pPyObject, pArgs);
	    
          Py_DECREF(pArgs);

          if(pResult == NULL){
            if(PyErr_Occurred())
              PyErr_Print();
            KRATOS_ERROR <<" pResult DO NOT exists "<<std::endl;
          }
          else{
            return PyFloat_AsDouble(pResult);
          }
	    
        }
        else{
          if (PyErr_Occurred())
            PyErr_Print();
          KRATOS_ERROR <<" pPyObject DO NOT exists or is NOT callable "<<std::endl;
        }
	
        Py_DECREF(pPyObject);

        KRATOS_CATCH("")
    }
    
    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    AssignVectorFieldToConditionsProcess& operator=(AssignVectorFieldToConditionsProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AssignVectorFieldToConditionsProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignVectorFieldToConditionsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignVectorFieldToConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ASSIGN_VECTOR_FIELD_TO_CONDITIONS_PROCESS_H_INCLUDED  defined
