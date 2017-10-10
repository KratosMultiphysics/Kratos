//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_APPLY_RIGID_BODY_ROTATION_FIELD_TO_NODES_PROCESS_H_INCLUDED )
#define  KRATOS_APPLY_RIGID_BODY_ROTATION_FIELD_TO_NODES_PROCESS_H_INCLUDED



// System includes

// External includes

// Project includes
#include "custom_processes/apply_rigid_body_rotation_to_nodes_process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for assigning a value to scalar variables or array_1d components processes in Kratos.
/** This function assigns a value to a variable belonging to all of the nodes in a given mesh
*/
class ApplyRigidBodyRotationFieldToNodesProcess : public ApplyRigidBodyRotationToNodesProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ApplyRigidBodyRotationFieldToNodesProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyRigidBodyRotationFieldToNodesProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    ApplyRigidBodyRotationFieldToNodesProcess(ModelPart& model_part,
					 PyObject* pPyObject,
					 const char* pPyMethodName,
					 const bool SpatialFieldFunction,
					 Parameters rParameters
				       ) : ApplyRigidBodyRotationToNodesProcess(model_part)
    {
        KRATOS_TRY
	
        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "variable_name": "VARIABLE_NAME",
                "direction" : [],
                "center" : []
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mvariable_name  = rParameters["variable_name"].GetString();

	mpPyObject      =  pPyObject;	
	mpPyMethodName  =  pPyMethodName;

	mIsSpatialField = SpatialFieldFunction;
	
	for( unsigned int i=0; i<3; i++)
	  {
	    mdirection[i] = rParameters["direction"][i].GetDouble();
	    mcenter[i] = rParameters["center"][i].GetDouble();
	  }

        KRATOS_CATCH("");
    }



    /// Destructor.
    virtual ~ApplyRigidBodyRotationFieldToNodesProcess() {}


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


    /// Execute method is used to execute the ApplyRigidBodyRotationFieldToNodesProcess algorithms.
    virtual void Execute() 
    {

        KRATOS_TRY;
	
	if( ! mIsSpatialField ){

	  const ProcessInfo& rCurrentProcessInfo = mr_model_part.GetProcessInfo();
	  const double& rCurrentTime  = rCurrentProcessInfo[TIME];
	  const ProcessInfo& rPreviousProcessInfo = rCurrentProcessInfo.GetPreviousTimeStepInfo();
	  const double& rPreviousTime = rPreviousProcessInfo[TIME];
	  
	  this->CallTimeFunction(rPreviousTime, mprevious_value);
	  this->CallTimeFunction(rCurrentTime, mvalue);
	  
	  ApplyRigidBodyRotationToNodesProcess::Execute();

	}
	else{
	  
	  ApplyRigidBodyRotation();
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
        return "ApplyRigidBodyRotationFieldToNodesProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ApplyRigidBodyRotationFieldToNodesProcess";
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
    ApplyRigidBodyRotationFieldToNodesProcess(ApplyRigidBodyRotationFieldToNodesProcess const& rOther);

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

    PyObject* mpPyObject;  
    const char* mpPyMethodName;
   
    bool mIsSpatialField;

    ///@}
    ///@name Private Operators
    ///@{

    void CallFunction(const Node<3>::Pointer& pNode, const double& time, double& rValue)
    {      
      KRATOS_TRY
	
      if( mIsSpatialField ){

	double x = pNode->X(), y = pNode->Y(), z = pNode->Z();
	   
	rValue = boost::python::call_method<double>(mpPyObject, mpPyMethodName, x, y, z, time);
	
      }
      else{
	
	rValue = boost::python::call_method<double>(mpPyObject, mpPyMethodName, 0.0, 0.0, 0.0, time);
	
      }
      
     KRATOS_CATCH( "" )
      
    }

    void CallTimeFunction(const double& time, double& rValue)
    {
      
      KRATOS_TRY
	
      rValue = boost::python::call_method<double>(mpPyObject, mpPyMethodName, 0.0, 0.0, 0.0, time);
	      
      KRATOS_CATCH( "" )
      
    }

    
    void ApplyRigidBodyRotation()
    {
      KRATOS_TRY
   
      const int nnodes = mr_model_part.GetMesh().Nodes().size();

      if(nnodes != 0)
        {
	  ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh().NodesBegin();

	  Matrix rotation_matrix;
	  Quaternion<double> total_quaternion;
	  array_1d<double,3> radius;
	  array_1d<double,3> distance;
	  array_1d<double,3> rotation;
	  array_1d<double,3> delta_rotation;

	  double norm = norm_2(mdirection);
	  if(norm!=0)
	    mdirection/=norm;


	  //Possible prescribed variables: ROTATION / ANGULAR_VELOCITY / ANGULAR_ACCELERATION
	  double value = 0.0;
	  bool dynamic_angular_velocity = false;
	  bool dynamic_angular_acceleration = false;

	  const ProcessInfo& rCurrentProcessInfo = mr_model_part.GetProcessInfo();
	  const double& rDeltaTime = rCurrentProcessInfo[DELTA_TIME];
	  const double& rCurrentTime = rCurrentProcessInfo[TIME];
	  const ProcessInfo& rPreviousProcessInfo = rCurrentProcessInfo.GetPreviousTimeStepInfo();
	  const double& rPreviousTime = rPreviousProcessInfo[TIME];
	  
	  array_1d<double,3> angular_velocity;
	  angular_velocity.clear();
	  array_1d<double,3> angular_acceleration;
	  angular_acceleration.clear();

	  double time_factor = 0.0;
	  if(mvariable_name == "ROTATION"){

	    time_factor = 1.0;
	    
	  }
	  else if(mvariable_name == "ANGULAR_VELOCITY"){
	    
	    dynamic_angular_velocity = true;
	    time_factor = rDeltaTime;
	      
	  }
	  else if(mvariable_name == "ANGULAR_ACCELERATION"){

	    dynamic_angular_velocity = true;
	    dynamic_angular_acceleration = true;
	    time_factor = rDeltaTime * rDeltaTime;
	  }

	  //#pragma omp parallel for  //it does not work in parallel
	  for(int i = 0; i<nnodes; i++)
            {
	      ModelPart::NodesContainerType::iterator it = it_begin + i;
	      
	      this->CallFunction(*(it.base()), rCurrentTime, value);
      
	      rotation = value * mdirection;
	      
	      rotation *= time_factor;
	      
	      if( dynamic_angular_velocity ){

		this->CallFunction(*(it.base()), rPreviousTime, value);
		delta_rotation  = rotation - time_factor * value * mdirection;	
		
		angular_velocity = delta_rotation / rDeltaTime;
		if( dynamic_angular_acceleration ){
		  angular_acceleration = angular_velocity / rDeltaTime;
		}
	      }
	      
	      //Get rotation matrix
	      total_quaternion = Quaternion<double>::FromRotationVector<array_1d<double,3> >(rotation);

	      distance = it->GetInitialPosition() - mcenter;
	      
	      total_quaternion.ToRotationMatrix(rotation_matrix);

	      noalias(radius) = prod(rotation_matrix, distance);

	      array_1d<double,3>& displacement = it->FastGetSolutionStepValue(DISPLACEMENT);
	      displacement =  radius - distance; //(mcenter + radius) - it->GetInitialPosition();
	      
	      if( dynamic_angular_velocity ){

		//compute the skewsymmmetric tensor of the angular velocity
		BeamMathUtils<double>::VectorToSkewSymmetricTensor(angular_velocity, rotation_matrix);

		//compute the contribution of the angular velocity to the velocity v = Wxr
		array_1d<double,3>& velocity = it->FastGetSolutionStepValue(VELOCITY);
		velocity = prod(rotation_matrix,radius);

		if( dynamic_angular_acceleration ){ 
		  //compute the contribution of the centripetal acceleration ac = Wx(Wxr)
		  array_1d<double,3>& acceleration = it->FastGetSolutionStepValue(ACCELERATION);
		  acceleration = prod(rotation_matrix,velocity);
		  
		  //compute the contribution of the angular acceleration to the acceleration a = Axr
		  BeamMathUtils<double>::VectorToSkewSymmetricTensor(angular_acceleration, rotation_matrix);
		  acceleration += prod(rotation_matrix,radius);	
		}
	      }
	    }
 
	}

      KRATOS_CATCH( "" )
    }


    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    ApplyRigidBodyRotationFieldToNodesProcess& operator=(ApplyRigidBodyRotationFieldToNodesProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class ApplyRigidBodyRotationFieldToNodesProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyRigidBodyRotationFieldToNodesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyRigidBodyRotationFieldToNodesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_APPLY_RIGID_BODY_ROTATION_FIELD_TO_NODES_PROCESS_H_INCLUDED  defined
