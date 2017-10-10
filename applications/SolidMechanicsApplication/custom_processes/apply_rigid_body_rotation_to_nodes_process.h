//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_APPLY_RIGID_BODY_ROTATION_TO_NODES_PROCESS_H_INCLUDED )
#define  KRATOS_APPLY_RIGID_BODY_ROTATION_TO_NODES_PROCESS_H_INCLUDED



// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "utilities/beam_math_utilities.hpp"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for assigning a value to scalar variables or array_1d components processes in Kratos.
/** This function assigns a value to a variable belonging to all of the nodes in a given mesh
*/
class ApplyRigidBodyRotationToNodesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ApplyRigidBodyRotationToNodesProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyRigidBodyRotationToNodesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

 
    ApplyRigidBodyRotationToNodesProcess(ModelPart& model_part) : Process(Flags()) , mr_model_part(model_part) {}


    
    ApplyRigidBodyRotationToNodesProcess(ModelPart& model_part,
				       Parameters rParameters
				       ) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY
			 
        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "variable_name": "VARIABLE_NAME",
                "modulus" : 1.0,
                "direction" : [],
                "center" : []
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mvariable_name = rParameters["variable_name"].GetString();


	if( KratosComponents< Variable<array_1d<double, 3> > >::Has( mvariable_name ) ) //case of array_1d variable
        {

	    mvalue = rParameters["modulus"].GetDouble();

	    for( unsigned int i=0; i<3; i++)
	      {
		mdirection[i] = rParameters["direction"][i].GetDouble();
		mcenter[i] = rParameters["center"][i].GetDouble();
	      }
	}
	else //case of bool variable
        {
	  KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is " << mvariable_name <<std::endl;
	  
	}

        KRATOS_CATCH("");
    }


    /// Destructor.
    virtual ~ApplyRigidBodyRotationToNodesProcess() {}


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


    /// Execute method is used to execute the ApplyRigidBodyRotationToNodesProcess algorithms.
    virtual void Execute() 
    {

        KRATOS_TRY;
		
	this->ApplyRigidBodyRotation();

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
	mprevious_value = mvalue;
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
        return "ApplyRigidBodyRotationToNodesProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ApplyRigidBodyRotationToNodesProcess";
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

    ModelPart& mr_model_part;
    std::string mvariable_name;
    double mvalue;
    double mprevious_value;
    array_1d<double,3> mdirection;
    array_1d<double,3> mcenter;
    
    ///@}
    ///@name Protected Operators
    ///@{

    /// Copy constructor.
    ApplyRigidBodyRotationToNodesProcess(ApplyRigidBodyRotationToNodesProcess const& rOther);

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


    ///@}
    ///@name Private Operators
    ///@{


    void ApplyRigidBodyRotation()
    {
      KRATOS_TRY
   
      const int nnodes = mr_model_part.GetMesh().Nodes().size();

      //std::cout<<" Apply rigid body rotation "<<std::endl;

      if(nnodes != 0)
        {
	  ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh().NodesBegin();

	  Matrix rotation_matrix;
	  array_1d<double,3> radius;
	  array_1d<double,3> distance;
	  
	  double norm = norm_2(mdirection);
	  if(norm!=0)
	    mdirection/=norm;


	  //Possible prescribed variables: ROTATION / ANGULAR_VELOCITY / ANGULAR_ACCELERATION
	  bool dynamic_angular_velocity = false;
	  bool dynamic_angular_acceleration = false;

	  const ProcessInfo& rCurrentProcessInfo = mr_model_part.GetProcessInfo();
	  const double& rDeltaTime = rCurrentProcessInfo[DELTA_TIME];

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
	   
	  array_1d<double,3> delta_rotation = time_factor * (mvalue - mprevious_value) * mdirection;
 
	  array_1d<double,3> rotation = time_factor * mvalue * mdirection;
	  
	  if( dynamic_angular_velocity ){
	    angular_velocity = delta_rotation / rDeltaTime;
	    if( dynamic_angular_acceleration ){
	      angular_acceleration = angular_velocity / rDeltaTime;
	    }
	  }
	  
	  //Get rotation matrix
	  Quaternion<double> total_quaternion = Quaternion<double>::FromRotationVector<array_1d<double,3> >(rotation);

          #pragma omp parallel for private(distance,radius,rotation_matrix)
	  for(int i = 0; i<nnodes; i++)
            {

	      ModelPart::NodesContainerType::iterator it = it_begin + i;
		
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
    ApplyRigidBodyRotationToNodesProcess& operator=(ApplyRigidBodyRotationToNodesProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class ApplyRigidBodyRotationToNodesProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyRigidBodyRotationToNodesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyRigidBodyRotationToNodesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_APPLY_RIGID_BODY_ROTATION_TO_NODES_PROCESS_H_INCLUDED  defined
