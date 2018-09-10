//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_ASSIGN_ROTATION_FIELD_ABOUT_AN_AXIS_TO_NODES_PROCESS_H_INCLUDED)
#define  KRATOS_ASSIGN_ROTATION_FIELD_ABOUT_AN_AXIS_TO_NODES_PROCESS_H_INCLUDED



// System includes

// External includes

// Project includes
#include "custom_processes/assign_rotation_about_an_axis_to_nodes_process.hpp"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for assigning a value to scalar variables or array_1d components processes in Kratos.
/** This function assigns a value to a variable belonging to all of the nodes in a given mesh
*/
class AssignRotationFieldAboutAnAxisToNodesProcess : public AssignRotationAboutAnAxisToNodesProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssignRotationFieldAboutAnAxisToNodesProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignRotationFieldAboutAnAxisToNodesProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    AssignRotationFieldAboutAnAxisToNodesProcess(ModelPart& model_part,
					 pybind11::object& rPyObject,
					 const std::string& rPyMethodName,
					 const bool SpatialFieldFunction,
					 Parameters rParameters
				       ) : AssignRotationAboutAnAxisToNodesProcess(model_part)
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

	mPyObject      =  rPyObject;
	mPyMethodName  =  rPyMethodName;

	mIsSpatialField = SpatialFieldFunction;

	for( unsigned int i=0; i<3; i++)
	  {
	    mdirection[i] = rParameters["direction"][i].GetDouble();
	    mcenter[i] = rParameters["center"][i].GetDouble();
	  }

	double norm = norm_2(mdirection);
	if(norm!=0)
	    mdirection/=norm;

        KRATOS_CATCH("");
    }



    /// Destructor.
    ~AssignRotationFieldAboutAnAxisToNodesProcess() override {}


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


    /// Execute method is used to execute the AssignRotationFieldAboutAnAxisToNodesProcess algorithms.
    void Execute()  override
    {

        KRATOS_TRY;

	if( ! mIsSpatialField ){

	  const ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();
	  const double& rCurrentTime  = rCurrentProcessInfo[TIME];
	  const ProcessInfo& rPreviousProcessInfo = rCurrentProcessInfo.GetPreviousTimeStepInfo();
	  const double& rPreviousTime = rPreviousProcessInfo[TIME];

	  this->CallTimeFunction(rPreviousTime, mprevious_value);
	  this->CallTimeFunction(rCurrentTime, mvalue);

	  AssignRotationAboutAnAxisToNodesProcess::Execute();

	}
	else{

	  AssignRotationAboutAnAxis();
	}

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
        return "AssignRotationFieldAboutAnAxisToNodesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignRotationFieldAboutAnAxisToNodesProcess";
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
    AssignRotationFieldAboutAnAxisToNodesProcess(AssignRotationFieldAboutAnAxisToNodesProcess const& rOther);

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

    pybind11::object mPyObject;
    std::string mPyMethodName;

    bool mIsSpatialField;

    ///@}
    ///@name Private Operators
    ///@{

    void CallFunction(const Node<3>::Pointer& pNode, const double& time, double& rValue)
    {
      KRATOS_TRY

      if( mIsSpatialField ){

	double x = pNode->X(), y = pNode->Y(), z = pNode->Z();
	rValue = mPyObject.attr(mPyMethodName.c_str())(x,y,z,time).cast<double>();

      }
      else{

        rValue = mPyObject.attr(mPyMethodName.c_str())(0.0,0.0,0.0,time).cast<double>();

      }

     KRATOS_CATCH( "" )

    }

    void CallTimeFunction(const double& time, double& rValue)
    {

      KRATOS_TRY

      rValue = mPyObject.attr(mPyMethodName.c_str())(0.0,0.0,0.0,time).cast<double>();

      KRATOS_CATCH( "" )

    }


    void AssignRotationAboutAnAxis()
    {
      KRATOS_TRY

      const int nnodes = mrModelPart.GetMesh().Nodes().size();

      if(nnodes != 0)
        {
	  ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh().NodesBegin();

	  Matrix rotation_matrix;
	  Quaternion<double> total_quaternion;
	  array_1d<double,3> radius;
	  array_1d<double,3> distance;
	  array_1d<double,3> rotation;
	  array_1d<double,3> delta_rotation;


	  //Possible prescribed variables: ROTATION / ANGULAR_VELOCITY / ANGULAR_ACCELERATION
	  double value = 0.0;
	  bool dynamic_angular_velocity = false;
	  bool dynamic_angular_acceleration = false;

	  const ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();
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
    AssignRotationFieldAboutAnAxisToNodesProcess& operator=(AssignRotationFieldAboutAnAxisToNodesProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AssignRotationFieldAboutAnAxisToNodesProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignRotationFieldAboutAnAxisToNodesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignRotationFieldAboutAnAxisToNodesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ASSIGN_ROTATION_FIELD_ABOUT_AN_AXIS_TO_NODES_PROCESS_H_INCLUDED  defined
