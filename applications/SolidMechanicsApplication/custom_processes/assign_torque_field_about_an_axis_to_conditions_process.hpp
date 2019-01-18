//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:             January 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_ASSIGN_TORQUE_FIELD_ABOUT_AN_AXIS_TO_CONDITIONS_PROCESS_H_INCLUDED)
#define  KRATOS_ASSIGN_TORQUE_FIELD_ABOUT_AN_AXIS_TO_CONDITIONS_PROCESS_H_INCLUDED



// System includes

// External includes

// Project includes
#include "custom_processes/assign_torque_about_an_axis_to_conditions_process.hpp"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for assigning a value to scalar variables or array_1d components processes in Kratos.
/** This function assigns a value to a variable belonging to all of the nodes in a given mesh
*/
class AssignTorqueFieldAboutAnAxisToConditionsProcess : public AssignTorqueAboutAnAxisToConditionsProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssignTorqueFieldAboutAnAxisToConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignTorqueFieldAboutAnAxisToConditionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    AssignTorqueFieldAboutAnAxisToConditionsProcess(ModelPart& model_part,
						    pybind11::object& pPyObject,
						    const std::string& pPyMethodName,
						    const bool SpatialFieldFunction,
						    Parameters rParameters
	                                         ) : AssignTorqueAboutAnAxisToConditionsProcess(model_part)
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

	    mPyObject      =  pPyObject;
	    mPyMethodName  =  pPyMethodName;

	    mIsSpatialField = SpatialFieldFunction;


	    for( unsigned int i=0; i<3; i++)
	    {
		mdirection[i] = rParameters["direction"][i].GetDouble();
		mcenter[i] = rParameters["center"][i].GetDouble();
	    }

	    double norm = norm_2(mdirection);
	    if(norm!=0)
	    mdirection/=norm;
	}
	else //case of other variable type
        {
	  KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is " << mvariable_name <<std::endl;
	}

        KRATOS_CATCH("");
    }


    /// Destructor.
    ~AssignTorqueFieldAboutAnAxisToConditionsProcess() override {}


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


    /// Execute method is used to execute the AssignTorqueFieldAboutAnAxisToConditionsProcess algorithms.
    void Execute()  override
    {

        KRATOS_TRY;


	if( ! mIsSpatialField ){

	  const ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();
	  const double& rCurrentTime  = rCurrentProcessInfo[TIME];

	  this->CallTimeFunction(rCurrentTime, mvalue);

	  AssignTorqueAboutAnAxisToConditionsProcess::Execute();

	}
	else //no spatial fields accepted (it have to be implemented if needed
        {
	  KRATOS_ERROR << "trying to set an spatial field....not implemented" << mvariable_name <<std::endl;
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
	AssignTorqueAboutAnAxisToConditionsProcess::ExecuteFinalize();
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
        return "AssignTorqueFieldAboutAnAxisToConditionsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignTorqueFieldAboutAnAxisToConditionsProcess";
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

    pybind11::object mPyObject;
    std::string mPyMethodName;

    bool mIsSpatialField;

    ///@}
    ///@name Protected Operators
    ///@{

    /// Copy constructor.
    AssignTorqueFieldAboutAnAxisToConditionsProcess(AssignTorqueFieldAboutAnAxisToConditionsProcess const& rOther);

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


    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    AssignTorqueFieldAboutAnAxisToConditionsProcess& operator=(AssignTorqueFieldAboutAnAxisToConditionsProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AssignTorqueFieldAboutAnAxisToConditionsProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignTorqueFieldAboutAnAxisToConditionsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignTorqueFieldAboutAnAxisToConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ASSIGN_TORQUE_FIELD_ABOUT_AN_AXIS_TO_CONDITIONS_PROCESS_H_INCLUDED  defined
