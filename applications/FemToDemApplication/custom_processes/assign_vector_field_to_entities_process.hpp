//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_ASSIGN_VECTOR_FIELD_TO_ENTITIES_PROCESS_H_INCLUDED)
#define  KRATOS_ASSIGN_VECTOR_FIELD_TO_ENTITIES_PROCESS_H_INCLUDED



// System includes

// External includes

// Project includes
#include "custom_processes/assign_scalar_field_to_entities_process.hpp"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for assigning a value to scalar variables or array_1d components processes in Kratos.
/** This function assigns a value to a variable belonging to all of the nodes in a given mesh
*/
class AssignVectorFieldToEntitiesProcess : public AssignScalarFieldToEntitiesProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssignVectorFieldToEntitiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignVectorFieldToEntitiesProcess);

    typedef AssignScalarFieldToEntitiesProcess   BaseType;

    ///@}
    ///@name Life Cycle
    ///@{
    AssignVectorFieldToEntitiesProcess(ModelPart& rModelPart,
                                       pybind11::object& rPyObject,
                                       const std::string& rPyMethodName,
                                       const bool SpatialFieldFunction,
                                       Parameters rParameters
                                       ) : BaseType(rModelPart, rPyObject, rPyMethodName, SpatialFieldFunction)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "variable_name": "VARIABLE_NAME",
                "entity_type": "NODES",
                "value" : [0.0, 0.0, 0.0],
                "local_axes" : {},
                "compound_assignment": "direct"
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        this->mvariable_name = rParameters["variable_name"].GetString();

	// Admissible values for local axes, are "empty" or
        //"local_axes" :{
        //    "origin" : [0.0, 0.0, 0.0]
        //    "axes"   : [ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ]
        //    }

	this->mHasLocalOrigin = false;
	if( rParameters["local_axes"].Has("origin") ){
	  this->mHasLocalOrigin = true;
	  this->mLocalOrigin.resize(3,false);
	  for( unsigned int i=0; i<3; i++)
	    this->mLocalOrigin[i] = rParameters["local_axes"]["origin"][i].GetDouble();
	}

	this->mHasLocalAxes = false;
	if( rParameters["local_axes"].Has("axes") ){
	  this->mHasLocalAxes = true;
	  this->mTransformationMatrix.resize(3,3,false);
	  for( unsigned int i=0; i<3; i++)
	    for( unsigned int j=0; j<3; j++)
	      this->mTransformationMatrix(i,j) = rParameters["local_axes"]["axes"][i][j].GetDouble();
	}

        if( rParameters["entity_type"].GetString() == "NODES" ){
          this->mEntity = EntityType::NODES;
        }
        else if(  rParameters["entity_type"].GetString() == "CONDITIONS" ){
          this->mEntity = EntityType::CONDITIONS;
        }
        else{
          KRATOS_ERROR <<" Entity type "<< rParameters["entity_type"].GetString() <<" is not supported "<<std::endl;
        }

        if( this->mEntity == EntityType::CONDITIONS ){

          if(KratosComponents< Variable<Vector> >::Has(this->mvariable_name) == false) //case of vector variable
          {
            KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is " << mvariable_name << std::endl;
          }
          else if( KratosComponents< Variable<array_1d<double,3> > >::Has( this->mvariable_name ) ) //case of array_1d variable
          {
            KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is " << mvariable_name << std::endl;
          }

        }
        else{
          KRATOS_ERROR << " Assignment to " << rParameters["entity_type"].GetString() << " not implemented "<< std::endl;
        }

        mvector_value[0] = rParameters["value"][0].GetDouble();
        mvector_value[1] = rParameters["value"][1].GetDouble();
        mvector_value[2] = rParameters["value"][2].GetDouble();


       this->SetAssignmentType(rParameters["compound_assignment"].GetString(), this->mAssignment);

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~AssignVectorFieldToEntitiesProcess() override {}


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


    /// Execute method is used to execute the AssignVectorFieldToEntitiesProcess algorithms.
    void Execute() override
    {

        KRATOS_TRY

	ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();

	const double& rCurrentTime = rCurrentProcessInfo[TIME];

 	if( KratosComponents< Variable<Vector> >::Has( this->mvariable_name ) ) //case of vector variable
        {

	  Vector Value;
	  //KratosComponents< Variable<Vector> >, Vector
	  AssignValueToConditions<>(KratosComponents< Variable<Vector> >::Get(this->mvariable_name), Value, rCurrentTime);

        }
	else if( KratosComponents< Variable<array_1d<double,3> > >::Has( this->mvariable_name ) ) //case of array_1d variable
        {

	  array_1d<double,3> Value;
	  //KratosComponents< Variable<array_1d<double,3> > >,array_1d<double,3>
	  AssignValueToConditions<>(KratosComponents< Variable<array_1d<double,3> > >::Get(this->mvariable_name), Value, rCurrentTime);

        }
	else
	{
	  KRATOS_ERROR << "Not able to set the variable. Attempting to set variable:" << this->mvariable_name << std::endl;
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

        KRATOS_TRY

        if( this->mEntity == EntityType::CONDITIONS ){

          mAssignment = AssignmentType::DIRECT;
          if( KratosComponents< Variable<Vector> >::Has( this->mvariable_name ) ) //case of vector variable
          {

            Vector Value;
            noalias(Value) = ZeroVector(3);
            BaseType::AssignValueToConditions<>(KratosComponents< Variable<Vector> >::Get(this->mvariable_name), Value);

          }
          else if( KratosComponents< Variable<array_1d<double,3> > >::Has( this->mvariable_name ) ) //case of array_1d variable
          {

            array_1d<double,3> Value;
            Value.clear();
            BaseType::AssignValueToConditions<>(KratosComponents< Variable<array_1d<double,3> > >::Get(this->mvariable_name), Value);

          }
          else
          {
            KRATOS_ERROR << "Not able to set the variable. Attempting to set variable:" << this->mvariable_name << std::endl;
          }

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
    std::string Info() const override
    {
        return "AssignVectorFieldToEntitiesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignVectorFieldToEntitiesProcess";
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
    ///@name Protected Operations
    ///@{

    /// Copy constructor.
    AssignVectorFieldToEntitiesProcess(AssignVectorFieldToEntitiesProcess const& rOther);

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

    array_1d<double,3> mvector_value;

    ///@}
    ///@name Private Operators
    ///@{

    template<class TDataType>
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
	    this->LocalAxesTransform(rConditionGeometry[i].X(), rConditionGeometry[i].Y(), rConditionGeometry[i].Z(), x, y, z);

            value = mPyObject.attr(this->mPyMethodName.c_str())(x,y,z,time).cast<double>();

	    for(unsigned int j=0; j<3; j++)
	      {
		rValue[counter] = value * mvector_value[j];
		counter++;
	      }
	  }

      }
      else{

        value = mPyObject.attr(this->mPyMethodName.c_str())(0.0,0.0,0.0,time).cast<double>();
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
    void AssignValueToConditions(TVarType& rVariable, TDataType& Value, const double& rTime )
    {

        typedef void (BaseType::*AssignmentMethodPointer) (ModelPart::ConditionType&, const TVarType&, const TDataType&);

        AssignmentMethodPointer AssignmentMethod = this->GetAssignmentMethod<AssignmentMethodPointer>();

        const int nconditions = mrModelPart.GetMesh().Conditions().size();

        if(nconditions != 0)
        {
            ModelPart::ConditionsContainerType::iterator it_begin = mrModelPart.GetMesh().ConditionsBegin();

            //#pragma omp parallel for  //it does not work in parallel
            for(int i = 0; i<nconditions; i++)
            {
                ModelPart::ConditionsContainerType::iterator it = it_begin + i;

		this->CallFunction<TDataType>(*(it.base()), rTime, Value);

                (this->*AssignmentMethod)(*it, rVariable, Value);
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
    AssignVectorFieldToEntitiesProcess& operator=(AssignVectorFieldToEntitiesProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AssignVectorFieldToEntitiesProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignVectorFieldToEntitiesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignVectorFieldToEntitiesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ASSIGN_VECTOR_FIELD_TO_ENTITIES_PROCESS_H_INCLUDED  defined
