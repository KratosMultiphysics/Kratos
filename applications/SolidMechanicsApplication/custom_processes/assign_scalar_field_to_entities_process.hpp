//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_ASSIGN_SCALAR_FIELD_TO_ENTITIES_PROCESS_H_INCLUDED)
#define  KRATOS_ASSIGN_SCALAR_FIELD_TO_ENTITIES_PROCESS_H_INCLUDED



// System includes

// External includes

// Project includes
#include "custom_processes/assign_scalar_variable_to_entities_process.hpp"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for assigning a value to scalar variables or array_1d components processes in Kratos.
/** This function assigns a value to a variable belonging to all of the nodes in a given mesh
*/
class AssignScalarFieldToEntitiesProcess : public AssignScalarVariableToEntitiesProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssignScalarFieldToEntitiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignScalarFieldToEntitiesProcess);

    typedef AssignScalarVariableToEntitiesProcess  BaseType;

    ///@}
    ///@name Life Cycle
    ///@{
    AssignScalarFieldToEntitiesProcess(ModelPart& rModelPart,
                                       pybind11::object& pPyObject,
                                       const std::string& pPyMethodName,
                                       const bool SpatialFieldFunction,
                                       Parameters rParameters
                                       ) : BaseType(rModelPart), mPyObject(pPyObject), mPyMethodName(pPyMethodName), mIsSpatialField(SpatialFieldFunction)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "variable_name": "VARIABLE_NAME",
                "entity_type": "NODES",
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

        if( rParameters["entity_type"].GetString() == "NODES" ){
          this->mEntity = EntityType::NODES;
        }
        else if(  rParameters["entity_type"].GetString() == "CONDITIONS" ){
          this->mEntity = EntityType::CONDITIONS;
        }
        else if(  rParameters["entity_type"].GetString() == "ELEMENTS" ){
          this->mEntity = EntityType::ELEMENTS;
        }
        else{
          KRATOS_ERROR <<" Entity type "<< rParameters["entity_type"].GetString() <<" is not supported "<<std::endl;
        }


        if( this->mEntity == EntityType::NODES ){

          if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(this->mvariable_name) ) //case of component variable
          {
            typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
            component_type var_component = KratosComponents< component_type >::Get(this->mvariable_name);

            if( rModelPart.GetNodalSolutionStepVariablesList().Has( var_component.GetSourceVariable() ) == false )
            {
              KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is " << this->mvariable_name << std::endl;
            }

          }
          else if( KratosComponents< Variable<double> >::Has( this->mvariable_name ) ) //case of double variable
          {
            if( rModelPart.GetNodalSolutionStepVariablesList().Has( KratosComponents< Variable<double> >::Get( this->mvariable_name ) ) == false )
            {
              KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is " << this->mvariable_name << std::endl;
            }

          }
          else
          {
            KRATOS_ERROR << "Not able to set the variable type/name. Attempting to set variable:" << this->mvariable_name << std::endl;
          }
        }
        else if( this->mEntity == EntityType::CONDITIONS || this->mEntity == EntityType::ELEMENTS ){

          if( KratosComponents< Variable<Vector> >::Has( this->mvariable_name ) == false ) //case of double variable
          {
            KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is " << this->mvariable_name << std::endl;
          }

        }
        else{
          KRATOS_ERROR << " Assignment to " << rParameters["entity_type"].GetString() << " not implemented "<< std::endl;
        }

        this->SetAssignmentType(rParameters["compound_assignment"].GetString(), mAssignment);

        KRATOS_CATCH("")
    }

    // Constructor.
    AssignScalarFieldToEntitiesProcess(ModelPart& rModelPart,
                                       pybind11::object& pPyObject,
                                       const std::string& pPyMethodName,
                                       const bool SpatialFieldFunction
                                       ) : BaseType(rModelPart), mPyObject(pPyObject), mPyMethodName(pPyMethodName), mIsSpatialField(SpatialFieldFunction)
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /// Destructor.
    ~AssignScalarFieldToEntitiesProcess() override {}


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


    /// Execute method is used to execute the AssignScalarFieldToEntitiesProcess algorithms.
    void Execute()  override
    {

        KRATOS_TRY

        ProcessInfo& rCurrentProcessInfo = this->mrModelPart.GetProcessInfo();

	const double& rCurrentTime = rCurrentProcessInfo[TIME];

        if( this->mEntity == EntityType::NODES || this->mEntity == EntityType::CONDITIONS ){

          if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(this->mvariable_name) ) //case of component variable
          {
            typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
            component_type var_component = KratosComponents< component_type >::Get(this->mvariable_name);
            AssignValueToNodes<component_type>(var_component, rCurrentTime);
          }
          else if( KratosComponents< Variable<double> >::Has( this->mvariable_name ) ) //case of double variable
          {
            AssignValueToNodes<>(KratosComponents< Variable<double> >::Get(this->mvariable_name), rCurrentTime);
          }
          else if( KratosComponents< Variable<Vector> >::Has( this->mvariable_name ) ) //case of vector variable
          {
            AssignValueToConditions<>(KratosComponents< Variable<Vector> >::Get(this->mvariable_name), rCurrentTime);
          }
          else
          {
            KRATOS_ERROR << "Not able to set the variable. Attempting to set variable:" << this->mvariable_name << std::endl;
          }

        }
        else if( this->mEntity == EntityType::ELEMENTS ){

          if( KratosComponents< Variable<Vector> >::Has( this->mvariable_name ) ) //case of vector variable
          {
            AssignValueToElements<>(KratosComponents< Variable<Vector> >::Get(this->mvariable_name), rCurrentTime);
          }
          else
          {
            KRATOS_ERROR << "Not able to set the variable. Attempting to set variable:" << this->mvariable_name << std::endl;
          }

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

      	if( KratosComponents< Variable<Vector> >::Has( this->mvariable_name ) ) //case of vector variable
        {
          mAssignment = AssignmentType::DIRECT;
	  Vector Value(3);
	  noalias(Value) = ZeroVector(3);
	  AssignValueToConditions<>(KratosComponents< Variable<Vector> >::Get(this->mvariable_name), Value);
        }
        else
        {
          KRATOS_ERROR << "Not able to set the variable. Attempting to set variable:" << this->mvariable_name << std::endl;
        }
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
    std::string Info() const override
    {
        return "AssignScalarFieldToEntitiesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignScalarFieldToEntitiesProcess";
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

    Vector mLocalOrigin;
    Matrix mTransformationMatrix;

    bool mIsSpatialField;

    bool mHasLocalOrigin;
    bool mHasLocalAxes;

    ///@}
    ///@name Protected Operators
    ///@{

    /// Copy constructor.
    AssignScalarFieldToEntitiesProcess(AssignScalarFieldToEntitiesProcess const& rOther);

    ///@}
    ///@name Protected Operations
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
 	rValue = mPyObject.attr(mPyMethodName.c_str())(x,y,z,time).cast<double>();

      }
      else{

       rValue = mPyObject.attr(mPyMethodName.c_str())(0.0,0.0,0.0,time).cast<double>();
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
            rValue[i] = mPyObject.attr(mPyMethodName.c_str())(x,y,z,time).cast<double>();
	  }

      }
      else{

        double value = mPyObject.attr(mPyMethodName.c_str())(0.0,0.0,0.0,time).cast<double>();
	for(unsigned int i=0; i<size; i++)
	  {
	    rValue[i] = value;
	  }

      }

    }


    void CallFunction(const Element::Pointer& pElement, const double& time, Vector& rValue)
    {

      Element::GeometryType& rElementGeometry = pElement->GetGeometry();
      unsigned int size = rElementGeometry.size();

      rValue.resize(size,false);

      if( mIsSpatialField ){

	double x = 0, y = 0, z = 0;

	for(unsigned int i=0; i<size; i++)
	  {
	    LocalAxesTransform(rElementGeometry[i].X(), rElementGeometry[i].Y(), rElementGeometry[i].Z(), x, y, z);
            rValue[i] = mPyObject.attr(mPyMethodName.c_str())(x,y,z,time).cast<double>();
	  }

      }
      else{

        double value = mPyObject.attr(mPyMethodName.c_str())(0.0,0.0,0.0,time).cast<double>();
	for(unsigned int i=0; i<size; i++)
	  {
	    rValue[i] = value;
	  }

      }

    }

    template< class TVarType >
    void AssignValueToNodes(TVarType& rVariable, const double& rTime)
    {
      if( this->mEntity == EntityType::NODES ){

        typedef void (BaseType::*AssignmentMethodPointer) (ModelPart::NodeType&, const TVarType&, const double&);

        AssignmentMethodPointer AssignmentMethod = this->GetAssignmentMethod<AssignmentMethodPointer>();

        const int nnodes = this->mrModelPart.GetMesh().Nodes().size();

	double Value = 0;

        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = this->mrModelPart.GetMesh().NodesBegin();

            //#pragma omp parallel for  //it does not work in parallel
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

		this->CallFunction(*(it.base()), rTime, Value);

		(this->*AssignmentMethod)(*it, rVariable, Value);
            }
        }

      }
    }

    template< class TVarType >
    void AssignValueToConditions(TVarType& rVariable, const double& rTime)
    {
      if( this->mEntity == EntityType::CONDITIONS ){

        typedef void (BaseType::*AssignmentMethodPointer) (ModelPart::ConditionType&, const Variable<Vector>&, const Vector&);

        AssignmentMethodPointer AssignmentMethod = this->GetAssignmentMethod<AssignmentMethodPointer>();

        const int nconditions = this->mrModelPart.GetMesh().Conditions().size();

	Vector Value;

        if(nconditions != 0)
        {
          ModelPart::ConditionsContainerType::iterator it_begin = this->mrModelPart.GetMesh().ConditionsBegin();

          //#pragma omp parallel for //it does not work in parallel
          for(int i = 0; i<nconditions; i++)
          {
            ModelPart::ConditionsContainerType::iterator it = it_begin + i;

            this->CallFunction(*(it.base()), rTime, Value);

            (this->*AssignmentMethod)(*it, rVariable, Value);
          }
        }

      }
    }

    template< class TVarType, class TDataType >
    void AssignValueToConditions(TVarType& rVariable, const TDataType Value)
    {

      if( this->mEntity == EntityType::CONDITIONS ){

        typedef void (BaseType::*AssignmentMethodPointer) (ModelPart::ConditionType&, const TVarType&, const TDataType&);
        AssignmentMethodPointer AssignmentMethod = this->GetAssignmentMethod<AssignmentMethodPointer>();

        const int nconditions = this->mrModelPart.GetMesh().Conditions().size();

        if(nconditions != 0)
        {
          ModelPart::ConditionsContainerType::iterator it_begin = this->mrModelPart.GetMesh().ConditionsBegin();

          #pragma omp parallel for
          for(int i = 0; i<nconditions; i++)
          {
            ModelPart::ConditionsContainerType::iterator it = it_begin + i;

            (this->*AssignmentMethod)(*it, rVariable, Value);
          }
        }

      }

    }


    template< class TVarType >
    void AssignValueToElements(TVarType& rVariable, const double& rTime)
    {
      if( this->mEntity == EntityType::ELEMENTS ){

        typedef void (BaseType::*AssignmentMethodPointer) (ModelPart::ElementType&, const Variable<Vector>&, const Vector&);

        AssignmentMethodPointer AssignmentMethod = this->GetAssignmentMethod<AssignmentMethodPointer>();

        const int nelements = this->mrModelPart.GetMesh().Elements().size();

	Vector Value;

        if(nelements != 0)
        {
          ModelPart::ElementsContainerType::iterator it_begin = this->mrModelPart.GetMesh().ElementsBegin();

          //#pragma omp parallel for //it does not work in parallel
          for(int i = 0; i<nelements; i++)
          {
            ModelPart::ElementsContainerType::iterator it = it_begin + i;

            this->CallFunction(*(it.base()), rTime, Value);

            (this->*AssignmentMethod)(*it, rVariable, Value);
          }
        }

      }
    }

    // template< class TVarType, class TDataType >
    // void AssignValueToElements(TVarType& rVariable, const TDataType Value)
    // {

    //   if( this->mEntity == ELEMENTS ){

    //     const int nelements = this->mrModelPart.GetMesh().Elements().size();

    //     if(nelements != 0)
    //     {
    //       ModelPart::ElementsContainerType::iterator it_begin = this->mrModelPart.GetMesh().ElementsBegin();

    //       #pragma omp parallel for
    //       for(int i = 0; i<nelements; i++)
    //       {
    //         ModelPart::ElementsContainerType::iterator it = it_begin + i;

    //         it->SetValue(rVariable, Value);
    //       }
    //     }

    //   }

    // }

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
    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    AssignScalarFieldToEntitiesProcess& operator=(AssignScalarFieldToEntitiesProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AssignScalarFieldToEntitiesProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignScalarFieldToEntitiesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignScalarFieldToEntitiesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ASSIGN_SCALAR_FIELD_TO_ENTITIES_PROCESS_H_INCLUDED  defined
