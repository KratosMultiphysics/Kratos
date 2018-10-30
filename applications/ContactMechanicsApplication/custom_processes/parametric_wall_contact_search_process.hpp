//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_PARAMETRIC_WALL_CONTACT_SEARCH_PROCESS_H_INCLUDED )
#define  KRATOS_PARAMETRIC_WALL_CONTACT_SEARCH_PROCESS_H_INCLUDED


// External includes

// System includes
#ifdef _OPENMP
#include <omp.h>
#endif

// Project includes
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

// Contact Point Conditions
#include "custom_conditions/rigid_contact/point_rigid_contact_penalty_2D_condition.hpp"
#include "custom_conditions/rigid_contact/axisym_point_rigid_contact_penalty_2D_condition.hpp"

#include "custom_conditions/rigid_contact/EP_point_rigid_contact_penalty_3D_condition.hpp"
#include "custom_conditions/rigid_contact/EP_point_rigid_contact_penalty_2D_condition.hpp"
#include "custom_conditions/rigid_contact/EP_point_rigid_contact_penalty_wP_3D_condition.hpp"
#include "custom_conditions/rigid_contact/EP_axisym_point_rigid_contact_penalty_2D_condition.hpp"


// #include "custom_conditions/rigid_contact/axisym_point_rigid_contact_penalty_water_2D_condition.hpp"
// #include "custom_conditions/beam_contact/beam_point_rigid_contact_penalty_3D_condition.hpp"
// #include "custom_conditions/beam_contact/beam_point_rigid_contact_LM_3D_condition.hpp"
// #include "custom_conditions/rigid_contact/rigid_body_point_rigid_contact_condition.hpp"

//#include "custom_friction/friction_law.hpp"

#include "contact_mechanics_application_variables.h"


namespace Kratos
{
  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{


  ///@name Kratos Classes
  ///@{

  /// The base class for all processes in Kratos.
  /** The process is the base class for all processes and defines a simple interface for them.
      Execute method is used to execute the Process algorithms. While the parameters of this method
      can be very different from one Process to other there is no way to create enough overridden
      versions of it. For this reason this method takes no argument and all Process parameters must
      be passed at construction time. The reason is that each constructor can take different set of
      argument without any dependency to other processes or the base Process class.
  */
  class ParametricWallContactSearchProcess
    : public Process
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( ParametricWallContactSearchProcess );

    typedef ModelPart::NodeType                   NodeType;
    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;
    typedef Point2D<ModelPart::NodeType>       Point2DType;
    typedef Point3D<ModelPart::NodeType>       Point3DType;
    typedef FrictionLaw::Pointer           FrictionLawType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ParametricWallContactSearchProcess(ModelPart& rMainModelPart): mrMainModelPart(rMainModelPart) {}


    ParametricWallContactSearchProcess( ModelPart& rMainModelPart,
					std::string rSubModelPartName,
					SpatialBoundingBox::Pointer pParametricWall,
					Parameters CustomParameters)
      : mrMainModelPart(rMainModelPart)
    {
      KRATOS_TRY

      mEchoLevel = 1;

      mpParametricWall = pParametricWall;

      Parameters DefaultParameters( R"(
            {
                   "contact_condition_type": "PointContactCondition2D1N",
                   "hydraulic_condition_type": "HydraulicPointContactCondition2D1N",
                   "kratos_module": "KratosMultiphysics.ContactMechanicsApplication",
                   "friction_law_type": "FrictionLaw",
                   "variables_of_properties":{
                     "FRICTION_ACTIVE": false,
                     "MU_STATIC": 0.3,
                     "MU_DYNAMIC": 0.2,
                     "PENALTY_PARAMETER": 1000,
                     "TANGENTIAL_PENALTY_RATIO": 0.1,
                     "TAU_STAB": 1
                   }

            }  )" );


      //validate against defaults -- this also ensures no type mismatch
      CustomParameters.ValidateAndAssignDefaults(DefaultParameters);

      //create condition prototype:
      mpConditionType = CreateConditionPrototype( CustomParameters );

      //std::cout<<" ConditionPointer "<<*mpConditionType<<std::endl;

      if(mpConditionType == NULL)
	std::cout<<" ERROR:: PROTOTYPE CONTACT WALL CONDITION NOT DEFINED PROPERLY "<<std::endl;

      //contact model part
      mContactModelPartName = rSubModelPartName;

      KRATOS_CATCH(" ")

    }

    /// Destructor.
    virtual ~ParametricWallContactSearchProcess() {}


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


    /// Execute method is used to execute the Process algorithms.
    void Execute() override
    {
      KRATOS_TRY

      if( mEchoLevel > 1 )
	std::cout<<"  [PARAMETRIC_CONTACT_SEARCH]:: -START- "<<std::endl;

      //update parametric wall position
      const ProcessInfo& rCurrentProcessInfo= mrMainModelPart.GetProcessInfo();
      const double Time = rCurrentProcessInfo[TIME];

      mpParametricWall->UpdateBoxPosition( Time );

      //reset CONTACT flag to all modelpart nodes
      ClearContactFlags();

      //search contact conditions
      SearchContactConditions();

      //create contact conditions
      this->CreateContactConditions();

      if( mEchoLevel > 1 )
	std::cout<<"  [PARAMETRIC_CONTACT_SEARCH]:: -END- "<<std::endl;

      KRATOS_CATCH( "" )
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
      KRATOS_TRY

      // set BOUNDARY flag to nodes from structural elements
      for(ModelPart::ElementsContainerType::iterator ie = mrMainModelPart.ElementsBegin(); ie!=mrMainModelPart.ElementsEnd(); ie++)
	{

	  if( ie->GetGeometry().size() == 2 ){
	    for( unsigned int i=0; i<ie->GetGeometry().size(); i++ )
	      {
		ie->GetGeometry()[i].Set(BOUNDARY,true);
	      }
	  }

	}

      // set BOUNDARY flag to nodes from RIGID sub model parts
      // for(ModelPart::SubModelPartIterator i_mp= mrMainModelPart.SubModelPartsBegin() ; i_mp!=mrMainModelPart.SubModelPartsEnd(); i_mp++)
      // 	{
      // 	  if( i_mp->Is(RIGID) ){

      // 	    for(ModelPart::ElementsContainerType::iterator ie = i_mp->ElementsBegin(); ie!=i_mp->ElementsEnd(); ie++)
      // 	      {

      // 		for( unsigned int i=0; i<ie->GetGeometry().size(); i++ )
      // 		  {
      // 		    ie->GetGeometry()[i].Set(BOUNDARY,true);
      // 		  }
      // 	      }
      // 	    for(ModelPart::ConditionsContainerType::iterator ic = i_mp->ConditionsBegin(); ic!=i_mp->ConditionsEnd(); ic++)
      // 	      {

      // 		for( unsigned int i=0; i<ic->GetGeometry().size(); i++ )
      // 		  {
      // 		    ic->GetGeometry()[i].Set(BOUNDARY,true);
      // 		  }
      // 	      }
      // 	  }

      // 	}

      KRATOS_CATCH( "" )
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
      return "ParametricWallContactSearchProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "ParametricWallContactSearchProcess";
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

    ModelPart&  mrMainModelPart;

    SpatialBoundingBox::Pointer  mpParametricWall;

    ConditionType::Pointer  mpConditionType;

    PropertiesType::Pointer mpProperties;

    std::string  mContactModelPartName;

    int  mEchoLevel;

    ///@}
    ///@name Protected Operators
    ///@{

    virtual void CreateContactConditions()
    {
      KRATOS_TRY

      ProcessInfo& rCurrentProcessInfo= mrMainModelPart.GetProcessInfo();
      double Dimension = rCurrentProcessInfo[SPACE_DIMENSION];

      ModelPart::ConditionsContainerType ContactConditions;

      ModelPart& rContactModelPart = mrMainModelPart.GetSubModelPart(mContactModelPartName);

      if( mEchoLevel > 1 ){
	std::cout<<"    ["<<rContactModelPart.Name()<<" :: CONDITIONS [OLD:"<<rContactModelPart.NumberOfConditions();
      }

      unsigned int id = mrMainModelPart.Conditions().back().Id() + 1;

      ModelPart::NodesContainerType& rNodes = mrMainModelPart.Nodes();

      // create contact condition for rigid and deformable bodies
      for(ModelPart::NodesContainerType::ptr_iterator nd = rNodes.ptr_begin(); nd != rNodes.ptr_end(); ++nd)
	{
	  if( (*nd)->Is(BOUNDARY) && (*nd)->Is(CONTACT) ){

	    ConditionType::Pointer pCondition;


	    if( (*nd)->Is(RIGID) ){  //rigid wall contacting with a rigid body

	      GeometryType::Pointer pGeometry;
	      if( Dimension == 2 )
		pGeometry = Kratos::make_shared<Point2DType>(*nd);
	      else if( Dimension == 3 )
		pGeometry = Kratos::make_shared<Point3DType>(*nd);

	      //pCondition= Kratos::make_shared<RigidBodyPointRigidContactCondition>(id, pGeometry, mpProperties, mpParametricWall);

	      ContactConditions.push_back(pCondition);

	    }
	    else{ //rigid wall contacting with a deformable body

	      Condition::NodesArrayType  pConditionNode;
	      pConditionNode.push_back( (*nd) );

	      ConditionType::Pointer pConditionType = FindPointCondition(rContactModelPart, (*nd) );

	      pCondition = pConditionType->Clone(id, pConditionNode);

	      pCondition->Set(CONTACT);

	      ContactConditions.push_back(pCondition);

	    }

	    id +=1;
	  }

	}


      rContactModelPart.Conditions().swap(ContactConditions);


      if( mEchoLevel > 1 ){
	std::cout<<" / NEW:"<<rContactModelPart.NumberOfConditions()<<"] "<<std::endl;
      }

      std::string ModelPartName;

      //Add contact conditions to computing domain
      for(ModelPart::SubModelPartIterator i_mp= mrMainModelPart.SubModelPartsBegin(); i_mp!=mrMainModelPart.SubModelPartsEnd(); i_mp++)
	{
	  if(i_mp->Is(SOLID) && i_mp->Is(ACTIVE))
	    ModelPartName = i_mp->Name();
	}

      AddContactConditions(rContactModelPart, mrMainModelPart.GetSubModelPart(ModelPartName));

      //Add contact conditions to  main domain( with AddCondition are already added )
      //AddContactConditions(rContactModelPart, mrMainModelPart);

      if( mEchoLevel >= 1 )
	std::cout<<"  [CONTACT CANDIDATES : "<<rContactModelPart.NumberOfConditions()<<"] ("<<mContactModelPartName<<") "<<std::endl;

      KRATOS_CATCH( "" )

    }


    //**************************************************************************
    //**************************************************************************

    void AddContactConditions(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart)
    {

      KRATOS_TRY

      //*******************************************************************
      //adding contact conditions
      //

      if( mEchoLevel > 1 ){
	std::cout<<"    ["<<rDestinationModelPart.Name()<<" :: CONDITIONS [OLD:"<<rDestinationModelPart.NumberOfConditions();
      }

      for(ModelPart::ConditionsContainerType::iterator ic = rOriginModelPart.ConditionsBegin(); ic!= rOriginModelPart.ConditionsEnd(); ic++)
	{

	  if(ic->Is(CONTACT))
	    rDestinationModelPart.AddCondition(*(ic.base()));

	}

      if( mEchoLevel > 1 ){
	std::cout<<" / NEW:"<<rDestinationModelPart.NumberOfConditions()<<"] "<<std::endl;
      }

      KRATOS_CATCH( "" )

    }

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

    //**************************************************************************
    //**************************************************************************
    ConditionType::Pointer CreateConditionPrototype( Parameters& CustomParameters )
    {
      KRATOS_TRY

      ProcessInfo& rCurrentProcessInfo= mrMainModelPart.GetProcessInfo();
      double Dimension = rCurrentProcessInfo[SPACE_DIMENSION];

      //create properties prototype for the contact conditions
      unsigned int NumberOfProperties = mrMainModelPart.NumberOfProperties();

      mpProperties = Kratos::make_shared<PropertiesType>(NumberOfProperties);

      // Friction Law is not a Kratos Component
      // std::string FrictionLawName = CustomParameters["friction_law_type"].GetString();
      // const FrictionLawType pCloneFrictionLaw = KratosComponents<FrictionLawType>::Get(FrictionLawName);
      // mpProperties->SetValue(FRICTION_LAW, pCloneFrictionLaw->Clone() );

      Parameters CustomProperties = CustomParameters["variables_of_properties"];

      mpProperties->SetValue(FRICTION_ACTIVE, CustomProperties["FRICTION_ACTIVE"].GetBool());
      mpProperties->SetValue(MU_STATIC, CustomProperties["MU_STATIC"].GetDouble());
      mpProperties->SetValue(MU_DYNAMIC, CustomProperties["MU_DYNAMIC"].GetDouble());
      mpProperties->SetValue(PENALTY_PARAMETER, CustomProperties["PENALTY_PARAMETER"].GetDouble());
      mpProperties->SetValue(TANGENTIAL_PENALTY_RATIO, CustomProperties["TANGENTIAL_PENALTY_RATIO"].GetDouble());
      mpProperties->SetValue(TAU_STAB, CustomProperties["TAU_STAB"].GetDouble());
      mpProperties->SetValue(THICKNESS, 1.0);
      mpProperties->SetValue(CONTACT_FRICTION_ANGLE, 0.0);

      mrMainModelPart.AddProperties(mpProperties);

      // create geometry prototype for the contact conditions
      GeometryType::Pointer pGeometry;
      if( Dimension == 2 )
	pGeometry = Kratos::make_shared<Point2DType>(*((mrMainModelPart.Nodes().begin()).base()));
      else if( Dimension == 3 )
	pGeometry = Kratos::make_shared<Point3DType>(*((mrMainModelPart.Nodes().begin()).base()));


      // create condition prototype
      std::string ConditionName = CustomParameters["contact_condition_type"].GetString();

      unsigned int LastConditionId = 1;
      if( mrMainModelPart.NumberOfConditions() != 0 )
	LastConditionId = mrMainModelPart.Conditions().back().Id() + 1;


      if(  ConditionName == "PointContactPenaltyCondition2D1N" ){
      	return Kratos::make_shared<PointRigidContactPenalty2DCondition>(LastConditionId, pGeometry, mpProperties, mpParametricWall);
      }
      else if(  ConditionName == "PointContactPenaltyCondition3D1N" ){
      	return Kratos::make_shared<PointRigidContactPenalty3DCondition>(LastConditionId, pGeometry, mpProperties, mpParametricWall);
      }
      else if(  ConditionName == "AxisymPointContactPenaltyCondition2D1N" ){
      	return Kratos::make_shared<AxisymPointRigidContactPenalty2DCondition>(LastConditionId, pGeometry, mpProperties, mpParametricWall);
      }
       else if(  ConditionName == "EPPointContactPenaltyCondition3D1N" ) {
         return Kratos::make_shared<EPPointRigidContactPenalty3DCondition>(LastConditionId, pGeometry, mpProperties, mpParametricWall);
     }
     else if(  ConditionName == "EPPointContactPenaltyCondition2D1N" ) {
       return Kratos::make_shared<EPPointRigidContactPenalty2DCondition>(LastConditionId, pGeometry, mpProperties, mpParametricWall);
     }
     else if(  ConditionName == "EPPointContactPenaltywPCondition3D1N" ) {
       return Kratos::make_shared<EPPointRigidContactPenaltywP3DCondition>(LastConditionId, pGeometry, mpProperties, mpParametricWall);
     }
     else if(  ConditionName == "EPAxisymPointContactPenaltyCondition2D1N" ) {
       return Kratos::make_shared<EPAxisymPointRigidContactPenalty2DCondition>(LastConditionId, pGeometry, mpProperties, mpParametricWall);
      } else {
        std::cout << ConditionName << std::endl;
        KRATOS_ERROR << "the specified contact condition does not exist " << std::endl;
      }
      // else if(  ConditionName == "AxisymPointWaterContactPenaltyCondition2D1N" ){
      // 	return Kratos::make_shared<AxisymPointRigidContactPenaltyWater2DCondition>(LastConditionId, pGeometry, mpProperties, mpParametricWall);
      // }
      // else if(  ConditionName == "BeamPointRigidContactPenalty3DCondition" ){
      // 	return Kratos::make_shared<BeamPointRigidContactPenalty3DCondition>(LastConditionId, pGeometry, mpProperties, mpParametricWall);
      // }


      //------------------//

      // When conditions are registered....
      // Condition::NodesArrayType ConditionNodes;
      // ConditionNodes.push_back( *((mrMainModelPart.Nodes().begin()).base()) );
      // ConditionType const& rCloneCondition = KratosComponents<ConditionType>::Get(ConditionName);
      // mpConditionType = rCloneCondition.Create(LastConditionId, ElementNodes, mpProperties);

      return NULL;

      KRATOS_CATCH( "" )
    }

    //**************************************************************************
    //**************************************************************************

    void SearchContactConditions()
    {
      KRATOS_TRY

      //Create Rigid Contact Conditions

      //update parametric wall position
      ProcessInfo& rCurrentProcessInfo= mrMainModelPart.GetProcessInfo();
      double Time = rCurrentProcessInfo[TIME];


#ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
#else
      int number_of_threads = 1;
#endif

      ModelPart::NodesContainerType& rNodes = mrMainModelPart.Nodes();

      vector<unsigned int> nodes_partition;
      OpenMPUtils::CreatePartition(number_of_threads, rNodes.size(), nodes_partition);


#pragma omp parallel
{

	int k = OpenMPUtils::ThisThread();
	ModelPart::NodesContainerType::iterator NodesBegin = rNodes.begin() + nodes_partition[k];
	ModelPart::NodesContainerType::iterator NodesEnd = rNodes.begin() + nodes_partition[k + 1];

	for(ModelPart::NodesContainerType::const_iterator nd = NodesBegin; nd != NodesEnd; nd++)
	  {
	    if( nd->Is(BOUNDARY) ){

	      if( nd->IsNot(RIGID) )
		{
		  //to perform contact with a tube radius must be set
		  double Radius = 0;

		  if( nd->IsNot(SLAVE) ){
		    //Radius = nd->GetValue(MEAN_RADIUS);
		  }
		  else{
		    Radius = 0;
		  }

		  Vector Point(3);
		  Point[0] = nd->X();
		  Point[1] = nd->Y();
		  Point[2] = nd->Z();

		  if( mpParametricWall->IsInside(Point,Time,Radius) ){
		    nd->Set(CONTACT);
		  }
		}

	    }

	  }



}

      // **************** Serial version of the same search:  **************** //

      // ModelPart::NodesContainerType& rNodes = mrMainModelPart.Nodes();

      // for(ModelPart::NodesContainerType::const_iterator nd = rNodes.begin(); nd != rNodes.end(); nd++)
      // 	{
      // 	  if( nd->Is(BOUNDARY) ){

      // 	    if( nd->IsNot(RIGID) )
      // 	      {
      // 		//to perform contact with a tube radius must be set
      // 		double Radius = 0;

      // 		if( nd->IsNot(SLAVE) ){
      // 		  //Radius = nd->GetValue(MEAN_RADIUS);
      // 		}
      // 		else{
      // 		  Radius = 0;
      // 		}

      // 		Vector Point(3);
      // 		Point[0] = nd->X();
      // 		Point[1] = nd->Y();
      // 		Point[2] = nd->Z();

      // 		if( mpParametricWall->IsInside(Point,Time,Radius) ){
      // 		  nd->Set(CONTACT);
      // 		}
      // 	      }

      // 	  }

      // 	}

      // **************** Serial version of the same search:  **************** //



      KRATOS_CATCH( "" )

    }


    //**************************************************************************
    //**************************************************************************

    Condition::Pointer FindPointCondition(ModelPart& rModelPart, Node<3>::Pointer pPoint)
    {

     KRATOS_TRY
      const ProcessInfo& rCurrentProcessInfo= mrMainModelPart.GetProcessInfo();
      if ( rCurrentProcessInfo.Has(IS_RESTARTED) && rCurrentProcessInfo.Has(LOAD_RESTART) ) {
         if ( rCurrentProcessInfo[IS_RESTARTED] == true) {
            if ( rCurrentProcessInfo[STEP] == rCurrentProcessInfo[LOAD_RESTART] ) {
std::cout << " doing my.... ";
               return mpConditionType;

            }
         }
      }

     for(ModelPart::ConditionsContainerType::iterator i_cond =rModelPart.ConditionsBegin(); i_cond!= rModelPart.ConditionsEnd(); i_cond++)
       {
	 if( i_cond->Is(CONTACT) && i_cond->GetGeometry().size() == 1 ){
	   if( i_cond->GetGeometry()[0].Id() == pPoint->Id() ){
	     return ( *(i_cond.base()) );
	   }
	 }
       }

     return  mpConditionType;

     KRATOS_CATCH( "" )

    }

    //**************************************************************************
    //**************************************************************************

    void ClearContactFlags ( )
    {
     KRATOS_TRY

      for(ModelPart::NodesContainerType::iterator i_node = mrMainModelPart.NodesBegin(); i_node!= mrMainModelPart.NodesEnd(); i_node++)
	{
	  if( i_node->Is(CONTACT) ){
	    i_node->Set(CONTACT,false);
	  }
	}

      KRATOS_CATCH( "" )
    }

    //**************************************************************************
    //**************************************************************************

    void RestoreContactFlags ( )
    {
      KRATOS_TRY

      for(ModelPart::ConditionsContainerType::iterator i_cond = mrMainModelPart.ConditionsBegin(); i_cond!= mrMainModelPart.ConditionsEnd(); i_cond++)
	{
	  if( i_cond->Is(CONTACT) ){
	    for(unsigned int i=0; i<i_cond->GetGeometry().size(); i++)
	      {
		i_cond->GetGeometry()[i].Set(CONTACT,true);
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


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ParametricWallContactSearchProcess& operator=(ParametricWallContactSearchProcess const& rOther);

    /// Copy constructor.
    //Process(Process const& rOther);


    ///@}

  }; // Class Process

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    ParametricWallContactSearchProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const ParametricWallContactSearchProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_PARAMETRIC_WALL_CONTACT_SEARCH_PROCESS_H_INCLUDED  defined
