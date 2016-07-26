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

// Contact Points
// #include "custom_conditions/point_rigid_contact_penalty_3D_condition.hpp"
// #include "custom_conditions/axisym_point_rigid_contact_penalty_2D_condition.hpp"
// #include "custom_conditions/axisym_point_rigid_contact_penalty_water_2D_condition.hpp"
// #include "custom_conditions/beam_point_rigid_contact_penalty_3D_condition.hpp"
// #include "custom_conditions/beam_point_rigid_contact_LM_3D_condition.hpp"
// #include "custom_conditions/rigid_body_point_rigid_contact_condition.hpp"

//#include "custom_conditions/custom_friction_laws/friction_law.hpp"

#include "contact_mechanics_application_variables.h"


namespace Kratos
{

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
    //typedef FrictionLaw::pointer           FrictionLawType;
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ParametricWallContactSearchProcess(ModelPart& rModelPart): mrModelPart(rModelPart) {}


    ParametricWallContactSearchProcess( ModelPart& rModelPart, 
					SpatialBoundingBox::Pointer pParametricWall, 
					Parameters CustomParameters) 
      : mrModelPart(rModelPart)
    {
      KRATOS_TRY

      mEchoLevel = 0;

      mpParametricWall = pParametricWall;
	   
      Parameters DefaultParameters( R"(
            {
                   "contact_condition_type": "PointContactCondition2D1N",
                   "friction_law_type": "FrictionLaw",
                   "implemented_in_module": "KratosMultiphysics.ContactMechanicsApplication",
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

      //create properties for the contact conditions
      unsigned int NumberOfProperties = rModelPart.NumberOfProperties();
      
      mpProperties = PropertiesType::Pointer(new PropertiesType(NumberOfProperties));

      // std::string FrictionLawName = CustomParameters["friction_law_type"].GetString();
      // FrictionLawType const& rCloneFrictionLaw = KratosComponents<FrictionLawType>::Get(FrictionLawName);
      // mpProperties->SetValue(FRICTION_LAW, rCloneFrictionLaw.Clone() );
      
      Parameters CustomProperties = CustomParameters["variables_of_properties"];

      mpProperties->SetValue(FRICTION_ACTIVE, CustomProperties["FRICTION_ACTIVE"].GetBool());
      mpProperties->SetValue(MU_STATIC, CustomProperties["MU_STATIC"].GetDouble());
      mpProperties->SetValue(MU_DYNAMIC, CustomProperties["MU_DYNAMIC"].GetDouble());
      mpProperties->SetValue(PENALTY_PARAMETER, CustomProperties["PENALTY_PARAMETER"].GetDouble());
      mpProperties->SetValue(TANGENTIAL_PENALTY_RATIO, CustomProperties["TANGENTIAL_PENALTY_RATIO"].GetDouble());
      mpProperties->SetValue(TAU_STAB, CustomProperties["TAU_STAB"].GetDouble());

      mrModelPart.AddProperties(mpProperties);

      // create condition prototype
      std::string ConditionName = CustomParameters[" PointContactCondition2D1N "].GetString();  

      //unsigned int LastConditionId  = rModelPart.Conditions().back().Id() + 1;

      ProcessInfo& rCurrentProcessInfo= mrModelPart.GetProcessInfo(); 
      double Dimension = rCurrentProcessInfo[DOMAIN_SIZE];

      GeometryType::Pointer pGeometry;
      if( Dimension == 2 )
	pGeometry = GeometryType::Pointer(new Point2DType( *((rModelPart.Nodes().begin()).base()) ));
      else if( Dimension == 3 )
	pGeometry = GeometryType::Pointer(new Point3DType( *((rModelPart.Nodes().begin()).base()) ));

      // if(  ElementName.compare("PointRigidContactPenalty2DCondition") == 0 ){
      // 	mpConditionType = ConditionType::Pointer(new PointRigidContactPenalty2DCondition(LastConditionId, pGeometry, mpProperties, mpParametricWall));
      // }
      // else if(  ElementName.compare("PointRigidContactPenalty3DCondition") == 0 ){
      // 	mpConditionType = ConditionType::Pointer(new PointRigidContactPenalty3DCondition(LastConditionId, pGeometry, mpProperties, mpParametricWall));
      // }
      // else if(  ElementName.compare("AxisymPointRigidContactPenalty2DCondition") == 0 ){
      // 	mpConditionType = ConditionType::Pointer(new AxisymPointRigidContactPenalty2DCondition(LastConditionId, pGeometry, mpProperties, mpParametricWall));
      // }
      // else if(  ElementName.compare("AxisymPointRigidContactPenaltyWater2DCondition") == 0 ){
      // 	mpConditionType = ConditionType::Pointer(new AxisymPointRigidContactPenaltyWater2DCondition(LastConditionId, pGeometry, mpProperties, mpParametricWall));
      // }
      // else if(  ElementName.compare("BeamPointRigidContactPenalty3DCondition") == 0 ){
      // 	mpConditionType = ConditionType::Pointer(new BeamPointRigidContactPenalty3DCondition(LastConditionId, pGeometry, mpProperties, mpParametricWall));
      // }



      // When conditions are registered....
      // Condition::NodesArrayType ConditionNodes;
      // ConditionNodes.push_back( *((rModelPart.Nodes().begin()).base()) );
      // ConditionType const& rCloneCondition = KratosComponents<ConditionType>::Get(ConditionName);
      // mpConditionType = rCloneCondition.Create(LastConditionId, ElementNodes, mpProperties);


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
    virtual void Execute()
    {

      KRATOS_TRY

      ProcessInfo& rCurrentProcessInfo= mrModelPart.GetProcessInfo(); 
      double Time      = rCurrentProcessInfo[TIME];
      double Dimension = rCurrentProcessInfo[DOMAIN_SIZE];

      if (Time == 0)
	KRATOS_THROW_ERROR( std::logic_error, "detected time = 0 in the Solution Scheme ... check if the time step is created correctly for the current model part", "" )


	  // update center position
	  mpParametricWall->UpdateBoxPosition( Time );
	    
	          

      // reset CONTACT flag to all mesh nodes
      ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();
	     
      for(ModelPart::NodesContainerType::iterator nd = rNodes.begin(); nd != rNodes.end(); ++nd)
	{
	  nd->Set(CONTACT,false);
	}
	     

      // set BOUNDARY flag to nodes from structural elements
      for(ModelPart::ElementsContainerType::iterator ie = mrModelPart.ElementsBegin(); ie!=mrModelPart.ElementsEnd(); ie++)
	{
	       
	  if( ie->GetGeometry().size() == 2 ){
	    for( unsigned int i=0; i<2; i++ )
	      {
		ie->GetGeometry()[i].Set(BOUNDARY,true);
	      }
	  }
	       
	}

      // set BOUNDARY flag to nodes from RIGID meshes
      ModelPart::MeshesContainerType& rMeshes = mrModelPart.GetMeshes();


      for(unsigned int m=0; m<rMeshes.size(); m++)
	{
	       
	  if( rMeshes[m].Is(RIGID) ){
		 		 
	    for(ModelPart::ElementsContainerType::iterator ie = rMeshes[m].ElementsBegin(); ie!=rMeshes[m].ElementsEnd(); ie++){
		   
	      for( unsigned int i=0; i<ie->GetGeometry().size(); i++ )
		{
		  ie->GetGeometry()[i].Set(BOUNDARY,true);
		}
	    }
	  }
	       
	}
	     


      //Create Rigid Contact Conditions
	     
      int id = mrModelPart.Conditions().back().Id() + 1;
     
      mConditionsNumber =  mrModelPart.Conditions().size();


#ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
#else
      int number_of_threads = 1;
#endif


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

      // Vector WallVelocity =  mpParametricWall->Velocity();
      // int  MovementLabel  =  mpParametricWall->GetMovementLabel();

      // //Check RIGID walls and search contacts
      // for(ModelPart::NodesContainerType::iterator nd = rNodes.begin(); nd != rNodes.end(); ++nd)
      // 	{
      // 	  //set point rigid wall condition : usually in non rigid_wall points
      // 	  if( nd->SolutionStepsDataHas(RIGID_WALL) ){

      // 	    if( nd->GetSolutionStepValue(RIGID_WALL) == MovementLabel ){

      // 	      //nd->Set(STRUCTURE);
      // 	      nd->Set(RIGID);

      // 	      //set new coordinates (MOVE MESH)
      // 	      nd->X() = nd->X0() + WallVelocity[0] * Time;
      // 	      nd->Y() = nd->Y0() + WallVelocity[1] * Time;
      // 	      nd->Z() = nd->Z0() + WallVelocity[2] * Time;

      // 	      //std::cout<<" node "<<(nd)->Id()<<" Position ("<<(nd)->X()<<", "<<(nd)->Y()<<" "<<(nd)->Z()<<") "<<std::endl;
      // 	    }
      // 	  }

      // 	  if( nd->Is(BOUNDARY) ){

      // 	    if( nd->IsNot(RIGID) )
      // 	      {
      // 		//to perform contact with a tube radius must be set
      // 		double Radius = 0;
      // 		if( !RigidBodyPresent ){
      // 		  Radius = nd->GetValue(MEAN_RADIUS);
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



      // create contact condition for rigid and deformable bodies

      for(ModelPart::NodesContainerType::ptr_iterator nd = rNodes.ptr_begin(); nd != rNodes.ptr_end(); ++nd)
	{

	  if( (*nd)->Is(BOUNDARY) && (*nd)->Is(CONTACT) ){

	    ConditionType::Pointer pCondition;
		   
		   
	    if( (*nd)->Is(RIGID) ){  //rigid wall contacting with a rigid body 
		     
	      GeometryType::Pointer pGeometry;
	      if( Dimension == 2 )
		pGeometry = GeometryType::Pointer(new Point2DType( (*nd) ));
	      else if( Dimension == 3 )
		pGeometry = GeometryType::Pointer(new Point3DType( (*nd) ));
	      
	      //pCondition= ModelPart::ConditionType::Pointer(new RigidBodyPointRigidContactCondition(id, pGeometry, mpProperties, mpParametricWall) );

	      mrModelPart.AddCondition(pCondition);
		     
	      id +=1;
		     
	    }
	    else{ //rigid wall contacting with a deformable body 

	      Condition::NodesArrayType  pConditionNode;
	      pConditionNode.push_back( (*nd) );

	      pCondition = mpConditionType->Clone(id, pConditionNode);

	      mrModelPart.AddCondition(pCondition);
		     
	      id +=1;
		     
	    }
		   
		   		
	  }       
	  else{ //clean nodal contact forces
	    //nd->FastGetSolutionStepValue(CONTACT_FORCE).clear();
	  }  
		     

	}

      //TODO: a comparison with the previous contact conditions to clone them and preserve the contact data. That implies that they have to be kept somewhere and not deleted in the FinalizeSolutionStep as is done now.


      if( mEchoLevel > 1 )
	std::cout<<"  [ Rigid Contacts : "<<mrModelPart.Conditions().size() - mConditionsNumber<<" ]"<<std::endl;
	     
      KRATOS_CATCH( "" )

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
   
      KRATOS_TRY

      //Clean Rigid Contact Conditions

      ModelPart::ConditionsContainerType NonRigidContactConditions;
	   
      unsigned int id = 0;
	   
      //std::cout<<" [ NUMBER OF CONDITIONS before rigid contact update: "<<mrModelPart.Conditions().size()<<" ]"<<std::endl;
	   
      for(ModelPart::ConditionsContainerType::iterator ic = mrModelPart.ConditionsBegin(); ic!= mrModelPart.ConditionsEnd(); ic++)
	{
	  if( id == mConditionsNumber )
	    break;
	       
	  NonRigidContactConditions.push_back(*(ic.base()));  
	       
	  id +=1;
	}
	   
      mrModelPart.Conditions().swap( NonRigidContactConditions );
	   
      //std::cout<<"  [ NUMBER OF CONDITIONS after  rigid contact update: "<<mrModelPart.Conditions().size()<<" ]"<<std::endl;
	   

      KRATOS_CATCH( "" )      
	
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
      return "ParametricWallContactSearchProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
      rOStream << "ParametricWallContactSearchProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}


  private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Static Member Variables
    ///@{
    ModelPart&  mrModelPart;

    SpatialBoundingBox::Pointer  mpParametricWall;

    ConditionType::Pointer  mpConditionType;

    PropertiesType::Pointer mpProperties;
     
    unsigned int  mConditionsNumber;

    int  mEchoLevel;

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


