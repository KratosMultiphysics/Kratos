//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                February 2016 $
//   Revision:            $Revision:                      0.0 $
//
//


#if !defined(KRATOS_REFINE_CONDITIONS_IN_CONTACT_MESHER_PROCESS_H_INCLUDED )
#define  KRATOS_REFINE_CONDITIONS_IN_CONTACT_MESHER_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "includes/model_part.h"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_processes/refine_conditions_mesher_process.hpp"

#include "pfem_solid_mechanics_application_variables.h"

///VARIABLES used:
//StepData: NODAL_H, NORMAL, CONTACT_FORCE, DISPLACEMENT
//Flags:    (checked) BOUNDARY, TO_SPLIT
//          (set)     BOUNDARY(nodes), TO_ERASE(conditions), NEW_ENTITY(conditions,nodes)(set), TO_SPLIT(conditions)->locally
//          (modified)
//          (reset)   TO_SPLIT
//(set):=(set in this process)

namespace Kratos
{

///@name Kratos Classes
///@{

/// Refine Mesh Boundary Process
/** The process labels the boundary conditions (TO_SPLIT)
    Dependencies: RemoveMeshNodesProcess.Execute()  is needed as a previous step

    Determines if new conditions must be inserted in boundary.
    If boundary must to be kept (CONSTRAINED),
    New conditions will be rebuild (spliting the old ones and inserting new nodes)
    Old conditions will be erased at the end.

*/

class KRATOS_API(PFEM_SOLID_MECHANICS_APPLICATION) RefineConditionsInContactMesherProcess
  : public RefineConditionsMesherProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( RefineConditionsInContactMesherProcess );

    typedef ModelPart::NodeType                   NodeType;
    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;
    typedef PointerVector<NodeType>        PointsArrayType;

    typedef PointerVectorSet<ConditionType, IndexedObject> ConditionsContainerType;
    typedef ConditionsContainerType::iterator                    ConditionIterator;
    typedef ConditionsContainerType::const_iterator      ConditionConstantIterator;


    struct RefineCounters
    {
    public:

     //counters:
     int number_of_contact_conditions;
     int number_of_contacts;
     int number_of_active_contacts;

     int number_contact_size_insertions;
     int number_contact_tip_insertions;

     int number_of_exterior_bounds;
     int number_of_tip_bounds;
     int number_of_energy_bounds;

      void Initialize()
      {
	//counters:
	number_of_contact_conditions = 0;

	number_of_contacts   = 0;
	number_of_active_contacts   = 0;

	number_contact_size_insertions = 0;
	number_contact_tip_insertions  = 0;

	number_of_exterior_bounds = 0;
	number_of_tip_bounds      = 0;
	number_of_energy_bounds   = 0;

      }
    };


    ///@}
    ///@name Life Cycle
    ///@{


    /// Default constructor
    RefineConditionsInContactMesherProcess(ModelPart& rModelPart,
				     std::vector<SpatialBoundingBox::Pointer> mRigidWalls,
				     MesherUtilities::MeshingParameters& rRemeshingParameters,
				     int EchoLevel)
      :RefineConditionsMesherProcess(rModelPart,rRemeshingParameters,EchoLevel)
    {


    }

    /// Destructor.
    virtual ~RefineConditionsInContactMesherProcess() {}


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
    void Execute()  override
    {

      KRATOS_TRY

      if( this->mEchoLevel > 0 ){
        std::cout<<" [ REFINE BOUNDARY : "<<std::endl;
	//std::cout<<"   Nodes and Conditions : "<<mrModelPart.Nodes().size()<<", "<<mrModelPart.Conditions().size()<<std::endl;
      }

      mrRemesh.Info->InsertedBoundaryConditions    = mrModelPart.NumberOfConditions();
      mrRemesh.Info->InsertedBoundaryNodes = mrModelPart.NumberOfNodes();

      RefineCounters LocalRefineInfo;
      LocalRefineInfo.Initialize();


      //if the insert switches are activated, we check if the boundaries got too coarse
      if( (mrRemesh.Refine->RefiningOptions.Is(MesherUtilities::REFINE_INSERT_NODES) || mrRemesh.Refine->RefiningOptions.Is(MesherUtilities::REFINE_ADD_NODES)) && mrRemesh.Refine->RefiningOptions.Is(MesherUtilities::REFINE_BOUNDARY) )
	{

	  std::vector<Point > list_of_points;
	  std::vector<ConditionType::Pointer> list_of_conditions;

	  unsigned int conditions_size = mrModelPart.Conditions().size();

	  if( mrModelPart.Name() == mrRemesh.SubModelPartName ){

	    list_of_points.reserve(conditions_size);
	    list_of_conditions.reserve(conditions_size);

	    RefineContactBoundary(mrModelPart, list_of_points, list_of_conditions, LocalRefineInfo);

	    RefineOtherBoundary(mrModelPart, list_of_points, list_of_conditions, LocalRefineInfo);

	    BuildNewConditions(mrModelPart, list_of_points, list_of_conditions, LocalRefineInfo);

	    CleanConditionsAndFlags(mrModelPart);

	  }
	  else{

	    ModelPart& rModelPart = mrModelPart.GetSubModelPart(mrRemesh.SubModelPartName);

	    conditions_size = rModelPart.Conditions().size();
	    list_of_points.reserve(conditions_size);
	    list_of_conditions.reserve(conditions_size);

	    RefineContactBoundary(rModelPart, list_of_points, list_of_conditions, LocalRefineInfo);

	    RefineOtherBoundary(rModelPart, list_of_points, list_of_conditions, LocalRefineInfo);

	    BuildNewConditions(rModelPart, list_of_points, list_of_conditions, LocalRefineInfo);

	    CleanConditionsAndFlags(rModelPart);
	  }



	} // REFINE END;



     mrRemesh.Info->InsertedBoundaryConditions    = mrModelPart.NumberOfConditions()-mrRemesh.Info->InsertedBoundaryConditions;
     mrRemesh.Info->InsertedBoundaryNodes = mrModelPart.NumberOfNodes()-mrRemesh.Info->InsertedBoundaryNodes;

     if( this->mEchoLevel > 0 ){
        std::cout<<"   [ CONDITIONS ( inserted : "<<mrRemesh.Info->InsertedBoundaryConditions<<" ) ]"<<std::endl;
        std::cout<<"   [ NODES      ( inserted : "<<mrRemesh.Info->InsertedBoundaryNodes<<" ) ]"<<std::endl;
        std::cout<<"   [ contact(TIP: "<<LocalRefineInfo.number_contact_tip_insertions<<", SIZE: "<<LocalRefineInfo.number_contact_size_insertions<<") -  bound(TIP: "<<LocalRefineInfo.number_of_tip_bounds<<", SIZE: "<<LocalRefineInfo.number_of_exterior_bounds<<")]"<<std::endl;


        std::cout<<"   REFINE BOUNDARY ]; "<<std::endl;
     }

     KRATOS_CATCH(" ")

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
        return "RefineConditionsInContactMesherProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "RefineConditionsInContactMesherProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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
    ///@name Member Variables
    ///@{

    std::vector<SpatialBoundingBox::Pointer> mRigidWalls;

    ///@}
    ///@name Un accessible methods
    ///@{

    //*******************************************************************************************
    //*******************************************************************************************

    void CleanConditionsAndFlags(ModelPart& rModelPart)
    {

      KRATOS_TRY

      //swap conditions for a temporary use
      ModelPart::ConditionsContainerType TemporaryConditions;
      TemporaryConditions.reserve(rModelPart.Conditions().size());
      TemporaryConditions.swap(rModelPart.Conditions());

      for(ModelPart::ConditionsContainerType::iterator ic = TemporaryConditions.begin(); ic!= TemporaryConditions.end(); ic++)
	{
	  Geometry< Node<3> > rGeometry =ic->GetGeometry();
	  for(unsigned int i=0; i<rGeometry.size(); ++i)
	    {
              rGeometry[i].Reset(TO_SPLIT);
	    }

           if(ic->IsNot(TO_ERASE)){
	     rModelPart.AddCondition(*(ic.base()));
           }
           // else{
           //   std::cout<<"   Condition RELEASED:"<<ic->Id()<<std::endl;
           // }
        }


      KRATOS_CATCH(" ")

    }


    //*******************************************************************************************
    //*******************************************************************************************

    void BuildNewConditions( ModelPart& rModelPart, std::vector<Point >& list_of_points, std::vector<ConditionType::Pointer>& list_of_conditions, RefineCounters& rLocalRefineInfo )
    {

      KRATOS_TRY

      //node to get the DOFs from
      Node<3>::DofsContainerType& reference_dofs = (rModelPart.NodesBegin())->GetDofs();
      //unsigned int step_data_size = rModelPart.GetNodalSolutionStepDataSize();
      double z = 0.0;

      unsigned int initial_node_size = mrModelPart.Nodes().size()+1; //total model part node size
      unsigned int initial_cond_size = mrModelPart.Conditions().size()+1; //total model part node size

      Node<3> new_point(0,0.0,0.0,0.0);

      int id=0;
      //if points were added, new nodes must be added to ModelPart
      for(unsigned int i = 0; i<list_of_points.size(); ++i)
        {
	  id   = initial_node_size+i;

	  double& x= list_of_points[i].X();
	  double& y= list_of_points[i].Y();

	  Node<3>::Pointer pnode = rModelPart.CreateNewNode(id,x,y,z);

	  pnode->SetBufferSize(rModelPart.NodesBegin()->GetBufferSize() );


	  //assign data to dofs

	  //2D edges:
	  Geometry< Node<3> >& rConditionGeometry  = list_of_conditions[i]->GetGeometry();


	  //generating the dofs
	  for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); ++iii)
	    {
              Node<3>::DofType& rDof = *iii;
              Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );

              if( rConditionGeometry[0].IsFixed(rDof.GetVariable()) && rConditionGeometry[1].IsFixed(rDof.GetVariable()) )
		(p_new_dof)->FixDof();
              else
		(p_new_dof)->FreeDof();
	    }

	  //assign data to dofs
	  VariablesList& variables_list = rModelPart.GetNodalSolutionStepVariablesList();

	  PointsArrayType  PointsArray;
	  PointsArray.push_back( rConditionGeometry(0) );
	  PointsArray.push_back( rConditionGeometry(1) );

	  Geometry<Node<3> > geom( PointsArray );

	  std::vector<double> N(2);
	  N[0] = 0.5;
	  N[1] = 0.5;

	  MeshDataTransferUtilities DataTransferUtilities;
	  DataTransferUtilities.Interpolate2Nodes( geom, N, variables_list, *pnode);

	  // //int cond_id = list_of_points[i].Id();
	  // //Geometry< Node<3> >& rConditionGeometry = (*(rModelPart.Conditions().find(cond_id).base()))->GetGeometry();

	  // unsigned int buffer_size = pnode->GetBufferSize();

	  // for(unsigned int step = 0; step<buffer_size; step++)
	  //   {
          //     //getting the data of the solution step
          //     double* step_data = (pnode)->SolutionStepData().Data(step);

          //     double* node0_data = rConditionGeometry[0].SolutionStepData().Data(step);
          //     double* node1_data = rConditionGeometry[1].SolutionStepData().Data(step);

          //     //copying this data in the position of the vector we are interested in
          //     for(unsigned int j= 0; j<step_data_size; j++)
	  // 	{
	  // 	  step_data[j] = 0.5*node0_data[j] + 0.5*node1_data[j];
	  // 	}
	  //   }


	  //set specific control values and flags:
	  pnode->Set(BOUNDARY);
	  pnode->Set(NEW_ENTITY);  //if boundary is rebuild, the flag INSERTED must be set to new conditions too
	  //std::cout<<"   Node ["<<pnode->Id()<<"] is a NEW_ENTITY "<<std::endl;

	  double& nodal_h = pnode->FastGetSolutionStepValue(NODAL_H);
	  //nodal_h = 0.5*(nodal_h+mrRemesh.Refine->CriticalSide); //modify nodal_h for security
	  nodal_h = mrRemesh.Refine->CriticalSide; //modify nodal_h for security

	  const array_1d<double,3> ZeroNormal(3,0.0);
	  //correct normal interpolation
	  noalias(pnode->GetSolutionStepValue(NORMAL)) = list_of_conditions[i]->GetValue(NORMAL);


	  //correct contact_normal interpolation (laplacian boundary projection uses it)
	  array_1d<double, 3 > & ContactForceNormal1  = rConditionGeometry[0].FastGetSolutionStepValue(CONTACT_FORCE);
	  array_1d<double, 3 > & ContactForceNormal2  = rConditionGeometry[1].FastGetSolutionStepValue(CONTACT_FORCE);
	  if(norm_2(ContactForceNormal1)==0 || norm_2(ContactForceNormal2)==0)
	    noalias(pnode->FastGetSolutionStepValue(CONTACT_FORCE)) = ZeroNormal;

	  //recover the original position of the node
	  const array_1d<double,3>& disp = pnode->FastGetSolutionStepValue(DISPLACEMENT);
	  pnode->X0() = pnode->X() - disp[0];
	  pnode->Y0() = pnode->Y() - disp[1];
	  pnode->Z0() = pnode->Z() - disp[2];

	  //Conditions must be also created with the add of a new node:
	  Condition::NodesArrayType face1;
	  Condition::NodesArrayType face2;
	  face1.reserve(2);
	  face2.reserve(2);

	  face1.push_back(rConditionGeometry(0));
	  face1.push_back(pnode);

	  face2.push_back(pnode);
	  face2.push_back(rConditionGeometry(1));

	  id   = initial_cond_size+(i*2);

	  ConditionType::Pointer pcond1      = list_of_conditions[i]->Clone(id, face1);
	  //std::cout<<" ID"<<id<<" 1s "<<pcond1->GetGeometry()[0].Id()<<" "<<pcond1->GetGeometry()[1].Id()<<std::endl;
	  id   = initial_cond_size+(i*2+1);
	  ConditionType::Pointer pcond2      = list_of_conditions[i]->Clone(id, face2);
	  //std::cout<<" ID"<<id<<" 2s "<<pcond2->GetGeometry()[0].Id()<<" "<<pcond2->GetGeometry()[1].Id()<<std::endl;

	  pcond1->Set(NEW_ENTITY);
	  pcond2->Set(NEW_ENTITY);

	  pcond1->SetValue(MASTER_NODES, list_of_conditions[i]->GetValue(MASTER_NODES) );
	  pcond1->SetValue(NORMAL, list_of_conditions[i]->GetValue(NORMAL) );
	  pcond1->SetValue(CAUCHY_STRESS_VECTOR,list_of_conditions[i]->GetValue(CAUCHY_STRESS_VECTOR));
	  pcond1->SetValue(DEFORMATION_GRADIENT,list_of_conditions[i]->GetValue(DEFORMATION_GRADIENT));

	  pcond2->SetValue(MASTER_NODES, list_of_conditions[i]->GetValue(MASTER_NODES) );
	  pcond2->SetValue(NORMAL, list_of_conditions[i]->GetValue(NORMAL) );
	  pcond2->SetValue(CAUCHY_STRESS_VECTOR,list_of_conditions[i]->GetValue(CAUCHY_STRESS_VECTOR));
	  pcond2->SetValue(DEFORMATION_GRADIENT,list_of_conditions[i]->GetValue(DEFORMATION_GRADIENT));


	  //std::cout<<" pcond0 "<<rModelPart.NumberOfConditions()<<" "<<mrRemesh.Info->InsertedBoundaryConditions<<std::endl;

	  rModelPart.AddCondition(pcond1);

	  //std::cout<<" pcond1 "<<rModelPart.NumberOfConditions()<<" "<<mrRemesh.Info->InsertedBoundaryConditions<<std::endl;

	  rModelPart.AddCondition(pcond2);

	  //std::cout<<" pcond2 "<<rModelPart.NumberOfConditions()<<" "<<mrRemesh.Info->InsertedBoundaryConditions<<std::endl;
        }


      KRATOS_CATCH(" ")

    }

    //*******************************************************************************************
    //*******************************************************************************************

    bool RefineContactBoundary(ModelPart& rModelPart, std::vector<Point >& list_of_points, std::vector<ConditionType::Pointer>& list_of_conditions, RefineCounters& rLocalRefineInfo )
    {

      KRATOS_TRY

	//ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

      //***SIZES :::: parameters do define the tolerance in mesh size:

      //DEFORMABLE CONTACT:
      double factor_for_tip_radius     = 0.2; //deformable contact tolerance in radius for detection tip sides to refine
      double factor_for_non_tip_side   = 3.0; // will be multiplied for nodal_h of the master node to compare with boundary nodes average nodal_h in a contact conditio which master node do not belongs to a tip

      double size_for_tip_contact_side      = 0.4 * mrRemesh.Refine->CriticalSide; // length size for the contact tip side
      double size_for_non_tip_contact_side  = 2.0 * mrRemesh.Refine->CriticalSide; // compared with contact size which master node do not belongs to a tip


      //*********************************************************************************
      // DETECTION OF NODES ON TIP CONTACTS START
      //*********************************************************************************

      unsigned int nodes_on_wall_tip = 0;
      double ContactFace = 0;
      for(ModelPart::ConditionsContainerType::iterator i_cond =rModelPart.ConditionsBegin(); i_cond!= rModelPart.ConditionsEnd(); i_cond++)
	{
	  if( i_cond->Is(CONTACT) && i_cond->GetGeometry().size() == 1 ){

	    // i_cond->GetValueOnIntegrationPoints(ContactFace, CONTACT_FACE, CurrentProcessInfo);

	    if( int(ContactFace) == 2 ){ // tip

	      i_cond->GetGeometry()[0].Set(TO_SPLIT);
	      nodes_on_wall_tip ++;

	    }
	  }
	}


      if( this->mEchoLevel > 0 )
        std::cout <<"   [ NODES ON WALL TIP: ( " <<nodes_on_wall_tip <<" ) ]"<<std::endl;

      //*********************************************************************************
      // DETECTION OF NODES ON TIP CONTACTS END
      //*********************************************************************************


      // std::vector<int> nodes_ids;
      // nodes_ids.resize(rModelPart.Conditions().size()); //mesh 0
      // std::fill( nodes_ids.begin(), nodes_ids.end(), 0 );

      //std::cout<<"   List of Conditions Reserved Size: "<<conditions_size<<std::endl;

      double tool_radius = 0;
      double side_length = 0;

      bool size_insert     = false;
      bool radius_insert   = false;
      bool contact_active  = false;

      bool contact_semi_active = false;
      bool tool_project        = false;

      std::vector<bool> semi_active_nodes;
      Node<3> new_point(0,0.0,0.0,0.0);


      for(ModelPart::ConditionsContainerType::iterator ic = mrModelPart.ConditionsBegin(); ic!= mrModelPart.ConditionsEnd(); ic++)
        {

	  size_insert    = false;
	  radius_insert  = false;
	  tool_project   = false;
	  contact_active = false;

	  tool_radius    = 0;
	  side_length    = 0;

	  Geometry< Node<3> > rConditionGeometry;
	  array_1d<double,3> tip_center;
	  tip_center.clear();


	  //LOOP TO CONSIDER ONLY CONTACT CONDITIONS
	  if( ic->Is(CONTACT) )   //Refine radius on the workpiece for the ContactDomain zone
	    {

              Node<3>  MasterNode;
              bool condition_found = false;
              ConditionType::Pointer MasterCondition  = mMesherUtilities.FindMasterCondition(*(ic.base()),MasterNode,rModelPart.Conditions(),condition_found);


              if(condition_found){

		if( MasterCondition->IsNot(TO_ERASE) ){


		  rConditionGeometry  = MasterCondition->GetGeometry();

		  //to recover TIP definition on conditions
		  if( MasterNode.SolutionStepsDataHas( WALL_TIP_RADIUS ) ) //master node in tool -->  refine workpiece  //
                    {

		      tool_radius = MasterNode.FastGetSolutionStepValue( WALL_TIP_RADIUS );
		      tip_center  = MasterNode.FastGetSolutionStepValue( WALL_REFERENCE_POINT );
		      // WARNING THE UPDATE OF THE TIP CENTER IS NEEDED !!!!

		      array_1d<double, 3 > radius;
		      radius[0]=rConditionGeometry[0].X()-tip_center[0];
		      radius[1]=rConditionGeometry[0].Y()-tip_center[1];
		      radius[2]=rConditionGeometry[0].Z()-tip_center[2];
		      double distance1=norm_2(radius);

		      radius[0]=rConditionGeometry[1].X()-tip_center[0];
		      radius[1]=rConditionGeometry[1].Y()-tip_center[1];
		      radius[2]=rConditionGeometry[1].Z()-tip_center[2];

		      double distance2=norm_2(radius);


		      // TO SPLIT DETECTION START
		      //If a node is detected in the wall tip is set TO_SPLIT
		      //the criteria to splitting will be applied later in the nodes marked as TO_SPLIT

		      if( (1-factor_for_tip_radius)*tool_radius < distance1 &&  distance1 < (1+factor_for_tip_radius)*tool_radius )
			rConditionGeometry[0].Set(TO_SPLIT);

		      if( (1-factor_for_tip_radius)*tool_radius < distance2 &&  distance2 < (1+factor_for_tip_radius)*tool_radius )
			rConditionGeometry[1].Set(TO_SPLIT);

		      // TO SPLIT DETECTION END


		      // ACTIVE CONTACT DETECTION START

		      contact_active = mMesherUtilities.CheckContactActive(rConditionGeometry, contact_semi_active, semi_active_nodes);
		      if(contact_active){
			rLocalRefineInfo.number_of_active_contacts ++;
		      }

		      // ACTIVE CONTACT DETECTION END


		      side_length = mMesherUtilities.CalculateSideLength (rConditionGeometry[0],rConditionGeometry[1]);

		      if( side_length > size_for_tip_contact_side ){


			if( ((1-factor_for_tip_radius)*tool_radius < distance1 &&  (1-factor_for_tip_radius)*tool_radius < distance2) &&
			    (distance1 < (1+factor_for_tip_radius)*tool_radius  &&  distance2 < (1+factor_for_tip_radius)*tool_radius) )
                          {
			    radius_insert = true;

                          }
			// else if( side_length > size_for_tip_contact_side &&
			// 	       ( distance1 < (1 + (side_size_factor+factor_for_tip_radius))*tool_radius && distance2 < (1 + (side_size_factor+factor_for_tip_radius))*tool_radius) ) {

			// 	size_insert = true;
			// 	std::cout<<" insert on radius-size "<<std::endl;
			// }

		      }


		      if(radius_insert){

			if(!contact_active){

			  radius_insert = false;
			  // std::cout<<" contact_not_active "<<std::endl;
			  // double& nodal_h1 = rConditionGeometry[0].FastGetSolutionStepValue(NODAL_H);
			  // double& nodal_h2 = rConditionGeometry[1].FastGetSolutionStepValue(NODAL_H);
			  // double& nodal_h0 = MasterNode.FastGetSolutionStepValue( NODAL_H );

			  // double side = norm_2(rConditionGeometry[0]-rConditionGeometry[1]);
			  // // double d1 = mMesherUtilities.FindBoundaryH (rConditionGeometry[0]);
			  // // double d2 = mMesherUtilities.FindBoundaryH (rConditionGeometry[1]);
			  // // double d0 = mMesherUtilities.FindBoundaryH (MasterNode);
			  // // double size_master = nodal_h0;

			  // bool candidate =false;
			  // if( ((nodal_h1+nodal_h2)*0.5) > factor_for_non_tip_side * nodal_h0 ){
			  //   candidate = true;
			  // }

			  // double crit_factor = 2;
			  // if( (side > size_for_non_tip_contact_side) && candidate ){
			  //   radius_insert = true;
			  // }

			}

		      }


                    }
		  else{ //refine boundary with nodal_h sizes to large

		    double& nodal_h1 = rConditionGeometry[0].FastGetSolutionStepValue(NODAL_H);
		    double& nodal_h2 = rConditionGeometry[1].FastGetSolutionStepValue(NODAL_H);
		    double& nodal_h0 = MasterNode.FastGetSolutionStepValue( NODAL_H );

		    double side = norm_2(rConditionGeometry[0]-rConditionGeometry[1]);
		    // double d1 = mMesherUtilities.FindBoundaryH (rConditionGeometry[0]);
		    // double d2 = mMesherUtilities.FindBoundaryH (rConditionGeometry[1]);
		    // double d0 = mMesherUtilities.FindBoundaryH (MasterNode);
		    // double size_master = nodal_h0;


		    bool candidate =false;
		    if( ((nodal_h1+nodal_h2)*0.5) > factor_for_non_tip_side * nodal_h0 ){
		      candidate = true;
		    }

		    if( (side > size_for_non_tip_contact_side ) && candidate ){
		      size_insert = true;
		    }
		  }



		  if( radius_insert || size_insert ) //Boundary must be rebuild
                    {

		      //std::cout<<"   CONTACT DOMAIN ELEMENT REFINED "<<ic->Id()<<std::endl;

		      new_point.X() = 0.5*( rConditionGeometry[1].X() + rConditionGeometry[0].X() );
		      new_point.Y() = 0.5*( rConditionGeometry[1].Y() + rConditionGeometry[0].Y() );
		      new_point.Z() = 0.5*( rConditionGeometry[1].Z() + rConditionGeometry[0].Z() );


		      new_point.SetId(ic->Id()); //set condition Id

		      ConditionType::Pointer ContactMasterCondition  = ic->GetValue(MASTER_CONDITION);


		      if( (rConditionGeometry[0].Is(TO_SPLIT) && rConditionGeometry[1].Is(TO_SPLIT)) )
			tool_project = true;

		      if( (rConditionGeometry[0].Is(TO_SPLIT) || rConditionGeometry[1].Is(TO_SPLIT)) && contact_active)
			tool_project = true;


		      if(tool_project) //master node in tool -->  refine workpiece  // (tool_radius ==0 in workpiece nodes)
			{

                          if(new_point.Y()<(tip_center[1]) && new_point.Y()>(tip_center[1]-tool_radius)){

			    // std::cout<<"   new_point  ("<<new_point.X()<<", "<<new_point.Y()<<") "<<std::endl;
			    // std::cout<<"   tip_center ("<<tip_center[0]<<", "<<tip_center[1]<<") radius "<<tool_radius<<std::endl;

			    array_1d<double,3> tip_normal = tip_center-new_point;

			    if(norm_2(tip_normal)<tool_radius*0.95){ //if is in the tool tip
			      tip_normal -= (tool_radius/norm_2(tip_normal)) * tip_normal;
			      if(norm_2(tip_normal)<tool_radius*0.05)
				new_point  += tip_normal;

			      // std::cout<<"   A: Tool Tip Correction COND ("<<ContactMasterCondition->Id()<<") "<<std::endl;
			      // std::cout<<"   new_point ("<<new_point.X()<<", "<<new_point.Y()<<") "<<std::endl;
			    }
                          }
			}

		      if(radius_insert)
			rLocalRefineInfo.number_contact_tip_insertions++;
		      if(size_insert)
			rLocalRefineInfo.number_contact_size_insertions++;

		      // std::cout<<"   MasterCondition RELEASED (Id: "<<ContactMasterCondition->Id()<<") "<<std::endl;
		      ContactMasterCondition->Set(TO_ERASE);
		      list_of_points.push_back(new_point);
		      list_of_conditions.push_back(ContactMasterCondition);
                    }
		}

		rLocalRefineInfo.number_of_contacts ++;
              }
              // else{

              //   std::cout<<"   Master Condition not found "<<std::endl;

              // }

              rLocalRefineInfo.number_of_contact_conditions ++;

	    }
        }

        // std::cout<<"   [ Contact Conditions : "<<rLocalRefineInfo.number_of_contact_conditions<<", (contacts in domain: "<<rLocalRefineInfo.number_of_contacts<<", of them active: "<<rLocalRefineInfo.number_of_active_contacts<<") ] "<<std::endl;
        // std::cout<<"   Contact Search End ["<<list_of_conditions.size()<<" : "<<list_of_points.size()<<"]"<<std::endl;

      return true;

      KRATOS_CATCH(" ")

    }

    //*******************************************************************************************
    //*******************************************************************************************

    bool RefineOtherBoundary(ModelPart& rModelPart, std::vector<Point >& list_of_points, std::vector<ConditionType::Pointer>& list_of_conditions, RefineCounters& rLocalRefineInfo )
    {

      KRATOS_TRY

      ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

      //***SIZES :::: parameters do define the tolerance in mesh size:

      //RIGID WALL CONTACT:
      double factor_for_tip_radius               = 0.20; //deformable contact tolerance in radius for detection tip sides to refine

      double size_for_wall_tip_contact_side      = 0.50 * mrRemesh.Refine->CriticalSide;
      double size_for_wall_semi_tip_contact_side = 0.75 * mrRemesh.Refine->CriticalSide; // semi contact or contact which a node in a tip
      double size_for_wall_non_tip_contact_side  = 1.50 * mrRemesh.Refine->CriticalSide; // semi contact or contact which no node in a tip

      //NON CONTACT:
      double size_for_energy_side                = 1.50 * mrRemesh.Refine->CriticalSide; // non contact side which dissipates energy
      double size_for_non_contact_side           = 3.50  * mrRemesh.Refine->CriticalSide;


      double tool_radius= 0;
      double side_length= 0;
      double plastic_power=0;

      bool radius_insert = false;
      bool energy_insert = false;
      bool mesh_size_insert = false;
      bool contact_active = false;
      bool contact_semi_active = false;
      bool tool_project = false;

      Node<3> new_point(0,0.0,0.0,0.0);

      std::vector<bool> semi_active_nodes;


      //LOOP TO CONSIDER ALL SUBDOMAIN CONDITIONS
      double cond_counter=0;
      for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(); ic!= rModelPart.ConditionsEnd(); ic++)
	{
	  cond_counter ++;
	  bool refine_candidate = false;
	  if( mrRemesh.Options.Is(MesherUtilities::CONSTRAINED) ){
	    if( ic->Is(BOUNDARY) ) //ONLY SET TO THE BOUNDARY SKIN CONDITIONS (CompositeCondition)
	      refine_candidate = true;
	    else
	      refine_candidate = false;
	  }
	  else{
	    refine_candidate = true;
	  }


	  if( refine_candidate ){
	    if (mrRemesh.Refine->RefiningBoxSetFlag == true ){
	      refine_candidate = mMesherUtilities.CheckConditionInBox(*(ic.base()),*(mrRemesh.Refine->RefiningBox),CurrentProcessInfo);
	    }
	  }


	  if( refine_candidate ){

	    radius_insert    = false;
	    energy_insert    = false;
	    mesh_size_insert = false;
	    tool_project     = false;
	    contact_active      = false;
	    contact_semi_active = false;

	    side_length   = 0;
	    tool_radius   = 0;
	    plastic_power = 0;

	    //double condition_radius = 0;
	    Geometry< Node<3> > rConditionGeometry;
	    array_1d<double,3> tip_center;

	    if( ic->IsNot(TO_ERASE) ){

	      //*********************************************************************************
	      // RIGID CONTACT CONDITIONS ON TIP START
	      //*********************************************************************************

	      // TOOL TIP INSERT;


	      // ACTIVE CONTACT DETECTION START

	      rConditionGeometry = ic->GetGeometry();
	      contact_active = mMesherUtilities.CheckContactActive(rConditionGeometry, contact_semi_active, semi_active_nodes);

	      // ACTIVE CONTACT DETECTION END


	      if( contact_active ){

		side_length = mMesherUtilities.CalculateSideLength (rConditionGeometry[0],rConditionGeometry[1]);

		if( side_length > size_for_wall_tip_contact_side ){

		  bool on_tip = false;
		  if(rConditionGeometry[0].Is(TO_SPLIT) && rConditionGeometry[1].Is(TO_SPLIT)){
		    on_tip = true;
		  }
		  else if (rConditionGeometry[0].Is(TO_SPLIT) || rConditionGeometry[1].Is(TO_SPLIT)){
		    if( side_length > size_for_wall_tip_contact_side ){
		      on_tip = true;
		    }
		  }

		  bool on_radius = false;

		  if( on_tip && mRigidWalls.size() ){

		    Vector Point(3);
		    if( rConditionGeometry[0].Is(TO_SPLIT) ){

		      Point[0] = rConditionGeometry[0].X();
		      Point[1] = rConditionGeometry[0].Y();
		      Point[2] = rConditionGeometry[0].Z();
		      on_radius = true;

		    }
		    else if( rConditionGeometry[1].Is(TO_SPLIT) ){

		      Point[0] = rConditionGeometry[1].X();
		      Point[1] = rConditionGeometry[1].Y();
		      Point[2] = rConditionGeometry[1].Z();
		      on_radius = true;

		    }
		    else{
		      on_radius = false;
		    }


		    if( on_radius ){

		      ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
		      double Time = CurrentProcessInfo[TIME];
		      for( unsigned int i = 0; i <mRigidWalls.size(); ++i )
			{
			  if(mRigidWalls[i]->IsInside( Point, Time ) ){
			    tool_radius =mRigidWalls[i]->GetRadius(Point);
			    tip_center  =mRigidWalls[i]->GetCenter(Point);
			    break;
			  }
			}
		    }

		  }

		  if( on_radius && on_tip ) //master node in tool -->  refine workpiece  // (tool_radius ==0 in workpiece nodes)
		    {
		      Node<3> center (0,tip_center[0],tip_center[1],tip_center[2]);
		      array_1d<double, 3 > radius;
		      radius[0]=rConditionGeometry[0].X()-center.X();
		      radius[1]=rConditionGeometry[0].Y()-center.Y();
		      radius[2]=rConditionGeometry[0].Z()-center.Z();

		      double distance1=norm_2(radius);

		      radius[0]=rConditionGeometry[1].X()-center.X();
		      radius[1]=rConditionGeometry[1].Y()-center.Y();
		      radius[2]=rConditionGeometry[1].Z()-center.Z();

		      double distance2=norm_2(radius);


		      if( ((1-factor_for_tip_radius)*tool_radius < distance1 &&  (1-factor_for_tip_radius)*tool_radius < distance2) &&
			  (distance1 < (1+factor_for_tip_radius)*tool_radius  &&  distance2 < (1+factor_for_tip_radius)*tool_radius) )
			{
			  radius_insert = true;
			}

		    }

		}

	      }

	      //*********************************************************************************
	      // RIGID CONTACT CONDITIONS ON TIP END
	      //*********************************************************************************


	      //*********************************************************************************
	      // FREE BOUNDARY CONDITIONS ENERGY INSERTION START
	      //*********************************************************************************


	      // ENERGY INSERT

	      unsigned int vsize=ic->GetValue(MASTER_ELEMENTS).size();

	      if (!radius_insert && mrRemesh.Refine->RefiningOptions.Is(MesherUtilities::REFINE_BOUNDARY_ON_THRESHOLD) && vsize>0){

		Element::ElementType& MasterElement = ic->GetValue(MASTER_ELEMENTS)[vsize-1];

		plastic_power=0;
		std::vector<double> Value(1);

		MasterElement.GetValueOnIntegrationPoints(mrRemesh.Refine->GetThresholdVariable(),Value,CurrentProcessInfo);


		Geometry<Node<3> >& pGeom = MasterElement.GetGeometry();
		plastic_power = Value[0] * pGeom.Area();

		// if( Value[0] > 0 )
		//   std::cout<<" plastic_power "<<plastic_power<<std::endl;

		//computation of the condition master element radius start:
		//PointsArrayType& vertices = pGeom.Points();

		// double average_side_length= mMesherUtilities.CalculateAverageSideLength (vertices[0].X(), vertices[0].Y(),
		// 							      vertices[1].X(), vertices[1].Y(),
		// 							      vertices[2].X(), vertices[2].Y());
		//condition_radius = pGeom.Area()/average_side_length;
		//computation of the condition master element radius end;

		//condition_radius is side_length
		side_length = mMesherUtilities.CalculateSideLength (rConditionGeometry[0],rConditionGeometry[1]);

		//condition_radius = mMesherUtilities.CalculateTriangleRadius (pGeom);

		//if( plastic_power > mrRemesh.Refine->CriticalDissipation && condition_radius > mrRemesh.Refine->CriticalRadius )
		if( plastic_power > mrRemesh.Refine->ReferenceThreshold * pGeom.Area() && side_length > size_for_energy_side )
		  {
		    energy_insert = true;
		  }
	      }


	      //*********************************************************************************
	      // FREE BOUNDARY CONDITIONS ENERGY INSERTION END
	      //*********************************************************************************


	      //*********************************************************************************
	      // FREE BOUNDARY CONDITIONS SIZE INSERTION
	      //*********************************************************************************


	      // BOUNDARY SIZE INSERT

	      if( (!radius_insert || !energy_insert) && vsize>0 ){

		Element::ElementType& MasterElement = ic->GetValue(MASTER_ELEMENTS)[vsize-1];

		//std::cout<<" MASTER_ELEMENT "<<MasterElement.Id()<<std::endl;


		Geometry<Node<3> >& vertices = MasterElement.GetGeometry();
		double Alpha =  mrRemesh.Refine->Alpha;

		bool accepted = mMesherUtilities.AlphaShape(Alpha,vertices,2);


		//condition_radius is side_length
		side_length = mMesherUtilities.CalculateSideLength (rConditionGeometry[0],rConditionGeometry[1]);

		//condition_radius = mMesherUtilities.CalculateTriangleRadius (pGeom);
		double critical_side_size = 0;

		bool on_tip = false;
		if( contact_semi_active ){

		  if (rConditionGeometry[0].Is(TO_SPLIT) || rConditionGeometry[1].Is(TO_SPLIT))
		    on_tip = true;

		  if( on_tip == true )
		    critical_side_size = size_for_wall_semi_tip_contact_side;
		  else
		    critical_side_size = size_for_wall_non_tip_contact_side;



		}
		else if( contact_active ){

		  if (rConditionGeometry[0].Is(TO_SPLIT) || rConditionGeometry[1].Is(TO_SPLIT))
		    on_tip = true;

		  if( on_tip == true )
		    critical_side_size = size_for_wall_semi_tip_contact_side;
		  else
		    critical_side_size = size_for_wall_non_tip_contact_side;

		}
		else{

		  critical_side_size  = size_for_non_contact_side;
		}


		//if(plastic_power > mrRemesh.Refine->CriticalDissipation && condition_radius > mrRemesh.Refine->CriticalRadius)
		if( !accepted && !contact_semi_active && !contact_active && side_length > critical_side_size )
		  {
		    mesh_size_insert = true;

		    // std::cout<<"   insert on mesh size "<<std::endl;
		  }
		else if( contact_semi_active && side_length > critical_side_size )
		  {
		    mesh_size_insert = true;

		    // std::cout<<"   insert on mesh size semi_contact "<<std::endl;

		  }
		else if( contact_active && side_length > critical_side_size )
		  {
		    mesh_size_insert = true;

		    // std::cout<<"   insert on mesh size on contact "<<std::endl;

		  }


	      }


	      //*********************************************************************************
	      //                   BOUNDARY REBUILD START                                      //
	      //*********************************************************************************


	      if( radius_insert || energy_insert || mesh_size_insert ) //Boundary must be rebuild
		{

		  // std::cout<<"   BOUNDARY DOMAIN ELEMENT REFINED "<<ic->Id()<<std::endl;

		  new_point.X() = 0.5*( rConditionGeometry[1].X() + rConditionGeometry[0].X() );
		  new_point.Y() = 0.5*( rConditionGeometry[1].Y() + rConditionGeometry[0].Y() );
		  new_point.Z() = 0.5*( rConditionGeometry[1].Z() + rConditionGeometry[0].Z() );

		  if( this->mEchoLevel > 0 ){
		    std::cout<<"   radius_insert "<<radius_insert<<" energy_insert "<<energy_insert<<" mesh_size_insert "<<mesh_size_insert<<std::endl;
		    std::cout<<"   NEW NODE  "<<new_point<<std::endl;
		  }

		  new_point.SetId(ic->Id()); //set condition Id

		  //it will be good if the node is detected in the tool tip using the rigid contact standards:

		  if( (rConditionGeometry[0].Is(TO_SPLIT) && rConditionGeometry[1].Is(TO_SPLIT)) )
		    tool_project = true;

		  if( (rConditionGeometry[0].Is(TO_SPLIT) || rConditionGeometry[1].Is(TO_SPLIT)) && contact_active)
		    tool_project = true;

		  if( (rConditionGeometry[0].Is(TO_SPLIT) || rConditionGeometry[1].Is(TO_SPLIT)) && contact_semi_active)
		    tool_project = true;

		  bool on_radius = false;

		  if( tool_project && mRigidWalls.size() ){

		    Vector Point(3);
		    if( rConditionGeometry[0].Is(TO_SPLIT) ){

		      Point[0] = rConditionGeometry[0].X();
		      Point[1] = rConditionGeometry[0].Y();
		      Point[2] = rConditionGeometry[0].Z();
		      on_radius = true;

		    }
		    else if( rConditionGeometry[1].Is(TO_SPLIT) ){

		      Point[0] = rConditionGeometry[1].X();
		      Point[1] = rConditionGeometry[1].Y();
		      Point[2] = rConditionGeometry[1].Z();
		      on_radius = true;

		    }
		    else{
		      on_radius = false;
		    }


		    if( on_radius ){

		      ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
		      double Time = CurrentProcessInfo[TIME];
		      for( unsigned int i = 0; i <mRigidWalls.size(); ++i )
			{
			  if(mRigidWalls[i]->IsInside( Point, Time ) ){
			    tool_radius =mRigidWalls[i]->GetRadius(Point);
			    tip_center  =mRigidWalls[i]->GetCenter(Point);
			    break;
			  }
			}
		    }




		    if(on_radius){

		      if(new_point.Y()<(tip_center[1]) && new_point.Y()>(tip_center[1]-tool_radius)){

			array_1d<double,3> tip_normal = tip_center-new_point;

			if(norm_2(tip_normal)<tool_radius){ //if is in the tool tip

			  tip_normal -= (tool_radius/norm_2(tip_normal)) * tip_normal;

			  if(norm_2(tip_normal)<tool_radius*0.08)
			    new_point  += tip_normal;
			}

		      }

		      if( this->mEchoLevel > 0 )
			std::cout<<"   TOOL PROJECT::on radius  "<<new_point<<std::endl;


		    }


		  }

		  if(radius_insert)
		    rLocalRefineInfo.number_of_tip_bounds ++;
		  if(energy_insert)
		    rLocalRefineInfo.number_of_energy_bounds ++;
		  if(mesh_size_insert)
		    rLocalRefineInfo.number_of_exterior_bounds++;

		  ic->Set(TO_ERASE);

		  if( this->mEchoLevel > 0 )
		    std::cout<<"   INSERTED NODE  "<<new_point<<std::endl;

		  list_of_points.push_back(new_point);
		  list_of_conditions.push_back(*(ic.base()));


		  // if( this->mEchoLevel > 0 ){
		  //   std::cout<<"   Refine Boundary  (Id:"<<ic->Id()<<"): ["<<rConditionGeometry[0].Id()<<", "<<rConditionGeometry[1].Id()<<"]"<<std::endl;
		  //   std::cout<<"   (x1:"<<rConditionGeometry[0].X()<<", y1: "<<rConditionGeometry[0].Y()<<") "<<" (x2:"<<rConditionGeometry[1].X()<<", y2: "<<rConditionGeometry[1].Y()<<") "<<std::endl;

		  //   //std::cout<<" Added Node [Rcrit:"<<condition_radius<<",Scrit:"<<side_length<<",PlasticPower:"<<plastic_power<<"]"<<std::endl;
		  //   std::cout<<"   Added Node [Scrit:"<<side_length<<",PlasticPower:"<<plastic_power<<"]"<<std::endl;
		  //   //std::cout<<"   Conditions [Rcrit:"<<mrRemesh.Refine->CriticalRadius<<",Scrit:"<<mrRemesh.Refine->CriticalSide<<",PlasticPower:"<<mrRemesh.Refine->ReferenceThreshold<<"]"<<std::endl;
		  //   std::cout<<"   Conditions [Scrit:"<<mrRemesh.Refine->CriticalSide<<",PlasticPower:"<<mrRemesh.Refine->ReferenceThreshold<<"]"<<std::endl;
		  // }

		}


	      //*********************************************************************************
	      //                   BOUNDARY REBUILD END                                        //
	      //*********************************************************************************


	    }
	    else{
	      if( this->mEchoLevel > 0 )
		std::cout<<" Condition "<<ic->Id()<<" Released "<<std::endl;
	    }

	  }
        }

      return true;

      KRATOS_CATCH(" ")

    }

    /// Assignment operator.
    RefineConditionsInContactMesherProcess& operator=(RefineConditionsInContactMesherProcess const& rOther);


    /// this function is a private function


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
                                  RefineConditionsInContactMesherProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RefineConditionsInContactMesherProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_REFINE_CONDITIONS_IN_CONTACT_MESHER_PROCESS_H_INCLUDED  defined
