//
//   Project Name:                       MachiningApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                     May 2014 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes


// Project includes
#include "custom_modelers/contact_domain_2D_modeler.hpp"

#include "pfem_solid_mechanics_application_variables.h"


namespace Kratos
{



  //*******************************************************************************************
  //*******************************************************************************************

  void ContactDomain2DModeler::TransferContactBoundaryData(ModelPart& rModelPart, bool initial)
  {
    KRATOS_TRY

    //Needed to compute effective gaps in the contact domain with lagrangian multipliers
    MeshDataTransferUtilities    MeshDataTransfer;

    //be careful: it must be done once only after the step solution: 
    if(initial)
      MeshDataTransfer.TransferBoundaryData(rModelPart,MeshDataTransferUtilities::INITIALIZATION);
    else
      MeshDataTransfer.TransferBoundaryData(rModelPart,MeshDataTransferUtilities::MASTER_ELEMENT_TO_NODE);
	              	
    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void ContactDomain2DModeler::GenerateContactMesh (ModelPart& rModelPart,
						    Element   const& rReferenceElement,
						    Condition const& rReferenceCondition,
						    bool   ConstrainedFlag,
						    double AlphaParameter,
						    double SizeFactor,
						    double OffsetFactor,
						    double PenaltyParameter,
						    double StabilityParameter,
						    bool   FrictionFlag,
						    double StaticFrictionCoefficient,
						    double DynamicFrictionCoefficient)
    
  {
    KRATOS_TRY

    MeshModeler::MeshingVariables MeshingVars;
    ContactVariables ContactVars;
    
    MeshingVars.SetReferenceElement   (rReferenceElement);
    MeshingVars.SetReferenceCondition (rReferenceCondition);

    MeshingVars.Refine.SizeFactor  = SizeFactor;
    MeshingVars.AlphaParameter     = AlphaParameter;
    MeshingVars.OffsetFactor       = OffsetFactor;

    ContactVars.OffsetFactor       = OffsetFactor;
    //std::cout<<" OffsetFactor "<<ContactVars.OffsetFactor<<std::endl;

    ContactVars.PenaltyParameter   = PenaltyParameter;
    ContactVars.StabilityParameter = StabilityParameter;

    if(FrictionFlag)
      ContactVars.FrictionFlag = 1;
    else
      ContactVars.FrictionFlag = 0;

    ContactVars.StaticFrictionCoefficient  = StaticFrictionCoefficient;
    ContactVars.DynamicFrictionCoefficient = DynamicFrictionCoefficient;

    //Sort Conditions
    unsigned int consecutive_index = 1;
    for(ModelPart::ConditionsContainerType::iterator it = rModelPart.ConditionsBegin(); it!=rModelPart.ConditionsEnd(); it++)
      it->SetId(consecutive_index++);
	
    //Update Boundary Normals before Contact Search
    BoundaryNormalsCalculationUtilities BoundaryComputation;
    BoundaryComputation.CalculateBoundaryNormals(rModelPart, mEchoLevel);

    rModelPart.Conditions().Sort();
    rModelPart.Conditions().Unique();

	
    std::cout<<" --------------                     -------------- "<<std::endl;
    std::cout<<" --------------      BOUND          -------------- "<<std::endl;
	
    //Generate Delaunay Triangulation for Contact Detection
    if(!ConstrainedFlag)
      GenerateContactDT(rModelPart,MeshingVars,ContactVars);
    else
      GenerateContactCDT(rModelPart,MeshingVars,ContactVars);
	
    std::cout<<" --------------                     -------------- "<<std::endl;
    std::cout<<" -------------- CONTACT SEARCH DONE -------------- "<<std::endl;
    std::cout<<" --------------                     -------------- "<<std::endl;


    //Renumerate conditions to add in the end of the Elements array (for writing purposes in the ID)
    //rModelPart.Elements().Sort();
    int LastElementId   = (rModelPart.Elements().end()-1)->Id();
    int LastConditionId = (rModelPart.Conditions().end()-1)->Id();
    if(LastElementId>LastConditionId){
      consecutive_index = LastElementId+1;
    }
    else{
      consecutive_index = LastConditionId+1;
    }
	
    for(ModelPart::ConditionsContainerType::iterator it = rModelPart.ConditionsBegin(); it!=rModelPart.ConditionsEnd(); it++){
      if(it->Is(CONTACT)){
	it->SetId(consecutive_index); 
	consecutive_index++;
      }
    }
    
    KRATOS_CATCH(" ")
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void ContactDomain2DModeler::GenerateContactDT(ModelPart& rModelPart,
					 MeshModeler::MeshingVariables& rMeshingVariables,
					 ContactVariables& rContactVariables)
  {

    KRATOS_TRY


    ClearContactConditions(rModelPart);

    //*********************************************************************
    struct triangulateio in;
    struct triangulateio out;

    ModelPart::NodesContainerType  BoundaryNodes;

    
    ////////////////////////////////////////////////////////////
    SetTriangulationShrankNodes(rModelPart,BoundaryNodes,rMeshingVariables,rContactVariables,in,out);
    ////////////////////////////////////////////////////////////
    
    
    //*********************************************************************

    rMeshingVariables.MeshingOptions.Set(MeshModeler::RECONNECT);
    ////////////////////////////////////////////////////////////
    this->GenerateTriangulation(rMeshingVariables.MeshingOptions,rMeshingVariables.RefiningOptions,in,out);
    ////////////////////////////////////////////////////////////
    rMeshingVariables.MeshingOptions.Reset(MeshModeler::RECONNECT);

    
    //*********************************************************************

    ////////////////////////////////////////////////////////////
    BuildContactConditions(rModelPart,BoundaryNodes,rMeshingVariables,rContactVariables,out);
    ////////////////////////////////////////////////////////////

    //free the rest of the memory
    DeletePointsList(in);
    DeleteTrianglesList(out);


    KRATOS_CATCH( "" )
      }

  //*******************************************************************************************
  //*******************************************************************************************

  void ContactDomain2DModeler::GenerateContactCDT(ModelPart& rModelPart,
						  MeshModeler::MeshingVariables& rMeshingVariables,
						  ContactVariables & rContactVariables)
  {

    KRATOS_TRY

    ClearContactConditions(rModelPart);

    //*********************************************************************
    struct triangulateio in;
    struct triangulateio out;

    ModelPart::NodesContainerType  BoundaryNodes;

    rMeshingVariables.MeshingOptions.Set(MeshModeler::CONSTRAINED_MESH);
    ////////////////////////////////////////////////////////////
    SetTriangulationShrankNodes(rModelPart,BoundaryNodes,rMeshingVariables,rContactVariables,in,out);
    ////////////////////////////////////////////////////////////
    
    //*********************************************************************

    rMeshingVariables.MeshingOptions.Set(MeshModeler::RECONNECT);

    int fail = this->GenerateTriangulation(rMeshingVariables.MeshingOptions,rMeshingVariables.RefiningOptions,in,out);

    
    if(fail){
      rMeshingVariables.MeshingOptions.Reset(MeshModeler::CONSTRAINED_MESH);
      ////////////////////////////////////////////////////////////
      fail = this->GenerateTriangulation(rMeshingVariables.MeshingOptions,rMeshingVariables.RefiningOptions,in, out);
      ////////////////////////////////////////////////////////////
      rMeshingVariables.MeshingOptions.Set(MeshModeler::CONSTRAINED_MESH);
    }

    rMeshingVariables.MeshingOptions.Reset(MeshModeler::RECONNECT);

    //*********************************************************************

    ////////////////////////////////////////////////////////////
    BuildContactConditions(rModelPart,BoundaryNodes,rMeshingVariables,rContactVariables,out);
    ////////////////////////////////////////////////////////////

    rMeshingVariables.MeshingOptions.Reset(MeshModeler::CONSTRAINED_MESH);

    //free the rest of the memory
    DeletePointsList(in);
    DeleteTrianglesList(out);


    KRATOS_CATCH( "" )
      }


  //*******************************************************************************************
  //*******************************************************************************************

  void ContactDomain2DModeler::SetTriangulationShrankNodes(ModelPart& rModelPart,
							   ModelPart::NodesContainerType& rBoundaryNodes,
							   MeshModeler::MeshingVariables& rMeshingVariables,
							   ContactVariables& rContactVariables,
							   struct triangulateio& in,
							   struct triangulateio& out)
  {

    KRATOS_TRY

    //writing the points coordinates in a vector and reordening the Id's from an initial id
    ModelPart::NodesContainerType::iterator nodes_begin = rModelPart.NodesBegin();
	
    const array_1d<double,3> ZeroOffset(3,0.0);

    for(unsigned int i = 0; i<rModelPart.Nodes().size(); i++)
      {
	if((nodes_begin + i)->Is(BOUNDARY) ){ //&& (nodes_begin + i)->IsNot(STRUCTURE)){
	  rBoundaryNodes.push_back( *((nodes_begin+i).base()) ); 
	}
	    
	noalias((nodes_begin + i)->FastGetSolutionStepValue(OFFSET)) = ZeroOffset;
      }


    //*********************************************************************
    ModelPart::NodesContainerType::iterator boundary_nodes_begin = rBoundaryNodes.begin();

    double nodal_h_min=1e20;
    double nodal_h=0;

    for(unsigned int i = 0; i<rBoundaryNodes.size(); i++)
      {
	nodal_h=(boundary_nodes_begin + i)->FastGetSolutionStepValue(NODAL_H); 
	
	if(nodal_h_min>nodal_h)  
	  nodal_h_min=nodal_h;
	  
      }

    double hnodal_offset_conversion = 0.35;
    if( rContactVariables.OffsetFactor > nodal_h_min*hnodal_offset_conversion || rContactVariables.OffsetFactor < nodal_h_min*0.01){
      rContactVariables.OffsetFactor = nodal_h_min*hnodal_offset_conversion;
      rMeshingVariables.OffsetFactor = rContactVariables.OffsetFactor;
    }

    std::cout<<" nodal_h_min "<<nodal_h_min<<" OffsetFactor "<<rContactVariables.OffsetFactor<<std::endl;

    //*********************************************************************

    ClearTrianglesList(in);
    ClearTrianglesList(out);

    //*********************************************************************
    if(rMeshingVariables.MeshingOptions.Is(MeshModeler::CONSTRAINED_MESH)){
      
      std::cout<<"  Constrained Contact Meshing "<<std::endl;
      
      //PART 1: node list
      double extra_radius = rContactVariables.OffsetFactor*4; 
      SpatialBoundingBox DomainBox (rModelPart, extra_radius);
      
      ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
      std::vector<Vector> BoxVertices = DomainBox.GetVertices( CurrentProcessInfo[TIME] );

      //input mesh: NODES
      in.numberofpoints = rBoundaryNodes.size() + BoxVertices.size();
      in.pointlist      = new REAL[in.numberofpoints * 2];

      rMeshingVariables.NodalPreIds.resize(in.numberofpoints+1);
      rMeshingVariables.NodalPreIds[0]=0;
        
      int base=0;
      int id  =0;
      for(unsigned int i = 0; i<rBoundaryNodes.size(); i++)
	{
	  //from now on it is consecutive
	  rMeshingVariables.NodalPreIds[id+1]=(boundary_nodes_begin + i)->Id();
	  (boundary_nodes_begin + i)->SetId(id+1);

	  array_1d<double, 3>&  Normal=(boundary_nodes_begin + i)->FastGetSolutionStepValue(NORMAL); //BOUNDARY_NORMAL must be set as nodal variable
	  double Shrink = (boundary_nodes_begin + i)->FastGetSolutionStepValue(SHRINK_FACTOR);   //SHRINK_FACTOR   must be set as nodal variable

	  //std::cout<<" Normal "<<Normal<<" Shrink "<<Shrink<<std::endl;

	  //if normal not normalized
	  //Shrink /=norm_2(Normal);
	    
	  array_1d<double, 3>& Offset = (boundary_nodes_begin + i)->FastGetSolutionStepValue(OFFSET);

	  if( norm_2(Normal)!= 0 ){
	    Normal /= norm_2(Normal);
	  }
	  else{
	    std::cout<<" Boundary NORMAL is Zero in node ["<<(boundary_nodes_begin + i)->Id()<<"] : something is wrong with the normals search "<<std::endl;
	  };

	  Offset[0] = ( (-1) * Normal[0] * Shrink * rContactVariables.OffsetFactor );
	  Offset[1] = ( (-1) * Normal[1] * Shrink * rContactVariables.OffsetFactor );

	  //std::cout<<" Id "<<rMeshingVariables.NodalPreIds[i+1]<<" real ID "<<(boundary_nodes_begin + i)->Id()<<" Shrink "<<Shrink<<" Normal "<<Normal<<" Offset "<<Offset<<std::endl;

	  in.pointlist[base]   = (boundary_nodes_begin + i)->X() + Offset[0];
	  in.pointlist[base+1] = (boundary_nodes_begin + i)->Y() + Offset[1];

	  //std::cout<<" vertex["<<id<<"] "<< in.pointlist[base] <<" "<<in.pointlist[base+1]<<std::endl; 

	  base+=2;
	  id  +=1;
	}

      
      std::vector<int> vertices_ids;
      for(unsigned int i = 0; i<BoxVertices.size(); i++)
	{
	  rMeshingVariables.NodalPreIds[id+1] = -1;
	  vertices_ids.push_back (id+1);

	  in.pointlist[base]   = BoxVertices[i][0];
	  in.pointlist[base+1] = BoxVertices[i][1];
	  
	  //std::cout<<"  BoxVertices ["<<i<<"]= ("<<BoxVertices[i][0]<<", "<<BoxVertices[i][1]<<"). Id = "<<vertices_ids[i]<<std::endl;

	  base+=2; 
	  id  +=1;
	}

     
      //PART 2: faced list (we can have holes in facets != area holes)
      
      in.numberofsegments           = rModelPart.NumberOfConditions() + BoxVertices.size();
      in.segmentmarkerlist          = new int[in.numberofsegments];
      in.segmentlist                = new int[in.numberofsegments*2];

      ModelPart::ConditionsContainerType::iterator conditions_begin = rModelPart.ConditionsBegin();
     
     
      base = 0;
      for(unsigned int i = 0; i<rModelPart.Conditions().size(); i++)
	{
	  Geometry< Node<3> >& rGeometry = (conditions_begin + i)->GetGeometry();
	  in.segmentlist[base]   = rGeometry[0].Id();
	  in.segmentlist[base+1] = rGeometry[1].Id();
	 
	  base+=2;
	}  

      for(unsigned int i = 0; i<BoxVertices.size()-1; i++) //2d (rectangular box of 4 sides)
	{
	  in.segmentlist[base]   = vertices_ids[i];
	  in.segmentlist[base+1] = vertices_ids[i+1];
	 
	  base+=2;
	}  

      in.segmentlist[base]   = vertices_ids[BoxVertices.size()-1];
      in.segmentlist[base+1] = vertices_ids[0];
      
      //PART 3: (area) hole list

      std::vector<Vector> Holes = DomainBox.GetHoles( rModelPart );

      //holes
      in.numberofholes              = Holes.size();
      in.holelist                   = new REAL[in.numberofholes * 2];

      for(unsigned int hl=0; hl<Holes.size();hl++)
	{
	  //inside point of the hole:
	  in.holelist[hl*2+0] = Holes[hl][0];
	  in.holelist[hl*2+1] = Holes[hl][1];
	}


      //PART 4: region attributes list
      in.numberofregions =0;
      in.regionlist =NULL;

    }
    else{

      //input mesh: NODES
      in.numberofpoints = rBoundaryNodes.size();
      in.pointlist      = new REAL[in.numberofpoints * 2];

      rMeshingVariables.NodalPreIds.resize(in.numberofpoints+1);
      rMeshingVariables.NodalPreIds[0]=0;
        
      int base=0;
      for(unsigned int i = 0; i<rBoundaryNodes.size(); i++)
	{
	  //from now on it is consecutive
	  rMeshingVariables.NodalPreIds[i+1]=(boundary_nodes_begin + i)->Id();
	  (boundary_nodes_begin + i)->SetId(i+1);

	  array_1d<double, 3>&  Normal=(boundary_nodes_begin + i)->FastGetSolutionStepValue(NORMAL); //BOUNDARY_NORMAL must be set as nodal variable
	  double Shrink = (boundary_nodes_begin + i)->FastGetSolutionStepValue(SHRINK_FACTOR);   //SHRINK_FACTOR   must be set as nodal variable


	  //if normal not normalized
	  //Shrink /=norm_2(Normal);
	    
	  array_1d<double, 3>& Offset = (boundary_nodes_begin + i)->FastGetSolutionStepValue(OFFSET);

	  Normal /= norm_2(Normal);
	  Offset[0] = ( (-1) * Normal[0] * Shrink * rContactVariables.OffsetFactor );
	  Offset[1] = ( (-1) * Normal[1] * Shrink * rContactVariables.OffsetFactor );

	  //std::cout<<" Id "<<rMeshingVariables.NodalPreIds[i+1]<<" real ID "<<(boundary_nodes_begin + i)->Id()<<" Shrink "<<Shrink<<" Normal "<<Normal<<" Offset "<<Offset<<std::endl;

	  in.pointlist[base]   = (boundary_nodes_begin + i)->X() + Offset[0];
	  in.pointlist[base+1] = (boundary_nodes_begin + i)->Y() + Offset[1];

	  base+=2;
	}

    }
    //*********************************************************************

    KRATOS_CATCH( "" )
  }


  //*******************************************************************************************
  //*******************************************************************************************

	
  //Set contact elements in model_part after the Delaunay Tesselation
  void ContactDomain2DModeler::BuildContactConditions(ModelPart& rModelPart,
						      ModelPart::NodesContainerType& rBoundaryNodes,
						      MeshModeler::MeshingVariables& rMeshingVariables,
						      ContactVariables& rContactVariables,
						      struct triangulateio& out)
  {
    KRATOS_TRY

    //*******************************************************************
    //selecting elements
    rMeshingVariables.RefiningOptions.Set(MeshModeler::SELECT_ELEMENTS);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::CONTACT_SEARCH);

    this->SelectMeshElements(rBoundaryNodes,rMeshingVariables,out);

    rMeshingVariables.RefiningOptions.Reset(MeshModeler::SELECT_ELEMENTS);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::CONTACT_SEARCH);

    //*******************************************************************
    //setting new elements
    //

    //properties to be used in the generation
    //int number_properties = rModelPart.NumberOfProperties();
    //std::cout<<" Mesh 0 Id: "<<rModelPart.GetMesh().pGetProperties(0)->Id()<<" number props "<<number_properties<<std::endl;


    PropertiesContainerType::ContainerType PropertiesArray = rModelPart.GetMesh().PropertiesArray();
    PropertiesArray[0] = PropertiesType::Pointer(new PropertiesType(0));
    PropertiesArray[0]->Data() =PropertiesArray[1]->Data();

    //Properties 0 in order to change the Id to 0 and then write contact elements in another layer
    Properties::Pointer properties = PropertiesArray[0];
	
    properties->SetValue(THICKNESS,PropertiesArray[0]->GetValue(THICKNESS));
    properties->SetValue(PENALTY_PARAMETER,rContactVariables.PenaltyParameter);
    properties->SetValue(TAU_STAB,rContactVariables.StabilityParameter);
    properties->SetValue(FRICTION_ACTIVE,rContactVariables.FrictionFlag);
    properties->SetValue(MU_STATIC,rContactVariables.StaticFrictionCoefficient);
    properties->SetValue(MU_DYNAMIC,rContactVariables.DynamicFrictionCoefficient);

	
    // for(unsigned int p=0; p<PropertiesArray.size(); p++)
    //   std::cout<<" Properties ID ["<<p<<"]: "<<PropertiesArray[p]->Id()<<std::endl;


    std::cout<<"   Properties have been SET: ["<<properties->Id()<<"] "<<std::endl;
	
	
    //node list can not be changed from now on
    //ModelPart::NodesContainerType& model_nodes = rModelPart.Nodes();
    //ModelPart::NodesContainerType& model_nodes = rBoundaryNodes;

    //restore global ID's
    int idr=0;
    for(ModelPart::NodesContainerType::iterator in = rBoundaryNodes.begin() ; in != rBoundaryNodes.end() ; in++)
      {
	idr= in->Id();

	// std::cout<<" current id "<<idr<<" total IDs "<<rMeshingVariables.NodalPreIds.size()<<std::endl;
	// std::cout<<" previous id "<<rMeshingVariables.NodalPreIds[ idr ]<<std::endl;
	in->SetId( rMeshingVariables.NodalPreIds[ idr ] );

	// double& Shrink=in->FastGetSolutionStepValue(SHRINK_FACTOR);   //SHRINK_FACTOR   must be set as nodal variable
	// array_1d<double, 3>&  Normal=in->FastGetSolutionStepValue(NORMAL); 
	// array_1d<double, 3>&  Offset=in->FastGetSolutionStepValue(OFFSET); 
	//std::cout<<" Id "<<in->Id()<<" real ID "<<idr<<" Shrink "<<Shrink<<" Normal "<<Normal<<" Offset "<<Offset<<std::endl;


      }

    Condition const & rReferenceCondition=rMeshingVariables.GetReferenceCondition(); //contact element
    //Condition const & rReferenceCondition=KratosComponents<Condition>::Get("ContactDomain2DCondition");
    
    std::cout<<"   [START contact Element Generation "<<std::endl;

    //generate kratos elements (conditions are not touched)
    int contact_Id = rModelPart.Conditions().size();
    int id = contact_Id;
    for(int el = 0; el< out.numberoftriangles; el++)
      {
	if(rMeshingVariables.PreservedElements[el])
	  {
	    Geometry<Node<3> > vertices;
	    for(int pn=0; pn<3; pn++)		   
	      {
		//note that out.trianglelist, starts from node 1, not from node 0, it can be directly assigned to rMeshingVariables.NodalPreIds.
		//vertices.push_back( *((model_nodes).find( rMeshingVariables.NodalPreIds[out.trianglelist[el*3+pn]] ).base() ) );
		vertices.push_back(rModelPart.pGetNode(rMeshingVariables.NodalPreIds[out.trianglelist[el*3+pn]]));
		vertices.back().Set(CONTACT);
	      }

	    id += 1;
		
	    Condition::Pointer p_contact_cond = rReferenceCondition.Create(id, vertices, properties);
		
	    //search rModelPart.Conditions() associated to this contact element
	    //assign the MASTER ELEMENT and the MASTER NODE
				
	    bool condition_found=false;
	    Condition::Pointer p_master_cond = this->mModelerUtilities.FindMasterCondition(p_contact_cond,rModelPart.Conditions(),condition_found);

	    if(!condition_found){
	      std::cout<<" MASTER CONDITION NOT FOUND: Contact Element Release "<<std::endl;
	      id -= 1;
	    }
	    else{

	      //std::cout<<" contact master "<<std::endl;
	      //p_master_cond->GetValue(MASTER_ELEMENTS)[0].PrintInfo(std::cout);
	      //p_master_cond->GetValue(MASTER_ELEMENTS)[0].PrintData(std::cout);
	      //p_master_cond->GetValue(MASTER_ELEMENTS)[0].GetProperties().PrintData(std::cout);
	      //std::cout<<std::endl;
	      p_contact_cond->SetValue(MASTER_CONDITION, p_master_cond );
	      p_contact_cond->SetValue(MASTER_ELEMENTS, p_master_cond->GetValue(MASTER_ELEMENTS) );
	      p_contact_cond->SetValue(MASTER_NODES, p_master_cond->GetValue(MASTER_NODES) );
	      p_contact_cond->SetValue(NORMAL, p_master_cond->GetValue(NORMAL) );
	      p_contact_cond->Set(CONTACT);

	      //setting new elements
	      (rModelPart.Conditions()).push_back(p_contact_cond);
	    }
	    // if( (rModelPart.Conditions()).back().Is(CONTACT) )
	    //     std::cout<<" IS a CONTACT CONDITION "<<std::endl;
	  }
      }

    std::cout<<"   [END   contact Elements Generation ["<<id-contact_Id<<"] ]"<<std::endl;

    std::cout<<"   Total Conditions AFTER: ["<<rModelPart.Conditions().size()<<"] ];"<<std::endl;

    KRATOS_CATCH( "" )
  }

  //*******************************************************************************************
  //*******************************************************************************************


 
  void ContactDomain2DModeler::ClearContactConditions(ModelPart& rModelPart)
  {

    KRATOS_TRY

    //*******************************************************************
    //clearing elements
    //

    std::cout<<" Total Conditions BEFORE: ["<<rModelPart.Conditions().size()<<"]"<<std::endl;

    ModelPart::ConditionsContainerType RemoveConditions;

    int contact_Id=0;
    for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(); ic!= rModelPart.ConditionsEnd(); ic++)
      {

	if(ic->IsNot(CONTACT)){
	  contact_Id+=1;
	  RemoveConditions.push_back(*(ic.base()));
	  RemoveConditions.back().SetId(contact_Id);
	}
      }
	
    rModelPart.Conditions().swap(RemoveConditions);
	
    // computationally slow:
    // std::vector<int> RemoveConditions;
    // for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(); ic!= rModelPart.ConditionsEnd(); ic++)
    // {
    //     if(ic->Is(CONTACT)){
    // 	//std::cout<<" Remove Condition "<<ic->Id()<<std::endl;
    // 	RemoveConditions.push_back(ic->Id());
    //     }
    // }
	
    // for(unsigned int i=0; i<RemoveConditions.size(); i++)
    //     rModelPart.RemoveCondition(RemoveConditions[i]);

    rModelPart.Conditions().Sort();
    rModelPart.Conditions().Unique();

    std::cout<<" Total Conditions CLEAN: ["<<rModelPart.Conditions().size()<<"]"<<std::endl;

    KRATOS_CATCH( "" )

  }



} // Namespace Kratos

