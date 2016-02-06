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
#include "custom_modelers/laplacian_smoothing.hpp"
#include "custom_modelers/mesh_modeler.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

  /**
   * Flags related to the meshing parameters
   */

  //meshing options
  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, REMESH,               0 );
  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, RECONNECT,            1 );
  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, REFINE_MESH,          2 );

  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, CONSTRAINED_MESH,     3 );
  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, BOUNDARIES_SEARCH,    4 );
  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, NEIGHBOURS_SEARCH,    5 );

  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, SET_DOF,              6 );

  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, CONTACT_SEARCH,       7 );


  //refining options
  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, SELECT_ELEMENTS,      0 );
  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, PASS_ALPHA_SHAPE,     1 );

  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, REFINE_INSERT_NODES,  2 );
  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, REFINE_ADD_NODES,     3 );

  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, REFINE_ELEMENTS,      4 );
  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, REFINE_BOUNDARY,      5 );

  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, REMOVE_NODES,         6 );
  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, REMOVE_ON_BOUNDARY,   7 );

  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, CRITERION_ERROR,      8 );
  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, CRITERION_ENERGY,     9 );
  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, CRITERION_DISTANCE,  10 );

  KRATOS_CREATE_LOCAL_FLAG ( MeshModeler, ENGAGED_NODES,       11 );




  //*******************************************************************************************
  //*******************************************************************************************

  void MeshModeler::Initialize ( int NumberOfDomains )
  {
    KRATOS_TRY
    
    if( mEchoLevel > 0 )
      std::cout<<" INITIALIZE MESH DATA: [ number of domains = "<<NumberOfDomains<<" ]"<<std::endl;
    
    MeshingVariables Variables;
    Variables.Initialize();
    mMeshingVariables.push_back(Variables);
    mMeshingVariables.back().Initialize();
   
    if(NumberOfDomains>1){
      for(int i=1; i<=NumberOfDomains; i++)
	{
	  MeshingVariables Variables;
	  Variables.Initialize();
	  mMeshingVariables.push_back(Variables);
	  mMeshingVariables.back().Initialize();
	} 
    }

    KRATOS_CATCH(" ")
  }
  
  //*******************************************************************************************
  //*******************************************************************************************

  void MeshModeler::SetRemeshData( Element   const& rReferenceElement,
				   Condition const& rReferenceCondition,
				   bool   RemeshFlag,
				   bool   ConstrainedFlag,
				   bool   MeshSmoothingFlag,
				   bool   JacobiSmoothingFlag,
				   bool   AvoidTipElementsFlag,
				   double AlphaParameter,
				   double OffsetFactor,
				   int    MeshId )
  {
    KRATOS_TRY
    
    if(mMeshingVariables.size() == 1 && MeshId != 0)
      std::cout<<" Something wrong with mesh ID "<<MeshId<<std::endl;

    mMeshingVariables[MeshId].SetReferenceElement   (rReferenceElement);

    mMeshingVariables[MeshId].SetReferenceCondition (rReferenceCondition);

    mMeshingVariables[MeshId].RemeshFlag            = RemeshFlag;

    mMeshingVariables[MeshId].ConstrainedFlag       = ConstrainedFlag;

    mMeshingVariables[MeshId].MeshSmoothingFlag     = MeshSmoothingFlag;

    mMeshingVariables[MeshId].JacobiSmoothingFlag   = JacobiSmoothingFlag;

    mMeshingVariables[MeshId].AvoidTipElementsFlag  = AvoidTipElementsFlag;

    mMeshingVariables[MeshId].AlphaParameter        = AlphaParameter;

    mMeshingVariables[MeshId].OffsetFactor          = OffsetFactor;

    mMeshingVariables[MeshId].RefineFlag            = false;

    mMeshingVariables[MeshId].BoundingBox.IsSetFlag = false;

    if( mEchoLevel > 0 )
      std::cout<<" SetRemeshData ["<<MeshId<<"]: [ RefineFlag: "<<mMeshingVariables[MeshId].RefineFlag<<" RemeshFlag : "<<mMeshingVariables[MeshId].RemeshFlag<<" ] "<<std::endl;

    KRATOS_CATCH(" ")
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void MeshModeler::SetRefineData( bool   RefineFlag,
				   double SizeFactor,
				   double Radius,
				   const Variable<double>& DissipationVariable,
				   double Dissipation,
				   const Variable<double>& ErrorVariable,
				   double Error,
				   int    MeshId )
  {    
    KRATOS_TRY

    if(mMeshingVariables.size() == 1 && MeshId != 0)
      std::cout<<" Something wrong with mesh ID "<<MeshId<<std::endl;

    //other reference variables
    double SideTolerance = 3; //6;

    mMeshingVariables[MeshId].RefineFlag                 = RefineFlag;

    mMeshingVariables[MeshId].Refine.SizeFactor          = SizeFactor; //nodal_h
    
    //std::cout<<" INPUT VARIABLE ["<<MeshId<<"]"<<DissipationVariable<<std::endl;    

    mMeshingVariables[MeshId].Refine.SetDissipationVariable(DissipationVariable); 
    //mMeshingVariables[MeshId].Refine.SetDissipationVariable(PLASTIC_DISSIPATION);
    //mMeshingVariables[MeshId].Refine.SetDissipationVariable(PLASTIC_STRAIN);

    //std::cout<<" VARIABLE ["<<MeshId<<"]"<<mMeshingVariables[MeshId].Refine.GetDissipationVariable()<<std::endl;

    mMeshingVariables[MeshId].Refine.CriticalDissipation = Dissipation; //40;  400;

    mMeshingVariables[MeshId].Refine.CriticalRadius      = Radius; //0.0004 m

    mMeshingVariables[MeshId].Refine.CriticalSide        = SideTolerance * Radius;  //0.02;

    //std::cout<<" INPUT VARIABLE ["<<MeshId<<"]"<<ErrorVariable<<std::endl;    

    mMeshingVariables[MeshId].Refine.SetErrorVariable(ErrorVariable);
    //mMeshingVariables[MeshId].Refine.SetErrorVariable(PLASTIC_STRAIN);
    //mMeshingVariables[MeshId].Refine.SetErrorVariable(NORM_ISOCHORIC_STRESS);

    //std::cout<<" VARIABLE ["<<MeshId<<"]"<<mMeshingVariables[MeshId].Refine.GetErrorVariable()<<std::endl;

    mMeshingVariables[MeshId].Refine.ReferenceError      = Error;  //2;

    if( mEchoLevel > 0 ){
      std::cout<<" CRITICAL RADIUS      : "<<mMeshingVariables[MeshId].Refine.CriticalRadius<<std::endl;
      std::cout<<" CRITICAL SIDE        : "<<mMeshingVariables[MeshId].Refine.CriticalSide<<std::endl;
      std::cout<<" CRITICAL DISSIPATION : "<<mMeshingVariables[MeshId].Refine.CriticalDissipation<<std::endl;
    }

    KRATOS_CATCH(" ")
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void MeshModeler::SetRigidWall ( SpatialBoundingBox::Pointer pRigidWall )
  {
    KRATOS_TRY

    for(unsigned int i=0; i<mMeshingVariables.size(); i++)
      {
	mMeshingVariables[i].RigidWallSetFlag =true;
	mMeshingVariables[i].RigidWalls.push_back(pRigidWall);
      }

    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void MeshModeler::SetRefiningBox ( double Radius,
				     Vector Center,
				     Vector Velocity )
    
  {
    KRATOS_TRY

    for(unsigned int i=0; i<mMeshingVariables.size(); i++)
      {
	mMeshingVariables[i].BoundingBox.IsSetFlag = true;
	mMeshingVariables[i].BoundingBox.Radius    = Radius;
	mMeshingVariables[i].BoundingBox.Center    = Center;
	mMeshingVariables[i].BoundingBox.Velocity  = Velocity;

	if( mEchoLevel > 1 )
	  std::cout<<" Bounding Box [ Radius: "<<Radius<<" Center: "<<Center<<" Velocity: "<<Velocity<<" ] "<<std::endl;

      }

    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void MeshModeler::SetRefiningBox ( double Radius,
				     Vector Center,
				     Vector Velocity,
				     int MeshId )
    
  {
    KRATOS_TRY

    mMeshingVariables[MeshId].BoundingBox.IsSetFlag = true;
    mMeshingVariables[MeshId].BoundingBox.Radius    = Radius;
    mMeshingVariables[MeshId].BoundingBox.Center    = Center;
    mMeshingVariables[MeshId].BoundingBox.Velocity  = Velocity;

    if( mEchoLevel > 1 )
      std::cout<<" Bounding Box [ Radius: "<<Radius<<" Center: "<<Center<<" Velocity: "<<Velocity<<" ] "<<std::endl;

    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************


  void MeshModeler::GenerateMesh ( ModelPart& rModelPart )
  {
    KRATOS_TRY

    unsigned int start=0;
    unsigned int NumberOfDomains = mMeshingVariables.size()-1;
    unsigned int NumberOfMeshes  = rModelPart.NumberOfMeshes();

    if( NumberOfDomains > NumberOfMeshes )
      std::cout<<" ERROR: MORE DOMAINS THAN MESHES "<<std::endl;

    if( NumberOfDomains <= NumberOfMeshes )
      NumberOfMeshes = NumberOfDomains;

    if(NumberOfMeshes>1) 
      start=1;

    // bool out_buffer_active = true;
    // std::streambuf* buffer = NULL;
    // if( mEchoLevel == 0 ){
    //   //std::cout<<" Deactivate cout "<<std::endl;
    //   buffer = std::cout.rdbuf();
    //   std::ofstream fout("/dev/null");
    //   std::cout.rdbuf(fout.rdbuf());
    //   //std::cout<<output(off,buffer);
    //   out_buffer_active = false;
    // }

    //By the way: set meshes options from bools
    for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
      {
	if( GetEchoLevel() > 0 ){
	  std::cout<<"  GetRemeshData ["<<MeshId<<"]: [ RefineFlag: "<<mMeshingVariables[MeshId].RefineFlag<<"; RemeshFlag : "<<mMeshingVariables[MeshId].RemeshFlag<<" ] "<<std::endl;
	}
	
	if(mMeshingVariables[MeshId].RemeshFlag)
	  rModelPart.GetMesh(MeshId).Set( REMESH );
	else
	  rModelPart.GetMesh(MeshId).Reset( REMESH );

	if(mMeshingVariables[MeshId].RefineFlag)
	  rModelPart.GetMesh(MeshId).Set( REFINE_MESH );
	
	if(mMeshingVariables[MeshId].ConstrainedFlag){
	  rModelPart.GetMesh(MeshId).Set( CONSTRAINED_MESH );
	  
	  //Update Boundary Normals before Constrained Meshing
	  // BoundaryNormalsCalculationUtilities BoundaryComputation;
	  // BoundaryComputation.CalculateBoundaryNormals(rModelPart);
	}
	
	// check mesh size introduced :: warning must be shown
	// if(!out_buffer_active)
	//   std::cout.rdbuf(buffer);
	
	mModelerUtilities.CheckCriticalRadius(rModelPart, mMeshingVariables[MeshId].Refine.CriticalRadius, MeshId);
		
	// if(!out_buffer_active){
	//   buffer = std::cout.rdbuf();
	//   std::ofstream fout("/dev/null");
	//   std::cout.rdbuf(fout.rdbuf());
	// }
	// check mesh size introduced :: warning must be shown

      }

    //Sort Conditions
    unsigned int consecutive_index = 1;
    for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(); ic!=rModelPart.ConditionsEnd(); ic++)
      ic->SetId(consecutive_index++);

    rModelPart.Conditions().Sort();
    rModelPart.Conditions().Unique();


    // //Sort Elements
    // consecutive_index = 1;
    // for(ModelPart::ElementsContainerType::iterator ie = rModelPart.ElementsBegin(); ie!=rModelPart.ElementsEnd(); ie++)
    //   ie->SetId(consecutive_index++);

    // rModelPart.Elements().Sort();
    // rModelPart.Elements().Unique();

    // //Sort Nodes, set STRUCTURE NODES (RIGID TOOL NODES) AT END
    // consecutive_index = 1;
    // unsigned int reverse_index = rModelPart.Nodes().size();
    // for(ModelPart::NodesContainerType::iterator in = rModelPart.NodesBegin(); in!=rModelPart.NodesEnd(); in++){
    //   if(in->IsNot(STRUCTURE) ){
    //     in->SetId(consecutive_index++);
    //   }
    //   else{
    //     in->SetId(reverse_index--);
    //   }
    // }

    // rModelPart.Nodes().Sort();
    // rModelPart.Nodes().Unique();


    LaplacianSmoothing   MeshGeometricSmoothing(rModelPart);
    MeshGeometricSmoothing.SetEchoLevel(GetEchoLevel());
	
    bool remesh_performed=false;
	
    if( GetEchoLevel() > 0 ){
      std::cout<<" --------------                     -------------- "<<std::endl;
      std::cout<<" --------------       DOMAIN        -------------- "<<std::endl;
    }
	
    ModelPart::MeshesContainerType Meshes = rModelPart.GetMeshes();
	  	  
    for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
      {
	
	//initialize info
	mMeshingVariables[MeshId].RemeshInfo.Initialize();

	mMeshingVariables[MeshId].RemeshInfo.NumberOfElements   = rModelPart.NumberOfElements(MeshId);
	mMeshingVariables[MeshId].RemeshInfo.NumberOfNodes      = rModelPart.NumberOfNodes(MeshId);
	mMeshingVariables[MeshId].RemeshInfo.NumberOfConditions = rModelPart.NumberOfConditions(MeshId);
	  
	//set properties in all meshes
	if(rModelPart.NumberOfProperties(MeshId)!=rModelPart.NumberOfProperties())
	  rModelPart.GetMesh(MeshId).SetProperties(rModelPart.GetMesh().pProperties());

	//transfer DETERMINANT_F FOR VARIABLES SMOOTHING to nodes
	if(mMeshingVariables[MeshId].JacobiSmoothingFlag){
	  mDataTransferUtilities.TransferElementalValuesToNodes(DETERMINANT_F,rModelPart,MeshId);
	}
	    
	if(Meshes[MeshId].Is( REMESH )){


	  if(Meshes[MeshId].Is( CONSTRAINED_MESH ))
	    {
	      if(Meshes[MeshId].Is( REFINE_MESH )){
		//Constrained Delaunay Triangulation
		if( GetEchoLevel() > 0 )
		  std::cout<<" [ MESH: "<<MeshId<<" REFINE RCDT START]:"<<std::endl;
		this->GenerateRCDT(rModelPart,mMeshingVariables[MeshId],MeshId);	
		if( GetEchoLevel() > 0 )
		  std::cout<<" [ MESH: "<<MeshId<<" REFINE RCDT END]"<<std::endl;	  
	      }
	      else{ 	
		//Generate Constrained Delaunay Triangulation
		if( GetEchoLevel() > 0 )
		  std::cout<<" [ MESH: "<<MeshId<<" REMESH CDT ]:"<<std::endl;
		this->GenerateCDT(rModelPart,mMeshingVariables[MeshId],MeshId);		
		if( GetEchoLevel() > 0 )
		  std::cout<<" [ MESH: "<<MeshId<<" REMESH END ]"<<std::endl;    		    
	      }
	    }
	  else
	    {
	      if(Meshes[MeshId].Is( REFINE_MESH )){ 
		//Constrained Delaunay Triangulation
		if( GetEchoLevel() > 0 )
		  std::cout<<" [ MESH: "<<MeshId<<" REFINE RDT START]:"<<std::endl;
		this->GenerateRDT(rModelPart,mMeshingVariables[MeshId],MeshId);	
		if( GetEchoLevel() > 0 )
		  std::cout<<" [ MESH: "<<MeshId<<" REFINE RDT END]"<<std::endl;	  
	      }
	      else{
		//Generate Delaunay Triangulation
		if( GetEchoLevel() > 0 )
		  std::cout<<" [ MESH: "<<MeshId<<" REMESH DT START ]:"<<std::endl;
		this->GenerateDT(rModelPart,mMeshingVariables[MeshId],MeshId);
		if( GetEchoLevel() > 0 )
		  std::cout<<" [ MESH: "<<MeshId<<" REMESH DT END ]"<<std::endl;
	      }
	    }
		
	  if( GetEchoLevel() > 0 ){
	    std::cout<<" --------------                     -------------- "<<std::endl;
	    std::cout<<" --------------  REMESH PERFORMED   -------------- "<<std::endl;
	    std::cout<<" --------------                     -------------- "<<std::endl;
	  }

	  remesh_performed=true;
	}
	else{

	  if(mMeshingVariables[MeshId].MeshSmoothingFlag){

	    if( GetEchoLevel() > 0 )
	      std::cout<<" [ MESH: "<<MeshId<<" TRANSFER START ]:"<<std::endl;

	    //transfer only is done if the remesh option is active
	    rModelPart.GetMesh(MeshId).Set( REMESH );  
	    //and if there is a minimum of inserted or removed nodes
	    mMeshingVariables[MeshId].RemeshInfo.InsertedNodes = mMeshingVariables[MeshId].RemeshInfo.NumberOfNewNodes;

	    this->PerformTransferOnly(rModelPart,mMeshingVariables[MeshId],MeshId);

	    remesh_performed=true;

	    if( GetEchoLevel() > 0 )
	      std::cout<<" [ MESH: "<<MeshId<<" TRANSFER END ]"<<std::endl;
	  }
	  else{
	    if( GetEchoLevel() > 0 )
	      std::cout<<" [ MESH: "<<MeshId<<" NO REMESH ]"<<std::endl;
	  }

	}


	//finalize info
	mMeshingVariables[MeshId].RemeshInfo.NumberOfNewElements   = rModelPart.NumberOfElements(MeshId) - mMeshingVariables[MeshId].RemeshInfo.NumberOfElements;
	mMeshingVariables[MeshId].RemeshInfo.NumberOfNewNodes      = rModelPart.NumberOfNodes(MeshId) - mMeshingVariables[MeshId].RemeshInfo.NumberOfNodes;
	mMeshingVariables[MeshId].RemeshInfo.NumberOfNewConditions = rModelPart.NumberOfConditions(MeshId) - mMeshingVariables[MeshId].RemeshInfo.NumberOfConditions;


      }


    //Once all meshes are build, the main mesh Id=0 must be reassigned
    if(NumberOfMeshes>1){
      mModelerUtilities.BuildTotalMesh(rModelPart, GetEchoLevel());
    }
    else{
      mModelerUtilities.CleanMeshFlags(rModelPart,0);
    }
	

    //remesh_performed = false; //deactivate searches
	   
    if(remesh_performed){


      //CHANGE ELEMENTS WITH ALL NODES IN BOUNDARY  (TIP ELEMENTS)
      for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
	{
		
	  if(Meshes[MeshId].Is( REMESH )){
	    ChangeTipElementsUtilities TipElements;
	    if(mMeshingVariables[MeshId].AvoidTipElementsFlag){
	      //TipElements.SwapDiagonals(rModelPart,MeshId);
	    }
	  }
	}		
	    
	    
      //NEIGHBOURS SEARCH
      NodalNeighboursSearchProcess FindNeighbours(rModelPart);
      FindNeighbours.Execute();

      //NODAL_H SEARCH	    
      //FindNodalHProcess FindNodalH(rModelPart);
      //FindNodalH.Execute();

      //CONDITIONS MASTER_ELEMENTS and MASTER_NODES SEARCH
      BoundarySkinBuildProcess BoundarySkinProcess(rModelPart,2);
      for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
	{
	  BoundarySkinProcess.SearchConditionMasters(MeshId);
	}


      //BOUNDARY NORMALS SEARCH // SHRINKAGE FACTOR
      //ComputeBoundaryNormals BoundUtils;
      BoundaryNormalsCalculationUtilities BoundaryComputation;
      for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
	{
	  BoundaryComputation.CalculateMeshBoundaryNormals(rModelPart, MeshId, mEchoLevel);
	}


      //LAPLACIAN SMOOTHING
      for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
	{

	  //Check Mesh Info to perform smoothing:
	  mMeshingVariables[MeshId].RemeshInfo.CheckGeometricalSmoothing();
	  mMeshingVariables[MeshId].RemeshInfo.CheckMechanicalSmoothing();

	  if( mMeshingVariables[MeshId].MeshSmoothingFlag && mMeshingVariables[MeshId].RemeshInfo.GeometricalSmoothingRequired ){
	    //MeshGeometricSmoothing.ApplyMeshSmoothing(rModelPart,LaplacianSmoothing::SMOOTH_ALL,MeshId);
	  }
		  
	  if( mMeshingVariables[MeshId].JacobiSmoothingFlag && mMeshingVariables[MeshId].RemeshInfo.MechanicalSmoothingRequired ){

	    //recover DETERMINANT_F FOR VARIABLES SMOOTHING from nodes
	    if( mMeshingVariables[MeshId].RemeshInfo.CriticalElements > 0 ){
	      
	      double critical_value = 0;

	      //Smoothing performed only in critical elements (based on Plastic Strain)	      
	      mDataTransferUtilities.TransferNodalValuesToElements(DETERMINANT_F,PLASTIC_STRAIN,critical_value,rModelPart,MeshId);

	      //Smoothing performed only in critical elements (based on Set Dissipation Variable)	      
	      //critical_value = mMeshingVariables[MeshId].Refine.CriticalDissipation;
	      //mDataTransferUtilities.TransferNodalValuesToElements(DETERMINANT_F,mMeshingVariables[MeshId].Refine.GetDissipationVariable(),critical_value,rModelPart,MeshId);

	      //Smoothing performed only in critical elements (based on Plastic Energy Dissipation)
	      //critical_value = mMeshingVariables[MeshId].Refine.CriticalDissipation;
	      //mDataTransferUtilities.TransferNodalValuesToElements(DETERMINANT_F,PLASTIC_DISSIPATION,critical_value,rModelPart,MeshId);

	    }
	    else{

	      //Smoothing performed to all mesh
	      mDataTransferUtilities.TransferNodalValuesToElements(DETERMINANT_F,rModelPart,MeshId);

	    }

	  }

	}

    }

    // if(!out_buffer_active)
    //   std::cout.rdbuf(buffer);


    KRATOS_CATCH(" ")
  }



} // Namespace Kratos

