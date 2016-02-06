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
#include "custom_modelers/triangular_mesh_2D_modeler.hpp"

#include "pfem_solid_mechanics_application_variables.h"


namespace Kratos
{
  
  //*******************************************************************************************
  //*******************************************************************************************
  void TriangularMesh2DModeler::PerformTransferOnly(ModelPart& rModelPart,
						    MeshingVariables& rMeshingVariables,
						    ModelPart::IndexType MeshId)
  {

    KRATOS_TRY

    if( this->GetEchoLevel() > 0 ){
      std::cout<<" [ [ [ ] ] ]"<<std::endl;
      std::cout<<" [ Trigen PFEM Transfer Only ]"<<std::endl;
      std::cout<<" [ PREVIOUS MESH (Elements: "<<rModelPart.NumberOfElements(MeshId)<<" Nodes: "<<rModelPart.NumberOfNodes(MeshId)<<" Conditions: "<<rModelPart.NumberOfConditions(MeshId)<<") ] MESH_ID: ["<<MeshId<<"]"<<std::endl;
    }

    //*********************************************************************
    struct triangulateio in;
    struct triangulateio out;

	
    rMeshingVariables.NodalIdsSetFlag=false;

    rMeshingVariables.MeshingOptions.Set(MeshModeler::SET_DOF);
    ////////////////////////////////////////////////////////////
    SetTriangulationNodes(rModelPart,rMeshingVariables,in,out,MeshId);
    ////////////////////////////////////////////////////////////
    rMeshingVariables.MeshingOptions.Reset(MeshModeler::SET_DOF);


    //*********************************************************************
    //input mesh: ELEMENTS
    in.numberoftriangles = rModelPart.Elements(MeshId).size();
    in.trianglelist = new int [in.numberoftriangles * 3];
    in.neighborlist = new int [in.numberoftriangles * 3];

    //copying the elements in the input file	    
    ModelPart::ElementsContainerType::iterator elem_begin = rModelPart.ElementsBegin(MeshId);

    int base=0;
    for(unsigned int el = 0; el<rModelPart.Elements(MeshId).size(); el++)
      {
	(elem_begin+el)->SetId(el+1);
	Geometry<Node<3> >& geom = (elem_begin+el)->GetGeometry();

	in.trianglelist[base]   = geom[0].Id();
	in.trianglelist[base+1] = geom[1].Id();
	in.trianglelist[base+2] = geom[2].Id();
	base+=3;
      }

	
    for(unsigned int el = 0; el<rModelPart.Elements(MeshId).size(); el++)
      {
	// Geometry<Node<3> >& geom = (elem_begin+el)->GetGeometry();
	WeakPointerVector<Element >& rE = (elem_begin+el)->GetValue(NEIGHBOUR_ELEMENTS);

	for(int pn=0; pn<3; pn++){
	  if( (elem_begin+el)->Id() == rE[pn].Id() )
	    in.neighborlist[el*3+pn] = 0;
	  else
	    in.neighborlist[el*3+pn] = rE[pn].Id();
	}
	      
      }
	  
    //*********************************************************************


    rMeshingVariables.RefiningOptions.Set(MeshModeler::SELECT_ELEMENTS);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::PASS_ALPHA_SHAPE);	  

    ////////////////////////////////////////////////////////////
    BuildMeshElements(rModelPart,rMeshingVariables,in,MeshId);
    ////////////////////////////////////////////////////////////

    rMeshingVariables.RefiningOptions.Reset(MeshModeler::PASS_ALPHA_SHAPE);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::SELECT_ELEMENTS);

    //*********************************************************************


    //*********************************************************************
	  
    ////////////////////////////////////////////////////////////
    BuildMeshBoundary(rModelPart,rMeshingVariables,MeshId);
    ////////////////////////////////////////////////////////////

    //*********************************************************************

    //free memory
    DeletePointsList(in);
    //delete [] in.trianglelist;

    if( this->GetEchoLevel() > 0 ){
      std::cout<<" [ NEW MESH (Elements: "<<rModelPart.Elements(MeshId).size()<<" Nodes: "<<rModelPart.Nodes(MeshId).size()<<" Conditions: "<<rModelPart.Conditions(MeshId).size()<<" ] "<<std::endl;
      std::cout<<" [ Finished Remeshing ] "<<std::endl;
      std::cout<<" [ [ [ ] ] ]"<<std::endl;
    }

    KRATOS_CATCH( "" )
  }

  
  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::GenerateDT(ModelPart& rModelPart,
					   MeshingVariables& rMeshingVariables,
					   ModelPart::IndexType MeshId)
  {

    KRATOS_TRY
 
    if( this->GetEchoLevel() > 0 ){
      std::cout<<" [ [ [ ] ] ]"<<std::endl;
      std::cout<<" [ Trigen PFEM DT Mesher ]"<<std::endl;
      std::cout<<" [ PREVIOUS MESH (Elements: "<<rModelPart.NumberOfElements(MeshId)<<" Nodes: "<<rModelPart.NumberOfNodes(MeshId)<<" Conditions: "<<rModelPart.NumberOfConditions(MeshId)<<") ] MESH_ID: ["<<MeshId<<"]"<<std::endl;
    }

    //int step_data_size = rModelPart.GetNodalSolutionStepDataSize();

    if(rMeshingVariables.RefineFlag){

      rMeshingVariables.RefiningOptions.Set(MeshModeler::REMOVE_NODES);
      rMeshingVariables.RefiningOptions.Set(MeshModeler::CRITERION_DISTANCE);
      ////////////////////////////////////////////////////////////
      RemoveCloseNodes (rModelPart,rMeshingVariables,MeshId);
      ////////////////////////////////////////////////////////////
      rMeshingVariables.RefiningOptions.Reset(MeshModeler::CRITERION_DISTANCE);
      rMeshingVariables.RefiningOptions.Reset(MeshModeler::REMOVE_NODES);
    }


    ////////////////////////////////////////////////////////////
    //Creating the containers for the input and output
    struct triangulateio in;
    struct triangulateio out;

    rMeshingVariables.NodalIdsSetFlag=false;

    rMeshingVariables.MeshingOptions.Set(MeshModeler::SET_DOF);
    ////////////////////////////////////////////////////////////
    SetTriangulationNodes(rModelPart,rMeshingVariables,in,out,MeshId);
    ////////////////////////////////////////////////////////////
    rMeshingVariables.MeshingOptions.Reset(MeshModeler::SET_DOF);


    ////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////
    boost::timer auxiliary;

    rMeshingVariables.MeshingOptions.Set(MeshModeler::NEIGHBOURS_SEARCH);
    GenerateTriangulation(rMeshingVariables.MeshingOptions,rMeshingVariables.RefiningOptions,in,out);
    rMeshingVariables.MeshingOptions.Reset(MeshModeler::NEIGHBOURS_SEARCH);

    //print out the mesh generation time
    if( this->GetEchoLevel() > 0 )
      std::cout<<" [ MESH GENERATION (TIME = "<<auxiliary.elapsed()<<") ] "<<std::endl;
    ////////////////////////////////////////////////////////////

    if(rMeshingVariables.RefineFlag){
	  
      ////////////////////////////////////////////////////////////
      //Select Elements to be preserved after passing Alpha-Shape
      rMeshingVariables.RefiningOptions.Set(MeshModeler::PASS_ALPHA_SHAPE);
      rMeshingVariables.RefiningOptions.Set(MeshModeler::SELECT_ELEMENTS);
	  	  
      SelectMeshElements(rModelPart.Nodes(MeshId),rMeshingVariables,out); //passing alpha shape and returning the Elements preserved  

      rMeshingVariables.RefiningOptions.Reset(MeshModeler::SELECT_ELEMENTS);
      rMeshingVariables.RefiningOptions.Reset(MeshModeler::PASS_ALPHA_SHAPE);
      ////////////////////////////////////////////////////////////

      //free the memory used in the first step, preserve out
      DeletePointsList(in);
      DeleteTrianglesList(in);
 	  
      ////////////////////////////////////////////////////////////
      //PERFORM ADAPTIVE REMESHING:
      //i.e. insert and remove nodes based upon mesh quality and prescribed sizes	
      struct triangulateio mid_out;

      rMeshingVariables.NodalIdsSetFlag=true; //second set must be true	  
      ////////////////////////////////////////////////////////////
      SetTriangulationNodes (rModelPart,rMeshingVariables,in,mid_out,MeshId);
      ////////////////////////////////////////////////////////////
	  

      rMeshingVariables.RefiningOptions.Set(MeshModeler::REFINE_ELEMENTS);
      ////////////////////////////////////////////////////////////
      RefineElements (rModelPart,rMeshingVariables,in,out,MeshId);
      ////////////////////////////////////////////////////////////
      rMeshingVariables.RefiningOptions.Reset(MeshModeler::REFINE_ELEMENTS);


      ////////////////////////////////////////////////////////////

      //free the memory used in the first step, free out
      ClearTrianglesList(out);

      ////////////////////////////////////////////////////////////
      rMeshingVariables.MeshingOptions.Set(MeshModeler::REFINE_MESH);
      rMeshingVariables.RefiningOptions.Set(MeshModeler::REFINE_ADD_NODES); //optional

      GenerateTriangulation(rMeshingVariables.MeshingOptions,rMeshingVariables.RefiningOptions,in,out);

      rMeshingVariables.MeshingOptions.Reset(MeshModeler::REFINE_MESH);
      rMeshingVariables.RefiningOptions.Reset(MeshModeler::REFINE_ADD_NODES); //optional
      ////////////////////////////////////////////////////////////

      //Building the entities for new nodes:
      GenerateNewParticles(rModelPart,rMeshingVariables,in,out,MeshId);

      ////////////////////////////////////////////////////////////
    }
    else{

      rMeshingVariables.RefiningOptions.Set(MeshModeler::SELECT_ELEMENTS);
      rMeshingVariables.RefiningOptions.Set(MeshModeler::PASS_ALPHA_SHAPE);	  
    }


    //*********************************************************************
    //input mesh: NODES //to not perturb the node position because of the meshing
    out.numberofpoints = rModelPart.Nodes(MeshId).size();

    //writing the points coordinates in a vector and reordening the Id's
    ModelPart::NodesContainerType::iterator nodes_begin = rModelPart.NodesBegin(MeshId);

    //std::cout<<"  [ SET NODES: ";
    int base=0;
    for(unsigned int i = 0; i<rModelPart.Nodes(MeshId).size(); i++)
      {
	out.pointlist[base]   = (nodes_begin + i)->X();
	out.pointlist[base+1] = (nodes_begin + i)->Y();
	    
	base+=2;
      }
    //*********************************************************************

    ////////////////////////////////////////////////////////////
    BuildMeshElements(rModelPart,rMeshingVariables,out,MeshId);
    ////////////////////////////////////////////////////////////

    rMeshingVariables.RefiningOptions.Reset(MeshModeler::PASS_ALPHA_SHAPE);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::SELECT_ELEMENTS);

    //*********************************************************************

    ////////////////////////////////////////////////////////////
    BuildMeshBoundary(rModelPart,rMeshingVariables, MeshId);
    ////////////////////////////////////////////////////////////

    //*********************************************************************

    //sort elements
    //rModelPart.Elements().Sort();

    //free memory
    DeletePointsList(in);
    delete [] in.trianglelist;
    DeleteTrianglesList(out);

    if( this->GetEchoLevel() > 0 ){
      std::cout<<" [ NEW MESH (Elements: "<<rModelPart.Elements(MeshId).size()<<" Nodes: "<<rModelPart.Nodes(MeshId).size()<<" Conditions: "<<rModelPart.Conditions(MeshId).size()<<" ] "<<std::endl;
      std::cout<<" [ Finished Remeshing ] "<<std::endl;
      std::cout<<" [ [ [ ] ] ]"<<std::endl;
    }

    KRATOS_CATCH( "" )
  }



  //*******************************************************************************************
  //*******************************************************************************************
  void TriangularMesh2DModeler::GenerateCDT(ModelPart& rModelPart,
					    MeshingVariables& rMeshingVariables,
					    ModelPart::IndexType MeshId)
  {
    KRATOS_TRY

    if( this->GetEchoLevel() > 0 ){
      std::cout<<" [ [ [ ] ] ]"<<std::endl;
      std::cout<<" [ Trigen PFEM Conforming Constrained Delaunay Mesher ]"<<std::endl;
      std::cout<<" [ PREVIOUS MESH (Elements: "<<rModelPart.NumberOfElements(MeshId)<<" Nodes: "<<rModelPart.NumberOfNodes(MeshId)<<" Conditions: "<<rModelPart.NumberOfConditions(MeshId)<<") ] MESH_ID: ["<<MeshId<<"]"<<std::endl;
    }

    //*********************************************************************
    struct triangulateio in;
    struct triangulateio out;

    rMeshingVariables.NodalIdsSetFlag=false;

    ////////////////////////////////////////////////////////////
    SetTriangulationNodes(rModelPart,rMeshingVariables,in,out,MeshId);
    ////////////////////////////////////////////////////////////

    //*********************************************************************
    //input mesh: ELEMENTS
    in.numberoftriangles = rModelPart.Elements(MeshId).size();
    in.trianglelist = new int [in.numberoftriangles * 3];

    //copying the elements in the input file	    


    ModelPart::ElementsContainerType::iterator elem_begin = rModelPart.ElementsBegin(MeshId);

    int base=0;
    for(unsigned int el = 0; el<rModelPart.Elements(MeshId).size(); el++)
      {
	Geometry<Node<3> >& geom = (elem_begin+el)->GetGeometry();

	in.trianglelist[base]   = geom[0].Id();
	in.trianglelist[base+1] = geom[1].Id();
	in.trianglelist[base+2] = geom[2].Id();
	base+=3;
      }
	
    ////////////////////////////////////////////////////////////

    //read and regenerate the existing mesh ... to generate the boundaries
    struct triangulateio mid;
    ClearTrianglesList(mid);

    rMeshingVariables.MeshingOptions.Set(MeshModeler::BOUNDARIES_SEARCH);
    GenerateTriangulation(rMeshingVariables.MeshingOptions,rMeshingVariables.RefiningOptions,in, mid);
    rMeshingVariables.MeshingOptions.Reset(MeshModeler::BOUNDARIES_SEARCH);

    // KRATOS_WATCH( in.numberofsegments )
    // KRATOS_WATCH( in.numberofpoints )
    // KRATOS_WATCH( in.numberoftriangles )
    // KRATOS_WATCH( in.numberofholes )

    //free the memory used in the first step
    DeletePointsList(in);
    DeleteTrianglesList(in);

    //uses the boundary list generated at the previous step to generate the "skin"
    mid.numberoftriangles = 0;
    delete [] mid.trianglelist; //delete only triangles

    rMeshingVariables.MeshingOptions.Set(MeshModeler::NEIGHBOURS_SEARCH);
    rMeshingVariables.MeshingOptions.Set(MeshModeler::CONSTRAINED_MESH);

    int fail = GenerateTriangulation(rMeshingVariables.MeshingOptions,rMeshingVariables.RefiningOptions,mid, out);

    if(fail){
      rMeshingVariables.MeshingOptions.Reset(MeshModeler::CONSTRAINED_MESH);
      fail = GenerateTriangulation(rMeshingVariables.MeshingOptions,rMeshingVariables.RefiningOptions,mid, out);
      rMeshingVariables.MeshingOptions.Set(MeshModeler::CONSTRAINED_MESH);
    }
    
    rMeshingVariables.MeshingOptions.Reset(MeshModeler::CONSTRAINED_MESH);
    rMeshingVariables.MeshingOptions.Reset(MeshModeler::NEIGHBOURS_SEARCH);

    // KRATOS_WATCH( out.numberofsegments )
    // KRATOS_WATCH( out.numberofpoints )
    // KRATOS_WATCH( out.numberoftriangles )
    // KRATOS_WATCH( out.numberofholes )
	
    ////////////////////////////////////////////////////////////
    BuildMeshElements(rModelPart,rMeshingVariables,out,MeshId);	  
    ////////////////////////////////////////////////////////////
	
    //*********************************************************************

    ////////////////////////////////////////////////////////////
    BuildMeshBoundary(rModelPart,rMeshingVariables,MeshId);
    ////////////////////////////////////////////////////////////

    //*********************************************************************

    //free the memory used in the intermediate step
    //DeleteTrianglesList(mid);
    //DeletePointsList(mid);
    ClearTrianglesList(mid);
    delete [] mid.trianglelist;

    //free the rest of the memory
    //DeleteTrianglesList(out);
    ClearTrianglesList(out);
	

    if( this->GetEchoLevel() > 0 ){
      std::cout<<" [ NEW MESH (Elements: "<<rModelPart.Elements(MeshId).size()<<" Nodes: "<<rModelPart.Nodes(MeshId).size()<<" Conditions: "<<rModelPart.Conditions(MeshId).size()<<" ] "<<std::endl;
      std::cout<<" [ Finished Remeshing ] "<<std::endl;
      std::cout<<" [ [ [ ] ] ]"<<std::endl;
    }

    KRATOS_CATCH( "" )
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::GenerateRDT(ModelPart& rModelPart,
					    MeshingVariables& rMeshingVariables,
					    ModelPart::IndexType MeshId)
  {
    KRATOS_TRY

    if( this->GetEchoLevel() > 0 ){
      std::cout<<" [ [ [ ] ] ]"<<std::endl;
      std::cout<<" [ Trigen PFEM DT Refine Mesher ]"<<std::endl;
      std::cout<<" [ PREVIOUS MESH (Elements: "<<rModelPart.NumberOfElements(MeshId)<<" Nodes: "<<rModelPart.NumberOfNodes(MeshId)<<" Conditions: "<<rModelPart.NumberOfConditions(MeshId)<<") ] MESH_ID: ["<<MeshId<<"]"<<std::endl;
    }
		    
    //*********************************************************************

    //Needed in RefineElements
    ////////////////////////////////////////////////////////////
    SetDissipativeElements (rModelPart,rMeshingVariables,MeshId);  
    ////////////////////////////////////////////////////////////


    //*********************************************************************

    rMeshingVariables.RefiningOptions.Set(MeshModeler::REFINE_INSERT_NODES);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::CRITERION_ENERGY);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::REFINE_BOUNDARY);

    ////////////////////////////////////////////////////////////
    RefineBoundary (rModelPart,rMeshingVariables,MeshId);
    ////////////////////////////////////////////////////////////

    rMeshingVariables.RefiningOptions.Reset(MeshModeler::REFINE_BOUNDARY);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::CRITERION_ENERGY);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::REFINE_INSERT_NODES);

    //********************************************************************

    //we need to redefine tool_tip boundaries after refining them !!

    rMeshingVariables.RefiningOptions.Set(MeshModeler::REMOVE_NODES);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::CRITERION_ERROR);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::CRITERION_DISTANCE);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::REMOVE_ON_BOUNDARY);

    ////////////////////////////////////////////////////////////
    RemoveCloseNodes (rModelPart,rMeshingVariables,MeshId);
    ////////////////////////////////////////////////////////////

    rMeshingVariables.RefiningOptions.Reset(MeshModeler::REMOVE_NODES);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::CRITERION_ERROR);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::CRITERION_DISTANCE);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::REMOVE_ON_BOUNDARY);


    ////////////////////////////////////////////////////////////
    //Creating the containers for the input and output
    struct triangulateio in;
    struct triangulateio out;

    rMeshingVariables.NodalIdsSetFlag=false;

    rMeshingVariables.MeshingOptions.Set(MeshModeler::SET_DOF);
    ////////////////////////////////////////////////////////////
    SetTriangulationNodes(rModelPart,rMeshingVariables,in,out,MeshId);
    ////////////////////////////////////////////////////////////
    rMeshingVariables.MeshingOptions.Reset(MeshModeler::SET_DOF);


    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    boost::timer auxiliary;

    rMeshingVariables.MeshingOptions.Set(MeshModeler::NEIGHBOURS_SEARCH);
    int fail = GenerateTriangulation(rMeshingVariables.MeshingOptions,rMeshingVariables.RefiningOptions,in,out);

    if(fail){
      std::cout<<" Mesher Failed first RCDT "<<std::endl;
    }

    rMeshingVariables.MeshingOptions.Reset(MeshModeler::NEIGHBOURS_SEARCH);

    if(in.numberofpoints!=out.numberofpoints){
      std::cout<<" [ MESH GENERATION FAILED: point insertion (initial = "<<in.numberofpoints<<" final = "<<out.numberofpoints<<") ] "<<std::endl;
    }

    //print out the mesh generation time
    if( this->GetEchoLevel() > 0 )
      std::cout<<" [ MESH GENERATION (TIME = "<<auxiliary.elapsed()<<") ] "<<std::endl;
    ////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////
    //Select Elements to be preserved after passing Alpha-Shape
    rMeshingVariables.RefiningOptions.Set(MeshModeler::PASS_ALPHA_SHAPE);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::SELECT_ELEMENTS);
		
    SelectMeshElements(rModelPart.Nodes(MeshId),rMeshingVariables,out); //passing alpha shape and returning the Elements preserved

    rMeshingVariables.RefiningOptions.Reset(MeshModeler::SELECT_ELEMENTS);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::PASS_ALPHA_SHAPE);
    ////////////////////////////////////////////////////////////

    //free the memory used in the first step, preserve out
    DeletePointsList(in);
    DeleteTrianglesList(in);

    ////////////////////////////////////////////////////////////
    //PERFORM ADAPTIVE REMESHING:
    //1.- Select Triangles to Refine via changing the nodal_h
    struct triangulateio mid_out;

    rMeshingVariables.NodalIdsSetFlag=true;  //second set must be true
    ////////////////////////////////////////////////////////////
    SetTriangulationNodes (rModelPart,rMeshingVariables,in,mid_out,MeshId);
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    rMeshingVariables.RefiningOptions.Set(MeshModeler::REFINE_ELEMENTS);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::CRITERION_ENERGY);
    //rMeshingVariables.RefiningOptions.Set(MeshModeler::REFINE_BOUNDARY);
		
    RefineElements (rModelPart,rMeshingVariables,in,out,MeshId);

    //rMeshingVariables.RefiningOptions.Reset(MeshModeler::REFINE_BOUNDARY);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::REFINE_ELEMENTS);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::CRITERION_ENERGY);
    ////////////////////////////////////////////////////////////


    //free the memory used in the first step, free out
    ClearTrianglesList(out);

    ////////////////////////////////////////////////////////////
    rMeshingVariables.MeshingOptions.Set(MeshModeler::REFINE_MESH);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::REFINE_ADD_NODES); //optional

    fail = GenerateTriangulation(rMeshingVariables.MeshingOptions,rMeshingVariables.RefiningOptions,in,out);

    if(fail){
      std::cout<<" Mesher Failed second RCDT "<<std::endl;
    }

    rMeshingVariables.MeshingOptions.Reset(MeshModeler::REFINE_MESH);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::REFINE_ADD_NODES); //optional
    ////////////////////////////////////////////////////////////

    //Building the entities for new nodes:
    GenerateNewParticles(rModelPart,rMeshingVariables,in,out,MeshId);

    ////////////////////////////////////////////////////////////
    BuildMeshElements(rModelPart,rMeshingVariables,out,MeshId);
    ////////////////////////////////////////////////////////////

    //std::cout<<" Check After BuildMeshElements "<<std::endl;
    //mModelerUtilities.CheckParticles (rModelPart,MeshId);

    //*********************************************************************

    ////////////////////////////////////////////////////////////
    BuildMeshBoundary(rModelPart,rMeshingVariables,MeshId);
    ////////////////////////////////////////////////////////////

    //*********************************************************************
    //sort elements
    //rModelPart.Elements().Sort();

    //free memory
    DeletePointsList(in);
    delete [] in.trianglelist;
    DeleteTrianglesList(out);

    if( this->GetEchoLevel() > 0 ){
      std::cout<<" [ NEW MESH (Elements: "<<rModelPart.Elements(MeshId).size()<<" Nodes: "<<rModelPart.Nodes(MeshId).size()<<" Conditions: "<<rModelPart.Conditions(MeshId).size()<<" ] "<<std::endl;
      std::cout<<" [ Finished Remeshing ] "<<std::endl;
      std::cout<<" [ [ [ ] ] ]"<<std::endl;
    }

    KRATOS_CATCH( "" )
  }


  //*******************************************************************************************
  //*******************************************************************************************
  
  void TriangularMesh2DModeler::GenerateRCDT(ModelPart& rModelPart,
					     MeshingVariables& rMeshingVariables,
					     ModelPart::IndexType MeshId)
  {
    KRATOS_TRY

    if( this->GetEchoLevel() > 0 ){
      std::cout<<" [ [ [ ] ] ]"<<std::endl;
      std::cout<<" [ Trigen PFEM CDT Refine Mesher ]"<<std::endl;
      std::cout<<" [ PREVIOUS MESH (Elements: "<<rModelPart.NumberOfElements(MeshId)<<" Nodes: "<<rModelPart.NumberOfNodes(MeshId)<<" Conditions: "<<rModelPart.NumberOfConditions(MeshId)<<") ] MESH_ID: ["<<MeshId<<"]"<<std::endl;
    }

    //*********************************************************************

    //Needed in RefineElements
    ////////////////////////////////////////////////////////////
    SetDissipativeElements (rModelPart,rMeshingVariables,MeshId);  
    ////////////////////////////////////////////////////////////

    //check initial mesh:
    // for(ModelPart::ElementsContainerType::const_iterator iii = rModelPart.ElementsBegin(MeshId);
    // 	iii != rModelPart.ElementsEnd(MeshId); iii++)
    //   {
    // 	std::cout<<" START TRIANGLE ["<<iii->Id()<<"] ["<<iii->GetGeometry()[0].Id()<<", "<<iii->GetGeometry()[1].Id()<<", "<<iii->GetGeometry()[2].Id()<<"] "<<std::endl;
    // 	std::cout<<" START PRESSURE ["<<iii->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)<<", "<<iii->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)<<",  "<<iii->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)<<"] "<<std::endl;
    //   }

    rMeshingVariables.MeshingOptions.Set(MeshModeler::CONSTRAINED_MESH);

    //*********************************************************************

    rMeshingVariables.RefiningOptions.Set(MeshModeler::REFINE_INSERT_NODES);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::CRITERION_ENERGY);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::REFINE_BOUNDARY);

    ////////////////////////////////////////////////////////////
    RefineBoundary (rModelPart,rMeshingVariables,MeshId);

    ////////////////////////////////////////////////////////////

    rMeshingVariables.RefiningOptions.Reset(MeshModeler::REFINE_BOUNDARY);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::CRITERION_ENERGY);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::REFINE_INSERT_NODES);

    //*********************************************************************

    //we need to redefine tool_tip boundaries after refining them !!

    rMeshingVariables.RefiningOptions.Set(MeshModeler::REMOVE_NODES);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::CRITERION_ERROR);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::CRITERION_DISTANCE);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::REMOVE_ON_BOUNDARY);

    ////////////////////////////////////////////////////////////
    RemoveCloseNodes (rModelPart,rMeshingVariables,MeshId);
    ////////////////////////////////////////////////////////////

    rMeshingVariables.RefiningOptions.Reset(MeshModeler::REMOVE_NODES);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::CRITERION_ERROR);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::CRITERION_DISTANCE);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::REMOVE_ON_BOUNDARY);


    ////////////////////////////////////////////////////////////
    //Creating the containers for the input and output
    struct triangulateio in;
    struct triangulateio out;

    rMeshingVariables.NodalIdsSetFlag=false;

    rMeshingVariables.MeshingOptions.Set(MeshModeler::SET_DOF);
    ////////////////////////////////////////////////////////////
    SetTriangulationNodes(rModelPart,rMeshingVariables,in,out,MeshId); 
    ////////////////////////////////////////////////////////////
    rMeshingVariables.MeshingOptions.Reset(MeshModeler::SET_DOF);


    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    boost::timer auxiliary;

    rMeshingVariables.MeshingOptions.Set(MeshModeler::NEIGHBOURS_SEARCH);

    int fail = GenerateTriangulation(rMeshingVariables.MeshingOptions,rMeshingVariables.RefiningOptions,in,out);

    if(fail){
      rMeshingVariables.MeshingOptions.Reset(MeshModeler::CONSTRAINED_MESH);
      fail = GenerateTriangulation(rMeshingVariables.MeshingOptions,rMeshingVariables.RefiningOptions,in, out);
      rMeshingVariables.MeshingOptions.Set(MeshModeler::CONSTRAINED_MESH);
    }

    rMeshingVariables.MeshingOptions.Reset(MeshModeler::NEIGHBOURS_SEARCH);

    if(in.numberofpoints!=out.numberofpoints){
      std::cout<<" [ MESH GENERATION FAILED: point insertion (initial = "<<in.numberofpoints<<" final = "<<out.numberofpoints<<") ] "<<std::endl;
    }

    if( rMeshingVariables.MeshingOptions.Is(MeshModeler::CONSTRAINED_MESH) )
      RecoverBoundaryPosition(rModelPart,rMeshingVariables,in,out,MeshId);

    //print out the mesh generation time
    if( this->GetEchoLevel() > 0 )
      std::cout<<" [ MESH GENERATION (TIME = "<<auxiliary.elapsed()<<") ] "<<std::endl;
    ////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////
    //Select Elements to be preserved after passing Alpha-Shape
    rMeshingVariables.RefiningOptions.Set(MeshModeler::PASS_ALPHA_SHAPE);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::SELECT_ELEMENTS);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::ENGAGED_NODES);

    SelectMeshElements(rModelPart.Nodes(MeshId),rMeshingVariables,out); //passing alpha shape and returning the Elements preserved

    rMeshingVariables.RefiningOptions.Reset(MeshModeler::ENGAGED_NODES);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::SELECT_ELEMENTS);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::PASS_ALPHA_SHAPE);
    ////////////////////////////////////////////////////////////

    //free the memory used in the first step, preserve out
    DeletePointsList(in);
    DeleteTrianglesList(in);

    ////////////////////////////////////////////////////////////
    //PERFORM ADAPTIVE REMESHING:
    //1.- Select Triangles to Refine via changing the nodal_h
    struct triangulateio mid_out;

    rMeshingVariables.NodalIdsSetFlag=true; //second set must be true
    ////////////////////////////////////////////////////////////
    SetTriangulationNodes (rModelPart,rMeshingVariables,in,mid_out,MeshId);
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    rMeshingVariables.RefiningOptions.Set(MeshModeler::REFINE_ELEMENTS);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::CRITERION_ENERGY);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::REFINE_BOUNDARY);

    RefineElements (rModelPart,rMeshingVariables,in,out,MeshId);

    rMeshingVariables.RefiningOptions.Reset(MeshModeler::REFINE_BOUNDARY);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::REFINE_ELEMENTS);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::CRITERION_ENERGY);
    ////////////////////////////////////////////////////////////


    //free the memory used in the first step, free out
    ClearTrianglesList(out);

    ////////////////////////////////////////////////////////////
    //to generate it constrained it is necessary to change the strategy:
    //a. pass a set of nodes to triangulate ok
    //b. pass a set of segments == conditions to apply the constraint ok
    //c. pass a set of holes if domains are not totally convex. ok
	
    rMeshingVariables.MeshingOptions.Set(MeshModeler::NEIGHBOURS_SEARCH);
    rMeshingVariables.MeshingOptions.Set(MeshModeler::REFINE_MESH);
    rMeshingVariables.RefiningOptions.Set(MeshModeler::REFINE_ADD_NODES); //optional

    fail = GenerateTriangulation(rMeshingVariables.MeshingOptions,rMeshingVariables.RefiningOptions,in,out);

    if(fail){
      rMeshingVariables.MeshingOptions.Reset(MeshModeler::CONSTRAINED_MESH);
      fail = GenerateTriangulation(rMeshingVariables.MeshingOptions,rMeshingVariables.RefiningOptions,in, out);
      rMeshingVariables.MeshingOptions.Set(MeshModeler::CONSTRAINED_MESH);
    }

    rMeshingVariables.MeshingOptions.Reset(MeshModeler::REFINE_MESH);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::REFINE_ADD_NODES); //optional
    rMeshingVariables.MeshingOptions.Reset(MeshModeler::NEIGHBOURS_SEARCH);


    if( rMeshingVariables.MeshingOptions.Is(MeshModeler::CONSTRAINED_MESH) )
      RecoverBoundaryPosition(rModelPart,rMeshingVariables,in,out,MeshId);


    //check if something changes:
    // for(ModelPart::ElementsContainerType::const_iterator iii = rModelPart.ElementsBegin(MeshId);
    // 	iii != rModelPart.ElementsEnd(MeshId); iii++)
    //   {
    // 	std::cout<<" SET TRIANGLE ["<<iii->Id()<<"] ["<<rMeshingVariables.NodalPreIds[iii->GetGeometry()[0].Id()]<<", "<<rMeshingVariables.NodalPreIds[iii->GetGeometry()[1].Id()]<<", "<<rMeshingVariables.NodalPreIds[iii->GetGeometry()[2].Id()]<<"] "<<std::endl;
    // 	std::cout<<" SET PRESSURE ["<<iii->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)<<", "<<iii->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)<<",  "<<iii->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)<<"] "<<std::endl;
    //   }

    ////////////////////////////////////////////////////////////

    //Building the entities for new nodes:
    GenerateNewParticles(rModelPart,rMeshingVariables,in,out,MeshId);

    ////////////////////////////////////////////////////////////
    BuildMeshElements(rModelPart,rMeshingVariables,out,MeshId);
    ////////////////////////////////////////////////////////////

    //std::cout<<" Check After BuildMeshElements "<<std::endl;
    //mModelerUtilities.CheckParticles (rModelPart,MeshId);

    //*********************************************************************

    ////////////////////////////////////////////////////////////
    BuildMeshBoundary(rModelPart,rMeshingVariables,MeshId);
    ////////////////////////////////////////////////////////////

    //*********************************************************************

    rMeshingVariables.MeshingOptions.Reset(MeshModeler::CONSTRAINED_MESH);

    //sort elements
    //rModelPart.Elements().Sort();

    //free memory
    DeletePointsList(in);
    delete [] in.trianglelist;
    DeleteTrianglesList(out);

    if( this->GetEchoLevel() > 0 ){
      std::cout<<" [ NEW MESH (Elements: "<<rModelPart.Elements(MeshId).size()<<" Nodes: "<<rModelPart.Nodes(MeshId).size()<<" Conditions: "<<rModelPart.Conditions(MeshId).size()<<" ] "<<std::endl;
      std::cout<<" [ Finished Remeshing ] "<<std::endl;
      std::cout<<" [ [ [ ] ] ]"<<std::endl;
    }

    KRATOS_CATCH( "" )
  }
    


  //*******************************************************************************************
  //METHODS CALLED BEFORE TESSELLATION
  //*******************************************************************************************

  void TriangularMesh2DModeler::SetTriangulationNodes(ModelPart& rModelPart,
						      MeshingVariables& rMeshingVariables,
						      struct triangulateio& in,
						      struct triangulateio& out,
						      ModelPart::IndexType MeshId)
  {
    KRATOS_TRY

    ClearTrianglesList(in);
    ClearTrianglesList(out);

    //*********************************************************************
    //input mesh: NODES
    in.numberofpoints = rModelPart.Nodes(MeshId).size();
    in.pointlist      = new REAL[in.numberofpoints * 2];

    if(!rMeshingVariables.NodalIdsSetFlag){
      rMeshingVariables.NodalPreIds.resize(in.numberofpoints+1);
      std::fill( rMeshingVariables.NodalPreIds.begin(), rMeshingVariables.NodalPreIds.end(), 0 );

      rMeshingVariables.NodalNewIds.resize(rModelPart.Nodes().size()+1);
      std::fill( rMeshingVariables.NodalNewIds.begin(), rMeshingVariables.NodalNewIds.end(), 0 );
    }
	 
    //writing the points coordinates in a vector and reordening the Id's
    ModelPart::NodesContainerType::iterator nodes_begin = rModelPart.NodesBegin(MeshId);

    //std::cout<<"  [ SET NODES: ";
    int base = 0;
    int direct = 1;
    for(unsigned int i = 0; i<rModelPart.Nodes(MeshId).size(); i++)
      {
	//if( (nodes_begin + i)->IsNot(STRUCTURE) ){
	//from now on it is consecutive
	if(!rMeshingVariables.NodalIdsSetFlag){
	  //std::cout<<" "<<(nodes_begin + i)->Id()<<"(";
	  rMeshingVariables.NodalNewIds[(nodes_begin + i)->Id()] = direct;
	  rMeshingVariables.NodalPreIds[direct]=(nodes_begin + i)->Id();
	  (nodes_begin + i)->SetId(direct);
	}
	   
	if(rMeshingVariables.MeshingOptions.Is(MeshModeler::CONSTRAINED_MESH)){

	  if( (nodes_begin + i)->Is(BOUNDARY) ){
	       
	    array_1d<double, 3>&  Normal=(nodes_begin + i)->FastGetSolutionStepValue(NORMAL); //BOUNDARY_NORMAL must be set as nodal variable
	    double Shrink = (nodes_begin + i)->FastGetSolutionStepValue(SHRINK_FACTOR);   //SHRINK_FACTOR   must be set as nodal variable
	       
	    //if normal not normalized
	    //Shrink /=norm_2(Normal);
	    
	    array_1d<double, 3> Offset;

	    //std::cout<<" Shrink "<<Shrink<<" offset "<<rMeshingVariables.OffsetFactor<<std::endl;
	    Normal /= norm_2(Normal);
	    Offset[0] = ( (-1) * Normal[0] * Shrink * rMeshingVariables.OffsetFactor * 0.25 );
	    Offset[1] = ( (-1) * Normal[1] * Shrink * rMeshingVariables.OffsetFactor * 0.25 );

	    //std::cout<<" off[0] "<<Offset[0]<<" off[1] "<<Offset[1]<<std::endl;

	    in.pointlist[base]   = (nodes_begin + i)->X() + Offset[0];
	    in.pointlist[base+1] = (nodes_begin + i)->Y() + Offset[1];
	  }
	  else{
	    in.pointlist[base]   = (nodes_begin + i)->X();
	    in.pointlist[base+1] = (nodes_begin + i)->Y();
	  }

	}
	else{
	  in.pointlist[base]   = (nodes_begin + i)->X();
	  in.pointlist[base+1] = (nodes_begin + i)->Y();
	}
	   
	//std::cout<<(nodes_begin + i)->X()<<", "<<(nodes_begin + i)->Y()<<") ";
	    
	if(rMeshingVariables.MeshingOptions.Is(MeshModeler::SET_DOF))
	  {
	    Node<3>::DofsContainerType& node_dofs = (nodes_begin + i)->GetDofs();
	    for(Node<3>::DofsContainerType::iterator iii = node_dofs.begin(); iii != node_dofs.end(); iii++)
	      {
		iii->SetId(direct);
	      }
	  }
	    
	base+=2;
	direct+=1;
	// }
	//std::cout<<" node: (local:"<<(nodes_begin + i)->Id()<<", global: "<<rMeshingVariables.NodalPreIds[(nodes_begin + i)->Id()]<<") "<<std::endl;
      }


    // std::cout<<"), "<<std::endl;
    //std::cout<<"    SET NODES ]; "<<std::endl;

    if(!rMeshingVariables.NodalIdsSetFlag){
      rMeshingVariables.NodalIdsSetFlag=true;
    }
    //*********************************************************************

    //SetFaces (segments)

    //PART 1: node list

    if(rMeshingVariables.MeshingOptions.Is(MeshModeler::CONSTRAINED_MESH)){

      //PART 2: faced list (we can have holes in facets != area holes)
      in.numberofsegments           = rModelPart.NumberOfConditions(MeshId);
      in.segmentmarkerlist          = new int[in.numberofsegments];
      in.segmentlist                = new int[in.numberofsegments*2];


      ModelPart::ConditionsContainerType::iterator conditions_begin = rModelPart.ConditionsBegin(MeshId);


      base = 0;
      for(unsigned int i = 0; i<rModelPart.Conditions(MeshId).size(); i++)
	{
	  if( (conditions_begin + i)->Is(TO_ERASE) )
	    std::cout<<" ERROR: condition to erase present "<<std::endl;

	  Geometry< Node<3> >& rGeometry = (conditions_begin + i)->GetGeometry();
	  in.segmentlist[base]   = rGeometry[0].Id();
	  in.segmentlist[base+1] = rGeometry[1].Id();
	      
	  base+=2;
	}  

      //PART 3: (area) hole list

      //holes
      in.numberofholes              = 0;
      in.holelist                   = (REAL*) NULL;

      //PART 4: region attributes list
      in.numberofregions            = 1;
      in.regionlist                 = new REAL[in.numberofregions * 4];

      //regions
      double inside_factor = 2;
      Geometry< Node<3> >& rGeometry = (conditions_begin)->GetGeometry();
      array_1d<double, 3>&  Normal   = rGeometry[0].FastGetSolutionStepValue(NORMAL); 
      double NormNormal = norm_2(Normal);
      if( NormNormal != 0)
	Normal /= NormNormal;

      //inside point of the region:
      in.regionlist[0] = rGeometry[0][0]-((-1)*Normal[0]*rMeshingVariables.OffsetFactor*inside_factor);
      in.regionlist[1] = rGeometry[0][1]-((-1)*Normal[1]*rMeshingVariables.OffsetFactor*inside_factor);
	    
      //region attribute (regional attribute or marker "A" must be switched)
      in.regionlist[2] = MeshId; 

      //region maximum volume attribute (maximum area attribute "a" (with no number following) must be switched)
      in.regionlist[3] = -1;


    }

    KRATOS_CATCH( "" )

  }


  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::RecoverBoundaryPosition(ModelPart& rModelPart,
							MeshingVariables& rMeshingVariables,
							struct triangulateio& in,
							struct triangulateio& out,
							ModelPart::IndexType MeshId)
  {
    KRATOS_TRY
    
    ModelPart::NodesContainerType::iterator nodes_begin = rModelPart.NodesBegin(MeshId);

    int base=0;
    for(unsigned int i = 0; i<rModelPart.Nodes(MeshId).size(); i++)
      { 
	   
	if( (nodes_begin + i)->Is(BOUNDARY) ){
	     
	  in.pointlist[base]   = (nodes_begin + i)->X();
	  in.pointlist[base+1] = (nodes_begin + i)->Y();
	     
	  out.pointlist[base]   = (nodes_begin + i)->X();
	  out.pointlist[base+1] = (nodes_begin + i)->Y();
	}
	   
	base+=2;
      }

    KRATOS_CATCH( "" )
  }


  //*******************************************************************************************
  //METHODS CALLED FOR THE TESSELLATION
  //*******************************************************************************************

  int TriangularMesh2DModeler::GenerateTriangulation(Flags& MeshingOptions,
						     Flags& RefiningOptions,
						     struct triangulateio& in,
						     struct triangulateio& out)
  {
    KRATOS_TRY

    int fail=0;

    struct triangulateio vorout;

    //initilize all to avoid memory problems
    ClearTrianglesList(vorout);

    //mesh options
    char  meshing_options[255];
    std::string meshing_info;
	
    //switches: https://www.cs.cmu.edu/~quake/triangle.switch.html


    if(MeshingOptions.Is(MeshModeler::CONSTRAINED_MESH) && MeshingOptions.Is(MeshModeler::NEIGHBOURS_SEARCH)){ //to mesh constrained delaunay
      strcpy (meshing_options, "pnBYYQ");
      meshing_info = "Constrained remeshing executed";
    }
    else if(MeshingOptions.Is(MeshModeler::BOUNDARIES_SEARCH) && MeshingOptions.Is(MeshModeler::NEIGHBOURS_SEARCH)){  //to get conectivities and boundaries only
      strcpy (meshing_options, "ncEBQ");
      meshing_info = "Boundaries-Neighbours remeshing executed";
    }
    else{

      if(MeshingOptions.Is(MeshModeler::BOUNDARIES_SEARCH)){  //to get conectivities and boundaries only
	strcpy (meshing_options, "rcEBQ");
	meshing_info = "Boundaries remeshing executed";
      }
	  
      if(MeshingOptions.Is(MeshModeler::NEIGHBOURS_SEARCH)){  //to get conectivities and neighbours only
	strcpy (meshing_options, "PneQ");
	meshing_info = "Neighbours remeshing executed";
      }

    }

    if(MeshingOptions.Is(MeshModeler::REFINE_MESH))
      {
	if(RefiningOptions.Is(MeshModeler::REFINE_INSERT_NODES)){ //to insert a set of given points and refine the mesh
	  strcpy (meshing_options, "riYYJQ");     // "riYYJQ"; //"riYYJQ" // "riJQ" //"riQ"
	  meshing_info = "Inserted remeshing executed";
	}
	else if(RefiningOptions.Is(MeshModeler::REFINE_ADD_NODES))  //to add_nodes automatically and refine the mesh ("q"-quality mesh and "a"-area constraint switches)
	  {
	    if( MeshingOptions.Is(MeshModeler::CONSTRAINED_MESH) )
	      {
		strcpy (meshing_options, "pYJq1.4arnCQ");  //"YYJaqrn" "YJq1.4arn" "Jq1.4arn"
		meshing_info = "Adaptive constrained remeshing executed";
	      }
	    else
	      {

		strcpy (meshing_options, "YJq1.4arnQ");  //"YYJaqrn" "YJq1.4arn" "Jq1.4arn"
		meshing_info = "Adaptive remeshing executed";
	      }


	  }
	else
	  //refine without adding nodes
	  {
	    strcpy (meshing_options, "YJrn");
	    meshing_info = "Non-Adaptive remeshing executed";
	  }

	    

      }

    if(MeshingOptions.Is(MeshModeler::RECONNECT)){  //to reconect a set of points only
      if(MeshingOptions.Is(MeshModeler::CONSTRAINED_MESH)){
	strcpy (meshing_options, "pBYYQ");
	meshing_info = "Constrained Reconnection remeshing executed";
      }
      else{
	strcpy (meshing_options, "QNP");
	meshing_info = "Reconnection remeshing executed";
      }
    }

    if( this->GetEchoLevel() > 0 )
      std::cout<<" [ REMESH: (in POINTS "<<in.numberofpoints<<") "<<std::endl;

    //perform the meshing
    try {
      triangulate (meshing_options,&in,&out,&vorout);
    }

    catch( int error_code ){

      switch(TriangleErrors(error_code))
	{
	case INPUT_MEMORY_ERROR:       fail=1;
	  break;
	case INTERNAL_ERROR:           fail=2;
	  break;
	case INVALID_GEOMETRY_ERROR:   fail=3;
	  break;
	default:                       fail=0;
	  //create new connections
	  if( this->GetEchoLevel() > 0 )
	    std::cout<<" triangulation done "<<std::endl;
	  break;
	}
    }
    
   

    if(MeshingOptions.IsNot(MeshModeler::REFINE_MESH) && in.numberofpoints<out.numberofpoints){
      fail=3;
      std::cout<<"  fail error: [NODES ADDED] something is wrong with the geometry "<<std::endl;
    }
    
    if( this->GetEchoLevel() > 0 ){
      std::cout<<"  -( "<<meshing_info<<" )- "<<std::endl;
      std::cout<<"  (out POINTS "<<out.numberofpoints<<") :  REMESH ]; "<<std::endl;
      std::cout<<std::endl;
    }

    return fail;
   
    KRATOS_CATCH( "" )

  }


  //*******************************************************************************************
  //METHODS CALLED AFTER TESSELLATION  
  //*******************************************************************************************
  
  //Set elements in model_part after the Delaunay Tesselation
  void TriangularMesh2DModeler::BuildMeshElements(ModelPart& rModelPart,
						  MeshingVariables& rMeshingVariables,
						  struct triangulateio &out,
						  ModelPart::IndexType MeshId)
  {
    KRATOS_TRY
    
    //*******************************************************************
    //clearing elements
    //rModelPart.Elements(MeshId).clear();

    //*******************************************************************
    //selecting elements
    rMeshingVariables.RefiningOptions.Set(MeshModeler::ENGAGED_NODES);
    SelectMeshElements(rModelPart.Nodes(MeshId),rMeshingVariables,out);
    rMeshingVariables.RefiningOptions.Reset(MeshModeler::ENGAGED_NODES);

    //*******************************************************************
    //setting new elements
    //(rModelPart.Elements(MeshId)).reserve(rMeshingVariables.Refine.NumberOfElements);

		
    //*******************************************************************
    //All nodes in boundary element change
    if(rMeshingVariables.AvoidTipElementsFlag){ //is not working correctly some dispositions not considered
      if( this->GetEchoLevel() > 0 )
	std::cout<<"[   AVOID TIP ELEMENTS START ]"<<std::endl;

      ChangeTipElementsUtilities TipElements;
      //TipElements.SwapDiagonals(rModelPart,out,rMeshingVariables.PreservedElements,MeshId);

      if( this->GetEchoLevel() > 0 )
	std::cout<<"[   AVOID TIP ELEMENTS END ]"<<std::endl;
    }
    //*******************************************************************


    //properties to be used in the generation
    int number_properties = rModelPart.NumberOfProperties();
    Properties::Pointer properties = rModelPart.GetMesh(MeshId).pGetProperties(number_properties-1);
    ModelPart::NodesContainerType::iterator nodes_begin = rModelPart.NodesBegin(MeshId);

    // properties->PrintData(std::cout);
    // std::cout<<std::endl;

    const Element & rReferenceElement=rMeshingVariables.GetReferenceElement();

    PointPointerVector list_of_element_centers;
    std::vector<Geometry<Node<3> > > list_of_element_vertices;
    //find the center and "radius" of the element
    double xc, yc, zc=0, radius;


    //generate kratos elements (conditions are not touched)
    int id = 0;
    std::vector<std::vector<int> > EmptyNeighList;
    rMeshingVariables.NeighbourList.swap(EmptyNeighList); 
    rMeshingVariables.NeighbourList.clear(); //destroy all elements

    int faces = 0;
    for(int el = 0; el<out.numberoftriangles; el++)
      {
	if(rMeshingVariables.PreservedElements[el])
	  {
	    Geometry<Node<3> > vertices;
	    std::vector<int >  neighbours (3);
	
	    for(int pn=0; pn<3; pn++)
	      {
		//note that out.trianglelist, starts from node 1, not from node 0, it can be directly assigned to rMeshingVariables.NodalPreIds.
		//vertices.push_back( *((model_nodes).find( rMeshingVariables.NodalPreIds[out.trianglelist[el*3+pn]] ).base() ) );
		vertices.push_back(*(nodes_begin + out.trianglelist[el*3+pn]-1).base());
		//vertices.push_back(rModelPart.pGetNode(out.trianglelist[el*3+pn],MeshId));
		   
		if(vertices.back().Is(TO_ERASE))
		  std::cout<<" WARNING:: mesh vertex RELEASED "<<vertices.back().Id()<<std::endl;
		  
		//std::cout<<" out.neighborlist "<<out.neighborlist[el*3+pn]<<std::endl;
		  		 
		if( out.neighborlist[el*3+pn]>0 )
		  {

		    if(rMeshingVariables.PreservedElements[ out.neighborlist[el*3+pn]-1 ])
		      {
			neighbours[pn]= out.neighborlist[el*3+pn];
		      }
		    else
		      {
			neighbours[pn]=-1;
			faces++;
		      }
		    
		  }
		else
		  {
		    neighbours[pn]=-1;
		    faces++;
		  }

		  

	      }


	    id += 1;

	    rMeshingVariables.PreservedElements[el] = id;
	    rMeshingVariables.NeighbourList.push_back(neighbours);
		
		
	    //std::cout<<" neigbours ["<<id<<"]: ("<<neighbours[0]<<", "<<neighbours[1]<<", "<<neighbours[2]<<") "<<std::endl;
		

	    //*******************************************************************
	    //1) Store Preserved elements in an array of vertices (Geometry<Node<3> > vertices;)

	    this->mDataTransferUtilities.CalculateCenterAndSearchRadius( vertices[0].X(), vertices[0].Y(),
								   vertices[1].X(), vertices[1].Y(),
								   vertices[2].X(), vertices[2].Y(),
								   xc, yc, radius );

	    //std::cout<<" XC ["<<id<<"]: ("<<xc<<" "<<yc<<") "<<std::endl;
	    //std::cout<<" vertices "<<vertices[0].X()<<" "<<vertices[2].X()<<std::endl;
	    //*******************************************************************
		
	    PointPointerType p_center = PointPointerType( new PointType(id,xc,yc,zc) );

	    //*******************************************************************
	    //2) Create list_of_centers 

	    list_of_element_centers.push_back( p_center );
	    list_of_element_vertices.push_back( vertices );
		
	    //*******************************************************************
		
	    // std::cout<<" list of centers "<<list_of_element_centers.back()->X()<<" "<<list_of_element_centers.back()->Y()<<std::endl;
	    // std::cout<<" list of vertices ";
	    // std::cout.flush();
	    // std::cout<<" vertices "<<list_of_element_vertices.back()[0].X()<<" "<<list_of_element_vertices.back()[2].X()<<std::endl;
	    // std::cout.flush();

		

	  }
	else{
	  rMeshingVariables.PreservedElements[el] = -1;
	}
    
	    
	// if(rMeshingVariables.PreservedElements[el])
	//   {
	// 	std::cout<<" Neighbours ["<<id-1<<"] :";
	// 	for(int pn=0; pn<3; pn++)
	// 	  {
	// 	    std::cout<<" neighborlist ("<<pn<<") "<<rMeshingVariables.NeighbourList[id-1][pn]<<std::endl;
	// 	  }
	//   }

	//std::cout<<" NodalPreIds ["<<el<<"] :"<<NodalPreIds[el]<<std::endl;
		
      }

    if( this->GetEchoLevel() > 0 )
      std::cout<<" [ FACES "<<faces<<"]"<<std::endl;
	
    //*******************************************************************
    //5) Laplacian Smoothing

    //Check Mesh Info to perform smoothing:
    rMeshingVariables.RemeshInfo.CheckGeometricalSmoothing();

    //if(rMeshingVariables.smoothing && rMeshingVariables.remesh && rMeshingVariables.RemeshInfo.GeometricalSmoothingRequired ){
    if( rMeshingVariables.MeshSmoothingFlag && rMeshingVariables.RemeshInfo.GeometricalSmoothingRequired ){
      LaplacianSmoothing  MeshGeometricSmoothing(rModelPart);
      MeshGeometricSmoothing.SetEchoLevel(this->GetEchoLevel());
      MeshGeometricSmoothing.ApplyMeshSmoothing(rModelPart,rMeshingVariables.PreservedElements,out,list_of_element_vertices,MeshId);
    }
    //*******************************************************************


    //*******************************************************************
    //6) Pass  rReferenceElement and transfer variables
    this->mDataTransferUtilities.TransferData (rModelPart,rReferenceElement,list_of_element_centers,list_of_element_vertices,MeshDataTransferUtilities::ELEMENT_TO_ELEMENT,MeshId);
    //*******************************************************************


    //*******************************************************************w
    //std::cout<<" Number of Nodes "<<rModelPart.Nodes(MeshId).size()<<" Number Of Ids "<<rMeshingVariables.NodalPreIds.size()<<std::endl;
	
    //7) Restore global ID's
    for(ModelPart::NodesContainerType::iterator in = rModelPart.NodesBegin(MeshId) ; in != rModelPart.NodesEnd(MeshId) ; in++)
      {
	//std::cout<<" node (local:"<<in->Id()<<", global:"<<rMeshingVariables.NodalPreIds[ in->Id() ]<<")"<<std::endl;
	in->SetId( rMeshingVariables.NodalPreIds[ in->Id() ] );
      }
    //*******************************************************************
	
	
    //*******************************************************************
    
    //8) Filling the neighbour list
    SetElementNeighbours(rModelPart,rMeshingVariables,MeshId);

    //*******************************************************************
	
    KRATOS_CATCH( "" )
  }



  //*******************************************************************************************
  //*******************************************************************************************

  //Select elements after the Delaunay Tesselation
  void TriangularMesh2DModeler::SelectMeshElements(ModelPart::NodesContainerType& rNodes,
						   MeshingVariables& rMeshingVariables,
						   struct triangulateio& out)
  {
    KRATOS_TRY

    if( this->GetEchoLevel() > 0 )
      std::cout<<" [ SELECT MESH ELEMENTS: ("<<(out.numberoftriangles)<<") "<<std::endl;

    rMeshingVariables.PreservedElements.clear();
    rMeshingVariables.PreservedElements.resize(out.numberoftriangles);
    std::fill( rMeshingVariables.PreservedElements.begin(), rMeshingVariables.PreservedElements.end(), 0 );
	
    rMeshingVariables.Refine.NumberOfElements=0;
    
    bool box_side_element = false;
    bool wrong_added_node = false;
    if(rMeshingVariables.RefiningOptions.IsNot(MeshModeler::SELECT_ELEMENTS))
      {

	for(int el=0; el<out.numberoftriangles; el++)
	  {
	    rMeshingVariables.PreservedElements[el]=1;
	    rMeshingVariables.Refine.NumberOfElements+=1;
	  }
      }
    else
      {
	if( this->GetEchoLevel() > 0 )
	  std::cout<<" Start Element Selection "<<out.numberoftriangles<<std::endl;
	int el;
	int number=0;
	//#pragma omp parallel for reduction(+:number) private(el)
	for(el=0; el<out.numberoftriangles; el++)
	  {
	    Geometry<Node<3> > vertices;
	    //double Alpha   = 0;
	    //double nodal_h = 0;
	    //double param   = 0.3333333;

	    // int  numflying=0;
	    // int  numlayer =0;
	    //int  numfixed =0;

	    int  numfreesurf =0;
	    int  numboundary =0;

	    // std::cout<<" num nodes "<<rNodes.size()<<std::endl;
	    // std::cout<<" selected vertices [ "<<out.trianglelist[el*3]<<", "<<out.trianglelist[el*3+1]<<", "<<out.trianglelist[el*3+2]<<"] "<<std::endl;
	    box_side_element = false;
	    for(int pn=0; pn<3; pn++)
	      {
		//set vertices
		if(rMeshingVariables.NodalPreIds[out.trianglelist[el*3+pn]]<0){
		  box_side_element = true;
		  break;
		}
		

		if( (unsigned int)out.trianglelist[el*3+pn] > rMeshingVariables.NodalPreIds.size() ){
		  wrong_added_node = true;
		  std::cout<<" ERROR: something is wrong: node out of bounds "<<std::endl;
		  break;
		}

		//vertices.push_back( *((rNodes).find( out.trianglelist[el*3+pn] ).base() ) );
		vertices.push_back(rNodes(out.trianglelist[el*3+pn]));

		//check flags on nodes
		if(vertices.back().Is(FREE_SURFACE))
		  numfreesurf++;

		if(vertices.back().Is(BOUNDARY))
		  numboundary++;

		// if(VertexPa[pn].match(_wall_))
		// 	numfixed++;

		// if(VertexPa[pn].match(_flying_))
		// 	numflying++;

		// if(VertexPa[pn].match(_layer_))
		// 	numlayer++;

		//nodal_h+=vertices.back().FastGetSolutionStepValue(NODAL_H);

	      }

	    

	    if(box_side_element || wrong_added_node){
	      //std::cout<<" Box_Side_Element "<<std::endl;
	      continue;
	    }

	    //1.- to not consider wall elements
	    // if(numfixed==3)
	    //   Alpha=0;

	    //2.- alpha shape:
	    //Alpha  = nodal_h * param;
	    //Alpha *= rMeshingVariables.AlphaParameter; //1.4; 1.35;

	    //2.1.- correction to avoid big elements on boundaries
	    // if(numflying>0){
	    //   Alpha*=0.8;
	    // }
	    // else{
	    //   if(numfixed+numsurf<=2){
	    //     //2.2.- correction to avoid voids in the fixed boundaries
	    //     if(numfixed>0)
	    // 	Alpha*=1.4;

	    //     //2.3.- correction to avoid voids on the free surface
	    //     if(numsurf>0)
	    // 	Alpha*=1.3;

	    //     //2.4.- correction to avoid voids in the next layer after fixed boundaries
	    //     if(numlayer>0 && !numsurf)
	    // 	Alpha*=1.2;
	    //   }

	    // }

	    //std::cout<<" ******** ELEMENT "<<el+1<<" ********** "<<std::endl;

	    double Alpha =  rMeshingVariables.AlphaParameter;
	    if(numboundary>=2)
	      Alpha*=1.8;
	
	    // std::cout<<" vertices for the contact element "<<std::endl;
	    // std::cout<<" (1): ["<<rMeshingVariables.NodalPreIds[vertices[0].Id()]<<"] "<<vertices[0]<<std::endl;
	    // std::cout<<" (2): ["<<rMeshingVariables.NodalPreIds[vertices[1].Id()]<<"] "<<vertices[1]<<std::endl;
	    // std::cout<<" (3): ["<<rMeshingVariables.NodalPreIds[vertices[2].Id()]<<"] "<<vertices[2]<<std::endl;

	    // std::cout<<" vertices for the subdomain element "<<std::endl;
	    // std::cout<<" (1): ["<<vertices[0].Id()<<"] "<<vertices[0]<<std::endl;
	    // std::cout<<" (2): ["<<vertices[1].Id()<<"] "<<vertices[1]<<std::endl;
	    // std::cout<<" (3): ["<<vertices[2].Id()<<"] "<<vertices[2]<<std::endl;

	    // std::cout<<" Element "<<el<<" with alpha "<<rMeshingVariables.AlphaParameter<<"("<<Alpha<<")"<<std::endl;

	    bool accepted=false;

	    if(rMeshingVariables.RefiningOptions.Is(MeshModeler::CONTACT_SEARCH))
	      {
		accepted=mModelerUtilities.ShrankAlphaShape(Alpha,vertices,rMeshingVariables.OffsetFactor,2);
	      }
	    else
	      {
		accepted=mModelerUtilities.AlphaShape(Alpha,vertices,2);
	      }
		
	   
	    //3.- to control all nodes from the same subdomain (problem, domain is not already set for new inserted particles on mesher)
	    // if(accepted)
	    // {
	    //   std::cout<<" Element passed Alpha Shape "<<std::endl;
	    //     if(rMeshingVariables.RefiningOptions.IsNot(MeshModeler::CONTACT_SEARCH))
	    //   	accepted=mModelerUtilities.CheckSubdomain(vertices);
	    // }

	    //3.1.-
	    bool self_contact = false;
	    if(rMeshingVariables.RefiningOptions.Is(MeshModeler::CONTACT_SEARCH))
	      self_contact = mModelerUtilities.CheckSubdomain(vertices);
	    	    
	    //4.- to control that the element is inside of the domain boundaries
	    if(accepted)
	      {
		if(rMeshingVariables.RefiningOptions.Is(MeshModeler::CONTACT_SEARCH))
		  {
		    accepted=mModelerUtilities.CheckOuterCentre(vertices,rMeshingVariables.OffsetFactor, self_contact);
		  }
		else
		  {
		    accepted=mModelerUtilities.CheckInnerCentre(vertices);
		  }
	      }
	    // else{
	    
	        
	    //   std::cout<<" Element DID NOT pass Alpha Shape ("<<Alpha<<") "<<std::endl;
	    // }
	      

	    if(accepted)
	      {
		//std::cout<<" Element ACCEPTED after cheking Center "<<number<<std::endl;
		rMeshingVariables.PreservedElements[el] = 1;
		number+=1;
	      }
	    // else{
	      
	    //   std::cout<<" Element DID NOT pass INNER/OUTER check "<<std::endl;
	    // }


	  }

	rMeshingVariables.Refine.NumberOfElements=number;

      }

    //std::cout<<" Number of Preserved Elements "<<rMeshingVariables.Refine.NumberOfElements<<std::endl;

    if(rMeshingVariables.RefiningOptions.Is(MeshModeler::ENGAGED_NODES)){

      //check engaged nodes
      for(int el=0; el<out.numberoftriangles; el++)
	{
	  if( rMeshingVariables.PreservedElements[el]){
	    for(int pn=0; pn<3; pn++)
	      {
		//set vertices
		rNodes[out.trianglelist[el*3+pn]].Set(MeshModeler::ENGAGED_NODES);
	      }
	  }
	    
	}

      int count_release = 0;
      for(ModelPart::NodesContainerType::iterator i_node = rNodes.begin() ; i_node != rNodes.end() ; i_node++)
	{
	  if( i_node->IsNot(MeshModeler::ENGAGED_NODES)  ){
	    i_node->Set(TO_ERASE);
	    if( this->GetEchoLevel() > 0 )
	      std::cout<<" NODE "<<i_node->Id()<<" RELEASE "<<std::endl;
	    if( i_node->IsNot(MeshModeler::ENGAGED_NODES) )
	      std::cout<<" ERROR: node "<<i_node->Id()<<" IS BOUNDARY RELEASE "<<std::endl;
	    count_release++;
	  }
	      
	  i_node->Reset(MeshModeler::ENGAGED_NODES);
	}
	  
      if( this->GetEchoLevel() > 0 )
	std::cout<<"   NUMBER OF RELEASED NODES "<<count_release<<std::endl;

    }

    if( this->GetEchoLevel() > 0 ){
      std::cout<<"   Generated_Elements :"<<out.numberoftriangles<<std::endl;
      std::cout<<"   Passed_AlphaShape  :"<<rMeshingVariables.Refine.NumberOfElements<<std::endl;
      std::cout<<"   SELECT MESH ELEMENTS ]; "<<std::endl;
    }

    KRATOS_CATCH( "" )

  }


  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::SetElementNeighbours(ModelPart& rModelPart,
						     MeshingVariables & rMeshingVariables,
						     ModelPart::IndexType MeshId)
  {
    KRATOS_TRY

    if( this->GetEchoLevel() > 0 ){
      std::cout<<" [ SET ELEMENT NEIGHBOURS : "<<std::endl;
      std::cout<<"   Initial Faces : "<<rModelPart.Conditions(MeshId).size()<<std::endl;
    }

    ModelPart::ElementsContainerType::const_iterator el_begin = rModelPart.ElementsBegin(MeshId);
	
    int facecounter=0;
    for(ModelPart::ElementsContainerType::const_iterator iii = rModelPart.ElementsBegin(MeshId);
	iii != rModelPart.ElementsEnd(MeshId); iii++)
      {

	int Id= iii->Id() - 1;
	//std::cout<<" Id ELNEIG "<<Id<<std::endl;

	(iii->GetValue(NEIGHBOUR_ELEMENTS)).resize(3);
	WeakPointerVector< Element >& neighb = iii->GetValue(NEIGHBOUR_ELEMENTS);

	for(int i = 0; i<3; i++)
	  {
	    int index = rMeshingVariables.NeighbourList[Id][i];
				
	    if(index > 0)
	      {
		//std::cout<<" Element "<<Id<<" size "<<rMeshingVariables.PreservedElements.size()<<std::endl;			    
		//std::cout<<" Index pre "<<index<<" size "<<rMeshingVariables.PreservedElements.size()<<std::endl;
		index = rMeshingVariables.PreservedElements[index-1];
		//std::cout<<" Index post "<<index<<std::endl;
	      }

	    if(index > 0)
	      {
		neighb(i) = *((el_begin + index -1 ).base());
	      }
	    else
	      {
		//neighb(i) = Element::WeakPointer();
		neighb(i) = *(iii.base());
		facecounter++;
	      }
	  }
      }
	
    if( this->GetEchoLevel() > 0 ){
      std::cout<<"   Final Faces : "<<facecounter<<std::endl;
      std::cout<<"   SET ELEMENT NEIGHBOURS ]; "<<std::endl;
    }

    KRATOS_CATCH( "" )

  }


  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::BuildMeshBoundary(ModelPart& rModelPart,
						  MeshingVariables& rMeshingVariables,
						  ModelPart::IndexType MeshId)
  {
    KRATOS_TRY

    if( this->GetEchoLevel() > 0 ){
      std::cout<<" [ SET BOUNDARY CONDITIONS : "<<std::endl;
      std::cout<<"   Initial Conditions : "<<rModelPart.Conditions(MeshId).size()<<std::endl;
    }

    //properties to be used in the generation
    int number_properties = rModelPart.NumberOfProperties();
    Properties::Pointer properties = rModelPart.GetMesh().pGetProperties(number_properties-1);

    //reset the boundary flag
    for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(MeshId); in!=rModelPart.NodesEnd(MeshId); in++)
      {
	in->Reset(BOUNDARY);
      }

    //filling the elemental neighbours list (from now on the elements list can not change)
    ModelPart::ElementsContainerType::const_iterator el_begin  = rModelPart.ElementsBegin(MeshId);
    ModelPart::ElementsContainerType::iterator elements_end = rModelPart.ElementsEnd(MeshId);


    //now the boundary faces
    int id=0;

    //set consecutive ids in the mesh conditions
    unsigned int condId=1;
    for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(MeshId); ic!= rModelPart.ConditionsEnd(MeshId); ic++)
      {
	Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();
	if( rConditionGeom[0].Is(TO_ERASE) || rConditionGeom[1].Is(TO_ERASE) )
	  ic->Set(TO_ERASE);
	    
	ic->SetId(condId);
	condId++;
      }
	

    //control the previous mesh conditions
    std::vector<int> PreservedConditions( rModelPart.Conditions(MeshId).size() );
    std::fill( PreservedConditions.begin(), PreservedConditions.end(), 0 );

    //swap previous conditions
    ModelPart::ConditionsContainerType temporal_conditions;
    temporal_conditions.reserve(rModelPart.Conditions(MeshId).size());
	    
    temporal_conditions.swap(rModelPart.Conditions(MeshId));


    // std::cout<<"   Preserved Conditions ( ";
    // for (unsigned int i=0; i<PreservedConditions.size(); i++)
    //   {
    //     std::cout<<" "<<PreservedConditions[i]<<"  ";
    //   }
    // std::cout<<" ) "<<std::endl;

    for(ModelPart::ElementsContainerType::const_iterator ie = rModelPart.ElementsBegin(MeshId); ie != rModelPart.ElementsEnd(MeshId); ie++)
      {
	int Id=ie->Id() -1 ;

	ModelPart::ElementsContainerType::iterator el_neighb;
	/*each face is opposite to the corresponding node number so in 2D
	  0 ----- 1 2
	  1 ----- 2 0
	  2 ----- 0 1
	*/

	//finding boundaries and creating the "skin"
	//
	//********************************************************************

	Geometry< Node<3> >& rGeom = ie->GetGeometry();
	boost::numeric::ublas::matrix<unsigned int> lpofa; //points that define the faces

	//Get the standard ReferenceCondition
	const Condition & rReferenceCondition=rMeshingVariables.GetReferenceCondition();

	for(unsigned int i = 0; i<rGeom.size(); i++)
	  {

	    int index = rMeshingVariables.NeighbourList[Id][i];
		
	    if(index > 0)
	      {
		index = rMeshingVariables.PreservedElements[index-1];
	      }
		
	    if( index > 0)
	      {
		el_neighb = (rModelPart.Elements(MeshId)).find( el_begin->Id() + index-1 ); //if not found-> returns the last element
	      }
	    else
	      {
		el_neighb = elements_end;
	      }

	      
	    if( el_neighb == elements_end )
	      {
		rGeom.NodesInFaces(lpofa);   

		//if no neighbour is present => the face is free surface
		rGeom[lpofa(1,i)].Set(BOUNDARY);
		rGeom[lpofa(2,i)].Set(BOUNDARY);
	
		//Get the correct ReferenceCondition
		Condition::Pointer pBoundaryCondition;
		bool condition_found = false;
		bool point_condition = false; 

		    
		for(ModelPart::ConditionsContainerType::iterator ic = temporal_conditions.begin(); ic!= temporal_conditions.end(); ic++)
		  {
		    if( ic->IsNot(TO_ERASE) ){

		      if( ic->IsNot(CONTACT) ){
			  
			Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();

			bool inserted = false;

			if( rConditionGeom[0].Is(NEW_ENTITY) || rConditionGeom[1].Is(NEW_ENTITY) ){
			  inserted = true;
			  // if(rConditionGeom[0].Is(NEW_ENTITY))
			  //   std::cout<<" node 0: Id-> "<<rConditionGeom[0].Id()<<" is new entity "<<std::endl;
			  // if(rConditionGeom[1].Is(NEW_ENTITY))
			  //   std::cout<<" node 1: Id-> "<<rConditionGeom[1].Id()<<" is new entity "<<std::endl;

			}

			if(ic->Is(NEW_ENTITY))
			  {
			    // std::cout<<" Inserted Cond "<<ic->Id()<<" nodes "<<rConditionGeom[0].Id()<<" "<<rConditionGeom[1].Id()<<std::endl;
			    inserted = false;
			  }


			if( !inserted ){
			      
			  if( PreservedConditions[ic->Id()-1] == 0 ){
			      
			    if( (rConditionGeom[0].Id() == rGeom[lpofa(1,i)].Id() 
				 && rConditionGeom[1].Id() == rGeom[lpofa(2,i)].Id() ) || 
				(rConditionGeom[0].Id() == rGeom[lpofa(2,i)].Id() 
				 && rConditionGeom[1].Id() == rGeom[lpofa(1,i)].Id() ) ){
				
			      pBoundaryCondition = (*(ic.base())); //accessing boost::shared_ptr  get() to obtain the raw pointer
			      PreservedConditions[ic->Id()-1] += 1; //add each time is used
			      if( rConditionGeom.PointsNumber() == 1 )
 			          point_condition = true;	

			      condition_found=true;
			    }
			  }
			    
			}
			else{
			    
			  if( PreservedConditions[ic->Id()-1] < 2 ){
			    
			    if( rConditionGeom[0].Id() == rGeom[lpofa(1,i)].Id() ||
				rConditionGeom[1].Id() == rGeom[lpofa(2,i)].Id() || 
				rConditionGeom[0].Id() == rGeom[lpofa(2,i)].Id() ||
				rConditionGeom[1].Id() == rGeom[lpofa(1,i)].Id() ){
				
			      pBoundaryCondition = (*(ic.base())); //accessing boost::shared_ptr  get() to obtain the raw pointer
			      PreservedConditions[ic->Id()-1] += 1; //add each time is used
			      if( rConditionGeom.PointsNumber() == 1 )
 			          point_condition = true;
	
			      condition_found=true;
			    }	

			    // std::cout<<" INSERTED COND "<<ic->Id()<<std::endl;
			  }			
			  
			}
			    			    			    			    
		      }
		      else{
			    
			PreservedConditions[ic->Id()-1] += 1;  //will not be restored
			//std::cout<<" Condition Contact "<<ic->Id()<<std::endl;

		      }
			
		    }
			  

		    if(condition_found==true){
		      // std::cout<<" Condition Found:  "<<ic->Id()<<" ("<<ic->GetGeometry()[0].Id()<<", "<<ic->GetGeometry()[1].Id()<<") == ("<<rGeom[lpofa(1,i)].Id()<<" "<<rGeom[lpofa(2,i)].Id()<<") ->  Used : "<<PreservedConditions[ic->Id()-1]<<" times "<<std::endl;
		      break;
		    }
		  }

		//Generate condition
		Condition::NodesArrayType face;
		face.reserve(2);
		face.push_back(rGeom(lpofa(1,i)));
		face.push_back(rGeom(lpofa(2,i)));
		id ++;

		//std::cout<<" id "<<id<<std::endl;

		Condition::Pointer p_cond;
		if(condition_found){
		  p_cond = pBoundaryCondition->Clone(id, face);
		      
		  p_cond->Data() =pBoundaryCondition->Data();

		  if( !point_condition ){
		      WeakPointerVector<Element > master_elems;
		      master_elems.push_back(Element::WeakPointer( *(ie.base()) ));
		      p_cond->SetValue(MASTER_ELEMENTS, master_elems );
		      
		      
		      WeakPointerVector<Node<3> > master_nodes;
		      master_nodes.push_back( Node<3>::WeakPointer( rGeom(lpofa(0,i)) ));
		      p_cond->SetValue(MASTER_NODES, master_nodes );
		  }    

		}
		else{
		  
		  if( this->GetEchoLevel() > 0 )
		    std::cout<<"   NOT FOUND CONDITION :: CREATED-> "<<id<<"("<<face[0].Id()<<","<<face[1].Id()<<")"<<std::endl;
		  p_cond = rReferenceCondition.Create(id, face, properties);
		      
		  //if a condition is created new nodes must be labeled TO_REFINE
		  face[0].Set(TO_REFINE);
		  face[1].Set(TO_REFINE);

		  Vector StressVector=ZeroVector(3);
		  Matrix DeformationGradient=identity_matrix<double>( 2 );

		  p_cond->SetValue(CAUCHY_STRESS_VECTOR,StressVector);
		  p_cond->SetValue(DEFORMATION_GRADIENT,DeformationGradient);
		  p_cond->GetValue(MASTER_ELEMENTS).push_back( Element::WeakPointer( *(ie.base()) ) );
		  p_cond->GetValue(MASTER_NODES).push_back( Node<3>::WeakPointer( rGeom(lpofa(0,i)) ) );
		}

		//usually one MasterElement and one MasterNode in 2D
		

		rModelPart.Conditions(MeshId).push_back(p_cond);

	      }

	  }

      }
    // std::cout<<"   Preserved Conditions ( ";
    // for (unsigned int i=0; i<PreservedConditions.size(); i++)
    //   {
    //     std::cout<<" "<<PreservedConditions[i]<<"  ";
    //   }
    // std::cout<<" ) "<<std::endl;

    if( this->GetEchoLevel() > 0 )
      std::cout<<"   Boundary Conditions LOCATED ["<<rModelPart.Conditions(MeshId).size()<<"]"<<std::endl;
    //all previous conditions have to be added
    for(ModelPart::ConditionsContainerType::iterator ic = temporal_conditions.begin(); ic!= temporal_conditions.end(); ic++)
      {
	bool node_not_preserved = false;
	bool condition_not_preserved = false;
	if( PreservedConditions[ic->Id()-1] == 0 ){

	  Geometry< Node<3> >& rGeom = ic->GetGeometry();
	  Condition::NodesArrayType face;

	  face.reserve(rGeom.size() );

	  for(unsigned int j=0; j<rGeom.size(); j++)
	    {
	      face.push_back(rGeom(j));
	    }

	  //if a condition is created new nodes must be labeled to REFINE
	  if( ic->Is(TO_ERASE) )
	    condition_not_preserved = true;
	  // if( face[0].IsNot(BOUNDARY) || face[1].IsNot(BOUNDARY) )
	  //   node_not_preserved = true;

	  if( face[0].Is(TO_ERASE) || face[1].Is(TO_ERASE) )
	    node_not_preserved = true;

	  if( face[0].Is(TO_REFINE) || face[1].Is(TO_REFINE) )
	    node_not_preserved = true;
	  
	  if(node_not_preserved == true || condition_not_preserved == true)
	    continue;

	  PreservedConditions[ic->Id()-1] += 1;

	  id +=1;

	  rModelPart.Conditions(MeshId).push_back(ic->Clone(id,face));

	  rModelPart.Conditions(MeshId).back().Data() = ic->Data();

	  if( this->GetEchoLevel() > 0 ){
	    std::cout<<" Temporal Condition Not Set "<<ic->Id()<<"("<<ic->GetGeometry()[0].Id()<<","<<ic->GetGeometry()[1].Id()<<")"<<std::endl;
	    std::cout<<" Push Back Not Set Conditions "<<id<<"("<<face[0].Id()<<","<<face[1].Id()<<")"<<std::endl;
	  }
	}
      }

    //control if previous conditions have been assigned
    int all_assigned = 0; 
    for(unsigned int i=0; i<PreservedConditions.size(); i++)
      {
	if( PreservedConditions[i] == 0 )
	  all_assigned ++;
      }

    if( this->GetEchoLevel() > 0 ){
      if(all_assigned == 0)
	std::cout<<"   Boundary Conditions RELOCATED ["<<all_assigned<<"]"<<std::endl;
      else
	std::cout<<"   Boundary Conditions NOT relocated ["<<all_assigned<<"]"<<std::endl;
	
      std::cout<<"   Final Conditions: "<<rModelPart.Conditions(MeshId).size()<<std::endl;
      std::cout<<"   SET BOUNDARY CONDITIONS ]; "<<std::endl;
    }
    //rModelPart.Conditions(MeshId).Sort();
    //rModelPart.Conditions(MeshId).Unique();
	
    KRATOS_CATCH( "" )

  }



  //*******************************************************************************************
  //METHODS CALLED BEFORE REFINING THE TESSELLATION
  //*******************************************************************************************

  void TriangularMesh2DModeler::SetDissipativeElements(ModelPart& rModelPart,
						       MeshingVariables& rMeshingVariables,
						       ModelPart::IndexType MeshId)
  {
    KRATOS_TRY
    
    ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
      
    //set label refine in elements that must be refined due to dissipation
    for(ModelPart::ElementsContainerType::const_iterator iii = rModelPart.ElementsBegin(MeshId);
	iii != rModelPart.ElementsEnd(MeshId); iii++)
      {
	double plastic_power=0;
	std::vector<double> Value(1);

	(iii)->GetValueOnIntegrationPoints(rMeshingVariables.Refine.GetDissipationVariable(),Value,CurrentProcessInfo);
	    
	plastic_power = Value[0] * iii->GetGeometry().Area();

	double critical_dissipation = rMeshingVariables.Refine.CriticalDissipation; // * iii->GetGeometry().Area();
	
	// if(plastic_power>0)
	//   std::cout<<" Element ["<<iii->Id()<<" plastic_power "<<plastic_power<<" CriticalDissipation "<<critical_dissipation<<" Area "<<iii->GetGeometry().Area()<<std::endl;

	if( plastic_power > critical_dissipation )
	  {
	    //std::cout<<" Refine element "<<std::endl;
	    Geometry< Node<3> >& rGeom = iii->GetGeometry();
	    for(unsigned int i = 0; i<rGeom.size(); i++)
	      {
		if(rGeom[i].IsNot(BOUNDARY))
		  rGeom[i].Set(TO_REFINE);
	      }
	  }
		    
      }
    
    KRATOS_CATCH( "" )
 
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::RefineBoundary(ModelPart& rModelPart,
        MeshingVariables& rMeshingVariables,
        ModelPart::IndexType MeshId)
  {

     KRATOS_TRY

        if( this->GetEchoLevel() > 0 ){
           std::cout<<" [ REFINE BOUNDARY : "<<std::endl;
           //std::cout<<"   Nodes and Conditions : "<<rModelPart.Nodes(MeshId).size()<<", "<<rModelPart.Conditions(MeshId).size()<<std::endl;
        }

     rMeshingVariables.RemeshInfo.InsertedConditions     = rModelPart.NumberOfConditions(MeshId);
     rMeshingVariables.RemeshInfo.InsertedBoundaryNodes = rModelPart.NumberOfNodes(MeshId);

     //***SIZES :::: parameters do define the tolerance in mesh size: 

     //DEFORMABLE CONTACT:
     double factor_for_tip_radius     = 0.2; //deformable contact tolerance in radius for detection tip sides to refine
     double factor_for_non_tip_side   = 3.0; // will be multiplied for nodal_h of the master node to compare with boundary nodes average nodal_h in a contact conditio which master node do not belongs to a tip

     double size_for_tip_contact_side      = 0.4 * rMeshingVariables.Refine.CriticalSide; // length size for the contact tip side
     double size_for_non_tip_contact_side  = 2.0 * rMeshingVariables.Refine.CriticalSide; //compared with contact size wich master node do not belongs to a tip

     //RIGID WALL CONTACT:
     double size_for_wall_tip_contact_side      = 0.50 * rMeshingVariables.Refine.CriticalSide; 
     double size_for_wall_semi_tip_contact_side = 0.75 * rMeshingVariables.Refine.CriticalSide; // semi contact or contact which a node in a tip
     double size_for_wall_non_tip_contact_side  = 1.25 * rMeshingVariables.Refine.CriticalSide; // semi contact or contact which no node in a tip

     //NON CONTACT:
     double size_for_energy_side                = 1.50 * rMeshingVariables.Refine.CriticalSide; // non contact side which dissipates energy
     double size_for_non_contact_side           = 3.0  * rMeshingVariables.Refine.CriticalSide;


     ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
     SpatialBoundingBox RefiningBox (rMeshingVariables.BoundingBox.Center,rMeshingVariables.BoundingBox.Radius,rMeshingVariables.BoundingBox.Velocity);

     unsigned int  conditions_size = 0;

     //counters:
     int total_contact_conditions = 0;
     int number_contacts_domain   = 0;
     int number_contacts_active   = 0;
     int contact_size   = 0;   
     int contact_tip    = 0;
     int exterior_bound = 0;
     int tip_bound      = 0;
     int energy_bound   = 0;



     //*********************************************************************************
     // DETECTION OF NODES ON TIP CONTACTS START
     //*********************************************************************************

     unsigned int nodes_on_wall_tip = 0;
     for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(MeshId); in!=rModelPart.NodesEnd(MeshId); in++)
     {
        if( mModelerUtilities.CheckNodeCloseWallTip(rMeshingVariables.RigidWalls,(*in),CurrentProcessInfo,factor_for_tip_radius) ){
           in->Set(TO_SPLIT);
           nodes_on_wall_tip ++;
        }

     }

     if( this->GetEchoLevel() > 0 )
        std::cout <<"   [ NODES ON WALL TIP: ( " <<nodes_on_wall_tip <<" ) ]"<<std::endl;

     //*********************************************************************************
     // DETECTION OF NODES ON TIP CONTACTS END
     //*********************************************************************************




     //if the insert switches are activated, we check if the boundaries got too coarse
     if (rMeshingVariables.RefiningOptions.Is(MeshModeler::REFINE_INSERT_NODES) && rMeshingVariables.RefiningOptions.Is(MeshModeler::REFINE_BOUNDARY) )
     {

        PointVector list_of_nodes;
        std::vector<Condition::Pointer> list_of_conditions;

        conditions_size = rModelPart.Conditions(MeshId).size();
        list_of_nodes.reserve(conditions_size);
        list_of_conditions.reserve(conditions_size);

        // std::vector<int> nodes_ids;
        // nodes_ids.resize(rModelPart.Conditions().size()); //mesh 0
        // std::fill( nodes_ids.begin(), nodes_ids.end(), 0 );

        //std::cout<<"   List of Conditions Reserved Size: "<<conditions_size<<std::endl;

        double tool_radius= 0;
        double side_length= 0;
        double plastic_power=0;
        bool size_insert = false;
        bool radius_insert = false;
        bool energy_insert = false;
        bool mesh_size_insert = false;
        bool contact_active = false;
        bool contact_semi_active = false;
        bool tool_project = false;

        std::vector<bool> semi_active_nodes;
        Node<3> new_point(0,0.0,0.0,0.0);

        //*********************************************************************************
        // DEFORMABLE CONTACT CONDITIONS START
        //*********************************************************************************

        for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(); ic!= rModelPart.ConditionsEnd(); ic++)
        {

           size_insert    = false;
           radius_insert  = false;
           energy_insert  = false;
           tool_project   = false;
           contact_active = false;
           side_length = 0;
           tool_radius = 0;
           plastic_power = 0;
           Geometry< Node<3> > rConditionGeom;
           array_1d<double,3> tip_center;
	   tip_center.clear();


           //LOOP TO CONSIDER ONLY CONTACT CONDITIONS
           if( ic->Is(CONTACT) )   //Refine radius on the workpiece for the ContactDomain zone
           {

              PointType  MasterNode;
              bool condition_found = false;
              Condition::Pointer MasterCondition  = mModelerUtilities.FindMasterCondition(*(ic.base()),MasterNode,rModelPart.Conditions(MeshId),condition_found);


              if(condition_found){

                 if( MasterCondition->IsNot(TO_ERASE) ){


                    rConditionGeom  = MasterCondition->GetGeometry(); 

                    //to recover TIP definition on conditions		  
                    if( MasterNode.SolutionStepsDataHas( WALL_TIP_RADIUS ) ) //master node in tool -->  refine workpiece  // 
                    {

                       tool_radius = MasterNode.FastGetSolutionStepValue( WALL_TIP_RADIUS );
                       tip_center  = MasterNode.FastGetSolutionStepValue( WALL_REFERENCE_POINT );
                       // WARNING THE UPDATED OF THE TIP CENTER IS NEEDED !!!!

                       array_1d<double, 3 > radius;
                       radius[0]=rConditionGeom[0].X()-tip_center[0];
                       radius[1]=rConditionGeom[0].Y()-tip_center[1];
                       radius[2]=rConditionGeom[0].Z()-tip_center[2];
                       double distance1=norm_2(radius);

                       radius[0]=rConditionGeom[1].X()-tip_center[0];
                       radius[1]=rConditionGeom[1].Y()-tip_center[1];
                       radius[2]=rConditionGeom[1].Z()-tip_center[2];

                       double distance2=norm_2(radius);


                       // TO SPLIT DETECTION START
                       //If a node is detected in the wall tip is set TO_SPLIT
                       //the criteria to splitting will be applied later in the nodes marked as TO_SPLIT

                       if( (1-factor_for_tip_radius)*tool_radius < distance1 &&  distance1 < (1+factor_for_tip_radius)*tool_radius )
                          rConditionGeom[0].Set(TO_SPLIT);

                       if( (1-factor_for_tip_radius)*tool_radius < distance2 &&  distance2 < (1+factor_for_tip_radius)*tool_radius )
                          rConditionGeom[1].Set(TO_SPLIT);

                       // TO SPLIT DETECTION END			  


                       // ACTIVE CONTACT DETECTION START

                       contact_active = mModelerUtilities.CheckContactActive(rConditionGeom, contact_semi_active, semi_active_nodes);
                       if(contact_active){
                          number_contacts_active ++;
                       }

                       // ACTIVE CONTACT DETECTION END


                       side_length = mModelerUtilities.CalculateSideLength (rConditionGeom[0],rConditionGeom[1]);		    	     

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
                             // double& nodal_h1 = rConditionGeom[0].FastGetSolutionStepValue(NODAL_H);
                             // double& nodal_h2 = rConditionGeom[1].FastGetSolutionStepValue(NODAL_H);
                             // double& nodal_h0 = MasterNode.FastGetSolutionStepValue( NODAL_H );

                             // double side = norm_2(rConditionGeom[0]-rConditionGeom[1]);
                             // // double d1 = mModelerUtilities.FindBoundaryH (rConditionGeom[0]);
                             // // double d2 = mModelerUtilities.FindBoundaryH (rConditionGeom[1]);
                             // // double d0 = mModelerUtilities.FindBoundaryH (MasterNode);
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

                       double& nodal_h1 = rConditionGeom[0].FastGetSolutionStepValue(NODAL_H);
                       double& nodal_h2 = rConditionGeom[1].FastGetSolutionStepValue(NODAL_H);
                       double& nodal_h0 = MasterNode.FastGetSolutionStepValue( NODAL_H );

                       double side = norm_2(rConditionGeom[0]-rConditionGeom[1]);
                       // double d1 = mModelerUtilities.FindBoundaryH (rConditionGeom[0]);
                       // double d2 = mModelerUtilities.FindBoundaryH (rConditionGeom[1]);
                       // double d0 = mModelerUtilities.FindBoundaryH (MasterNode);
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

                       new_point.X() = 0.5*( rConditionGeom[1].X() + rConditionGeom[0].X() );
                       new_point.Y() = 0.5*( rConditionGeom[1].Y() + rConditionGeom[0].Y() );
                       new_point.Z() = 0.5*( rConditionGeom[1].Z() + rConditionGeom[0].Z() );


                       new_point.SetId(ic->Id()); //set condition Id

                       Condition::Pointer ContactMasterCondition  = ic->GetValue(MASTER_CONDITION);


                       if( (rConditionGeom[0].Is(TO_SPLIT) && rConditionGeom[1].Is(TO_SPLIT)) )
                          tool_project = true;

                       if( (rConditionGeom[0].Is(TO_SPLIT) || rConditionGeom[1].Is(TO_SPLIT)) && contact_active)
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
                          contact_tip++;
                       if(size_insert)
                          contact_size++;

                       // std::cout<<"   MasterCondition RELEASED (Id: "<<ContactMasterCondition->Id()<<") "<<std::endl;
                       ContactMasterCondition->Set(TO_ERASE);
                       list_of_nodes.push_back(new_point);
                       list_of_conditions.push_back(ContactMasterCondition);
                    }		    
                 }

                 number_contacts_domain ++;
              }
              // else{

              //   std::cout<<"   Master Condition not found "<<std::endl;

              // }

              total_contact_conditions ++;

           }
        }

        // std::cout<<"   [ Contact Conditions : "<<total_contact_conditions<<", (contacts in domain: "<<number_contacts_domain<<", of them active: "<<number_contacts_active<<") ] "<<std::endl;
        // std::cout<<"   Contact Search End ["<<list_of_conditions.size()<<" : "<<list_of_nodes.size()<<"]"<<std::endl;

        //*********************************************************************************
        // DEFORMABLE CONTACT CONDITIONS END
        //*********************************************************************************



        //*********************************************************************************
        // RIGID CONTACT CONDITIONS AND OTHER BOUNDARY CONDITIONS START
        //*********************************************************************************


        //LOOP TO CONSIDER ALL SUBDOMAIN CONDITIONS
        double cond_counter=0;
        for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(MeshId); ic!= rModelPart.ConditionsEnd(MeshId); ic++)
        {
           cond_counter ++;
           bool refine_candidate = false;
           if( rMeshingVariables.MeshingOptions.Is(MeshModeler::CONSTRAINED_MESH) ){
              if( ic->Is(BOUNDARY) ) //ONLY SET TO THE BOUNDARY SKIN CONDITIONS (CompositeCondition)
                 refine_candidate = true;
              else
                 refine_candidate = false;
           }
           else{
              refine_candidate = true; 
           }


           if( refine_candidate ){
              if (rMeshingVariables.BoundingBox.IsSetFlag == true ){
                 refine_candidate = mModelerUtilities.CheckConditionInBox(*(ic.base()),RefiningBox,CurrentProcessInfo);
              }
           }


           if( refine_candidate ){

              radius_insert = false;
              energy_insert = false;
              mesh_size_insert = false;
              tool_project = false;
              contact_active = false;
              contact_semi_active = false;
              side_length = 0;
              tool_radius = 0;
              plastic_power = 0;
              //double condition_radius = 0;
              Geometry< Node<3> > rConditionGeom;
              array_1d<double,3> tip_center;

              if( ic->IsNot(TO_ERASE) ){

                 //*********************************************************************************
                 // RIGID CONTACT CONDITIONS ON TIP START
                 //*********************************************************************************

                 // TOOL TIP INSERT;


                 // ACTIVE CONTACT DETECTION START

                 rConditionGeom = ic->GetGeometry();
                 contact_active = mModelerUtilities.CheckContactActive(rConditionGeom, contact_semi_active, semi_active_nodes);

                 // ACTIVE CONTACT DETECTION END


                 if( contact_active ){

                    side_length = mModelerUtilities.CalculateSideLength (rConditionGeom[0],rConditionGeom[1]);		    	     

                    if( side_length > size_for_wall_tip_contact_side ){

                       bool on_tip = false;
                       if(rConditionGeom[0].Is(TO_SPLIT) && rConditionGeom[1].Is(TO_SPLIT)){
                          on_tip = true;
                       }
                       else if (rConditionGeom[0].Is(TO_SPLIT) || rConditionGeom[1].Is(TO_SPLIT)){
                          if( side_length > size_for_wall_tip_contact_side ){
                             on_tip = true;
                          }
                       }		  

                       bool on_radius = false;

                       if( on_tip && rMeshingVariables.RigidWallSetFlag ){

                          Vector Point(3);
                          if( rConditionGeom[0].Is(TO_SPLIT) ){

                             Point[0] = rConditionGeom[0].X();
                             Point[1] = rConditionGeom[0].Y();
                             Point[2] = rConditionGeom[0].Z();
                             on_radius = true;

                          }
                          else if( rConditionGeom[1].Is(TO_SPLIT) ){

                             Point[0] = rConditionGeom[1].X();
                             Point[1] = rConditionGeom[1].Y();
                             Point[2] = rConditionGeom[1].Z();
                             on_radius = true;

                          }
                          else{
                             on_radius = false;
                          }


                          if( on_radius ){

                             ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
                             double Time = CurrentProcessInfo[TIME];  
                             for( unsigned int i = 0; i < rMeshingVariables.RigidWalls.size(); i++ )
                             {
                                if( rMeshingVariables.RigidWalls[i]->IsInside( Point, Time ) ){
                                   tool_radius = rMeshingVariables.RigidWalls[i]->GetRadius(Point);
                                   tip_center  = rMeshingVariables.RigidWalls[i]->GetCenter(Point);
                                   break;
                                }
                             }
                          }

                       }

                       if( on_radius && on_tip ) //master node in tool -->  refine workpiece  // (tool_radius ==0 in workpiece nodes)
                       {
                          PointType center (0,tip_center[0],tip_center[1],tip_center[2]);
                          array_1d<double, 3 > radius;
                          radius[0]=rConditionGeom[0].X()-center.X();
                          radius[1]=rConditionGeom[0].Y()-center.Y();
                          radius[2]=rConditionGeom[0].Z()-center.Z();

                          double distance1=norm_2(radius);

                          radius[0]=rConditionGeom[1].X()-center.X();
                          radius[1]=rConditionGeom[1].Y()-center.Y();
                          radius[2]=rConditionGeom[1].Z()-center.Z();

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

                 if (!radius_insert && rMeshingVariables.RefiningOptions.Is(MeshModeler::CRITERION_ENERGY) && vsize>0){

                    Element::ElementType& MasterElement = ic->GetValue(MASTER_ELEMENTS)[vsize-1];

                    plastic_power=0;
                    std::vector<double> Value(1);

                    MasterElement.GetValueOnIntegrationPoints(rMeshingVariables.Refine.GetDissipationVariable(),Value,CurrentProcessInfo);


                    Geometry<Node<3> >& pGeom = MasterElement.GetGeometry();
                    plastic_power = Value[0] * pGeom.Area();

                    //computation of the condition master element radius start: 
                    //PointsArrayType& vertices = pGeom.Points();

                    // double average_side_length= mModelerUtilities.CalculateAverageSideLength (vertices[0].X(), vertices[0].Y(),
                    // 							      vertices[1].X(), vertices[1].Y(),
                    // 							      vertices[2].X(), vertices[2].Y());
                    //condition_radius = pGeom.Area()/average_side_length;
                    //computation of the condition master element radius end;

                    //condition_radius is side_length
                    side_length = mModelerUtilities.CalculateSideLength (rConditionGeom[0],rConditionGeom[1]);

                    //condition_radius = mModelerUtilities.CalculateTriangleRadius (pGeom);

                    //if( plastic_power > rMeshingVariables.Refine.CriticalDissipation && condition_radius > rMeshingVariables.Refine.CriticalRadius )
                    if( plastic_power > rMeshingVariables.Refine.CriticalDissipation * pGeom.Area() && side_length > size_for_energy_side )
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
                    double Alpha =  rMeshingVariables.AlphaParameter;

                    bool accepted = mModelerUtilities.AlphaShape(Alpha,vertices,2);


                    //condition_radius is side_length
                    side_length = mModelerUtilities.CalculateSideLength (rConditionGeom[0],rConditionGeom[1]);

                    //condition_radius = mModelerUtilities.CalculateTriangleRadius (pGeom);
                    double critical_side_size = 0;

                    bool on_tip = false;
                    if( contact_semi_active ){

                       if (rConditionGeom[0].Is(TO_SPLIT) || rConditionGeom[1].Is(TO_SPLIT))
                          on_tip = true;

                       if( on_tip == true )
                          critical_side_size = size_for_wall_semi_tip_contact_side;
                       else
                          critical_side_size = size_for_wall_non_tip_contact_side;



                    }
                    else if( contact_active ){

                       if (rConditionGeom[0].Is(TO_SPLIT) || rConditionGeom[1].Is(TO_SPLIT))
                          on_tip = true;

                       if( on_tip == true )
                          critical_side_size = size_for_wall_semi_tip_contact_side;
                       else
                          critical_side_size = size_for_wall_non_tip_contact_side;

                    }
                    else{

                       critical_side_size  = size_for_non_contact_side;
                    }


                    //if(plastic_power > rMeshingVariables.Refine.CriticalDissipation && condition_radius > rMeshingVariables.Refine.CriticalRadius)
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
                 // RIGID CONTACT CONDITIONS AND OTHER BOUNDARY CONDITIONS END
                 //*********************************************************************************

                 //*********************************************************************************
                 //                   BOUNDARY REBUILD START                                      //
                 //*********************************************************************************


                 if( radius_insert || energy_insert || mesh_size_insert ) //Boundary must be rebuild 
                 {

                    // std::cout<<"   BOUNDARY DOMAIN ELEMENT REFINED "<<ic->Id()<<std::endl;

                    new_point.X() = 0.5*( rConditionGeom[1].X() + rConditionGeom[0].X() );
                    new_point.Y() = 0.5*( rConditionGeom[1].Y() + rConditionGeom[0].Y() );
                    new_point.Z() = 0.5*( rConditionGeom[1].Z() + rConditionGeom[0].Z() );

                    if( this->GetEchoLevel() > 0 )
                       std::cout<<"   NEW NODE  "<<new_point<<std::endl;

                    new_point.SetId(ic->Id()); //set condition Id

                    //it will be good if the node is detected in the tool tip using the rigid contact standards:

                    if( (rConditionGeom[0].Is(TO_SPLIT) && rConditionGeom[1].Is(TO_SPLIT)) )
                       tool_project = true;

                    if( (rConditionGeom[0].Is(TO_SPLIT) || rConditionGeom[1].Is(TO_SPLIT)) && contact_active)
                       tool_project = true;

                    if( (rConditionGeom[0].Is(TO_SPLIT) || rConditionGeom[1].Is(TO_SPLIT)) && contact_semi_active)
                       tool_project = true;

                    bool on_radius = false;

                    if( tool_project && rMeshingVariables.RigidWallSetFlag ){

                       Vector Point(3);
                       if( rConditionGeom[0].Is(TO_SPLIT) ){

                          Point[0] = rConditionGeom[0].X();
                          Point[1] = rConditionGeom[0].Y();
                          Point[2] = rConditionGeom[0].Z();
                          on_radius = true;

                       }
                       else if( rConditionGeom[1].Is(TO_SPLIT) ){

                          Point[0] = rConditionGeom[1].X();
                          Point[1] = rConditionGeom[1].Y();
                          Point[2] = rConditionGeom[1].Z();
                          on_radius = true;

                       }
                       else{
                          on_radius = false;
                       }


                       if( on_radius ){

                          ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
                          double Time = CurrentProcessInfo[TIME];  
                          for( unsigned int i = 0; i < rMeshingVariables.RigidWalls.size(); i++ )
                          {
                             if( rMeshingVariables.RigidWalls[i]->IsInside( Point, Time ) ){
                                tool_radius = rMeshingVariables.RigidWalls[i]->GetRadius(Point);
                                tip_center  = rMeshingVariables.RigidWalls[i]->GetCenter(Point);
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

                          if( this->GetEchoLevel() > 0 )
                             std::cout<<"   TOOL PROJECT::on radius  "<<new_point<<std::endl;


                       }


                    }

                    if(radius_insert)
                       tip_bound ++;
                    if(energy_insert)
                       energy_bound ++;
                    if(mesh_size_insert)
                       exterior_bound++;

                    ic->Set(TO_ERASE);

                    if( this->GetEchoLevel() > 0 )
                       std::cout<<"   INSERTED NODE  "<<new_point<<std::endl;

                    list_of_nodes.push_back(new_point);
                    list_of_conditions.push_back(*(ic.base()));



                    // std::cout<<"   Refine Boundary  (Id:"<<ic->Id()<<"): ["<<rConditionGeom[0].Id()<<", "<<rConditionGeom[1].Id()<<"]"<<std::endl;
                    // std::cout<<"   (x1:"<<rConditionGeom[0].X()<<", y1: "<<rConditionGeom[0].Y()<<") "<<" (x2:"<<rConditionGeom[1].X()<<", y2: "<<rConditionGeom[1].Y()<<") "<<std::endl;

                    //std::cout<<" Added Node [Rcrit:"<<condition_radius<<",Scrit:"<<side_length<<",PlasticPower:"<<plastic_power<<"]"<<std::endl;
                    //std::cout<<" Conditions [Rcrit:"<<rMeshingVariables.Refine.CriticalRadius<<",Scrit:"<<rMeshingVariables.Refine.CriticalSide<<",PlasticPower:"<<rMeshingVariables.Refine.CriticalDissipation<<"]"<<std::endl;

                 }


                 //*********************************************************************************
                 //                   BOUNDARY REBUILD END                                        //
                 //*********************************************************************************


              }
              else{
                 if( this->GetEchoLevel() > 0 )
                    std::cout<<" Condition "<<ic->Id()<<" Released "<<std::endl;
              }

           }
        }	  


        //*********************************************************************************
        //                   DOFS AND NEW CONDITIONS REBUILD START                       //
        //*********************************************************************************

        //node to get the DOFs from
        Node<3>::DofsContainerType& reference_dofs = (rModelPart.NodesBegin(MeshId))->GetDofs();
        unsigned int step_data_size = rModelPart.GetNodalSolutionStepDataSize();
        double z = 0.0;

        unsigned int initial_node_size = rModelPart.Nodes().size()+1; //total model part node size
        unsigned int initial_cond_size = rModelPart.Conditions().size()+1; //total model part node size
        int id=0;

        //if points were added, new nodes must be added to ModelPart
        for(unsigned int i = 0; i<list_of_nodes.size(); i++)
        {
           id   = initial_node_size+i;

           double& x= list_of_nodes[i].X();
           double& y= list_of_nodes[i].Y();

           Node<3>::Pointer pnode = rModelPart.CreateNewNode(id,x,y,z);

           //set to the main mesh (Mesh 0) to avoid problems in the NodalPreIds (number of nodes: change) in other methods
           if(MeshId!=0)
              (rModelPart.Nodes(MeshId)).push_back(pnode);


           pnode->SetBufferSize(rModelPart.NodesBegin(MeshId)->GetBufferSize() );


           //assign data to dofs
           unsigned int buffer_size = pnode->GetBufferSize();

           //2D edges:
           Geometry< Node<3> >& rConditionGeom  = list_of_conditions[i]->GetGeometry(); 


           //generating the dofs
           for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
           {
              Node<3>::DofType& rDof = *iii;
              Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );

              if( rConditionGeom[0].IsFixed(rDof.GetVariable()) && rConditionGeom[1].IsFixed(rDof.GetVariable()) )
                 (p_new_dof)->FixDof();
              else
                 (p_new_dof)->FreeDof();
           }



           //int cond_id = list_of_nodes[i].Id();
           //Geometry< Node<3> >& rConditionGeom = (*(rModelPart.Conditions(MeshId).find(cond_id).base()))->GetGeometry();

           for(unsigned int step = 0; step<buffer_size; step++)
           {
              //getting the data of the solution step
              double* step_data = (pnode)->SolutionStepData().Data(step);

              double* node0_data = rConditionGeom[0].SolutionStepData().Data(step);
              double* node1_data = rConditionGeom[1].SolutionStepData().Data(step);

              //copying this data in the position of the vector we are interested in
              for(unsigned int j= 0; j<step_data_size; j++)
              {
                 step_data[j] = 0.5*node0_data[j] + 0.5*node1_data[j];
              }
           }


           //set specific control values and flags:
           pnode->Set(BOUNDARY);
           pnode->Set(NEW_ENTITY);  //if boundary is rebuild, the flag INSERTED must be set to new conditions too
           //std::cout<<"   Node ["<<pnode->Id()<<"] is a NEW_ENTITY "<<std::endl;

           pnode->SetValue(DOMAIN_LABEL,MeshId);
           double& nodal_h = pnode->FastGetSolutionStepValue(NODAL_H);
           //nodal_h = 0.5*(nodal_h+rMeshingVariables.Refine.CriticalSide); //modify nodal_h for security
           nodal_h = rMeshingVariables.Refine.CriticalSide; //modify nodal_h for security

           const array_1d<double,3> ZeroNormal(3,0.0);
           //correct normal interpolation
           noalias(pnode->GetSolutionStepValue(NORMAL)) = list_of_conditions[i]->GetValue(NORMAL);


           //correct contact_normal interpolation (laplacian boundary projection uses it)
           array_1d<double, 3 > & ContactForceNormal1  = rConditionGeom[0].FastGetSolutionStepValue(CONTACT_FORCE);
           array_1d<double, 3 > & ContactForceNormal2  = rConditionGeom[1].FastGetSolutionStepValue(CONTACT_FORCE);
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

           face1.push_back(rConditionGeom(0));
           face1.push_back(pnode);

           face2.push_back(pnode);
           face2.push_back(rConditionGeom(1));

           id   = initial_cond_size+(i*2);

           Condition::Pointer pcond1      = list_of_conditions[i]->Clone(id, face1);
           // std::cout<<" ID"<<id<<" 1s "<<pcond1->GetGeometry()[0].Id()<<" "<<pcond1->GetGeometry()[1].Id()<<std::endl;
           id   = initial_cond_size+(i*2+1);
           Condition::Pointer pcond2      = list_of_conditions[i]->Clone(id, face2);
           // std::cout<<" ID"<<id<<" 2s "<<pcond2->GetGeometry()[0].Id()<<" "<<pcond2->GetGeometry()[1].Id()<<std::endl;

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

           (rModelPart.Conditions(MeshId)).push_back(pcond1);
           (rModelPart.Conditions(MeshId)).push_back(pcond2);

           // (rModelPart.Conditions()).push_back(pcond1);
           // (rModelPart.Conditions()).push_back(pcond2);
        }


        //*********************************************************************************
        //                   DOFS AND NEW CONDITIONS REBUILD END                         //
        //*********************************************************************************




        //*********************************************************************************
        //                   CLEAN CONDITIONS AND FLAGS START                            //
        //*********************************************************************************


        //Clean Conditions
        ModelPart::ConditionsContainerType RemoveConditions;

        //id = 0;
        for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(MeshId); ic!= rModelPart.ConditionsEnd(MeshId); ic++)
        {

           Geometry< Node<3> > rGeom =ic->GetGeometry();
           for(unsigned int i=0; i<rGeom.size(); i++)
           {
              rGeom[i].Reset(TO_SPLIT);
           }

           if(ic->IsNot(TO_ERASE)){
              //id+=1;
              RemoveConditions.push_back(*(ic.base()));
              //RemoveConditions.back().SetId(id);
           }
           // else{
           //   std::cout<<"   Condition RELEASED:"<<ic->Id()<<std::endl;
           // }
        }

        rModelPart.Conditions(MeshId).swap(RemoveConditions);


        //*********************************************************************************
        //                   CLEAN CONDITIONS AND FLAGS END                              //
        //*********************************************************************************




     } // REFINE END;



     rMeshingVariables.RemeshInfo.InsertedConditions     = rModelPart.NumberOfConditions(MeshId)-rMeshingVariables.RemeshInfo.InsertedConditions;
     rMeshingVariables.RemeshInfo.InsertedBoundaryNodes = rModelPart.NumberOfNodes(MeshId)-rMeshingVariables.RemeshInfo.InsertedBoundaryNodes;

     if( this->GetEchoLevel() > 0 ){
        std::cout<<"   [ CONDITIONS ( inserted : "<<rMeshingVariables.RemeshInfo.InsertedConditions<<" ) ]"<<std::endl;
        std::cout<<"   [ NODES      ( inserted : "<<rMeshingVariables.RemeshInfo.InsertedBoundaryNodes<<" ) ]"<<std::endl;
        std::cout<<"   [ contact(TIP: "<<contact_tip<<", SIZE: "<<contact_size<<") -  bound(TIP: "<<tip_bound<<", SIZE: "<<exterior_bound<<")]"<<std::endl;


        std::cout<<"   REFINE BOUNDARY ]; "<<std::endl;
     }

     KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************


  void TriangularMesh2DModeler::RemoveCloseNodes(ModelPart& rModelPart,
						 MeshingVariables& rMeshingVariables,
						 ModelPart::IndexType MeshId)
  {
    KRATOS_TRY

    if( this->GetEchoLevel() > 0 ){
      std::cout<<" [ REMOVE CLOSE NODES: "<<std::endl;
      //std::cout<<"   Nodes before erasing : "<<rModelPart.Nodes(MeshId).size()<<std::endl;
    }

    double RemovedConditions = rModelPart.NumberOfConditions(MeshId);
    double NumberOfNodes = rModelPart.NumberOfNodes(MeshId);

    bool any_node_removed = false;
    bool any_condition_removed = false;

    int error_remove = 0;
    int distance_remove = 0;
	
    double inside_nodes_removed   = 0;
    double boundary_nodes_removed = 0;


    //***SIZES :::: parameters do define the tolerance in mesh size: 
    double size_for_criterion_error   = 2.0 * rMeshingVariables.Refine.CriticalRadius; //compared with mean node radius
    double size_for_distance_inside   = 1.0 * rMeshingVariables.Refine.CriticalRadius; //compared with element radius
    double size_for_distance_boundary = 1.5 * size_for_distance_inside; //compared with element radius
    double size_for_wall_tip_contact_side = 0.5* rMeshingVariables.Refine.CriticalSide;
    size_for_wall_tip_contact_side *= 0.3;
    bool derefine_wall_tip_contact = false;

    //if the remove_node switch is activated, we check if the nodes got too close
    if (rMeshingVariables.RefiningOptions.Is(MeshModeler::REMOVE_NODES))
    {

       ////////////////////////////////////////////////////////////
       if (rMeshingVariables.RefiningOptions.Is(MeshModeler::CRITERION_ERROR))	      
       {
          MeshErrorCalculationUtilities MeshErrorDistribution;
          MeshErrorDistribution.SetEchoLevel(this->GetEchoLevel());

          std::vector<double> NodalError;
          std::vector<int>    nodes_ids;


          MeshErrorDistribution.NodalErrorCalculation(rModelPart,NodalError,nodes_ids,MeshId,rMeshingVariables.Refine.GetErrorVariable());

          for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(MeshId); in != rModelPart.NodesEnd(MeshId); in++)
          {

             WeakPointerVector<Node<3> >& rN = in->GetValue(NEIGHBOUR_NODES);
             int erased_nodes =0;
             for(unsigned int i = 0; i < rN.size(); i++)
             {
                if(rN[i].Is(TO_ERASE))
                   erased_nodes += 1;
             }


             if( in->IsNot(BOUNDARY) &&  in->IsNot(STRUCTURE) && erased_nodes < 1 )
             {
                double& MeanError = in->FastGetSolutionStepValue(MEAN_ERROR);
                MeanError = NodalError[nodes_ids[in->Id()]];

                WeakPointerVector<Element >& neighb_elems = in->GetValue(NEIGHBOUR_ELEMENTS);
                double mean_node_radius = 0;
                for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ne++)
                {
                   mean_node_radius+= mModelerUtilities.CalculateTriangleRadius(ne->GetGeometry());
                }

                mean_node_radius /= double(neighb_elems.size());

                if(NodalError[nodes_ids[in->Id()]] < rMeshingVariables.Refine.ReferenceError && mean_node_radius < size_for_criterion_error)
                {
                   //std::cout<<"   Energy : node remove ["<<in->Id()<<"] : "<<NodalError[nodes_ids[in->Id()]]<<std::endl;
                   //std::cout<<"   mean_node_radius "<<mean_node_radius<<" < "<<size_for_criterion_error<<" size_for_criterion_error"<<std::endl;
                   in->Set(TO_ERASE);
                   any_node_removed = true;
                   error_remove++;
                }
             }
          }

       }
       //////////////////////////////////////////////////////////// 


       ////////////////////////////////////////////////////////////
       if (rMeshingVariables.RefiningOptions.Is(MeshModeler::REMOVE_ON_BOUNDARY))
       {
          any_condition_removed = RemoveNonConvexBoundary(rModelPart,rMeshingVariables,MeshId);
       }
       //////////////////////////////////////////////////////////// 


       ////////////////////////////////////////////////////////////
       if (rMeshingVariables.RefiningOptions.Is(MeshModeler::CRITERION_DISTANCE))	      
       {

          //bucket size definition:
          unsigned int bucket_size = 20;


          //create the list of the nodes to be check during the search
          PointPointerVector list_of_nodes;
          list_of_nodes.reserve(rModelPart.NumberOfNodes(MeshId));
          for(ModelPart::NodesContainerType::iterator i_node = rModelPart.NodesBegin(MeshId) ; i_node != rModelPart.NodesEnd(MeshId) ; i_node++)
          {
             (list_of_nodes).push_back(*(i_node.base()));
          }

          KdtreeType nodes_tree(list_of_nodes.begin(),list_of_nodes.end(), bucket_size);

          ////////////////////////////////////////////////////////////

          //all of the nodes in this list will be preserved
          unsigned int num_neighbours = 100;

          PointPointerVector neighbours         (num_neighbours);
          DistanceVector     neighbour_distances(num_neighbours);


          //radius means the distance, if the distance between two nodes is closer to radius -> mark for removing
          double radius=0;
          Node<3> work_point(0,0.0,0.0,0.0);
          unsigned int n_points_in_radius;


          for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(MeshId); in != rModelPart.NodesEnd(MeshId); in++)
          {
             bool on_contact_tip = false;
             array_1d<double, 3 > & ContactForceNormal  = in->FastGetSolutionStepValue(CONTACT_FORCE);

		if(norm_2(ContactForceNormal)>0 || in->Is(TO_SPLIT) || in->Is(CONTACT) )
                on_contact_tip = true;				  

             if( in->IsNot(NEW_ENTITY) )
             {
                //radius=rMeshingVariables.Refine.SizeFactor*in->FastGetSolutionStepValue(NODAL_H);
                radius = size_for_distance_inside;

                work_point[0]=in->X();
                work_point[1]=in->Y();
                work_point[2]=in->Z();

                n_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, neighbours.begin(),neighbour_distances.begin(), num_neighbours);

                if (n_points_in_radius>1)
                {
                   //std::cout<<"     Points in Radius "<< n_points_in_radius<<" radius "<<radius<<std::endl;

                   //if( in->IsNot(STRUCTURE) ) {//MEANS DOFS FIXED

                   if ( in->IsNot(BOUNDARY) )
                   {

                      //look if we are already erasing any of the other nodes
                      unsigned int contact_nodes = 0;
                      unsigned int erased_nodes = 0;
                      for(PointPointerVectorIterator nn=neighbours.begin(); nn!=neighbours.begin() + n_points_in_radius ; nn++)
                      {
                         if( (*nn)->Is(BOUNDARY) && (*nn)->Is(CONTACT) )
                            contact_nodes += 1;

                         if( (*nn)->Is(TO_ERASE) )
                            erased_nodes += 1;
                      }

                      if( erased_nodes < 1 && contact_nodes < 1){ //we release the node if no other nodes neighbours are being erased
                         in->Set(TO_ERASE);
                         //std::cout<<"     Distance Criterion Node ["<<in->Id()<<"] TO_ERASE "<<std::endl;
                         any_node_removed = true;
                         inside_nodes_removed++;
                         //distance_remove++;
                      }

                   }
                   else if ( rMeshingVariables.RefiningOptions.Is(MeshModeler::REMOVE_ON_BOUNDARY) && (in)->IsNot(TO_ERASE)) //boundary nodes will be removed if they get REALLY close to another boundary node (0.2(=extra_factor) * h_factor)
                   {

                      //std::cout<<"  Remove close boundary nodes: Candidate ["<<in->Id()<<"]"<<std::endl;

                      //here we loop over the neighbouring nodes and if there are nodes
                      //with BOUNDARY flag and closer than 0.2*nodal_h from our node, we remove the node we are considering
                      unsigned int k = 0;
                      unsigned int counter = 0;
                      for(PointPointerVectorIterator nn=neighbours.begin(); nn!=neighbours.begin() + n_points_in_radius ; nn++)
                      {
                         bool nn_on_contact_tip = false;
                         array_1d<double, 3 > & ContactForceNormal  = (*nn)->FastGetSolutionStepValue(CONTACT_FORCE);

                         if(norm_2(ContactForceNormal)>0 || (*nn)->Is(TO_SPLIT) || (*nn)->Is(CONTACT) )
                            nn_on_contact_tip = true;				  

                         //std::cout<<" radius * extra_factor "<<(extra_factor*radius)<<" >? "<<neighbour_distances[k]<<std::endl;
                         if ( (*nn)->Is(BOUNDARY) && !nn_on_contact_tip && neighbour_distances[k] < size_for_distance_boundary && neighbour_distances[k] > 0.0 )
                         {
                            //KRATOS_WATCH( neighbours_distances[k] )
                            if((*nn)->IsNot(TO_ERASE)){
                               counter += 1;
                            }
                         }

                         if ( (*nn)->Is(BOUNDARY) && nn_on_contact_tip && neighbour_distances[k] < size_for_wall_tip_contact_side ) {
                            if ( (*nn)->IsNot(TO_ERASE)) { 
                               counter += 1;
                            }
                         }


                         k++;
                      }

                      if(counter > 1 && in->IsNot(NEW_ENTITY) && !on_contact_tip ){ //Can be inserted in the boundary refine
                         in->Set(TO_ERASE);
                         //std::cout<<"     Removed Boundary Node ["<<in->Id()<<"] on Distance "<<std::endl;
                         any_node_removed = true;
                         boundary_nodes_removed++;
                         //distance_remove ++;
                      }
                      else if ( counter > 2 && in->IsNot(NEW_ENTITY) && on_contact_tip && derefine_wall_tip_contact) {
                         in->Set(TO_ERASE);
                         std::cout << "     Removing a TIP POINT due to that criterion [" << in->Id() << "]" << std::endl;
                         any_node_removed = true;
                         boundary_nodes_removed++;
                      }

                   }

                   //}

                }

             }	
          }

          //Build boundary after removing boundary nodes due distance criterion
          if(boundary_nodes_removed){


             std::vector<std::vector<Condition::Pointer> > node_shared_conditions(rModelPart.NumberOfNodes()+1); //all domain nodes

             for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(MeshId); ic!= rModelPart.ConditionsEnd(MeshId); ic++)
             {	 
                if(ic->IsNot(NEW_ENTITY) && ic->IsNot(TO_ERASE)){
                   Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();
                   for(unsigned int i=0; i<rConditionGeom.size(); i++){
                      //std::cout<<"["<<ic->Id()<<"] i "<<i<<" condition "<<rConditionGeom[i].Id()<<std::endl;
                      if(rConditionGeom[i].Is(TO_ERASE)){
                         if( this->GetEchoLevel() > 0 )
                            std::cout<<"     Released node condition ["<<rConditionGeom[i].Id()<<"]: WARNING "<<std::endl;
                      }

                      node_shared_conditions[rConditionGeom[i].Id()].push_back(*(ic.base()));	  
                   }
                }

             }


             //nodes
             int i=0,j=0;
             unsigned int initial_cond_size = rModelPart.Conditions().size()+1; //total model part node size
             unsigned int id = 1;
             unsigned int new_id = 0;

             for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(MeshId); in != rModelPart.NodesEnd(MeshId); in++) 
             {

                if( in->Is(BOUNDARY) && in->IsNot(BLOCKED) && in->IsNot(NEW_ENTITY) && in->Is(TO_ERASE) ){

                   unsigned int nodeId = in->Id();

                   if(node_shared_conditions[nodeId].size()>=2){

                      // std::cout<<"     nodeId "<<nodeId<<std::endl;
                      if(node_shared_conditions[nodeId][0]->IsNot(TO_ERASE) && node_shared_conditions[nodeId][1]->IsNot(TO_ERASE)){

                         if(node_shared_conditions[nodeId][0]->GetGeometry()[0].Id() == in->Id()){
                            i = 1;
                            j = 0;
                         }
                         else{
                            i = 0;
                            j = 1;
                         }


                         Geometry< Node<3> >& rConditionGeom1 = node_shared_conditions[nodeId][i]->GetGeometry();
                         Geometry< Node<3> >& rConditionGeom2 = node_shared_conditions[nodeId][j]->GetGeometry();

                         //node in id Node1;

                         Node<3> & Node0 = rConditionGeom1[0]; // other node in condition [1]
                         Node<3> & Node2 = rConditionGeom2[1]; // other node in condition [2]

                         node_shared_conditions[nodeId][i]->Set(TO_ERASE); //release condition [1]
                         node_shared_conditions[nodeId][j]->Set(TO_ERASE); //release condition [2]

                         any_condition_removed = true;

                         Condition::Pointer NewCond = node_shared_conditions[nodeId][i];

                         Node0.Set(BLOCKED);
                         Node0.Set(MeshModeler::ENGAGED_NODES);

                         Node2.Set(BLOCKED);
                         Node2.Set(MeshModeler::ENGAGED_NODES);

                         //create new condition Node0-NodeB
                         Condition::NodesArrayType face;
                         face.reserve(2);

                         face.push_back(rConditionGeom1(0));
                         face.push_back(rConditionGeom2(1));

                         new_id = initial_cond_size + id;
                         //properties to be used in the generation
                         Condition::Pointer pcond       = NewCond->Clone(new_id, face);
                         // std::cout<<"     ID"<<id<<" 1s "<<pcond1->GetGeometry()[0].Id()<<" "<<pcond1->GetGeometry()[1].Id()<<std::endl;

                         pcond->Set(NEW_ENTITY);

                         //std::cout<<"     Condition INSERTED (Id: "<<new_id<<") ["<<rConditionGeom1[0].Id()<<", "<<rConditionGeom2[1].Id()<<"] "<<std::endl;

                         pcond->SetValue(NORMAL, NewCond->GetValue(NORMAL) );

                         pcond->SetValue(MASTER_NODES, NewCond->GetValue(MASTER_NODES) );
                         pcond->SetValue(CAUCHY_STRESS_VECTOR, NewCond->GetValue(CAUCHY_STRESS_VECTOR));
                         pcond->SetValue(DEFORMATION_GRADIENT, NewCond->GetValue(DEFORMATION_GRADIENT));

                         (rModelPart.Conditions(MeshId)).push_back(pcond);

                         id +=1;
                      }

                   }
                }

             }

             for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(MeshId); in != rModelPart.NodesEnd(MeshId); in++)
             {
                in->Reset(BLOCKED);
             }

          }
          //Build boundary after removing boundary nodes due distance criterion



       }
       // REMOVE ON DISTANCE
       ////////////////////////////////////////////////////////////


       if(any_node_removed)
          mModelerUtilities.CleanRemovedNodes(rModelPart,MeshId);

       if(any_condition_removed){
          //Clean Conditions
          ModelPart::ConditionsContainerType RemoveConditions;

          //id = 0;
          for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(MeshId); ic!= rModelPart.ConditionsEnd(MeshId); ic++)
          {

             if(ic->IsNot(TO_ERASE)){
                //id+=1;
                RemoveConditions.push_back(*(ic.base()));
                //RemoveConditions.back().SetId(id);
             }
	     else{
	       std::cout<<"   Condition RELEASED:"<<ic->Id()<<std::endl;
	     }
          }

          rModelPart.Conditions(MeshId).swap(RemoveConditions);

       }


    }


    // number of removed nodes:
    rMeshingVariables.RemeshInfo.RemovedNodes = NumberOfNodes - rModelPart.NumberOfNodes(MeshId);
    distance_remove =  inside_nodes_removed + boundary_nodes_removed;

    RemovedConditions -= rModelPart.NumberOfConditions(MeshId);

    if( this->GetEchoLevel() > 0 ){
       std::cout<<"   [ CONDITIONS ( removed : "<<RemovedConditions<<" ) ]"<<std::endl;
       std::cout<<"   [ NODES      ( removed : "<<rMeshingVariables.RemeshInfo.RemovedNodes<<" ) ]"<<std::endl;
       std::cout<<"   [ Error(removed: "<<error_remove<<"); Distance(removed: "<<distance_remove<<"; inside: "<<inside_nodes_removed<<"; boundary: "<<boundary_nodes_removed<<") ]"<<std::endl;


       //std::cout<<"   Nodes after  erasing : "<<rModelPart.Nodes(MeshId).size()<<std::endl;
       std::cout<<"   REMOVE CLOSE NODES ]; "<<std::endl;
    }

    KRATOS_CATCH(" ")
  }

  //*******************************************************************************************
  //*******************************************************************************************

  bool TriangularMesh2DModeler::RemoveNonConvexBoundary(ModelPart& rModelPart,
							MeshingVariables& rMeshingVariables,
							ModelPart::IndexType MeshId)
  {
    KRATOS_TRY

    if( this->GetEchoLevel() > 0 ){
      std::cout<<"   [ REMOVE NON CONVEX BOUNDARY : "<<std::endl;
      //std::cout<<"     Starting Conditions : "<<rModelPart.Conditions(MeshId).size()<<std::endl;
    }

    double RemovedConditions = rModelPart.NumberOfConditions(MeshId);

    //***SIZES :::: parameters do define the tolerance in mesh size: 
    double critical_angle        = -120;
    double size_for_side_normal  =  rMeshingVariables.Refine.CriticalRadius;


    std::vector<std::vector<Condition::Pointer> > node_shared_conditions(rModelPart.NumberOfNodes()+1); //all domain nodes
      
    //std::cout<<"     Shared Conditions Size "<<node_shared_conditions.size()<<std::endl;

    for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(MeshId); ic!= rModelPart.ConditionsEnd(MeshId); ic++)
      {	 
	if(ic->IsNot(NEW_ENTITY) && ic->IsNot(TO_ERASE)){
	  Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();
	  for(unsigned int i=0; i<rConditionGeom.size(); i++){
	    //std::cout<<"["<<ic->Id()<<"] i "<<i<<" condition "<<rConditionGeom[i].Id()<<std::endl;
	    if(rConditionGeom[i].Is(TO_ERASE))
	      std::cout<<"     WARNING: Released node condition "<<std::endl;

	    node_shared_conditions[rConditionGeom[i].Id()].push_back(*(ic.base()));	  
	  }
	}
      }

    //std::cout<<"     Node Shared Conditions (Pair of Condition Nodes) is now set "<<std::endl;

    //angles 
    double condition_angle = 0;

    //vector of the neighbour conditions
    array_1d<double,3> S1;
    array_1d<double,3> S2;
    S1.clear();
    S2.clear();

    //normals of the neighbour conditions
    array_1d<double,3> N1;
    array_1d<double,3> N2;
    N1.clear();
    N2.clear();

    //nodes
    int i=0,j=0;
      
    //condition id and size
    unsigned int initial_cond_size = rModelPart.Conditions().size()+1; //total model part node size
    unsigned int id = 1;
    unsigned int new_id = 0;
    int RemovedNodes =0;

    for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(MeshId); in != rModelPart.NodesEnd(MeshId); in++)
      {

	if( in->Is(BOUNDARY) && in->IsNot(BLOCKED) && in->IsNot(NEW_ENTITY) )
	  {
	    unsigned int nodeId = in->Id();

	    if(node_shared_conditions[nodeId].size()>=2){

	      // std::cout<<"     nodeId "<<nodeId<<std::endl;
	      if(node_shared_conditions[nodeId][0]->IsNot(TO_ERASE) && node_shared_conditions[nodeId][1]->IsNot(TO_ERASE)){
		
		if(node_shared_conditions[nodeId][0]->GetGeometry()[0].Id() == in->Id()){
		  i = 1;
		  j = 0;
		}
		else{
		  i = 0;
		  j = 1;
		}
	      
		  
		//Node1*  neighbour conditions in 2D:   (Node0) ---[1]--- (Node1*)

		//normal condition [1]
		N1 = node_shared_conditions[nodeId][i]->GetValue(NORMAL);
		//normal condition [2]
		N2 = node_shared_conditions[nodeId][j]->GetValue(NORMAL);
	      
		// std::cout<<"     N1 "<<N1<<std::endl;
		// std::cout<<"     N2 "<<N2<<std::endl;
		

		Geometry< Node<3> >& rConditionGeom1 = node_shared_conditions[nodeId][i]->GetGeometry();
		Geometry< Node<3> >& rConditionGeom2 = node_shared_conditions[nodeId][j]->GetGeometry();
	      
		//node in id Node1;

		Node<3> & Node0 = rConditionGeom1[0]; // other node in condition [1]
		Node<3> & Node2 = rConditionGeom2[1]; // other node in condition [2]


		// std::cout<<"     Node0: "<<rConditionGeom1[0].Id()<<" Node 1: "<<rConditionGeom1[1].Id()<<std::endl;
		// std::cout<<"     Node1: "<<rConditionGeom2[0].Id()<<" Node 2: "<<rConditionGeom2[1].Id()<<std::endl;
		//segment condition [1]
		S1[0] = rConditionGeom1[1].X() - rConditionGeom1[0].X();
		S1[1] = rConditionGeom1[1].Y() - rConditionGeom1[0].Y();
	      
		if(norm_2(S1)!=0)
		  S1/=norm_2(S1);

		//segment condition [2]
		S2[0] = rConditionGeom2[1].X() - rConditionGeom2[0].X();
		S2[1] = rConditionGeom2[1].Y() - rConditionGeom2[0].Y();

		if(norm_2(S2)!=0)
		  S2/=norm_2(S2);
		  
		// std::cout<<"     S1 "<<S1<<std::endl;
		// std::cout<<"     S2 "<<S2<<std::endl;


		bool remove_S1 = false;
		if(norm_2(S1)<size_for_side_normal)
		  remove_S1 = true;

		bool remove_S2 = false;
		if(norm_2(S2)<size_for_side_normal)
		  remove_S2 = true;

		if(remove_S1 || remove_S2){
		    
		  node_shared_conditions[nodeId][i]->Set(TO_ERASE); //release condition [1]
		  node_shared_conditions[nodeId][j]->Set(TO_ERASE); //release condition [2]
		  in->Set(TO_ERASE);    //release Node1*

		  Condition::Pointer NewCond = node_shared_conditions[nodeId][i];
		    
		  Node0.Set(BLOCKED);
		  Node2.Set(BLOCKED);

		  //create new condition Node0-NodeB
		  Condition::NodesArrayType face;
		  face.reserve(2);

		  face.push_back(rConditionGeom1(0));
		  face.push_back(rConditionGeom2(1));
		
		  new_id = initial_cond_size + id;
		  //properties to be used in the generation
		  Condition::Pointer pcond       = NewCond->Clone(new_id, face);
		  // std::cout<<"     ID"<<id<<" 1s "<<pcond1->GetGeometry()[0].Id()<<" "<<pcond1->GetGeometry()[1].Id()<<std::endl;

		  pcond->Set(NEW_ENTITY);

		  //std::cout<<"     Condition INSERTED (Id: "<<new_id<<") ["<<rConditionGeom1[0].Id()<<", "<<rConditionGeom2[1].Id()<<"] "<<std::endl;

		  pcond->SetValue(NORMAL, NewCond->GetValue(NORMAL) );

		  pcond->SetValue(MASTER_NODES, NewCond->GetValue(MASTER_NODES) );
		  pcond->SetValue(CAUCHY_STRESS_VECTOR, NewCond->GetValue(CAUCHY_STRESS_VECTOR));
		  pcond->SetValue(DEFORMATION_GRADIENT, NewCond->GetValue(DEFORMATION_GRADIENT));

		  (rModelPart.Conditions(MeshId)).push_back(pcond);

		  RemovedNodes += 1;
		  id +=1;

		   
		}
		else{

		  double projection_sides   = inner_prod(S1,S2);
		  double projection_normals = inner_prod(N1,N2);
		  double relative_angle = 0;

		  if(projection_normals!=0)
		    relative_angle = projection_sides/projection_normals;
		  
		  if(relative_angle<=1 && relative_angle>=-1 )
		    condition_angle = (180.0/3.14159) * std::acos(relative_angle);
	    
		  if(inner_prod(S1,N2)<0) 
		    condition_angle *=(-1);

		   // std::cout<<"     projection_sides "<<projection_sides<<std::endl;
		   // std::cout<<"     projection_normals "<<projection_normals<<std::endl;
		   // std::cout<<"     relative_angle "<<relative_angle<<std::endl;
		   // std::cout<<"     condition_angle "<<condition_angle<<" critical_angle "<<critical_angle<<std::endl;
		  

		  if( condition_angle < -40 ){		    
		    // std::cout<<"     B NODE "<<in->Id()<<std::endl;
		    // std::cout<<"     projection_sides "<<projection_sides<<std::endl;
		    // std::cout<<"     projection_normals "<<projection_normals<<std::endl;
		    // std::cout<<"     relative_angle "<<relative_angle<<std::endl;
		    // std::cout<<"     condition_angle "<<condition_angle<<" critical_angle "<<critical_angle<<std::endl;
		    in->Set(VISITED);

		    Node0.Set(VISITED);
		    Node2.Set(VISITED);
		  
		  }

		  if(condition_angle<critical_angle){
		

		    //Path of neighbour conditions in 2D:   (NodeA) ---[0]--- (Node0) ---[1]--- (Node1*) ---[2]--- (Node2) ---[2]--- (Node2) ---[3]--- (NodeB)
		    
		    //realease positions:
		    node_shared_conditions[nodeId][i]->Set(TO_ERASE); //release condition [1]
		    node_shared_conditions[nodeId][j]->Set(TO_ERASE); //release condition [2]

		    in->Set(TO_ERASE);    //release Node1*
		    Node2.Set(TO_ERASE);  //release Node2

		    if( this->GetEchoLevel() > 0 ){
		      std::cout<<"     Node Release/Modify  i "<<in->Id()<<std::endl;
		      std::cout<<"     Node Release/Modify  j "<<Node2.Id()<<std::endl;
		    }

		    //set Node0 to a new position (between 0 and 2)
		    Node0.X() = 0.5 * ( Node0.X() + Node2.X() );
		    Node0.Y() = 0.5 * ( Node0.Y() + Node2.Y() );
		    Node0.Z() = 0.5 * ( Node0.Z() + Node2.Z() );

		    //assign data to dofs
		    unsigned int buffer_size = Node0.GetBufferSize();
		    unsigned int step_data_size = rModelPart.GetNodalSolutionStepDataSize();

		    for(unsigned int step = 0; step<buffer_size; step++)
		      {
			//getting the data of the solution step
			double* step_data = Node0.SolutionStepData().Data(step);

			double* node0_data = Node0.SolutionStepData().Data(step);
			double* node1_data = Node0.SolutionStepData().Data(step);

			//copying this data in the position of the vector we are interested in
			for(unsigned int j= 0; j<step_data_size; j++)
			  {
			    step_data[j] = 0.5*node0_data[j] + 0.5*node1_data[j];
			  }
		      }
			
		    //recover the original position of the node
		    const array_1d<double,3>& disp = Node0.FastGetSolutionStepValue(DISPLACEMENT);
		    Node0.X0() = Node0.X() - disp[0];
		    Node0.Y0() = Node0.Y() - disp[1];
		    Node0.Z0() = Node0.Z() - disp[2];
	
		    //search shared condition of Node0 and Node A
		    if(node_shared_conditions[Node0.Id()][0]->Id() == Node0.Id()){
		      i = 1;
		    }
		    else{
		      i = 0;
		    }
		
		    Geometry< Node<3> >& rConditionGeom0 = node_shared_conditions[Node0.Id()][i]->GetGeometry();
		    Node<3> & NodeA = rConditionGeom0[0];

		    //search shared condition of Node2 and Node B
		    if(node_shared_conditions[Node2.Id()][0]->Id() == Node2.Id()){
		      i = 0;
		    }
		    else{
		      i = 1;
		    }
		
		    //New conditions profile in 2D:  (NodeA) ---[0]--- (Node0**) ---[3]--- (NodeB)   where (Node0**) is (Node0) in another position

		    Condition::Pointer NewCond = node_shared_conditions[Node2.Id()][i];
		    NewCond->Set(TO_ERASE);
		    Geometry< Node<3> >& rConditionGeom3 = NewCond->GetGeometry();
		    Node<3> & NodeB = rConditionGeom3[1];

		    NodeA.Set(MeshModeler::ENGAGED_NODES);
		    NodeB.Set(MeshModeler::ENGAGED_NODES);

		    Node0.Set(MeshModeler::ENGAGED_NODES);
		
		
		    //create new condition Node0-NodeB
		    Condition::NodesArrayType face;
		    face.reserve(2);

		    face.push_back(rConditionGeom1(0));
		    face.push_back(rConditionGeom3(1));
		
		    new_id = initial_cond_size + id;
		    //properties to be used in the generation
		    Condition::Pointer pcond       = NewCond->Clone(new_id, face);
		    // std::cout<<" ID"<<id<<" 1s "<<pcond1->GetGeometry()[0].Id()<<" "<<pcond1->GetGeometry()[1].Id()<<std::endl;

		    pcond->Set(NEW_ENTITY);

		    if( this->GetEchoLevel() > 0 ){
		      std::cout<<"     Condition INSERTED (Id: "<<new_id<<") ["<<rConditionGeom1[0].Id()<<", "<<rConditionGeom3[1].Id()<<"] "<<std::endl;
		    }

		    rConditionGeom1[0].Set(TO_ERASE,false);  // do not release Node1
		    rConditionGeom3[1].Set(TO_ERASE,false);  // do not release Node2


		    pcond->SetValue(NORMAL, NewCond->GetValue(NORMAL) );

		    pcond->SetValue(MASTER_NODES, NewCond->GetValue(MASTER_NODES) );
		    pcond->SetValue(CAUCHY_STRESS_VECTOR, NewCond->GetValue(CAUCHY_STRESS_VECTOR));
		    pcond->SetValue(DEFORMATION_GRADIENT, NewCond->GetValue(DEFORMATION_GRADIENT));

		    (rModelPart.Conditions(MeshId)).push_back(pcond);

		    RemovedNodes += 1;
		    id +=1;
		
		  }
		}
	      }
	    }
	  }
      }

    for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(MeshId); in != rModelPart.NodesEnd(MeshId); in++)
      {
	in->Reset(BLOCKED);
      }
	      
    
    RemovedConditions = rModelPart.Conditions(MeshId).size() - RemovedConditions;

    if( this->GetEchoLevel() > 0 ){
      std::cout<<"     [ CONDITIONS ( removed : "<<RemovedConditions<<" ) ]"<<std::endl;
      std::cout<<"     [ NODES      ( removed : "<<RemovedNodes<<" ) ]"<<std::endl;
    
      std::cout<<"     Ending   Conditions : "<<rModelPart.Conditions(MeshId).size()<<"  (Removed nodes: "<< RemovedNodes<<" ) "<<std::endl;
      std::cout<<"     REMOVE NON CONVEX BOUNDARY ]; "<<std::endl;
    }

    if(RemovedNodes)
      return true;
    else
      return false;
    
    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //METHODS CALLED AFTER REFINING THE TESSELLATION 
  //*******************************************************************************************

  void TriangularMesh2DModeler::RefineElements(ModelPart& rModelPart,
					       MeshingVariables& rMeshingVariables,
					       struct triangulateio& in,
					       struct triangulateio& out,
					       ModelPart::IndexType MeshId)
  {
    KRATOS_TRY

    if( this->GetEchoLevel() > 0 ){
      std::cout<<" [ SELECT ELEMENTS TO REFINE : "<<std::endl;
      //std::cout<<"   refine selection "<<std::endl;
    }

    //***SIZES :::: parameters do define the tolerance in mesh size: 
    double size_for_inside_elements   = 0.75 * rMeshingVariables.Refine.CriticalRadius;
    double size_for_boundary_elements = 1.50 * rMeshingVariables.Refine.CriticalRadius; 

    double nodal_h_refining_factor     = 0.75;
    double nodal_h_non_refining_factor = 2.00;

    ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

    SpatialBoundingBox RefiningBox (rMeshingVariables.BoundingBox.Center,rMeshingVariables.BoundingBox.Radius,rMeshingVariables.BoundingBox.Velocity);


    if(rMeshingVariables.RefiningOptions.Is(MeshModeler::REFINE_ELEMENTS))
      {

	in.numberoftriangles=rMeshingVariables.Refine.NumberOfElements;

	in.trianglelist     = new int  [in.numberoftriangles * 3];
	in.trianglearealist = new REAL [in.numberoftriangles];

	ModelPart::NodesContainerType::iterator nodes_begin = rModelPart.NodesBegin(MeshId);

	//PREPARE THE NODAL_H as a variable to control the automatic point insertion
	//**************************************************************************

	if(rMeshingVariables.RefiningOptions.IsNot(MeshModeler::REFINE_BOUNDARY)){

	  for(unsigned int i = 0; i<rModelPart.Nodes(MeshId).size(); i++)
	    {
	      ////Assign a huge NODAL_H to the free surface nodes, so that there no nodes will be added
	      // if ( (nodes_begin + i)->Is(FREE_SURFACE))
	      // {
	      // 	double & nodal_h=(nodes_begin + i)->FastGetSolutionStepValue(NODAL_H);
	      // 	nodal_h*=nodal_h_non_refining_factor;
	      // }

	      //Assign a huge NODAL_H to the Boundary nodes, so that there no nodes will be added
	      if ( (nodes_begin + i)->Is(BOUNDARY))
		{
		  double & nodal_h=(nodes_begin + i)->FastGetSolutionStepValue(NODAL_H);
		  nodal_h*=nodal_h_non_refining_factor;
		}


	    }

	}

	//SET THE REFINED ELEMENTS AND THE AREA (NODAL_H)
	//*********************************************************************

	int id = 0; 
	    
	for(int el = 0; el< out.numberoftriangles; el++)
	  {
	    if(rMeshingVariables.PreservedElements[el])
	      {

		double prescribed_h      = 0;
		bool   dissipative       = false;
		bool   refine_size       = false;

		int    count_dissipative = 0;
		int    count_boundary_inserted = 0;
		int    count_boundary = 0;
		int    count_contact_boundary = 0;

		Geometry<Node<3> > vertices;

		for(int pn=0; pn<3; pn++)
		  {
		    in.trianglelist[id*3+pn]= out.trianglelist[el*3+pn];
		      
		    vertices.push_back(*(nodes_begin + out.trianglelist[el*3+pn]-1).base());

		    prescribed_h += (nodes_begin + out.trianglelist[el*3+pn]-1)->FastGetSolutionStepValue(NODAL_H);
		      
		    if((nodes_begin + out.trianglelist[el*3+pn]-1)->Is(TO_REFINE))
		      count_dissipative+=1;

		    if((nodes_begin + out.trianglelist[el*3+pn]-1)->Is(NEW_ENTITY))
		      count_boundary_inserted+=1;	      

		    if((nodes_begin + out.trianglelist[el*3+pn]-1)->Is(BOUNDARY)){
		      count_boundary+=1;
		      array_1d<double, 3 > & ContactForceNormal = (nodes_begin + out.trianglelist[el*3+pn]-1)->FastGetSolutionStepValue(CONTACT_FORCE);
		      if( norm_2(ContactForceNormal) )
			count_contact_boundary+=1;
		    }
		    
		  }
		    

		// if(count_dissipative>0)
		//   std::cout<<" Count REFINE nodes "<<count_dissipative<<std::endl;
		// std::cout<<"   prescribed_h (el:"<<el<<") = "<<prescribed_h<<std::endl;
		
		bool refine_candidate = true;
		if (rMeshingVariables.BoundingBox.IsSetFlag == true ){
		  refine_candidate = mModelerUtilities.CheckVerticesInBox(vertices,RefiningBox,CurrentProcessInfo);
		}
		  

		  		  
		double element_area = 0;
		double element_radius = mModelerUtilities.CalculateTriangleRadius (vertices[0].X(),vertices[0].Y(),
										   vertices[1].X(),vertices[1].Y(),
										   vertices[2].X(),vertices[2].Y(),
										   element_area);
		  
	
		//calculate the prescribed h
		prescribed_h *= 0.3333;

		double h = rMeshingVariables.AlphaParameter * prescribed_h;
		//if h is the height of a equilateral triangle, the area is sqrt(3)*h*h/4
		double element_ideal_radius = sqrt(3.0) * 0.25 * ( h * h );
		
		//std::cout<<"   prescribed_h (el:"<<el<<") = "<<prescribed_h<<std::endl;

		if( refine_candidate ){


		  //********* PLASTIC POWER ENERGY REFINEMENT CRITERION (A)
		  if(count_dissipative>=2){
		    dissipative = true;
		    //Set Critical Elements
		    rMeshingVariables.RemeshInfo.CriticalElements += 1;
		    //std::cout<<" Dissipative True "<<std::endl;
		  }


		  //********* SIZE REFINEMENT CRITERION (B)

		  // if(dissipative)
		  //   std::cout<<" element_radius "<<element_radius<<" CriticalRadius "<<size_for_inside_elements<<std::endl;

		  double critical_size = size_for_inside_elements;
		  if( count_boundary >= 2 )
		    critical_size = size_for_boundary_elements;


		  if( element_radius > critical_size ){
		    refine_size = true;
		  }

		    
		  //Also a criteria for the CriticalDissipation (set in nodes)
		  if(rMeshingVariables.RefiningOptions.Is(MeshModeler::CRITERION_ENERGY))
		    {
		      // if( dissipative )
		      // 	std::cout<<" [ Refine Criteria ["<<id<<"] : (dissipative: "<<dissipative<<"; size: "<<refine_size<<"; prescribed h: "<<prescribed_h<<") ]"<<std::endl;


		      //********* PLASTIC POWER ENERGY REFINEMENT CRITERION (A)
		      if( (dissipative == true && refine_size == true) )
			{
			  in.trianglearealist[id] = nodal_h_refining_factor * element_area;

			  //std::cout<<" Area Factor Refine DISSIPATIVE :"<<in.trianglearealist[id]<<std::endl;
			}
		      else if (refine_size == true)
			{


			  in.trianglearealist[id] = element_ideal_radius;

			  if( count_boundary_inserted && count_contact_boundary){
			    in.trianglearealist[id] = nodal_h_refining_factor * element_ideal_radius;
			    //std::cout<<" count boundary inserted-contact on "<<std::endl;
			  }


			  //std::cout<<" Area Factor Refine SIZE :"<<in.trianglearealist[id]<<std::endl;

			}
		      else{
			
			  in.trianglearealist[id] = nodal_h_non_refining_factor * element_area;
			  //std::cout<<" Area Factor Refine NO :"<<in.trianglearealist[id]<<std::endl;
		      }


		      //std::cout<<" [ AREA: "<<Area<<"; NEW AREA: "<<in.trianglearealist[id]<<" ]"<<std::endl;
		    }
		  else
		    {
					    
		      //********* SIZE REFINEMENT CRITERION (B)
		      if( refine_size == true ){

			in.trianglearealist[id] = nodal_h_refining_factor * element_ideal_radius;
		      }
		      else{

			//in.trianglearealist[id] = element_area;
			in.trianglearealist[id] = element_ideal_radius;
		      }

		    }		    

		  //std::cout<<"   mod_prescribed_h (el:"<<el<<") = "<<prescribed_h<<" [ Triangle Area: "<<in.trianglearealist[id]<<" ]"<<std::endl;

		}
		else{

		  //in.trianglearealist[id] = element_ideal_radius;
		  in.trianglearealist[id] = nodal_h_non_refining_factor * element_area;

		}


		id += 1;					
		  	  
	      }



	  }



	//RESTORE THE NODAL_H ON BOUNDARY
	//*********************************************************************
	if(rMeshingVariables.RefiningOptions.IsNot(MeshModeler::REFINE_BOUNDARY)){

	  for(unsigned int i = 0; i<rModelPart.Nodes(MeshId).size(); i++)
	    {
	      // //Unassign the NODAL_H of the free surface nodes
	      // if ( (nodes_begin + i)->Is(FREE_SURFACE))
	      // {
	      // 	double & nodal_h=(nodes_begin + i)->FastGetSolutionStepValue(NODAL_H);
	      // 	nodal_h/=nodal_h_non_refining_factor;
	      // }

	      //Unassign the NODAL_H of the Boundary nodes
	      if ( (nodes_begin + i)->Is(BOUNDARY))
		{
		  double & nodal_h=(nodes_begin + i)->FastGetSolutionStepValue(NODAL_H);
		  nodal_h/=nodal_h_non_refining_factor;
		}


	    }
	}
      }

	
    //Check struct "in"
    //WriteTriangles(in);
    //WritePoints(in);

    if( this->GetEchoLevel() > 0 )
      std::cout<<"   SELECT ELEMENTS TO REFINE ]; "<<std::endl;

    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::GenerateNewParticles(ModelPart& rModelPart, 
						     MeshingVariables& rMeshingVariables,
						     struct triangulateio& in,
						     struct triangulateio& out,
						     ModelPart::IndexType MeshId)
  {
    KRATOS_TRY

    if( this->GetEchoLevel() > 0 )
      std::cout<<" [ GENERATE NEW NODES: "<<std::endl;

    //Find out where the new nodes belong to:

    //creating an auxiliary list for the new nodes
    PointPointerVector list_of_new_nodes;
    //std::vector<int> local_ids;

    //node to get the DOFs from
    Node<3>::DofsContainerType& reference_dofs = (rModelPart.NodesBegin(MeshId))->GetDofs();

    double z = 0.0;

    unsigned int initial_node_size = rModelPart.Nodes().size()+1; //total model part node size

    //std::cout<<" LAST ID "<<(*( (rModelPart.NodesBegin(MeshId) +  in.numberofpoints-1).base()))->Id()<<" "<<rMeshingVariables.NodalPreIds[(*( (rModelPart.NodesBegin(MeshId) +  in.numberofpoints-1).base()))->Id()]<<std::endl;

    //if points were added, new nodes must be added to ModelPart
    int j = 0;
    if (out.numberofpoints > in.numberofpoints)
      {
	for(int i = in.numberofpoints; i<out.numberofpoints; i++)
	  {
	    unsigned int id = initial_node_size + j ;
	    int base = i*2;
	    double& x= out.pointlist[base];
	    double& y= out.pointlist[base+1];

	    //std::cout<<" domain node id "<<id<<" local id "<<i+1<<std::endl;
	    //std::cout<<" node creation position ("<<x<<", "<<y<<")"<<std::endl;
	    Node<3>::Pointer pnode = rModelPart.CreateNewNode(id,x,y,z);
		
	    //set to the main mesh (Mesh 0) to avoid problems in the NodalPreIds (number of nodes: change) in other methods
	    pnode->SetBufferSize(rModelPart.NodesBegin(MeshId)->GetBufferSize() );

	    list_of_new_nodes.push_back( pnode );

	    if(rMeshingVariables.NodalIdsSetFlag){
	      rMeshingVariables.NodalPreIds.push_back( pnode->Id() );
	      pnode->SetId(i+1);
	    }
	    //local_ids.push_back(i+1);

	    //set to the main mesh (Mesh 0) to avoid problems in the NodalPreIds (number of nodes: change) in other methods
	    if(MeshId!=0)
	      (rModelPart.Nodes(MeshId)).push_back(pnode);

	    //generating the dofs
	    for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
	      {
		Node<3>::DofType& rDof = *iii;
		Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );

		(p_new_dof)->FreeDof();
	
	      }
		
	    j++;

	  }

      }


    //Inserted nodes
    rMeshingVariables.RemeshInfo.InsertedNodes = out.numberofpoints-in.numberofpoints;

    if( this->GetEchoLevel() > 0 )
      std::cout <<"   [ GENERATED NODES: ( added: " << rMeshingVariables.RemeshInfo.InsertedNodes <<" ) ]"<<std::endl;
		
    //Set new NodalPreIds in the rMeshingVariables.NodalPreIds
    // j=0;
    // if(rMeshingVariables.NodalIdsSetFlag){
    //   for( PointPointerVectorIterator it =  list_of_new_nodes.begin(); it!=list_of_new_nodes.end(); it++)
    //     {
    //       //set NodalPreIds
    //       rMeshingVariables.NodalPreIds.push_back((*it)->Id());
    //       //recover local ids
    //       (*it)->SetId(local_ids[j]);
    //       //std::cout<<" Add rMeshingVariables.NodalPreIds [ "<<(*it)->Id()<<" ]"<<std::endl;
    //       j++;
    //     }
    // }

    array_1d<double,3> N;
    array_1d<double,3> x1,x2,x3,xc;


    int point_base;
    unsigned int       MaximumNumberOfResults = list_of_new_nodes.size();
    PointPointerVector Results            (MaximumNumberOfResults);
    DistanceVector     ResultsDistances   (MaximumNumberOfResults);

    int step_data_size = rModelPart.GetNodalSolutionStepDataSize();

    //if points were added
    if(out.numberofpoints-in.numberofpoints > 0)
      {

	unsigned int   bucket_size = 20;
	Node<3> work_point(0,0.0,0.0,0.0);

	ModelPart::NodesContainerType::iterator nodes_begin = rModelPart.NodesBegin(MeshId);
	KdtreeType  nodes_tree(list_of_new_nodes.begin(),list_of_new_nodes.end(),bucket_size);

	for(int el = 0; el< in.numberoftriangles; el++)
	  {
	    int base = el * 3;
	    //coordinates
	    point_base = (in.trianglelist[base] - 1)*2;
	    x1[0] = in.pointlist[point_base];
	    x1[1] = in.pointlist[point_base+1];

	    point_base = (in.trianglelist[base+1] - 1)*2;
	    x2[0] = in.pointlist[point_base];
	    x2[1] = in.pointlist[point_base+1];

	    point_base = (in.trianglelist[base+2] - 1)*2;
	    x3[0] = in.pointlist[point_base];
	    x3[1] = in.pointlist[point_base+1];

	    //find the center and "radius" of the element
	    double xc,  yc, radius;
	    this->mDataTransferUtilities.CalculateCenterAndSearchRadius( x1[0], x1[1],
								   x2[0], x2[1],
								   x3[0], x3[1],
								   xc,yc,radius);

	    //find all of the new nodes within the radius
	    work_point.X() = xc;
	    work_point.Y() = yc;
	    work_point.Z() = 0.0;

	    int number_of_points_in_radius = nodes_tree.SearchInRadius (work_point, radius*1.01, Results.begin(), ResultsDistances.begin(),  MaximumNumberOfResults);

	    Triangle2D3<Node<3> > geom(*( (nodes_begin +  in.trianglelist[base]-1).base() 	),
				       *( (nodes_begin +  in.trianglelist[base+1]-1).base() ),
				       *( (nodes_begin +  in.trianglelist[base+2]-1).base() ) );

	    //check if inside and eventually interpolate
	    for( PointPointerVectorIterator it_found = Results.begin(); it_found != Results.begin() + number_of_points_in_radius; it_found++)
	      {
		//if((*it_found)->IsNot(STRUCTURE)){
		bool is_inside = false;
		is_inside = this->mDataTransferUtilities.CalculatePosition( x1[0], x1[1],
								      x2[0], x2[1],
								      x3[0], x3[1],
								      (*it_found)->X(), (*it_found)->Y(), N );


		if(is_inside == true)
		  {
		    // check:
		    // std::cout<<" REFINE TRIANGLE ["<<rMeshingVariables.NodalPreIds[(*( (nodes_begin +  in.trianglelist[base]-1).base()))->Id()]<<", "<<rMeshingVariables.NodalPreIds[(*( (nodes_begin +  in.trianglelist[base+1]-1).base()))->Id()]<<", "<<rMeshingVariables.NodalPreIds[(*( (nodes_begin +  in.trianglelist[base+2]-1).base()))->Id()]<<"] "<<std::endl;
		    // std::cout<<" REFINE PRESSURE ["<<(*( (nodes_begin +  in.trianglelist[base]-1).base()))->FastGetSolutionStepValue(PRESSURE)<<", "<<(*( (nodes_begin +  in.trianglelist[base+1]-1).base()))->FastGetSolutionStepValue(PRESSURE)<<", "<<(*( (nodes_begin +  in.trianglelist[base+2]-1).base()))->FastGetSolutionStepValue(PRESSURE)<<"] "<<std::endl;

		    // std::cout<<" PRESSURE Prev "<<(*it_found)->FastGetSolutionStepValue(PRESSURE)<<std::endl;

		    this->mDataTransferUtilities.Interpolate( geom, N, step_data_size, *(it_found ) );

		    // std::cout<<" PRESSURE Prev "<<(*it_found)->FastGetSolutionStepValue(PRESSURE)<<std::endl;
		  }
		//}
	      }
	  }
      }

    const array_1d<double,3> ZeroNormal(3,0.0);
    //set the coordinates to the original value
    for( PointPointerVectorIterator it =  list_of_new_nodes.begin(); it!=list_of_new_nodes.end(); it++)
      {
	const array_1d<double,3>& displacement = (*it)->FastGetSolutionStepValue(DISPLACEMENT);
	(*it)->X0() = (*it)->X() - displacement[0];
	(*it)->Y0() = (*it)->Y() - displacement[1];
	(*it)->Z0() = 0.0;

	//correct contact_normal interpolation
	noalias((*it)->GetSolutionStepValue(CONTACT_FORCE)) = ZeroNormal;
	noalias((*it)->GetSolutionStepValue(CONTACT_FORCE)) = ZeroNormal;
		    
	(*it)->SetValue(DOMAIN_LABEL,MeshId);

      }

    if( this->GetEchoLevel() > 0 )
      std::cout<<"   GENERATE NEW NODES ]; "<<std::endl;

    KRATOS_CATCH(" ")
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::WriteTriangles(struct triangulateio& tr)
  { 
    KRATOS_TRY

    for(int el = 0; el< tr.numberoftriangles; el++)
      {
	std::cout<<"   TRIANGLE "<<el<<" : [ _";
	for(int pn=0; pn<3; pn++)
	  {
	    std::cout<<tr.trianglelist[el*3+pn]<<"_";
	  }
	std::cout<<" ]   Area: "<<tr.trianglearealist[el]<<std::endl;
      }   

    KRATOS_CATCH(" ")
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::WritePoints(struct triangulateio& tr)
  { 
    KRATOS_TRY

    int base=0;
    for(int nd = 0; nd< tr.numberofpoints; nd++)
      {
	std::cout<<"   Point "<<nd+1<<" : [ ";
	std::cout<<tr.pointlist[base]<<" "<<tr.pointlist[base+1]<<" ]"<<std::endl;;
	base+=2;
      }

    KRATOS_CATCH(" ")
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::ClearTrianglesList(struct triangulateio& tr)
  {
    KRATOS_TRY

    tr.pointlist                  = (REAL*) NULL;
    tr.pointattributelist         = (REAL*) NULL;
    tr.pointmarkerlist            = (int*) NULL;
    tr.numberofpoints             = 0;
    tr.numberofpointattributes    = 0;

    tr.trianglelist               = (int*) NULL;
    tr.triangleattributelist      = (REAL*) NULL;
    tr.trianglearealist           = (REAL*) NULL;
    tr.neighborlist               = (int*) NULL;
    tr.numberoftriangles          = 0;
    tr.numberofcorners            = 3; //for three node triangles
    tr.numberoftriangleattributes = 0;

    tr.segmentlist                = (int*) NULL;
    tr.segmentmarkerlist          = (int*) NULL;
    tr.numberofsegments           = 0;

    tr.holelist                   = (REAL*) NULL;
    tr.numberofholes              = 0;

    tr.regionlist                 = (REAL*) NULL;
    tr.numberofregions            = 0;

    tr.edgelist                   = (int*) NULL;
    tr.edgemarkerlist             = (int*) NULL;
    tr.normlist                   = (REAL*) NULL;
    tr.numberofedges              = 0;

    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::DeleteTrianglesList(struct triangulateio& tr)
  {
    KRATOS_TRY

    //always for "out":
    trifree (tr.trianglelist);

    delete [] tr.triangleattributelist;
    delete [] tr.trianglearealist;

    //in case of n switch not used
    //delete [] tr.neighborlist;

    //if p is switched then in and out are pointed:(free only once)
    //delete [] tr.segmentlist;
    //delete [] tr.segmentmarkerlist;

    delete [] tr.holelist;
    delete [] tr.regionlist;

    //delete [] tr.edgelist;
    //delete [] tr.edgemarkerlist;

    KRATOS_CATCH(" ")
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::DeletePointsList (struct triangulateio& tr)
  {
    KRATOS_TRY

    delete [] tr.pointlist;
    delete [] tr.pointmarkerlist;
    delete [] tr.pointattributelist;

    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::FreeTrianglesList(struct triangulateio& tr)
  {
    KRATOS_TRY

    if(tr.pointlist != NULL) free(tr.pointlist );
    if(tr.pointattributelist != NULL) free(tr.pointattributelist );
    if(tr.pointmarkerlist != NULL) free(tr.pointmarkerlist   );
     
    if(tr.trianglelist != NULL) free(tr.trianglelist  );
    if(tr.triangleattributelist != NULL) free(tr.triangleattributelist );
    if(tr.trianglearealist != NULL) free(tr.trianglearealist );
    if(tr.neighborlist != NULL) free(tr.neighborlist   );

    if(tr.segmentlist != NULL) free(tr.segmentlist    );
    if(tr.segmentmarkerlist != NULL) free(tr.segmentmarkerlist   );

    if(tr.holelist != NULL) free(tr.holelist      );

    if(tr.regionlist != NULL) free(tr.regionlist  );

    if(tr.edgelist != NULL) free(tr.edgelist   );
    if(tr.edgemarkerlist != NULL) free(tr.edgemarkerlist   );
    if(tr.normlist != NULL) free(tr.normlist  );

    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************


} // Namespace Kratos

