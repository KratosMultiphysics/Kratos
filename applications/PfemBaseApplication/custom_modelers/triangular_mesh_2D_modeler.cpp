//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_modelers/laplacian_smoothing.hpp"
#include "custom_modelers/triangular_mesh_2D_modeler.hpp"

#include "pfem_base_application_variables.h"


namespace Kratos
{
  

  //*******************************************************************************************
  //*******************************************************************************************
  void TriangularMesh2DModeler::PerformTransferOnly(ModelPart& rModelPart,
						    MeshingParametersType& rMeshingVariables,
						    ModelPart::IndexType MeshId)
  {

    KRATOS_TRY

    //configuration Options false: REMESH, REFINE, CONSTRAINED
    //configuration Options true : MESH_SMOOTHING

    //execution Options: SET_DOF, SELECT_ELEMENTS, PASS_ALPHA_SHAPE

    this->StartEcho(rModelPart,"Trigen PFEM Transfer Only",MeshId);

    //*********************************************************************
    struct triangulateio in;
    struct triangulateio out;

	
    rMeshingVariables.NodalIdsSetFlag=false;

    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::SET_DOF);
    ////////////////////////////////////////////////////////////
    SetTriangulationNodes(rModelPart,rMeshingVariables,in,out,MeshId);
    ////////////////////////////////////////////////////////////
    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::SET_DOF);


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
	WeakPointerVector<Element >& rE = (elem_begin+el)->GetValue(NEIGHBOUR_ELEMENTS);

	for(int pn=0; pn<3; pn++){
	  if( (elem_begin+el)->Id() == rE[pn].Id() )
	    in.neighborlist[el*3+pn] = 0;
	  else
	    in.neighborlist[el*3+pn] = rE[pn].Id();
	}
	      
      }
	  
    //*********************************************************************


    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::SELECT_ELEMENTS);
    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::PASS_ALPHA_SHAPE);	  

    ////////////////////////////////////////////////////////////
    BuildMeshElements(rModelPart,rMeshingVariables,in,MeshId);
    ////////////////////////////////////////////////////////////

    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::PASS_ALPHA_SHAPE);
    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::SELECT_ELEMENTS);

    //*********************************************************************


    //*********************************************************************
	  
    ////////////////////////////////////////////////////////////
    this->FinalizeMeshGeneration(rModelPart,rMeshingVariables,MeshId);
    ////////////////////////////////////////////////////////////

    //*********************************************************************

    //free memory
    DeletePointsList(in);
    //delete [] in.trianglelist;

    this->EndEcho(rModelPart,"Trigen PFEM Transfer Only",MeshId);

    KRATOS_CATCH( "" )
  }

  
  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::GenerateDT(ModelPart& rModelPart,
					   MeshingParametersType& rMeshingVariables,
					   ModelPart::IndexType MeshId)
  {

    KRATOS_TRY
 
    //configuration Options false: REFINE, CONSTRAINED
    //configuration Options true : REMESH

    //execution Options: SET_DOF, NEIGHBOURS_SEARCH, SELECT_ELEMENTS, PASS_ALPHA_SHAPE

    //if refine true: REFINE_ELEMENTS, REFINE_ADD_NODES


    this->StartEcho(rModelPart,"Trigen PFEM DT Mesher",MeshId);

    //int step_data_size = rModelPart.GetNodalSolutionStepDataSize();

    if(rMeshingVariables.Options.Is(ModelerUtilities::REFINE)){

      //*********************************************************************
      
      ////////////////////////////////////////////////////////////
      this->InitializeMeshGeneration(rModelPart,rMeshingVariables,MeshId);
      ////////////////////////////////////////////////////////////
      
      //*********************************************************************

    }

    ////////////////////////////////////////////////////////////
    //Creating the containers for the input and output
    struct triangulateio in;
    struct triangulateio out;

    rMeshingVariables.NodalIdsSetFlag=false;

    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::SET_DOF);
    ////////////////////////////////////////////////////////////
    SetTriangulationNodes(rModelPart,rMeshingVariables,in,out,MeshId);
    ////////////////////////////////////////////////////////////
    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::SET_DOF);


    ////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////
    boost::timer auxiliary;

    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::NEIGHBOURS_SEARCH);

    GenerateTriangulation(rMeshingVariables.ExecutionOptions,rMeshingVariables.Refine->RefiningOptions,in,out);

    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::NEIGHBOURS_SEARCH);

    //print out the mesh generation time
    if( this->GetEchoLevel() > 0 )
      std::cout<<" [ MESH GENERATION (TIME = "<<auxiliary.elapsed()<<") ] "<<std::endl;
    ////////////////////////////////////////////////////////////

    if(rMeshingVariables.Options.Is(ModelerUtilities::REFINE)){
	  
      ////////////////////////////////////////////////////////////
      //Select Elements to be preserved after passing Alpha-Shape
      rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::PASS_ALPHA_SHAPE);
      rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::SELECT_ELEMENTS);
	  	  
      SelectMeshElements(rModelPart.Nodes(MeshId),rMeshingVariables,out); //passing alpha shape and returning the Elements preserved  

      rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::SELECT_ELEMENTS);
      rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::PASS_ALPHA_SHAPE);
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
	  
      ////////////////////////////////////////////////////////////
      RefineElements (rModelPart,rMeshingVariables,in,out,MeshId);
      ////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////

      //free the memory used in the first step, free out
      ClearTrianglesList(out);

      ////////////////////////////////////////////////////////////
      rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::REFINE);

      GenerateTriangulation(rMeshingVariables.ExecutionOptions,rMeshingVariables.Refine->RefiningOptions,in,out);

      rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::REFINE);
      ////////////////////////////////////////////////////////////

      //Building the entities for new nodes:
      GenerateNewParticles(rModelPart,rMeshingVariables,in,out,MeshId);

      ////////////////////////////////////////////////////////////
    }
    else{

      rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::SELECT_ELEMENTS);
      rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::PASS_ALPHA_SHAPE);	  
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

    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::PASS_ALPHA_SHAPE);
    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::SELECT_ELEMENTS);

    //*********************************************************************
	  
    ////////////////////////////////////////////////////////////
    this->FinalizeMeshGeneration(rModelPart,rMeshingVariables,MeshId);
    ////////////////////////////////////////////////////////////

    //*********************************************************************

    //sort elements
    //rModelPart.Elements().Sort();

    //free memory
    DeletePointsList(in);
    delete [] in.trianglelist;
    DeleteTrianglesList(out);

    this->EndEcho(rModelPart,"Trigen PFEM DT Mesher",MeshId);

    KRATOS_CATCH( "" )
  }



  //*******************************************************************************************
  //*******************************************************************************************
  void TriangularMesh2DModeler::GenerateCDT(ModelPart& rModelPart,
					    MeshingParametersType& rMeshingVariables,
					    ModelPart::IndexType MeshId)
  {
    KRATOS_TRY

    //configuration Options false: REFINE
    //configuration Options true : REMESH, CONSTRAINED

    //execution Options: BOUNDARIES_SEARCH(temp), NEIGHBOURS_SEARCH(temp), SELECT_ELEMENTS, PASS_ALPHA_SHAPE


    this->StartEcho(rModelPart,"Trigen PFEM CDT Mesher",MeshId);

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

    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::BOUNDARIES_SEARCH);

    GenerateTriangulation(rMeshingVariables.ExecutionOptions,rMeshingVariables.Refine->RefiningOptions,in, mid);

    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::BOUNDARIES_SEARCH);

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

    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::CONSTRAINED);
    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::NEIGHBOURS_SEARCH);

    int fail = GenerateTriangulation(rMeshingVariables.ExecutionOptions,rMeshingVariables.Refine->RefiningOptions,mid, out);

    if(fail){
      rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::CONSTRAINED);
      fail = GenerateTriangulation(rMeshingVariables.ExecutionOptions,rMeshingVariables.Refine->RefiningOptions,mid, out);
      rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::CONSTRAINED);
    }
    
    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::NEIGHBOURS_SEARCH);
    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::CONSTRAINED);

    // KRATOS_WATCH( out.numberofsegments )
    // KRATOS_WATCH( out.numberofpoints )
    // KRATOS_WATCH( out.numberoftriangles )
    // KRATOS_WATCH( out.numberofholes )
	
    ////////////////////////////////////////////////////////////
    BuildMeshElements(rModelPart,rMeshingVariables,out,MeshId);	  
    ////////////////////////////////////////////////////////////
	

    //*********************************************************************
	  
    ////////////////////////////////////////////////////////////
    this->FinalizeMeshGeneration(rModelPart,rMeshingVariables,MeshId);
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
	
    this->EndEcho(rModelPart,"Trigen PFEM CDT Mesher",MeshId);

    KRATOS_CATCH( "" )
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::GenerateRDT(ModelPart& rModelPart,
					    MeshingParametersType& rMeshingVariables,
					    ModelPart::IndexType MeshId)
  {
    KRATOS_TRY

    //configuration Options false: CONSTRAINED
    //configuration Options true : REMESH, REFINE

    //execution Options: SET_DOF, NEIGHBOURS_SEARCH(temp), SELECT_ELEMENTS, PASS_ALPHA_SHAPE
    //execution Refinine Options: REFINE_ADD_NODES, REFINE_ELEMENTS, CRITERION_ENERGY

    this->StartEcho(rModelPart,"Trigen PFEM RDT Mesher",MeshId);
	    
    //remove_nodes:
    //REMOVE_NODES, CRITERION_ERROR, REMOVE_ON_BOUNDARY, CRITERION_DISTANCE, ENGAGED_NODES(nodes)

    //refine_boundary:
    //REFINE_INSERT_NODES, REFINE_BOUNDARY // REBUILD_BOUNDARY(CONSTRAINED), CRITERION_ENERGY


    //*********************************************************************
    
    ////////////////////////////////////////////////////////////
    this->InitializeMeshGeneration(rModelPart,rMeshingVariables,MeshId);
    ////////////////////////////////////////////////////////////
    
    //*********************************************************************
    

    ////////////////////////////////////////////////////////////
    //Creating the containers for the input and output
    struct triangulateio in;
    struct triangulateio out;

    rMeshingVariables.NodalIdsSetFlag=false;

    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::SET_DOF);
    ////////////////////////////////////////////////////////////
    SetTriangulationNodes(rModelPart,rMeshingVariables,in,out,MeshId);
    ////////////////////////////////////////////////////////////
    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::SET_DOF);


    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    boost::timer auxiliary;

    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::NEIGHBOURS_SEARCH);
    int fail = GenerateTriangulation(rMeshingVariables.ExecutionOptions,rMeshingVariables.Refine->RefiningOptions,in,out);

    if(fail){
      std::cout<<" Mesher Failed first RDT "<<std::endl;
    }

    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::NEIGHBOURS_SEARCH);

    if(in.numberofpoints!=out.numberofpoints){
      std::cout<<" [ MESH GENERATION FAILED: point insertion (initial = "<<in.numberofpoints<<" final = "<<out.numberofpoints<<") ] "<<std::endl;
    }

    //print out the mesh generation time
    if( this->GetEchoLevel() > 0 )
      std::cout<<" [ MESH GENERATION (TIME = "<<auxiliary.elapsed()<<") ] "<<std::endl;
    ////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////
    //Select Elements to be preserved after passing Alpha-Shape
    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::PASS_ALPHA_SHAPE);
    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::SELECT_ELEMENTS);
		
    SelectMeshElements(rModelPart.Nodes(MeshId),rMeshingVariables,out); //passing alpha shape and returning the Elements preserved

    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::SELECT_ELEMENTS);
    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::PASS_ALPHA_SHAPE);
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
    RefineElements (rModelPart,rMeshingVariables,in,out,MeshId);
    ////////////////////////////////////////////////////////////


    //free the memory used in the first step, free out
    ClearTrianglesList(out);

    ////////////////////////////////////////////////////////////   
    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::REFINE);

    fail = GenerateTriangulation(rMeshingVariables.ExecutionOptions,rMeshingVariables.Refine->RefiningOptions,in,out);

    if(fail){
      std::cout<<" Mesher Failed second RDT "<<std::endl;
    }

    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::REFINE);
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
    this->FinalizeMeshGeneration(rModelPart,rMeshingVariables,MeshId);
    ////////////////////////////////////////////////////////////

    //*********************************************************************

    //sort elements
    //rModelPart.Elements().Sort();

    //free memory
    DeletePointsList(in);
    delete [] in.trianglelist;
    DeleteTrianglesList(out);

    this->EndEcho(rModelPart,"Trigen PFEM RDT Mesher",MeshId);

    KRATOS_CATCH( "" )
  }


  //*******************************************************************************************
  //*******************************************************************************************
  
  void TriangularMesh2DModeler::GenerateRCDT(ModelPart& rModelPart,
					     MeshingParametersType& rMeshingVariables,
					     ModelPart::IndexType MeshId)
  {
    KRATOS_TRY

    //configuration Options false: 
    //configuration Options true : REMESH, REFINE, CONSTRAINED

    //execution Options: SET_DOF, NEIGHBOURS_SEARCH(temp), SELECT_ELEMENTS, PASS_ALPHA_SHAPE, ENGAGED_NODES
    //execution Refinine Options: REFINE_ADD_NODES, REFINE_ELEMENTS, CRITERION_ENERGY, REFINE_BOUNDARY


    this->StartEcho(rModelPart,"Trigen PFEM RCDT Mesher",MeshId);


    //*********************************************************************
    
    ////////////////////////////////////////////////////////////
    this->InitializeMeshGeneration(rModelPart,rMeshingVariables,MeshId);
    ////////////////////////////////////////////////////////////
    
    //*********************************************************************
    

    ////////////////////////////////////////////////////////////
    //Creating the containers for the input and output
    struct triangulateio in;
    struct triangulateio out;

    rMeshingVariables.NodalIdsSetFlag=false;

    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::CONSTRAINED);
    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::SET_DOF);
    ////////////////////////////////////////////////////////////
    SetTriangulationNodes(rModelPart,rMeshingVariables,in,out,MeshId); 
    ////////////////////////////////////////////////////////////
    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::SET_DOF);
    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::CONSTRAINED);

    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    boost::timer auxiliary;

    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::CONSTRAINED);
    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::NEIGHBOURS_SEARCH);

    int fail = GenerateTriangulation(rMeshingVariables.ExecutionOptions,rMeshingVariables.Refine->RefiningOptions,in,out);

    if(fail){
      rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::CONSTRAINED);
      fail = GenerateTriangulation(rMeshingVariables.ExecutionOptions,rMeshingVariables.Refine->RefiningOptions,in, out);
      rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::CONSTRAINED);
    }

    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::NEIGHBOURS_SEARCH);
    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::CONSTRAINED);

    if(in.numberofpoints!=out.numberofpoints){
      std::cout<<" [ MESH GENERATION FAILED: point insertion (initial = "<<in.numberofpoints<<" final = "<<out.numberofpoints<<") ] "<<std::endl;
    }

    if( rMeshingVariables.Options.Is(ModelerUtilities::CONSTRAINED) )
      RecoverBoundaryPosition(rModelPart,rMeshingVariables,in,out,MeshId);

    //print out the mesh generation time
    if( this->GetEchoLevel() > 0 )
      std::cout<<" [ MESH GENERATION (TIME = "<<auxiliary.elapsed()<<") ] "<<std::endl;
    ////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////
    //Select Elements to be preserved after passing Alpha-Shape
    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::PASS_ALPHA_SHAPE);
    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::SELECT_ELEMENTS);
    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::ENGAGED_NODES);

    SelectMeshElements(rModelPart.Nodes(MeshId),rMeshingVariables,out); //passing alpha shape and returning the Elements preserved

    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::ENGAGED_NODES);
    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::SELECT_ELEMENTS);
    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::PASS_ALPHA_SHAPE);
    ////////////////////////////////////////////////////////////

    //free the memory used in the first step, preserve out
    DeletePointsList(in);
    DeleteTrianglesList(in);

    ////////////////////////////////////////////////////////////
    //PERFORM ADAPTIVE REMESHING:
    //1.- Select Triangles to Refine via changing the nodal_h
    struct triangulateio mid_out;

    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::CONSTRAINED);

    rMeshingVariables.NodalIdsSetFlag=true; //second set must be true
    ////////////////////////////////////////////////////////////
    SetTriangulationNodes (rModelPart,rMeshingVariables,in,mid_out,MeshId);
    ////////////////////////////////////////////////////////////

    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::CONSTRAINED);

    ////////////////////////////////////////////////////////////
    RefineElements (rModelPart,rMeshingVariables,in,out,MeshId);
    ////////////////////////////////////////////////////////////


    //free the memory used in the first step, free out
    ClearTrianglesList(out);

    ////////////////////////////////////////////////////////////
    //to generate it constrained it is necessary to change the strategy:
    //a. pass a set of nodes to triangulate ok
    //b. pass a set of segments == conditions to apply the constraint ok
    //c. pass a set of holes if domains are not totally convex. ok
	
    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::CONSTRAINED);
    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::NEIGHBOURS_SEARCH);
    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::REFINE);

    fail = GenerateTriangulation(rMeshingVariables.ExecutionOptions,rMeshingVariables.Refine->RefiningOptions,in,out);

    if(fail){
      rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::CONSTRAINED);
      fail = GenerateTriangulation(rMeshingVariables.ExecutionOptions,rMeshingVariables.Refine->RefiningOptions,in, out);
      rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::CONSTRAINED);
    }

    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::REFINE);
    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::NEIGHBOURS_SEARCH);
    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::CONSTRAINED);

    if( rMeshingVariables.Options.Is(ModelerUtilities::CONSTRAINED) )
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
    //ModelerUtilities::CheckParticles (rModelPart,MeshId);

    //*********************************************************************
	  
    ////////////////////////////////////////////////////////////
    this->FinalizeMeshGeneration(rModelPart,rMeshingVariables,MeshId);
    ////////////////////////////////////////////////////////////

    //*********************************************************************

    //sort elements
    //rModelPart.Elements().Sort();

    //free memory
    DeletePointsList(in);
    delete [] in.trianglelist;
    DeleteTrianglesList(out);

    this->EndEcho(rModelPart,"Trigen PFEM RCDT Mesher",MeshId);

    KRATOS_CATCH( "" )
  }
    


  //*******************************************************************************************
  //METHODS CALLED BEFORE TESSELLATION
  //*******************************************************************************************

  void TriangularMesh2DModeler::SetTriangulationNodes(ModelPart& rModelPart,
						      MeshingParametersType& rMeshingVariables,
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
	   
	if(rMeshingVariables.ExecutionOptions.Is(ModelerUtilities::CONSTRAINED)){

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
	    
	if(rMeshingVariables.ExecutionOptions.Is(ModelerUtilities::SET_DOF))
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

    if(rMeshingVariables.ExecutionOptions.Is(ModelerUtilities::CONSTRAINED)){

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
							MeshingParametersType& rMeshingVariables,
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

    if(MeshingOptions.Is(ModelerUtilities::REFINE)){

      if(RefiningOptions.Is(ModelerUtilities::REFINE_ADD_NODES)){  //to add_nodes automatically and refine the mesh ("q"-quality mesh and "a"-area constraint switches)
	  
	if( MeshingOptions.Is(ModelerUtilities::CONSTRAINED) ){
	  strcpy (meshing_options, "pYJq1.4arnCQ");  //"YYJaqrn" "YJq1.4arn" "Jq1.4arn"
	  meshing_info = "Adaptive constrained remeshing executed";
	}
	else{
	  strcpy (meshing_options, "YJq1.4arnQ");  //"YYJaqrn" "YJq1.4arn" "Jq1.4arn"
	  meshing_info = "Adaptive remeshing executed";
	}


      }
      else if(RefiningOptions.Is(ModelerUtilities::REFINE_INSERT_NODES)){ //to insert a set of given points and refine the mesh
	  
	if( MeshingOptions.Is(ModelerUtilities::CONSTRAINED) ){
	  strcpy (meshing_options, "rinYYJQ");     // "riYYJQ"; //"riYYJQ" // "riJQ" //"riQ"
	  meshing_info = "Inserted constrained remeshing executed";
	}
	else{ //to insert a set of given points and refine the mesh
	  strcpy (meshing_options, "rinJQ");
	  meshing_info = "Inserted remeshing executed";
	}  
	  
      }
      else{ //refine without adding nodes
	strcpy (meshing_options, "YJrn");
	meshing_info = "Non-Adaptive remeshing executed";
      }
	
    }
    else if(MeshingOptions.Is(ModelerUtilities::RECONNECT)){   //to reconect a set of points only

      if(MeshingOptions.Is(ModelerUtilities::CONSTRAINED)){
	strcpy (meshing_options, "pBYYQ");
	meshing_info = "Constrained Reconnection remeshing executed";
      }
      else{
	strcpy (meshing_options, "QNP");
	meshing_info = "Reconnection remeshing executed";
      }

    }
    else{ //to remesh (it is implicit if there is a call to the mesher)

      if( MeshingOptions.Is(ModelerUtilities::NEIGHBOURS_SEARCH)){
	
	if(MeshingOptions.Is(ModelerUtilities::CONSTRAINED)){ //to mesh constrained delaunay
	  strcpy (meshing_options, "pnBYYQ");
	  meshing_info = "Constrained remeshing executed";
	}
	else{ 

	  if(MeshingOptions.Is(ModelerUtilities::BOUNDARIES_SEARCH)){  //to get conectivities, boundaries and neighbours only
	    strcpy (meshing_options, "ncEBQ");
	    meshing_info = "Boundaries-Neighbours remeshing executed";
	  }
	  else{
	    //to get conectivities and neighbours only
	    strcpy (meshing_options, "PneQ");
	    meshing_info = "Neighbours remeshing executed";
	  }
	}
      }
      
      if(MeshingOptions.Is(ModelerUtilities::BOUNDARIES_SEARCH)){  //to get conectivities and boundaries only
	strcpy (meshing_options, "rcEBQ");
	meshing_info = "Boundaries remeshing executed";
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
    
   

    if(MeshingOptions.IsNot(ModelerUtilities::REFINE) && in.numberofpoints<out.numberofpoints){
      fail=3;
      std::cout<<"  fail error: [NODES ADDED] something is wrong with the geometry "<<std::endl;
    }
    
    if( this->GetEchoLevel() > 0 ){
      std::cout<<"  -( "<<meshing_info<<" )- "<<std::endl;
      std::cout<<"  (out ELEMENTS "<<out.numberoftriangles<<") "<<std::endl;
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
						  MeshingParametersType& rMeshingVariables,
						  struct triangulateio &out,
						  ModelPart::IndexType MeshId)
  {
    KRATOS_TRY
    
    //*******************************************************************
    //clearing elements
    //rModelPart.Elements(MeshId).clear();

    //*******************************************************************
    //selecting elements
    rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::ENGAGED_NODES);
    SelectMeshElements(rModelPart.Nodes(MeshId),rMeshingVariables,out);
    rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::ENGAGED_NODES);

    //*******************************************************************
    //setting new elements
    //(rModelPart.Elements(MeshId)).reserve(rMeshingVariables.Info->NumberOfElements);

		
    //*******************************************************************
    // mine 2016 TIP
    // //All nodes in boundary element change
    // if(rMeshingVariables.AvoidTipElementsFlag){ //is not working correctly some dispositions not considered
    //   if( this->GetEchoLevel() > 0 )
    // 	std::cout<<"[   AVOID TIP ELEMENTS START ]"<<std::endl;

    //   ChangeTipElementsUtilities TipElements;
    //   //TipElements.SwapDiagonals(rModelPart,out,rMeshingVariables.PreservedElements,MeshId);

    //   if( this->GetEchoLevel() > 0 )
    // 	std::cout<<"[   AVOID TIP ELEMENTS END ]"<<std::endl;
    // }
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

	    mpDataTransferUtilities->CalculateCenterAndSearchRadius( vertices[0].X(), vertices[0].Y(),
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
    rMeshingVariables.Info->CheckGeometricalSmoothing();

    //if(rMeshingVariables.smoothing && rMeshingVariables.remesh && rMeshingVariables.Info->GeometricalSmoothingRequired ){
    if( rMeshingVariables.Options.Is(ModelerUtilities::MESH_SMOOTHING) && rMeshingVariables.Info->GeometricalSmoothingRequired ){
      LaplacianSmoothing  MeshGeometricSmoothing(rModelPart);
      MeshGeometricSmoothing.SetEchoLevel(this->GetEchoLevel());
      MeshGeometricSmoothing.ApplyMeshSmoothing(rModelPart,rMeshingVariables.PreservedElements,out,list_of_element_vertices,MeshId);
    }
    //*******************************************************************


    //*******************************************************************
    //6) Pass  rReferenceElement and transfer variables
    mpDataTransferUtilities->TransferData(rModelPart,rReferenceElement,list_of_element_centers,list_of_element_vertices,MeshDataTransferUtilities::ELEMENT_TO_ELEMENT,MeshId);
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
						   MeshingParametersType& rMeshingVariables,
						   struct triangulateio& out)
  {
    KRATOS_TRY

    if( this->GetEchoLevel() > 0 )
      std::cout<<" [ SELECT MESH ELEMENTS: ("<<(out.numberoftriangles)<<") "<<std::endl;

    rMeshingVariables.PreservedElements.clear();
    rMeshingVariables.PreservedElements.resize(out.numberoftriangles);
    std::fill( rMeshingVariables.PreservedElements.begin(), rMeshingVariables.PreservedElements.end(), 0 );
	
    rMeshingVariables.Info->NumberOfElements=0;
    
    bool box_side_element = false;
    bool wrong_added_node = false;
    if(rMeshingVariables.ExecutionOptions.IsNot(ModelerUtilities::SELECT_ELEMENTS))
      {

	for(int el=0; el<out.numberoftriangles; el++)
	  {
	    rMeshingVariables.PreservedElements[el]=1;
	    rMeshingVariables.Info->NumberOfElements+=1;
	  }
      }
    else
      {
	if( this->GetEchoLevel() > 0 )
	  std::cout<<"   Start Element Selection "<<out.numberoftriangles<<std::endl;

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
		  if(rMeshingVariables.ExecutionOptions.IsNot(ModelerUtilities::CONTACT_SEARCH))
		    std::cout<<" ERROR: something is wrong: nodal id < 0 "<<std::endl;
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
	    
	    ModelerUtilities ModelerUtils;
	    if(rMeshingVariables.ExecutionOptions.Is(ModelerUtilities::PASS_ALPHA_SHAPE)){

	      if(rMeshingVariables.ExecutionOptions.Is(ModelerUtilities::CONTACT_SEARCH))
		{
		  accepted=ModelerUtils.ShrankAlphaShape(Alpha,vertices,rMeshingVariables.OffsetFactor,2);
		}
	      else
		{
		  accepted=ModelerUtils.AlphaShape(Alpha,vertices,2);
		}

	    }
	    else{

	      accepted = true;

	    }
	      	   
	    //3.- to control all nodes from the same subdomain (problem, domain is not already set for new inserted particles on mesher)
	    // if(accepted)
	    // {
	    //   std::cout<<" Element passed Alpha Shape "<<std::endl;
	    //     if(rMeshingVariables.Refine->Options.IsNot(ModelerUtilities::CONTACT_SEARCH))
	    //   	accepted=ModelerUtilities::CheckSubdomain(vertices);
	    // }

	    //3.1.-
	    bool self_contact = false;
	    if(rMeshingVariables.ExecutionOptions.Is(ModelerUtilities::CONTACT_SEARCH))
	      self_contact = ModelerUtils.CheckSubdomain(vertices);
	    	    
	    //4.- to control that the element is inside of the domain boundaries
	    if(accepted)
	      {
		if(rMeshingVariables.ExecutionOptions.Is(ModelerUtilities::CONTACT_SEARCH))
		  {
		    accepted=ModelerUtils.CheckOuterCentre(vertices,rMeshingVariables.OffsetFactor, self_contact);
		  }
		else
		  {
		    accepted=ModelerUtils.CheckInnerCentre(vertices);
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

	rMeshingVariables.Info->NumberOfElements=number;

      }

    //std::cout<<"   Number of Preserved Elements "<<rMeshingVariables.Info->NumberOfElements<<std::endl;

    if(rMeshingVariables.ExecutionOptions.Is(ModelerUtilities::ENGAGED_NODES)){

      //check engaged nodes
      for(int el=0; el<out.numberoftriangles; el++)
	{
	  if( rMeshingVariables.PreservedElements[el]){
	    for(int pn=0; pn<3; pn++)
	      {
		//set vertices
		rNodes[out.trianglelist[el*3+pn]].Set(ModelerUtilities::ENGAGED_NODES);
	      }
	  }
	    
	}

      int count_release = 0;
      for(ModelPart::NodesContainerType::iterator i_node = rNodes.begin() ; i_node != rNodes.end() ; i_node++)
	{
	  if( i_node->IsNot(ModelerUtilities::ENGAGED_NODES)  ){
	    i_node->Set(TO_ERASE);
	    if( this->GetEchoLevel() > 0 )
	      std::cout<<" NODE "<<i_node->Id()<<" RELEASE "<<std::endl;
	    if( i_node->IsNot(ModelerUtilities::ENGAGED_NODES) )
	      std::cout<<" ERROR: node "<<i_node->Id()<<" IS BOUNDARY RELEASE "<<std::endl;
	    count_release++;
	  }
	      
	  i_node->Reset(ModelerUtilities::ENGAGED_NODES);
	}
	  
      if( this->GetEchoLevel() > 0 )
	std::cout<<"   NUMBER OF RELEASED NODES "<<count_release<<std::endl;

    }

    if( this->GetEchoLevel() > 0 ){
      std::cout<<"   Generated_Elements :"<<out.numberoftriangles<<std::endl;
      std::cout<<"   Passed_AlphaShape  :"<<rMeshingVariables.Info->NumberOfElements<<std::endl;
      std::cout<<"   SELECT MESH ELEMENTS ]; "<<std::endl;
    }

    KRATOS_CATCH( "" )

  }


  //*******************************************************************************************
  //METHODS CALLED BEFORE REFINING THE TESSELLATION
  //*******************************************************************************************

  //...

  //*******************************************************************************************
  //METHODS CALLED AFTER REFINING THE TESSELLATION 
  //*******************************************************************************************

  void TriangularMesh2DModeler::RefineElements(ModelPart& rModelPart,
					       MeshingParametersType& rMeshingVariables,
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
    double size_for_inside_elements   = 0.75 * rMeshingVariables.Refine->CriticalRadius;
    double size_for_boundary_elements = 1.50 * rMeshingVariables.Refine->CriticalRadius; 

    double nodal_h_refining_factor     = 0.75;
    double nodal_h_non_refining_factor = 2.00;

    ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

    if(rMeshingVariables.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_ELEMENTS)
       && rMeshingVariables.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_ADD_NODES) )
      {

	in.numberoftriangles=rMeshingVariables.Info->NumberOfElements;

	in.trianglelist     = new int  [in.numberoftriangles * 3];
	in.trianglearealist = new REAL [in.numberoftriangles];

	ModelPart::NodesContainerType::iterator nodes_begin = rModelPart.NodesBegin(MeshId);

	//PREPARE THE NODAL_H as a variable to control the automatic point insertion
	//**************************************************************************

	if(rMeshingVariables.Refine->RefiningOptions.IsNot(ModelerUtilities::REFINE_BOUNDARY)){

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
		bool   dissipative       = false;  //dissipative means reference threshold is overwhelmed
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
		ModelerUtilities ModelerUtils;
		if (rMeshingVariables.Refine->RefiningBoxSetFlag == true ){
		  refine_candidate = ModelerUtils.CheckVerticesInBox(vertices,*(rMeshingVariables.Refine->RefiningBox),CurrentProcessInfo);
		}
		  

		  		  
		double element_area = 0;
		double element_radius = ModelerUtils.CalculateTriangleRadius (vertices[0].X(),vertices[0].Y(),
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
		    rMeshingVariables.Info->CriticalElements += 1;
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
		  if(rMeshingVariables.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_ELEMENTS_ON_THRESHOLD)
		     && rMeshingVariables.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_ELEMENTS_ON_DISTANCE))
		    {
		      // if( dissipative )
		      // 	std::cout<<" [ Refine Criteria ["<<id<<"] : (dissipative: "<<dissipative<<"; size: "<<refine_size<<"; prescribed h: "<<prescribed_h<<") ]"<<std::endl;


		      //********* THRESHOLD REFINEMENT CRITERION (A)
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
		  else if(rMeshingVariables.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_ELEMENTS_ON_DISTANCE))
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
		  else{
		    
		    //in.trianglearealist[id] = element_ideal_radius;
		    in.trianglearealist[id] = nodal_h_non_refining_factor * element_area;
		    
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
	if(rMeshingVariables.Refine->RefiningOptions.IsNot(ModelerUtilities::REFINE_BOUNDARY)){

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
						     MeshingParametersType& rMeshingVariables,
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
    rMeshingVariables.Info->InsertedNodes = out.numberofpoints-in.numberofpoints;

    if( this->GetEchoLevel() > 0 )
      std::cout <<"   [ GENERATED NODES: ( added: " << rMeshingVariables.Info->InsertedNodes <<" ) ]"<<std::endl;
		
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
	    mpDataTransferUtilities->CalculateCenterAndSearchRadius( x1[0], x1[1],
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
		is_inside = mpDataTransferUtilities->CalculatePosition( x1[0], x1[1],
								      x2[0], x2[1],
								      x3[0], x3[1],
								      (*it_found)->X(), (*it_found)->Y(), N );


		if(is_inside == true)
		  {
		    // check:
		    // std::cout<<" REFINE TRIANGLE ["<<rMeshingVariables.NodalPreIds[(*( (nodes_begin +  in.trianglelist[base]-1).base()))->Id()]<<", "<<rMeshingVariables.NodalPreIds[(*( (nodes_begin +  in.trianglelist[base+1]-1).base()))->Id()]<<", "<<rMeshingVariables.NodalPreIds[(*( (nodes_begin +  in.trianglelist[base+2]-1).base()))->Id()]<<"] "<<std::endl;
		    // std::cout<<" REFINE PRESSURE ["<<(*( (nodes_begin +  in.trianglelist[base]-1).base()))->FastGetSolutionStepValue(PRESSURE)<<", "<<(*( (nodes_begin +  in.trianglelist[base+1]-1).base()))->FastGetSolutionStepValue(PRESSURE)<<", "<<(*( (nodes_begin +  in.trianglelist[base+2]-1).base()))->FastGetSolutionStepValue(PRESSURE)<<"] "<<std::endl;

		    // std::cout<<" PRESSURE Prev "<<(*it_found)->FastGetSolutionStepValue(PRESSURE)<<std::endl;

		    mpDataTransferUtilities->Interpolate( geom, N, step_data_size, *(it_found ) );

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

