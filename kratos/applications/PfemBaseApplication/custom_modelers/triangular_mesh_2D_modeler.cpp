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

  void TriangularMesh2DModeler::Generate(ModelPart& rModelPart,
					 MeshingParametersType& rMeshingVariables)
  {

    KRATOS_TRY
 
    unsigned int& MeshId = rMeshingVariables.MeshId;

    this->StartEcho(rModelPart,"PFEM Base Remesh",MeshId);
    
    //*********************************************************************

    ////////////////////////////////////////////////////////////
    this->ExecutePreMeshingProcesses();   
    ////////////////////////////////////////////////////////////

    //*********************************************************************      

    //Creating the containers for the input and output
    struct triangulateio in;
    struct triangulateio out;
    ClearTrianglesList(out);

    BuildInput(rModelPart,rMeshingVariables,in);
    
    //*********************************************************************

    boost::timer auxiliary;

    //Generate Mesh
    ////////////////////////////////////////////////////////////
    int fail = GenerateTessellation(rMeshingVariables,in,out);
    ////////////////////////////////////////////////////////////

    // if(fail){
    //   if( rMeshingVariables.ExecutionOptions.Is(ModelerUtilities::CONSTRAINED) ){	
    // 	rMeshingVariables.ExecutionOptions.Reset(ModelerUtilities::CONSTRAINED);
    // 	////////////////////////////////////////////////////////////
    // 	fail = GenerateTessellation(rMeshingVariables,in, out);
    // 	////////////////////////////////////////////////////////////
    // 	rMeshingVariables.ExecutionOptions.Set(ModelerUtilities::CONSTRAINED);
    //   }
    // }

    if(fail || in.numberofpoints!=out.numberofpoints){
      std::cout<<" [ MESH GENERATION FAILED: point insertion (initial = "<<in.numberofpoints<<" final = "<<out.numberofpoints<<") ] "<<std::endl;
    }

    //Print out the mesh generation time
    if( this->GetEchoLevel() > 0 )
      std::cout<<" [ MESH GENERATION (TIME = "<<auxiliary.elapsed()<<") ] "<<std::endl;


    //*********************************************************************

    //GetOutput
    SetToContainer(rMeshingVariables.OutMesh,out);

    //*********************************************************************
    
    ////////////////////////////////////////////////////////////
    this->ExecutePostMeshingProcesses();
    ////////////////////////////////////////////////////////////
    

    //*********************************************************************

    //Free input memory or keep it to transfer it for next mesh generation
    if( rMeshingVariables.ExecutionOptions.Is(ModelerUtilities::FINALIZE_MESHER_INPUT) ){
      DeleteInContainer(rMeshingVariables.InMesh,in);
      rMeshingVariables.InputInitializedFlag = false;
    }

    //*********************************************************************

    //Free output memory
    if(rMeshingVariables.Options.Is(ModelerUtilities::REMESH))
       DeleteOutContainer(rMeshingVariables.OutMesh,out);

    
    this->EndEcho(rModelPart,"PFEM Base Remesh",MeshId);

    KRATOS_CATCH( "" )

  }



  //*******************************************************************************************
  //*******************************************************************************************

  int TriangularMesh2DModeler::GenerateTessellation(MeshingParametersType& rMeshingVariables,
						    struct triangulateio& in,
						    struct triangulateio& out)
  {
    KRATOS_TRY

    //if not remesh return 0 (no fail)
    if(rMeshingVariables.Options.IsNot(ModelerUtilities::REMESH)){
      // out.pointlist         = in.pointlist;
      // out.numberofpoints    = in.numberofpoints;
      // out.trianglelist      = in.trianglelist;
      // out.numberoftriangles = in.numberoftriangles;
      // out.trianglearealist  = in.trianglearealist;
      // out.neighborlist      = in.neighborlist;
      out = in;
      return 0;
    }

    int fail=0;

    struct triangulateio vorout;

    //initilize all to avoid memory problems
    ClearTrianglesList(vorout);

    //switches: https://www.cs.cmu.edu/~quake/triangle.switch.html

    if( this->GetEchoLevel() > 0 )
      std::cout<<" [ REMESH: (in POINTS "<<in.numberofpoints<<") "<<std::endl;

    //this->WritePoints(in);

    std::string str = rMeshingVariables.TessellationFlags;
    char *meshing_options = new char[str.length() + 1];
    strcpy(meshing_options, str.c_str());

    //perform the meshing
    try {
      triangulate(meshing_options,&in,&out,&vorout);
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
      
    delete [] meshing_options;

    if(rMeshingVariables.Options.IsNot(ModelerUtilities::REFINE) && in.numberofpoints<out.numberofpoints){
      fail=3;
      std::cout<<"  fail error: [NODES ADDED] something is wrong with the geometry "<<std::endl;
    }
    
    if( this->GetEchoLevel() > 0 ){
      std::cout<<"  -( "<<rMeshingVariables.TessellationInfo<<" )- "<<std::endl;
      std::cout<<"  (out ELEMENTS "<<out.numberoftriangles<<") "<<std::endl;
      std::cout<<"  (out POINTS "<<out.numberofpoints<<") :  REMESH ]; "<<std::endl;
      std::cout<<std::endl;
    }

    return fail;
   
    KRATOS_CATCH( "" )

  }


  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::BuildInput(ModelPart& rModelPart,
					   MeshingParametersType& rMeshingVariables,
					   struct triangulateio& in)
					   
  {
    KRATOS_TRY

    if( rMeshingVariables.ExecutionOptions.Is(ModelerUtilities::INITIALIZE_MESHER_INPUT) ){

      //Set Nodes
      if( rMeshingVariables.ExecutionOptions.Is(ModelerUtilities::TRANSFER_KRATOS_NODES_TO_MESHER) )
	this->SetNodes(rModelPart,rMeshingVariables);
      
      //Set Elements
      if( rMeshingVariables.ExecutionOptions.Is(ModelerUtilities::TRANSFER_KRATOS_ELEMENTS_TO_MESHER) )
	this->SetElements(rModelPart,rMeshingVariables);
      
      //Set Neighbours
      if( rMeshingVariables.ExecutionOptions.Is(ModelerUtilities::TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER) )
	this->SetNeighbours(rModelPart,rMeshingVariables);

      rMeshingVariables.InputInitializedFlag = true;
    }
    
    ClearTrianglesList(in);
    GetFromContainer(rMeshingVariables.InMesh,in);

    // if( rMeshingVariables.ExecutionOptions.IsNot(ModelerUtilities::INITIALIZE_MESHER_INPUT) ){
    //   WritePoints(in);
    //   WriteTriangles(in);
    // }

    //Set Faces
    if( rMeshingVariables.ExecutionOptions.Is(ModelerUtilities::TRANSFER_KRATOS_FACES_TO_MESHER) )
      this->SetFaces(rModelPart,rMeshingVariables, in);

   

    KRATOS_CATCH( "" )
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::SetFaces(ModelPart& rModelPart,
					 MeshingParametersType& rMeshingVariables,
					 struct triangulateio& in)
  {
     KRATOS_TRY

     unsigned int& MeshId = rMeshingVariables.MeshId;

     //*********************************************************************

     if( in.segmentlist != NULL ){
       delete [] in.segmentlist;
       in.numberofsegments = 0;
     }

     if( in.segmentmarkerlist != NULL )
       delete [] in.segmentmarkerlist;

     if( in.holelist != NULL ){
      delete [] in.holelist;
      in.numberofholes = 0;
     }

     if( in.regionlist != NULL ){
      delete [] in.regionlist;
      in.numberofregions = 0;
     }


     //PART 2: faced list (we can have holes in facets != area holes)
     in.numberofsegments           = rModelPart.NumberOfConditions(MeshId);
     in.segmentmarkerlist          = new int[in.numberofsegments];
     in.segmentlist                = new int[in.numberofsegments*2];
     
     
      ModelPart::ConditionsContainerType::iterator conditions_begin = rModelPart.ConditionsBegin(MeshId);
      
      
      int base = 0;
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


      KRATOS_CATCH( "" )

  }

  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::GetFromContainer(ModelerUtilities::MeshContainer& rMesh, struct triangulateio& tr)
  {

    KRATOS_TRY
      
    //get pointers
    tr.pointlist        = rMesh.GetPointList();
    tr.trianglelist     = rMesh.GetElementList();
    tr.trianglearealist = rMesh.GetElementSizeList();
    tr.neighborlist     = rMesh.GetElementNeighbourList();      
    
    if( rMesh.GetNumberOfPoints() != 0 )
      tr.numberofpoints = rMesh.GetNumberOfPoints();
      
    if( rMesh.GetNumberOfElements() != 0 )
      tr.numberoftriangles = rMesh.GetNumberOfElements();


    KRATOS_CATCH( "" )
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::SetToContainer(ModelerUtilities::MeshContainer& rMesh, struct triangulateio& tr)
  {

    KRATOS_TRY

    //set pointers
    rMesh.SetPointList(tr.pointlist);
    rMesh.SetElementList(tr.trianglelist);
    rMesh.SetElementSizeList(tr.trianglearealist);
    rMesh.SetElementNeighbourList(tr.neighborlist);

    // copy the numbers 
    if( tr.numberofpoints != 0 ){
      rMesh.SetNumberOfPoints(tr.numberofpoints);
    }

    if( tr.numberoftriangles != 0 ){
      rMesh.SetNumberOfElements(tr.numberoftriangles);
    }
    
    KRATOS_CATCH( "" )

  }

  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::DeleteInContainer(ModelerUtilities::MeshContainer& rMesh, struct triangulateio& tr)
  {
    KRATOS_TRY

    //delete triangle other structures
    //DeleteTrianglesList(tr);
    ClearTrianglesList(tr);

    //delete modeler container
    rMesh.Finalize();    

    KRATOS_CATCH( "" )
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void TriangularMesh2DModeler::DeleteOutContainer(ModelerUtilities::MeshContainer& rMesh, struct triangulateio& tr)
  {
    KRATOS_TRY

    //delete triangle other structures
    DeleteTrianglesList(tr);
    ClearTrianglesList(tr);

    //delete modeler container
    rMesh.Finalize();    

    KRATOS_CATCH( "" )
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
	//std::cout<<" ]   Area: "<<tr.trianglearealist[el]<<std::endl;
	std::cout<<" ] "<<std::endl;
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
    if(tr.numberoftriangles)
      trifree (tr.trianglelist);

    if( tr.triangleattributelist != NULL ) 
      delete [] tr.triangleattributelist;

    if( tr.holelist != NULL )
      delete [] tr.holelist;
    
    // if( tr.regionlist != NULL )
    //   delete [] tr.regionlist;      

    //if p is switched then in and out are pointed:(free only once)
    //delete [] tr.segmentlist;
    //delete [] tr.segmentmarkerlist;

    //in case of n switch not used
    //delete [] tr.neighborlist;
    //delete [] tr.trianglearealist;

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

