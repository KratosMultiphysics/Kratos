//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_meshers/laplacian_smoothing.hpp"
#include "custom_meshers/tetrahedral_mesh_3D_mesher.hpp"
#include "geometries/tetrahedra_3d_4.h"

#include "delaunay_meshing_application_variables.h"
//#include "custom_processes/print_output_mesh_process.hpp"

namespace Kratos
{

  //*******************************************************************************************
  //*******************************************************************************************

  void TetrahedralMesh3DMesher::Generate(ModelPart& rModelPart,
					  MeshingParametersType& rMeshingVariables)
  {

    KRATOS_TRY

    this->StartEcho(rModelPart,"DELAUNAY Remesh");

    //*********************************************************************

    ////////////////////////////////////////////////////////////
    this->ExecutePreMeshingProcesses();
    ////////////////////////////////////////////////////////////

    //*********************************************************************

    //Creating the containers for the input and output
    tetgenio in;
    tetgenio out;

    BuildInput(rModelPart,rMeshingVariables,in);

    //*********************************************************************

    double begin_time = OpenMPUtils::GetCurrentTime();

    //Generate Mesh
    ////////////////////////////////////////////////////////////
    int fail = GenerateTessellation(rMeshingVariables,in,out);
    ////////////////////////////////////////////////////////////

    // if(fail){
    //   if( rMeshingVariables.ExecutionOptions.Is(MesherUtilities::CONSTRAINED) ){
    // 	rMeshingVariables.ExecutionOptions.Reset(MesherUtilities::CONSTRAINED);
    // 	////////////////////////////////////////////////////////////
    // 	fail = GenerateTessellation(rMeshingVariables,in, out);
    // 	////////////////////////////////////////////////////////////
    // 	rMeshingVariables.ExecutionOptions.Set(MesherUtilities::CONSTRAINED);
    //   }
    // }

    if(fail || in.numberofpoints!=out.numberofpoints){
      std::cout<<" [ MESH GENERATION FAILED: point insertion (initial = "<<in.numberofpoints<<" final = "<<out.numberofpoints<<") ] "<<std::endl;
    }

    //print out the mesh generation time
    if( this->GetEchoLevel() > 0 ){
      double end_time = OpenMPUtils::GetCurrentTime();
      std::cout<<" [ MESH GENERATION (TIME = "<<end_time-begin_time<<") ] "<<std::endl;
    }

    //*********************************************************************

    //GetOutput
    SetToContainer(rMeshingVariables.OutMesh,out);

    //*********************************************************************

    ////////////////////////////////////////////////////////////
    this->ExecutePostMeshingProcesses();
    ////////////////////////////////////////////////////////////


    //*********************************************************************

    //Free input memory or keep it to transfer it for next mesh generation
    if( rMeshingVariables.ExecutionOptions.Is(MesherUtilities::FINALIZE_MESHER_INPUT) ){
      DeleteInContainer(rMeshingVariables.InMesh,in);
      rMeshingVariables.InputInitializedFlag = false;
    }

    //*********************************************************************

    //Free output memory
    if(rMeshingVariables.Options.Is(MesherUtilities::REMESH))
      DeleteOutContainer(rMeshingVariables.OutMesh,out);


    this->EndEcho(rModelPart,"DELAUNAY Remesh");

    KRATOS_CATCH( "" )

  }


  //*******************************************************************************************
  //*******************************************************************************************

  int TetrahedralMesh3DMesher::GenerateTessellation(MeshingParametersType& rMeshingVariables,
						     tetgenio& in,
						     tetgenio& out)
  {
    KRATOS_TRY

    //if not remesh return true
    if(rMeshingVariables.Options.IsNot(MesherUtilities::REMESH)){
      // out.pointlist             = in.pointlist;
      // out.numberofpoints        = in.numberofpoints;
      // out.tetrahedronlist       = in.tetrahedronlist;
      // out.numberoftetrahedra    = in.numberoftetrahedra;
      // out.tetrahedronvolumelist = in.tetrahedronvolumelist;
      // out.neighborlist          = in.neighborlist;
      out = in;
      return 1;
    }

    int fail=0;

    tetgenio vorout;

    if( this->GetEchoLevel() > 0 )
      std::cout<<" [ REMESH: (in POINTS "<<in.numberofpoints<<") "<<std::endl;

    //this->WritePoints(in);

    std::string str = rMeshingVariables.TessellationFlags;
    char *meshing_options = new char[str.length() + 1];
    strcpy(meshing_options, str.c_str());

    //perform the meshing
    try {
      tetrahedralize(meshing_options,&in,&out,&vorout);
    }

    catch( int error_code ){

      switch(TetgenErrors(error_code))
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
	    std::cout<<" tetrahedralization done "<<std::endl;
	  break;
	}
    }

    delete [] meshing_options;

    //this->CheckInOutPoints(in, out);
    //this->WritePoints(out);
    //this->WriteTetrahedra(out);

    if(rMeshingVariables.Options.IsNot(MesherUtilities::REFINE) && in.numberofpoints<out.numberofpoints){
      fail=3;
      std::cout<<"  warning : [NODES ADDED] check mesh geometry "<<std::endl;
    }

    if( this->GetEchoLevel() > 0 ){
      std::cout<<"  -( "<<rMeshingVariables.TessellationInfo<<" )- "<<std::endl;
      std::cout<<"  (out ELEMENTS "<<out.numberoftetrahedra<<") "<<std::endl;
      std::cout<<"  (out POINTS "<<out.numberofpoints<<") :  REMESH ]; "<<std::endl;
      std::cout<<std::endl;
    }

    return fail;

    KRATOS_CATCH( "" )

  }

  //*******************************************************************************************
  //*******************************************************************************************

  void TetrahedralMesh3DMesher::BuildInput(ModelPart& rModelPart,
					    MeshingParametersType& rMeshingVariables,
					    tetgenio& in)

  {
    KRATOS_TRY


    if( rMeshingVariables.ExecutionOptions.Is(MesherUtilities::INITIALIZE_MESHER_INPUT) ){

      //Set Nodes
      if( rMeshingVariables.ExecutionOptions.Is(MesherUtilities::TRANSFER_KRATOS_NODES_TO_MESHER) )
	this->SetNodes(rModelPart,rMeshingVariables);

      //Set Elements
      if( rMeshingVariables.ExecutionOptions.Is(MesherUtilities::TRANSFER_KRATOS_ELEMENTS_TO_MESHER) )
	this->SetElements(rModelPart,rMeshingVariables);

      //Set Neighbours
      if( rMeshingVariables.ExecutionOptions.Is(MesherUtilities::TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER) )
	this->SetNeighbours(rModelPart,rMeshingVariables);

      rMeshingVariables.InputInitializedFlag = true;
    }

    //input mesh: NODES
    in.firstnumber    = 1;
    in.mesh_dim       = 3;
    GetFromContainer(rMeshingVariables.InMesh,in);

    //Set Faces
    if( rMeshingVariables.ExecutionOptions.Is(MesherUtilities::TRANSFER_KRATOS_FACES_TO_MESHER) )
      this->SetFaces(rModelPart,rMeshingVariables, in);

    //std::cout<<" InMesh [nodes:"<<rMeshingVariables.InMesh.GetNumberOfPoints()<<", elements:"<<rMeshingVariables.InMesh.GetNumberOfElements()<<"]"<<std::endl;

    //Print input mesh
    //PrintOutputMeshProcess print_output = PrintOutputMeshProcess(rModelPart,rMeshingVariables,"input",1);
    //print_output.Execute();

    KRATOS_CATCH( "" )
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void TetrahedralMesh3DMesher::SetFaces(ModelPart& rModelPart,
					  MeshingParametersType& rMeshingVariables,
					  tetgenio& in)
  {
    KRATOS_TRY

    //*********************************************************************

    if(in.facetlist){
      delete [] in.facetlist;
      in.numberoffacets = 0;
    }

    if(in.facetmarkerlist){
      delete [] in.facetmarkerlist;
    }

    if(in.holelist){
      delete [] in.holelist;
      in.numberofholes = 0;
    }

    if(in.regionlist){
      delete [] in.regionlist;
      in.numberofregions = 0;
    }


    //PART 2: facet list (we can have holes in facets != volume holes)

    in.numberoffacets           = rModelPart.NumberOfConditions();
    in.facetmarkerlist          = new int[in.numberoffacets];
    in.facetlist                = new tetgenio::facet[in.numberoffacets];


    ModelPart::ConditionsContainerType::iterator conditions_begin = rModelPart.ConditionsBegin();

    //facets
    tetgenio::facet   *f;
    tetgenio::polygon *p;

    for(int fc=0; fc<in.numberoffacets; fc++)
      {
	f = &in.facetlist[fc];

	f->numberofpolygons = 1;
	f->polygonlist      = new tetgenio::polygon[f->numberofpolygons];
	f->numberofholes    = 0;
	f->holelist         = NULL;
	p = &f->polygonlist[0];
	p->numberofvertices = 3; //face is a triangle
	p->vertexlist       = new int[p->numberofvertices];


	if( (conditions_begin + fc)->Is(TO_ERASE) )
	  std::cout<<" ERROR: condition to erase present "<<std::endl;

	Geometry< Node<3> >& rGeometry = (conditions_begin + fc)->GetGeometry();

	for (int nd=0;nd<3;nd++)
	  {
	    p->vertexlist[nd] = rGeometry[nd].Id();
	  }

	in.facetmarkerlist[fc] = 0;

      }

    //PART 3: (volume) hole list

    //holes
    in.numberofholes            = 0;
    in.holelist                 = (REAL*) NULL;

    //PART 4: region attributes list

    //regions
    in.numberofregions          = 1;
    in.regionlist               = new REAL[in.numberofregions * 5];


    double inside_factor = 2;
    Geometry< Node<3> >& rGeometry = (conditions_begin)->GetGeometry();
    array_1d<double, 3>&  Normal   = rGeometry[0].FastGetSolutionStepValue(NORMAL);

    std::cout<<" Normal [NodeId= "<<rGeometry[0].Id()<<"] "<<Normal<<std::endl;

    double NormNormal = norm_2(Normal);
    if( NormNormal != 0)
      Normal /= NormNormal;

    //inside point of the region:
    in.regionlist[0] = rGeometry[0][0]-(Normal[0]*rMeshingVariables.OffsetFactor*inside_factor);
    in.regionlist[1] = rGeometry[0][1]-(Normal[1]*rMeshingVariables.OffsetFactor*inside_factor);
    in.regionlist[2] = rGeometry[0][2]-(Normal[2]*rMeshingVariables.OffsetFactor*inside_factor);

    //region attribute (regional attribute or marker "A" must be switched)
    in.regionlist[3] = 0;

    //region maximum volume attribute (maximum volume attribute "a" (with no number following) must be switched)
    in.regionlist[4] = -1;

    std::cout<<" Number of facets "<<in.numberoffacets<<" region ("<<in.regionlist[0]<<", "<<in.regionlist[1]<<", "<<in.regionlist[2]<<") normal:"<<Normal<<" Offset "<<rMeshingVariables.OffsetFactor*inside_factor<<std::endl;


    KRATOS_CATCH( "" )

  }

  //*******************************************************************************************
  //*******************************************************************************************
  void TetrahedralMesh3DMesher::GetFromContainer(MesherUtilities::MeshContainer& rMesh, tetgenio& tr)
  {

    KRATOS_TRY

    //get pointers
    tr.pointlist             = rMesh.GetPointList();
    tr.tetrahedronlist       = rMesh.GetElementList();
    tr.tetrahedronvolumelist = rMesh.GetElementSizeList();
    tr.neighborlist          = rMesh.GetElementNeighbourList();

    if( rMesh.GetNumberOfPoints() != 0 )
      tr.numberofpoints = rMesh.GetNumberOfPoints();

    if( rMesh.GetNumberOfElements() != 0 )
      tr.numberoftetrahedra = rMesh.GetNumberOfElements();

    KRATOS_CATCH( "" )
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void TetrahedralMesh3DMesher::SetToContainer(MesherUtilities::MeshContainer& rMesh, tetgenio& tr)
  {

    KRATOS_TRY

    //set pointers
    rMesh.SetPointList(tr.pointlist);
    rMesh.SetElementList(tr.tetrahedronlist);
    rMesh.SetElementSizeList(tr.tetrahedronvolumelist);
    rMesh.SetElementNeighbourList(tr.neighborlist);

    // copy the numbers
    if( tr.numberofpoints != 0 ){
      rMesh.SetNumberOfPoints(tr.numberofpoints);
    }

    if( tr.numberoftetrahedra != 0 ){
      rMesh.SetNumberOfElements(tr.numberoftetrahedra);
    }

    KRATOS_CATCH( "" )

  }

  //*******************************************************************************************
  //*******************************************************************************************

  void TetrahedralMesh3DMesher::DeleteInContainer(MesherUtilities::MeshContainer& rMesh, tetgenio& tr)
  {
    KRATOS_TRY

    //delete mesher container
    rMesh.Finalize();
    ClearTetgenIO(tr); // blocks tetgen automatic destructor deletetion of a NULL pointer []


    KRATOS_CATCH( "" )
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void TetrahedralMesh3DMesher::DeleteOutContainer(MesherUtilities::MeshContainer& rMesh, tetgenio& tr)
  {
    KRATOS_TRY

    //delete mesher container
    // rMesh.Finalize();
    // ClearTetgenIO(tr); // blocks tetgen automatic destructor deletetion of a NULL pointer []
    tr.deinitialize();
    tr.initialize();

    KRATOS_CATCH( "" )
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void TetrahedralMesh3DMesher::WriteTetrahedra(tetgenio& tr)
  {
    KRATOS_TRY

    std::cout<<" Write Tetrahedra "<<std::endl;
    for(int el = 0; el< tr.numberoftetrahedra; el++)
      {
	std::cout<<"   TETRAHEDRON "<<el<<" : [ _";
	for(int pn=0; pn<4; pn++)
	  {
	    std::cout<<tr.tetrahedronlist[el*4+pn]<<"_";
	  }
	std::cout<<" ] "<<std::endl; //  Volume: "<<tr.tetrahedronvolumelist[el]<<std::endl;
      }

    KRATOS_CATCH(" ")
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void TetrahedralMesh3DMesher::WritePoints(tetgenio& tr)
  {
    KRATOS_TRY

    int base=0;
    std::cout<<" numberofpoints "<<tr.numberofpoints<<" dimension "<<tr.mesh_dim<<std::endl;
    for(int nd = 0; nd< tr.numberofpoints; nd++)
      {
	std::cout<<"   Point "<<nd+1<<" : [ ";
	std::cout<<tr.pointlist[base]<<" "<<tr.pointlist[base+1]<<" "<<tr.pointlist[base+2]<<" ]"<<std::endl;;
	base+=3;
      }

    KRATOS_CATCH(" ")
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void TetrahedralMesh3DMesher::CheckInOutPoints(tetgenio& in,tetgenio& out)
  {
    KRATOS_TRY

    if( in.numberofpoints != out.numberofpoints )
      std::cout<<"  Input and Output points amount is not the same : [in:"<<in.numberofpoints<<",out:"<<out.numberofpoints<<"]"<<std::endl;

    int base=0;
    bool coincide = true;
    for(int nd = 0; nd< in.numberofpoints; nd++)
      {
	// std::cout<<"   Point "<<nd+1<<" : [ ";
	// std::cout<<in.pointlist[base]<<" "<<in.pointlist[base+1]<<" "<<in.pointlist[base+2]<<" ]"<<std::endl;
	// std::cout<<" ["<<out.pointlist[base]<<" "<<out.pointlist[base+1]<<" "<<out.pointlist[base+2]<<" ]"<<std::endl;

	if( fabs(in.pointlist[base]-out.pointlist[base])>1e-8 )
	  coincide = false;
	if( fabs(in.pointlist[base+1]-out.pointlist[base+1])>1e-8 )
	  coincide = false;
	if( fabs(in.pointlist[base+2]-out.pointlist[base+2])>1e-8 )
	  coincide = false;

	base+=3;
      }

    if(coincide == false)
      std::cout<<"  Input and Output points not coincide "<<std::endl;

    KRATOS_CATCH(" ")
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void TetrahedralMesh3DMesher::ClearTetgenIO(tetgenio& tr)
  {
    KRATOS_TRY

    tr.pointlist                     = (REAL*) NULL;
    tr.numberofpoints                = 0;
    tr.numberofpointattributes       = 0;

    tr.tetrahedronlist               = (int*) NULL;
    tr.tetrahedronvolumelist         = (REAL*) NULL;
    tr.neighborlist                  = (int*) NULL;
    tr.numberoftetrahedra            = 0;

    KRATOS_CATCH(" ")
  }



} // Namespace Kratos

