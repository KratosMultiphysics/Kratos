//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_modelers/contact_domain_3D_modeler.hpp"

#include "contact_mechanics_application_variables.h"

namespace Kratos
{

  //*******************************************************************************************
  //*******************************************************************************************

  void ContactDomain3DModeler::SetNodes(ModelPart& rModelPart,
					MeshingParametersType& rMeshingVariables)
  {

    KRATOS_TRY

    unsigned int& MeshId = rMeshingVariables.MeshId;

    const unsigned int dimension = rModelPart.ConditionsBegin(MeshId)->GetGeometry().WorkingSpaceDimension();

    //*********************************************************************

    //writing the points coordinates in a vector and reordening the Id's from an initial id
    ModelPart::NodesContainerType::iterator nodes_begin = rModelPart.NodesBegin();
	
    const array_1d<double,3> ZeroOffset(3,0.0);

    double nodal_h_min = std::numeric_limits<double>::max();
    double nodal_h     = 0;

    for(unsigned int i = 0; i<rModelPart.Nodes().size(); i++)
      {
	noalias((nodes_begin + i)->FastGetSolutionStepValue(OFFSET)) = ZeroOffset;

	nodal_h = (nodes_begin + i)->FastGetSolutionStepValue(NODAL_H); 
	if(nodal_h_min>nodal_h)  
	  nodal_h_min=nodal_h;
      }


    double hnodal_offset_conversion = 0.35;
    if( rMeshingVariables.OffsetFactor > nodal_h_min*hnodal_offset_conversion || rMeshingVariables.OffsetFactor < nodal_h_min*0.01){
      rMeshingVariables.OffsetFactor = nodal_h_min*hnodal_offset_conversion;
    }

    std::cout<<"   Minimum Nodal_h "<<nodal_h_min<<" OffsetFactor "<<rMeshingVariables.OffsetFactor<<std::endl;

    std::vector<PointType> BoxVertices;
    BoxVertices.resize(0);
    //*********************************************************************
    if(rMeshingVariables.Options.Is(ModelerUtilities::CONSTRAINED)){
      
      std::cout<<"   Constrained Contact Meshing "<<std::endl;
      
      //PART 1: node list
      double extra_radius = rMeshingVariables.OffsetFactor*4; 
      SpatialBoundingBox DomainBox (rModelPart, extra_radius);
      
      ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
      BoxVertices = DomainBox.GetVertices( CurrentProcessInfo[TIME], dimension );
    }

    //input mesh: NODES
    ModelerUtilities::MeshContainer& InMesh = rMeshingVariables.InMesh;
    
    InMesh.CreatePointList(rModelPart.NumberOfNodes() + BoxVertices.size(), dimension);
    
    double* PointList     = InMesh.GetPointList();
    int& NumberOfPoints   = InMesh.GetNumberOfPoints();
    
    if(!rMeshingVariables.InputInitializedFlag){
     
      rMeshingVariables.NodeMaxId = 0; 
      if((int)rMeshingVariables.NodalPreIds.size() != NumberOfPoints+1)
	rMeshingVariables.NodalPreIds.resize(NumberOfPoints+1);
      
      std::fill( rMeshingVariables.NodalPreIds.begin(), rMeshingVariables.NodalPreIds.end(), 0 );
    }
    
    //writing the points coordinates in a vector and reordening the Id's   
    int base   = 0;
    int direct = 1;
    
    double Shrink = 0;
    array_1d<double, 3> Offset;
    
    for(unsigned int i = 0; i<rModelPart.NumberOfNodes(); i++)
      {
	//from now on it is consecutive
	if(!rMeshingVariables.InputInitializedFlag){
	  rMeshingVariables.NodalPreIds[direct]=(nodes_begin + i)->Id();
	  (nodes_begin + i)->SetId(direct);
	  if( rMeshingVariables.NodalPreIds[direct] > (int)rMeshingVariables.NodeMaxId )
	    rMeshingVariables.NodeMaxId = rMeshingVariables.NodalPreIds[direct];

	}
	
	array_1d<double, 3>& Coordinates = (nodes_begin + i)->Coordinates();
	array_1d<double, 3>& Normal      = (nodes_begin + i)->FastGetSolutionStepValue(NORMAL); //BOUNDARY_NORMAL must be set as nodal variable
	Shrink = (nodes_begin + i)->FastGetSolutionStepValue(SHRINK_FACTOR);   //SHRINK_FACTOR   must be set as nodal variable	          
	
	Normal /= norm_2(Normal);
	for(unsigned int j=0; j<dimension; j++){
	  Offset[j] = ( (-1) * Normal[j] * Shrink * rMeshingVariables.OffsetFactor);
	}
	
	for(unsigned int j=0; j<dimension; j++){
	  PointList[base+j]   = Coordinates[j] + Offset[j];
	}
	
	base+=dimension;
	direct+=1;
      }

    if(BoxVertices.size() !=0 ){
	    
      std::vector<int> vertices_ids;
      for(unsigned int i = 0; i<BoxVertices.size(); i++)
	{
	  rMeshingVariables.NodalPreIds[direct] = -1;
	  vertices_ids.push_back(direct);

	  for(unsigned int j=0; j<dimension; j++){
	    PointList[base+j]   = BoxVertices[i][j];
	  }
 
	  //std::cout<<"   BoxVertices ["<<i<<"]= ("<<BoxVertices[i][0]<<", "<<BoxVertices[i][1]<<"). Id = "<<vertices_ids[i]<<std::endl;

	  base+=dimension; 
	  direct+=1;
	}
    }
     
    //InMesh.SetPointList(PointList);

    KRATOS_CATCH( "" )
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void ContactDomain3DModeler::SetFaces(ModelPart& rModelPart,
					MeshingParametersType& rMeshingVariables,
					tetgenio& in)
  {
    KRATOS_TRY
      
    unsigned int& MeshId = rMeshingVariables.MeshId;

    const unsigned int dimension = rModelPart.ConditionsBegin(MeshId)->GetGeometry().WorkingSpaceDimension();

    //*********************************************************************
    if( in.facetlist != NULL ){
      delete [] in.facetlist;
      in.numberoffacets = 0;
    }

    if( in.facetmarkerlist != NULL )
      delete [] in.facetmarkerlist;
    
    if( in.holelist != NULL ){
      delete [] in.holelist;
      in.numberofholes = 0;
    }
    
    if( in.regionlist != NULL ){
      delete [] in.regionlist;
      in.numberofregions = 0;
    }
    
    //PART 2: facet list (we can have holes in facets != volume holes)
    //*********************************************************************
    std::vector<PointType> BoxVertices;
    BoxVertices.resize(0);
    double extra_radius = rMeshingVariables.OffsetFactor*4; 
    SpatialBoundingBox DomainBox (rModelPart, extra_radius);
    
    ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
    BoxVertices = DomainBox.GetVertices( CurrentProcessInfo[TIME], dimension );
    
    
    in.numberoffacets       = rModelPart.NumberOfConditions() + BoxVertices.size();
    in.facetmarkerlist      = new int[in.numberoffacets];
    in.facetlist            = new tetgenio::facet[in.numberoffacets];
    
    ModelPart::ConditionsContainerType::iterator conditions_begin = rModelPart.ConditionsBegin();

    //facets
    tetgenio::facet   *f;
    tetgenio::polygon *p;
    
    for(unsigned int fc=0; fc<rModelPart.Conditions(MeshId).size(); fc++)
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

	in.facetmarkerlist[fc] = MeshId+1;
	
      }

    //BoundaryBox facets
    ModelerUtilities::MeshContainer& InMesh = rMeshingVariables.InMesh;
    int& NumberOfPoints   = InMesh.GetNumberOfPoints();

    int ids = NumberOfPoints - BoxVertices.size() + 1;
    for(int fc=rModelPart.Conditions(MeshId).size(); fc<in.numberoffacets; fc++) //3d (prismatic box of 6 sides)
      {
	f = &in.facetlist[fc];

	f->numberofpolygons = 1;
	f->polygonlist      = new tetgenio::polygon[f->numberofpolygons];

	f->numberofholes    = 0;
	f->holelist         = NULL;
            
	p = &f->polygonlist[0];

	p->numberofvertices = 4; //faces of a quadrilateral
	p->vertexlist       = new int[p->numberofvertices];

	if( (conditions_begin + fc)->Is(TO_ERASE) )
	  std::cout<<" ERROR: condition to erase present "<<std::endl;

	for (int nd=0;nd<4;nd++)
	  {
	    p->vertexlist[nd] = ids+nd;
	  }      

	in.facetmarkerlist[fc] = MeshId+1;	
      }  

    
    //PART 3: (volume) hole list    
    std::vector<bounded_vector<double, 3> >& Holes = rMeshingVariables.GetHoles();

    //holes
    in.numberofholes              = Holes.size();
    in.holelist                   = new REAL[in.numberofholes * 3];
    
     for(unsigned int hl=0; hl<Holes.size();hl++)
       {
	 //std::cout<<"   BoxHoles ["<<hl<<"]= ("<<Holes[hl][0]<<", "<<Holes[hl][1]<<", "<<Holes[hl][2]<<")"<<std::endl;

	 //inside point of the hole:
	 in.holelist[hl*3+0] = Holes[hl][0];
	 in.holelist[hl*3+1] = Holes[hl][1];
	 in.holelist[hl*3+2] = Holes[hl][2];
       }
     
     
     //PART 4: region attributes list
     in.numberofregions =0;
     in.regionlist =NULL;
     
     KRATOS_CATCH( "" )
  }
  

} // Namespace Kratos

