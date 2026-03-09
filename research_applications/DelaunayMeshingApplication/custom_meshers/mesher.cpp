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
#include "custom_meshers/mesher.hpp"

#include "delaunay_meshing_application_variables.h"

namespace Kratos
{


  //*******************************************************************************************
  //*******************************************************************************************

  void Mesher::Initialize()
  {
  }


  //*******************************************************************************************
  //*******************************************************************************************


  void Mesher::InitializeMesher( ModelPart& rModelPart )
  {
  }



  //*******************************************************************************************
  //*******************************************************************************************


  void Mesher::FinalizeMesher( ModelPart& rModelPart )
  {
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void Mesher::SetMeshingParameters( MeshingParametersType::Pointer& rMeshingParameters )
  {
    KRATOS_TRY

    mpMeshingVariables = rMeshingParameters;

    bool Remesh = false;
    if( mpMeshingVariables->Options.Is(MesherUtilities::REMESH) )
      Remesh = true;

    bool Refine = false;
    if( mpMeshingVariables->Options.Is(MesherUtilities::REFINE) )
      Refine = true;

    bool Transfer = false;
    if( mpMeshingVariables->Options.Is(MesherUtilities::TRANSFER) )
      Transfer = true;

    if( mEchoLevel > 0 )
      std::cout<<"  SetRemeshData: [ RefineFlag: "<<Refine<<" RemeshFlag: "<<Remesh<<" TransferFlag: "<<Transfer<<" ] "<<std::endl;

    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void Mesher::SetPreMeshingProcess( MesherProcess::Pointer pPreMeshingProcess )
  {
     KRATOS_TRY

     mPreMeshingProcesses.push_back(pPreMeshingProcess); //NOTE: order set = order of execution

     KRATOS_CATCH(" ")
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void Mesher::SetPostMeshingProcess( MesherProcess::Pointer pPostMeshingProcess )
  {
     KRATOS_TRY

     mPostMeshingProcesses.push_back(pPostMeshingProcess); //NOTE: order set = order of execution

     KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void Mesher::SetPreMeshingProcessVector( std::vector<MesherProcess::Pointer>& rPreMeshingProcessVector )
  {
     KRATOS_TRY

     mPreMeshingProcesses = rPreMeshingProcessVector;

     KRATOS_CATCH(" ")
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void Mesher::SetPostMeshingProcessVector( std::vector<MesherProcess::Pointer>& rPostMeshingProcessVector )
  {
     KRATOS_TRY

     mPostMeshingProcesses = rPostMeshingProcessVector;

     KRATOS_CATCH(" ")
  }



  //*******************************************************************************************
  //*******************************************************************************************

  void Mesher::ExecutePreMeshingProcesses()
  {
    KRATOS_TRY

    //Refine and Remove nodes processes
    ////////////////////////////////////////////////////////////
    if( mPreMeshingProcesses.size() )
      for(unsigned int i=0; i<mPreMeshingProcesses.size(); ++i)
	mPreMeshingProcesses[i]->Execute();
    ////////////////////////////////////////////////////////////

    KRATOS_CATCH( "" )
  }



  //*******************************************************************************************
  //*******************************************************************************************

  void Mesher::ExecutePostMeshingProcesses()
  {
    KRATOS_TRY

    //Rebuild Boundary processes
    ////////////////////////////////////////////////////////////
    if( mPostMeshingProcesses.size() )
      for(unsigned int i=0; i<mPostMeshingProcesses.size(); ++i)
	mPostMeshingProcesses[i]->Execute();
    ////////////////////////////////////////////////////////////


    KRATOS_CATCH( "" )

  }


  //*******************************************************************************************
  //*******************************************************************************************

  void Mesher::SetMesherUtilities(MesherUtilities::Pointer rMesherUtilities )
  {
    KRATOS_TRY

    mpMesherUtilities = rMesherUtilities;

    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void Mesher::SetDataTransferUtilities(MeshDataTransferUtilities::Pointer rDataTransferUtilities )
  {
    KRATOS_TRY

    mpDataTransferUtilities = rDataTransferUtilities;

    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************


  void Mesher::ExecuteMeshing ( ModelPart& rModelPart )
  {
    KRATOS_TRY


    if( GetEchoLevel() > 0 ){
      std::cout<<" [ GetRemeshData: [ RefineFlag: "<<mpMeshingVariables->Options.Is(MesherUtilities::REFINE)<<"; RemeshFlag : "<<mpMeshingVariables->Options.Is(MesherUtilities::REMESH)<<" ] ]"<<std::endl;
    }

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

    // Located in the beginning of the assignation:

    // check mesh size introduced :: warning must be shown
    // if(!out_buffer_active)
    //   std::cout.rdbuf(buffer);

    if(mpMeshingVariables->Options.Is( MesherUtilities::REFINE )){
      MesherUtilities MesherUtils;
      MesherUtils.CheckCriticalRadius(rModelPart, mpMeshingVariables->Refine->CriticalRadius);
    }

    // if(!out_buffer_active){
    //   buffer = std::cout.rdbuf();
    //   std::ofstream fout("/dev/null");
    //   std::cout.rdbuf(fout.rdbuf());
    // }
    // check mesh size introduced :: warning must be shown


    //bool remesh_performed=false;

    if( GetEchoLevel() > 0 ){
      std::cout<<" --------------                     -------------- "<<std::endl;
      std::cout<<" --------------       DOMAIN        -------------- "<<std::endl;
    }

    // Located in the beginning of the assignation:

    //generate mesh
    this->Generate(rModelPart,*(mpMeshingVariables));


    KRATOS_CATCH(" ")
  }

  //*******************************************************************************************
  //*******************************************************************************************
  void Mesher::StartEcho(ModelPart& rSubModelPart,
			      std::string GenerationMessage)
  {

    if( this->GetEchoLevel() > 0 ){
      std::cout<<" [ [ [ ] ] ]"<<std::endl;
      std::cout<<" [ "<<GenerationMessage <<" ]"<<std::endl;
      std::cout<<" [ PREVIOUS MESH ["<<rSubModelPart.Name()<<"] (Elements: "<<rSubModelPart.NumberOfElements()<<" Nodes: "<<rSubModelPart.NumberOfNodes()<<" Conditions: "<<rSubModelPart.NumberOfConditions()<<") ]"<<std::endl;
    }
  }

  //*******************************************************************************************
  //*******************************************************************************************
  void Mesher::EndEcho(ModelPart& rSubModelPart,
			    std::string GenerationMessage)
  {

    if( this->GetEchoLevel() > 0 ){
      std::cout<<" [ NEW MESH ["<<rSubModelPart.Name()<<"] (Elements: "<<rSubModelPart.Elements().size()<<" Nodes: "<<rSubModelPart.Nodes().size()<<" Conditions: "<<rSubModelPart.Conditions().size()<<") ]"<<std::endl;
      std::cout<<" [ "<<GenerationMessage <<" ]"<<std::endl;
      std::cout<<" [ Finished Remeshing ] "<<std::endl;
      std::cout<<" [ [ [ ] ] ]"<<std::endl;
    }

  }


  //*******************************************************************************************
  //*******************************************************************************************

  void Mesher::SetNodes(ModelPart& rModelPart,
			     MeshingParametersType& rMeshingVariables)
  {
    KRATOS_TRY

    MesherUtilities MesherUtils;
    MesherUtils.SetNodes(rModelPart,rMeshingVariables);

    KRATOS_CATCH( "" )

  }


  //*******************************************************************************************
  //*******************************************************************************************

  void Mesher::SetElements(ModelPart& rModelPart,
				MeshingParametersType& rMeshingVariables)
  {
    KRATOS_TRY

    MesherUtilities MesherUtils;
    MesherUtils.SetElements(rModelPart,rMeshingVariables);

    KRATOS_CATCH( "" )

  }


  //*******************************************************************************************
  //*******************************************************************************************

  void Mesher::SetNeighbours(ModelPart& rModelPart,
                             MeshingParametersType& rMeshingVariables)
  {
    KRATOS_TRY

    //*********************************************************************
    //input mesh: NEIGHBOURELEMENTS
    ModelPart::ElementsContainerType::iterator element_begin = rModelPart.ElementsBegin();
    const unsigned int nds          = element_begin->GetGeometry().size();

    MesherUtilities::MeshContainer& InMesh = rMeshingVariables.InMesh;

    InMesh.CreateElementNeighbourList(rModelPart.Elements().size(), nds);

    int* ElementNeighbourList      = InMesh.GetElementNeighbourList();

    for(unsigned int el = 0; el<rModelPart.Elements().size(); el++)
      {
	ElementWeakPtrVectorType& nElements = (element_begin+el)->GetValue(NEIGHBOUR_ELEMENTS);

        int counter = 0;
        for(auto& i_nelem : nElements)
        {
	  if( (element_begin+el)->Id() == i_nelem.Id() )
	    ElementNeighbourList[el*nds+counter] = 0;
	  else
	    ElementNeighbourList[el*nds+counter] = i_nelem.Id();
          ++counter;
	}

      }

    KRATOS_CATCH( "" )
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void Mesher::SetElementNeighbours(ModelPart& rModelPart,
                                    MeshingParametersType & rMeshingVariables)
  {
    KRATOS_TRY

    if( this->GetEchoLevel() > 0 ){
      std::cout<<" [ SET ELEMENT NEIGHBORS : "<<std::endl;
      std::cout<<"   Initial Faces : "<<rModelPart.Conditions().size()<<std::endl;
    }

    ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();

    int facecounter=0;
    for(ModelPart::ElementsContainerType::iterator ci_elem = rModelPart.ElementsBegin();
	ci_elem != rModelPart.ElementsEnd(); ++ci_elem)
      {
	int Id= ci_elem->Id() - 1;

	int number_of_faces = ci_elem->GetGeometry().FacesNumber(); //defined for triangles and tetrahedra
	ElementWeakPtrVectorType& nElements = ci_elem->GetValue(NEIGHBOUR_ELEMENTS);
        nElements.resize(number_of_faces);

	for(int i = 0; i<number_of_faces; i++)
	  {
	    int index = rMeshingVariables.NeighbourList[Id][i];

	    if(index > 0)
	      {
		index = rMeshingVariables.PreservedElements[index-1];
	      }

	    if(index > 0)
	      {
		nElements(i) = GlobalPointer<Element>(*(el_begin + index -1).base());
	      }
	    else
	      {
		nElements(i) = GlobalPointer<Element>(*ci_elem.base());
		facecounter++;
	      }
	  }
      }

    if( this->GetEchoLevel() > 0 ){
      std::cout<<"   Final Faces : "<<facecounter<<std::endl;
      std::cout<<"   SET ELEMENT NEIGHBORS ]; "<<std::endl;
    }

    KRATOS_CATCH( "" )

  }


  //*******************************************************************************************
  //*******************************************************************************************

  void Mesher::RecoverBoundaryPosition(ModelPart& rModelPart,
					    MeshingParametersType& rMeshingVariables)
  {
    KRATOS_TRY

    const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

    //*********************************************************************
    //input mesh: ELEMENTS

    MesherUtilities::MeshContainer& InMesh = rMeshingVariables.InMesh;
    double* InPointList  = InMesh.GetPointList();

    MesherUtilities::MeshContainer& OutMesh = rMeshingVariables.OutMesh;
    double* OutPointList = OutMesh.GetPointList();

    ModelPart::NodesContainerType::iterator nodes_begin = rModelPart.NodesBegin();

    int base=0;
    for(unsigned int i = 0; i<rModelPart.Nodes().size(); i++)
      {

	if( (nodes_begin + i)->Is(BOUNDARY) ){

	  array_1d<double, 3>& Position = (nodes_begin + i)->Coordinates();
	  for( unsigned int j=0; j<dimension; j++)
	    {
	      InPointList[base+j]    = Position[j];
	      OutPointList[base+j]   = Position[j];
	    }
	}

	base+=dimension;
      }

    KRATOS_CATCH( "" )
  }


} // Namespace Kratos
