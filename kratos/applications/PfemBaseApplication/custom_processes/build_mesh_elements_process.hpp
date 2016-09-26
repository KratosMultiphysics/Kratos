//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_BUILD_MESH_ELEMENTS_PROCESS_H_INCLUDED )
#define  KRATOS_BUILD_MESH_ELEMENTS_PROCESS_H_INCLUDED


// System includes
#include <string>
#include <iostream>

// External includes
#include <boost/timer.hpp>

// Project includes
#include "includes/model_part.h"
#include "custom_modelers/laplacian_smoothing.hpp"
#include "custom_utilities/modeler_utilities.hpp"

///VARIABLES used:
//Data:     MASTER_ELEMENTS(set), MASTER_NODES(set)
//StepData: 
//Flags:    (checked) TO_ERASE, TO_REFINE, CONTACT, NEW_ENTITY
//          (set)     BOUNDARY(set),  [TO_REFINE(nodes), TO_ERASE(condition)]->locally to not preserve condition
//          (modified)  
//          (reset)   
// (set):=(set in this process)

namespace Kratos
{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}
  ///@name Kratos Classes
  ///@{

  /// Short class definition.
  /** Detail class definition.
   */
  class BuildMeshElementsProcess
    : public Process
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BuildMeshElementsProcess
    KRATOS_CLASS_POINTER_DEFINITION( BuildMeshElementsProcess );

    typedef ModelPart::NodeType                   NodeType;
    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BuildMeshElementsProcess(ModelPart& rModelPart,
			     ModelerUtilities::MeshingParameters& rRemeshingParameters,
			     int EchoLevel)
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    { 
      mMeshId = mrRemesh.MeshId;
      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~BuildMeshElementsProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
      Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    virtual void Execute()
    {
      KRATOS_TRY

      if( mEchoLevel > 0 )
	std::cout<<" [ GENERATE NEW ELEMENTS: "<<std::endl;

      if( mrModelPart.Name() != mrRemesh.SubModelPartName )
	std::cout<<" ModelPart Supplied do not corresponds to the Meshing Domain: ("<<mrModelPart.Name()<<" != "<<mrRemesh.SubModelPartName<<")"<<std::endl;
    
      //*******************************************************************
      //selecting elements
      
	if( !mrRemesh.MeshElementsSelectedFlag )  //Select Mesh Elements not performed  ... is needed to be done before building new elements
	{  
	  std::cout<<" ERROR : no selection of elements performed before building the elements "<<std::endl;
	  SelectMeshElementsProcess SelectElements(mrModelPart,mrRemesh,mEchoLevel);
	  SelectElements.Execute();
	}


      //*******************************************************************
      //setting new elements
      //(mrModelPart.Elements(MeshId)).reserve(mrRemesh.Info->NumberOfElements);

		
      //*******************************************************************
      // mine 2016 TIP
      // //All nodes in boundary element change
      // if(mrRemesh.AvoidTipElementsFlag){ //is not working correctly some dispositions not considered
      //   if( mEchoLevel > 0 )
      // 	std::cout<<"[   AVOID TIP ELEMENTS START ]"<<std::endl;

      //   ChangeTipElementsUtilities TipElements;
      //   //TipElements.SwapDiagonals(mrModelPart,out,mrRemesh.PreservedElements,MeshId);
      
      //   if( mEchoLevel > 0 )
      // 	std::cout<<"[   AVOID TIP ELEMENTS END ]"<<std::endl;
      // }
      //*******************************************************************


      //properties to be used in the generation
      int number_properties = mrModelPart.GetParentModelPart()->NumberOfProperties();
      Properties::Pointer properties = mrModelPart.GetParentModelPart()->GetMesh(mMeshId).pGetProperties(number_properties-1);
      ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin(mMeshId);	  
      
      ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin(mMeshId);

      // properties->PrintData(std::cout);
      // std::cout<<std::endl;
      
      MeshDataTransferUtilities DataTransferUtilities;

      const Element & rReferenceElement = mrRemesh.GetReferenceElement();
      
      const unsigned int nds = element_begin->GetGeometry().size();

      int& OutNumberOfElements = mrRemesh.OutMesh.GetNumberOfElements();
      int* OutElementList      = mrRemesh.OutMesh.GetElementList();

      std::vector<NodeType::Pointer>    list_of_element_centers;
      std::vector<Geometry<NodeType > > list_of_element_vertices; //is this list needed?
      //find the center and "radius" of the element
      double xc=0;
      double yc=0;
      double zc=0;
      double radius=0;    
      
      //generate kratos elements (conditions are not touched)
      int id = 0;
      // std::vector<std::vector<int> > EmptyNeighList;
      // mrRemesh.NeighbourList.swap(EmptyNeighList); 
      // mrRemesh.NeighbourList.clear(); //destroy all elements
      
      for(int el = 0; el<OutNumberOfElements; el++)
	{
	  if(mrRemesh.PreservedElements[el])
	    {
	      Geometry<NodeType> vertices;
	      std::vector<int>   neighbours (nds);
	      
	      for(unsigned int i=0; i<nds; i++)
		{
		  //note that OutElementList, starts from node 1, not from node 0, it can be directly assigned to mrRemesh.NodalPreIds.
		  vertices.push_back(*(nodes_begin + OutElementList[el*nds+i]-1).base());
		  //vertices.push_back(mrModelPart.pGetNode(OutElementList[el*3+pn],MeshId));
		  
		  if(vertices.back().Is(TO_ERASE))
		    std::cout<<" WARNING:: mesh vertex RELEASED "<<vertices.back().Id()<<std::endl;		  
		}
	      
	      id += 1;
	      
	      mrRemesh.PreservedElements[el] = id;
	      
	      
	      //*******************************************************************
	      //1) Store Preserved elements in an array of vertices (Geometry<NodeType > vertices;)
	      
	      
	      if( vertices.size() == 3 ){
		DataTransferUtilities.CalculateCenterAndSearchRadius( vertices[0].X(), vertices[0].Y(),
									vertices[1].X(), vertices[1].Y(),
									vertices[2].X(), vertices[2].Y(),
									xc,yc,radius );
	      }
	      else if( vertices.size() == 4 ){
		DataTransferUtilities.CalculateCenterAndSearchRadius( vertices[0].X(), vertices[0].Y(), vertices[0].Z(),
									vertices[1].X(), vertices[1].Y(), vertices[1].Z(),
									vertices[2].X(), vertices[2].Y(), vertices[2].Z(),
									vertices[3].X(), vertices[3].Y(), vertices[3].Z(),
									xc,yc,zc,radius );	      
	      }
	      
	      //std::cout<<" XC ["<<id<<"]: ("<<xc<<" "<<yc<<") "<<std::endl;
	      //std::cout<<" vertices "<<vertices[0].X()<<" "<<vertices[2].X()<<std::endl;
	      //*******************************************************************

	      NodeType::Pointer p_center = boost::make_shared< NodeType >( id, xc, yc, zc );
	      
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
	    mrRemesh.PreservedElements[el] = -1;
	  }
    
	    
	  // if(mrRemesh.PreservedElements[el])
	  //   {
	  // 	std::cout<<" Neighbours ["<<id-1<<"] :";
	  // 	for(int pn=0; pn<3; pn++)
	  // 	  {
	  // 	    std::cout<<" neighborlist ("<<pn<<") "<<mrRemesh.NeighbourList[id-1][pn]<<std::endl;
	  // 	  }
	  //   }
	  
	  //std::cout<<" NodalPreIds ["<<el<<"] :"<<NodalPreIds[el]<<std::endl;
	  
	}
      
      
      //*******************************************************************
      //5) Laplacian Smoothing
      
      //Check Mesh Info to perform smoothing:
      if( mrRemesh.Options.Is(ModelerUtilities::REFINE) )
	mrRemesh.Info->CheckGeometricalSmoothing();
      else
	mrRemesh.Info->GeometricalSmoothingRequired = true;
      
      //if(mrRemesh.smoothing && mrRemesh.remesh && mrRemesh.Info->GeometricalSmoothingRequired ){
      if( mrRemesh.Options.Is(ModelerUtilities::MESH_SMOOTHING) && mrRemesh.Info->GeometricalSmoothingRequired ){
	int& OutNumberOfPoints = mrRemesh.OutMesh.GetNumberOfPoints();
	LaplacianSmoothing  MeshGeometricSmoothing(mrModelPart);
	MeshGeometricSmoothing.SetEchoLevel(mEchoLevel);
	MeshGeometricSmoothing.ApplyMeshSmoothing(mrModelPart,mrRemesh.PreservedElements,OutElementList,OutNumberOfPoints,mMeshId);
      }
      //*******************************************************************
      
      
      //*******************************************************************
      //6) Pass  rReferenceElement and transfer variables
      DataTransferUtilities.TransferData(mrModelPart,rReferenceElement,list_of_element_centers,list_of_element_vertices,MeshDataTransferUtilities::ELEMENT_TO_ELEMENT,mMeshId);
      //*******************************************************************
      
      
      //*******************************************************************w
      //std::cout<<" Number of Nodes "<<mrModelPart.Nodes(MeshId).size()<<" Number Of Ids "<<mrRemesh.NodalPreIds.size()<<std::endl;
      
      //7) Restore global ID's
      for(ModelPart::NodesContainerType::iterator in = mrModelPart.NodesBegin(mMeshId) ; in != mrModelPart.NodesEnd(mMeshId) ; in++)
	{
	  //std::cout<<" node (local:"<<in->Id()<<", global:"<<mrRemesh.NodalPreIds[ in->Id() ]<<")"<<std::endl;
	  in->SetId( mrRemesh.NodalPreIds[ in->Id() ] );
	}
      //*******************************************************************
      
      
      //*******************************************************************
      
      //8) Filling the neighbour list
      SetElementNeighbours(mrModelPart);

      //*******************************************************************
      
      
      if( mEchoLevel > 0 )
	std::cout<<"   GENERATE NEW ELEMENTS ]; "<<std::endl;


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
    virtual std::string Info() const
    {
      return "BuildMeshElementsProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
      rOStream << "BuildMeshElementsProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

  protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

  private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    ModelPart& mrModelPart;

    ModelerUtilities::MeshingParameters& mrRemesh;

    int mMeshId;
    int mEchoLevel;

    ///@}
    ///@name Private Operators
    ///@{


    //**************************************************************************
    //**************************************************************************

    void SetElementNeighbours(ModelPart& rModelPart)
    {
      KRATOS_TRY

	if( mEchoLevel > 0 ){
	  std::cout<<" [ SET ELEMENT NEIGHBOURS : "<<std::endl;
	  std::cout<<"   Initial Faces : "<<rModelPart.Conditions(mMeshId).size()<<std::endl;
	}

      ModelPart::ElementsContainerType::iterator element_begin = rModelPart.ElementsBegin(mMeshId);	  

      const unsigned int nds = element_begin->GetGeometry().size();

      int* OutElementNeighbourList = mrRemesh.OutMesh.GetElementNeighbourList();

      int facecounter=0;
      int Id = 0;
      for(ModelPart::ElementsContainerType::const_iterator iii = rModelPart.ElementsBegin(mMeshId);
	  iii != rModelPart.ElementsEnd(mMeshId); iii++)
	{
	  
	  for(unsigned int i= 0; i<mrRemesh.PreservedElements.size(); i++)
	    {
	      if( mrRemesh.PreservedElements[Id] == -1)
		Id++;
	      else
		break;
	    }
	  
	  int number_of_faces = iii->GetGeometry().FacesNumber(); //defined for triangles and tetrahedra
	  (iii->GetValue(NEIGHBOUR_ELEMENTS)).resize(number_of_faces);
	  WeakPointerVector< Element >& neighb = iii->GetValue(NEIGHBOUR_ELEMENTS);

	  int index = 0;
	  for(int i = 0; i<number_of_faces; i++)
	    {
	      index = OutElementNeighbourList[Id*nds+i];  //mrRemesh.NeighbourList[Id][i];
			
	      if(index > 0)
		{
		  //std::cout<<" Element "<<Id<<" size "<<mrRemesh.PreservedElements.size()<<std::endl;			    
		  //std::cout<<" Index pre "<<index<<" size "<<mrRemesh.PreservedElements.size()<<std::endl;
		  index =  mrRemesh.PreservedElements[index-1]; 
		  //std::cout<<" Index post "<<index<<std::endl;
		}

	      if(index > 0)
		{
		  neighb(i) = *((element_begin + index -1 ).base());
		}
	      else
		{
		  //neighb(i) = Element::WeakPointer();
		  neighb(i) = *(iii.base());
		  facecounter++;
		}
	    }

	  Id++;
	}
	
      if( mEchoLevel > 0 ){
	std::cout<<"   Final Faces : "<<facecounter<<std::endl;
	std::cout<<"   SET ELEMENT NEIGHBOURS ]; "<<std::endl;
      }

      KRATOS_CATCH( "" )

    }

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    BuildMeshElementsProcess& operator=(BuildMeshElementsProcess const& rOther);

    /// Copy constructor.
    //BuildMeshElementsProcess(BuildMeshElementsProcess const& rOther);


    ///@}

  }; // Class BuildMeshElementsProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    BuildMeshElementsProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const BuildMeshElementsProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_BUILD_MESH_ELEMENTS_PROCESS_H_INCLUDED  defined 
