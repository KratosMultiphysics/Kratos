//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_GENERATE_NEW_ELEMENTS_MESHER_PROCESS_H_INCLUDED )
#define  KRATOS_GENERATE_NEW_ELEMENTS_MESHER_PROCESS_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_meshers/laplacian_smoothing.hpp"
#include "custom_processes/mesher_process.hpp"

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
  class GenerateNewElementsMesherProcess
    : public MesherProcess
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GenerateNewElementsMesherProcess
    KRATOS_CLASS_POINTER_DEFINITION( GenerateNewElementsMesherProcess );

    typedef ModelPart::NodeType                   NodeType;
    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GenerateNewElementsMesherProcess(ModelPart& rModelPart,
			     MesherUtilities::MeshingParameters& rRemeshingParameters,
			     int EchoLevel)
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    {
      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~GenerateNewElementsMesherProcess()
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

    void Execute() override
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
	  SelectElementsMesherProcess SelectElements(mrModelPart,mrRemesh,mEchoLevel);
	  SelectElements.Execute();
	}


      //*******************************************************************
      //setting new elements
      //(mrModelPart.Elements()).reserve(mrRemesh.Info->NumberOfElements);


      //*******************************************************************
      // mine 2016 TIP
      // //All nodes in boundary element change
      // if(mrRemesh.AvoidTipElementsFlag){ //is not working correctly some dispositions not considered
      //   if( mEchoLevel > 0 )
      // 	std::cout<<"[   AVOID TIP ELEMENTS START ]"<<std::endl;

      //   ChangeTipElementsUtilities TipElements;
      //   //TipElements.SwapDiagonals(mrModelPart,out,mrRemesh.PreservedElements);

      //   if( mEchoLevel > 0 )
      // 	std::cout<<"[   AVOID TIP ELEMENTS END ]"<<std::endl;
      // }
      //*******************************************************************


      //properties to be used in the generation
      int number_properties = mrModelPart.GetParentModelPart()->NumberOfProperties();
      Properties::Pointer properties = mrModelPart.GetParentModelPart()->GetMesh().pGetProperties(number_properties-1);
      ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin();

      ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin();

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

      for(int el = 0; el<OutNumberOfElements; ++el)
	{
	  if(mrRemesh.PreservedElements[el])
	    {
	      Geometry<NodeType> vertices;
	      std::vector<int>   neighbours (nds);

	      for(unsigned int i=0; i<nds; ++i)
		{
		  //note that OutElementList, starts from node 1, not from node 0, it can be directly assigned to mrRemesh.NodalPreIds.
		  vertices.push_back(*(nodes_begin + OutElementList[el*nds+i]-1).base());
		  //vertices.push_back(mrModelPart.pGetNode(OutElementList[el*3+pn]));

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

	      NodeType::Pointer p_center = Kratos::make_shared< NodeType >( id, xc, yc, zc );

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
	  // 	for(int pn=0; pn<3; ++pn)
	  // 	  {
	  // 	    std::cout<<" neighborlist ("<<pn<<") "<<mrRemesh.NeighbourList[id-1][pn]<<std::endl;
	  // 	  }
	  //   }

	  //std::cout<<" NodalPreIds ["<<el<<"] :"<<NodalPreIds[el]<<std::endl;

	}


      //*******************************************************************
      //5) Laplacian Smoothing

      //Check Mesh Info to perform smoothing:
      if( mrRemesh.Options.Is(MesherUtilities::REFINE) )
	mrRemesh.Info->CheckGeometricalSmoothing();
      else
	mrRemesh.Info->GeometricalSmoothingRequired = true;

      // check geometrical smoothing is better for solid cutting problems commented July 2018 for fluid testing
      //if( mrRemesh.Options.Is(MesherUtilities::MESH_SMOOTHING) && mrRemesh.Info->GeometricalSmoothingRequired ){
      if( mrRemesh.Options.Is(MesherUtilities::MESH_SMOOTHING) ){
      	int& OutNumberOfPoints = mrRemesh.OutMesh.GetNumberOfPoints();
	LaplacianSmoothing  MeshGeometricSmoothing(mrModelPart);
	MeshGeometricSmoothing.SetEchoLevel(mEchoLevel);
	MeshGeometricSmoothing.ApplyMeshSmoothing(mrModelPart,mrRemesh.PreservedElements,OutElementList,OutNumberOfPoints);
      }
      //*******************************************************************


      //*******************************************************************
      //6) Pass  rReferenceElement and transfer variables
      DataTransferUtilities.TransferData(mrModelPart,rReferenceElement,list_of_element_centers,list_of_element_vertices,MeshDataTransferUtilities::ELEMENT_TO_ELEMENT);
      //*******************************************************************


      //*******************************************************************w
      //std::cout<<" Number of Nodes "<<mrModelPart.Nodes().size()<<" Number Of Ids "<<mrRemesh.NodalPreIds.size()<<std::endl;

      //7) Restore global ID's
      for(ModelPart::NodesContainerType::iterator in = mrModelPart.NodesBegin() ; in != mrModelPart.NodesEnd() ; ++in)
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
    std::string Info() const override
    {
      return "GenerateNewElementsMesherProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "GenerateNewElementsMesherProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

    MesherUtilities::MeshingParameters& mrRemesh;

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
	  std::cout<<"   Initial Faces : "<<rModelPart.Conditions().size()<<std::endl;
	}

      ModelPart::ElementsContainerType::iterator element_begin = rModelPart.ElementsBegin();

      const unsigned int nds = element_begin->GetGeometry().size();

      int* OutElementNeighbourList = mrRemesh.OutMesh.GetElementNeighbourList();

      int facecounter=0;
      int Id = 0;
      for(ModelPart::ElementsContainerType::const_iterator ie = rModelPart.ElementsBegin();
	  ie != rModelPart.ElementsEnd(); ++ie)
	{

	  for(unsigned int i= 0; i<mrRemesh.PreservedElements.size(); ++i)
	    {
	      if( mrRemesh.PreservedElements[Id] == -1)
	  	Id++;
	      else
	  	break;
	    }

	  //Id = ie->Id()-1;

	  unsigned int number_of_faces = ie->GetGeometry().FacesNumber(); //defined for triangles and tetrahedra
	  (ie->GetValue(NEIGHBOUR_ELEMENTS)).resize(number_of_faces);
	  WeakPointerVector< Element >& neighb = ie->GetValue(NEIGHBOUR_ELEMENTS);

	  int index = 0;
	  for(unsigned int iface = 0; iface<number_of_faces; ++iface)
	    {

	      index = OutElementNeighbourList[Id*nds+iface];

	      if(index > 0)
		{
		  //std::cout<<" Element "<<Id<<" size "<<mrRemesh.PreservedElements.size()<<std::endl;
		  //std::cout<<" Index pre "<<index<<" size "<<mrRemesh.PreservedElements.size()<<std::endl;
		  index =  mrRemesh.PreservedElements[index-1];
		  //std::cout<<" Index post "<<index<<std::endl;
		}

	      if(index > 0)
		{
		  //neighb(iface) = (mrModelPart.Elements()).find( elements_begin->Id() + index -1 );
		  neighb(iface) = *((element_begin + index -1 ).base());
		}
	      else
		{
		  //neighb(iface) = Element::WeakPointer();
		  neighb(iface) = *(ie.base());
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
    GenerateNewElementsMesherProcess& operator=(GenerateNewElementsMesherProcess const& rOther);

    /// Copy constructor.
    //GenerateNewElementsMesherProcess(GenerateNewElementsMesherProcess const& rOther);


    ///@}

  }; // Class GenerateNewElementsMesherProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    GenerateNewElementsMesherProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const GenerateNewElementsMesherProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_GENERATE_NEW_ELEMENTS_MESHER_PROCESS_H_INCLUDED  defined
