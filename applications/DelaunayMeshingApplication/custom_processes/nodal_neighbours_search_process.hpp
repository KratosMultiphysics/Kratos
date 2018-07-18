//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_NODAL_NEIGHBOURS_SEARCH_PROCESS_H_INCLUDED )
#define  KRATOS_NODAL_NEIGHBOURS_SEARCH_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "custom_processes/mesher_process.hpp"

///VARIABLES used:
//Data:     NEIGHBOUR_ELEMENTS(set), NEIGHBOUR_NODES(set)
//StepData:
//Flags:    (checked)
//          (set)     STRUCTURE(set)->when a dof is fixed
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
  typedef  ModelPart::NodesContainerType NodesContainerType;
  typedef  ModelPart::ElementsContainerType ElementsContainerType;


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
  class NodalNeighboursSearchProcess
    : public MesherProcess
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NodalNeighboursSearchProcess
    KRATOS_CLASS_POINTER_DEFINITION( NodalNeighboursSearchProcess );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /// avg_elems ------ expected number of neighbour elements per node.,
    /// avg_nodes ------ expected number of neighbour Nodes
    /// the better the guess for the quantities above the less memory occupied and the fastest the algorithm
    NodalNeighboursSearchProcess(ModelPart& rModelPart,
				 int EchoLevel = 0,
				 int AverageElements = 10,
				 int AverageNodes = 10)
      : mrModelPart(rModelPart)
    {
      mAverageElements = AverageNodes;
      mAverageNodes = AverageElements;
      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~NodalNeighboursSearchProcess()
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
      bool success=false;

      int method = 0; //Kratos or Lohner method

      double begin_time = OpenMPUtils::GetCurrentTime();

      if(method==0)
	{
	  success=KratosSearch();
	}
      else
	{
	  success=LohnerSearch(); //seems to be worse (needs to be optimized)
	}



      if(!success)
	{
	  std::cout<<" ERROR:  Nodal Neighbours Search FAILED !!! "<<std::endl;
	}
      else
	{
	  //print out the mesh generation time
	  if( mEchoLevel > 1 ){
            double end_time = OpenMPUtils::GetCurrentTime();
	    std::cout<<"  Neighbour Nodes Search time = "<<end_time-begin_time<<std::endl;
          }
	  //PrintNodeNeighbours();
	}

    };

    void ClearNeighbours()
    {
      NodesContainerType& rNodes = mrModelPart.Nodes();
      for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); ++in)
	{
	  WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);
	  rE.erase(rE.begin(),rE.end());

	  WeakPointerVector<Node<3> >& rN = in->GetValue(NEIGHBOUR_NODES);
	  rN.erase(rN.begin(),rN.end() );
	}

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
      return "NodalNeighboursSearchProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "NodalNeighboursSearchProcess";
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
    int mAverageElements;
    int mAverageNodes;
    int mEchoLevel;

    ///@}
    ///@name Private Operators
    ///@{

    //******************************************************************************************
    //******************************************************************************************
    template< class TDataType > void  AddUniqueWeakPointer
    (WeakPointerVector< TDataType >& v, const typename TDataType::WeakPointer candidate)
    {
      typename WeakPointerVector< TDataType >::iterator i = v.begin();
      typename WeakPointerVector< TDataType >::iterator endit = v.end();
      while ( i != endit && (i)->Id() != (candidate.lock())->Id())
	{
	  ++i;
	}
      if( i == endit )
	{
	  v.push_back(candidate);
	}

    }

    void CleanNodeNeighbours()
    {
      NodesContainerType&    rNodes = mrModelPart.Nodes();
      //*************  Erase old node neighbours  *************//
      for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); ++in)
	{
	  (in->GetValue(NEIGHBOUR_NODES)).reserve(mAverageNodes);
	  WeakPointerVector<Node<3> >& rN = in->GetValue(NEIGHBOUR_NODES);
	  rN.erase(rN.begin(),rN.end() );

	  (in->GetValue(NEIGHBOUR_ELEMENTS)).reserve(mAverageElements);
	  WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);
	  rE.erase(rE.begin(),rE.end() );

	  //set fixed nodes as Nodes<3>::STRUCTURE  to not be removed in the meshing
	  Node<3>::DofsContainerType& node_dofs = in->GetDofs();
	  for( Node<3>::DofsContainerType::const_iterator i_dof = node_dofs.begin() ; i_dof != node_dofs.end() ; ++i_dof)
	    {
	      if(i_dof->IsFixed()){
		in->Set(STRUCTURE);
		break;
	      }
	    }

	}

      //std::cout<<"  [ Node Neighbours CLEAN ] "<<std::endl;
    }


    void PrintNodeNeighbours()
    {

      NodesContainerType& rNodes = mrModelPart.Nodes();
      std::cout<<" NODES: neighbour elems: "<<std::endl;
      for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); ++in)
	{
	  std::cout<<"["<<in->Id()<<"]:"<<std::endl;
	  std::cout<<"( ";
	  WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);
	  for(unsigned int i = 0; i < rE.size(); ++i)
	    {
	      std::cout<< rE[i].Id()<<", ";
	    }
	  std::cout<<" )"<<std::endl;
	}

      std::cout<<std::endl;

      std::cout<<" NODES: neighbour nodes: "<<std::endl;

      for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); ++in)
	{
	  std::cout<<"["<<in->Id()<<"]:"<<std::endl;
	  std::cout<<"( ";
	  WeakPointerVector<Node<3> >& rN = in->GetValue(NEIGHBOUR_NODES);
	  for(unsigned int i = 0; i < rN.size(); ++i)
	    {
	      std::cout<< rN[i].Id()<<", ";
	    }
	  std::cout<<" )"<<std::endl;
	}

      std::cout<<std::endl;

    }
    ///@}
    ///@name Private Operations
    ///@{

    bool KratosSearch()
    {
      NodesContainerType&    rNodes = mrModelPart.Nodes();
      ElementsContainerType& rElems = mrModelPart.Elements();

      //*************  Erase old node neighbours  *************//
      CleanNodeNeighbours();

      //*************  Neighbours of nodes  ************//

      //add the neighbour elements to all the nodes in the mesh
      for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ++ie)
	{
	  Element::GeometryType& pGeom = ie->GetGeometry();
	  for(unsigned int i = 0; i < pGeom.size(); ++i)
	    {
	      //KRATOS_WATCH( pGeom[i] )
	      (pGeom[i].GetValue(NEIGHBOUR_ELEMENTS)).push_back( Element::WeakPointer( *(ie.base()) ) );
	      //KRATOS_WATCH( (pGeom[i].GetValue(NEIGHBOUR_ELEMENTS)).size() )
	    }
	}


      //adding the neighbouring nodes to all nodes in the mesh
      for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); ++in)
	{
	  WeakPointerVector< Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);
	  //KRATOS_WATCH ( *in )
	  //KRATOS_WATCH ( rE.size() )
	  for(unsigned int ie = 0; ie < rE.size(); ++ie)
	    {
	      Element::GeometryType& pGeom = rE[ie].GetGeometry();
	      for(unsigned int i = 0; i < pGeom.size(); ++i)
		{
		  //std::cout<<" inside pgeom loop {"<<i<<"} rE.size() : "<<rE.size()<<std::endl;
		  if( pGeom[i].Id() != in->Id() )
		    {
		      Element::NodeType::WeakPointer temp = pGeom(i);
		      WeakPointerVector< Node<3> >& rN = in->GetValue(NEIGHBOUR_NODES);
		      AddUniqueWeakPointer< Node<3> >(rN, temp);
		      //std::cout<<" inside add unique {"<<i<<"} rE.size() : "<<rE.size()<<std::endl;
		    }

		}
	    }
	}

      return true;
    }


    bool LohnerSearch()
    {
      NodesContainerType&    rNodes = mrModelPart.Nodes();
      ElementsContainerType& rElems = mrModelPart.Elements();

      //*************  Erase old node neighbours  *************//
      CleanNodeNeighbours();


      //*************  Neighbours of nodes  ************//
      unsigned int Ne=rElems.size();
      unsigned int Np=rNodes.size();

      if(Ne==1) return false; //there are no elements

      //Elements  Position vs Nodes
      std::vector<unsigned int> PSharedE (Np+1);   //vector counting shared Elements
      PSharedE.clear();
      //Particles Position vs Particles
      std::vector<unsigned int> PSharedN (Np+1);   //vector counting shared Particles
      PSharedN.clear();

      std::vector<std::vector<unsigned int> >   PSurroundN (Np+1); //General matrix to account Neighbour Particles


      //ELEMENTS SURROUNDING A POINT
      //1.- Count the number of elements connected to each point
      int ipoi1=0;
      PSharedE[0]=1;
      PSharedN[0]=1;

      for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ++ie)
	{
	  Element::GeometryType& pGeom = ie->GetGeometry();
	  for(unsigned int i = 0; i < pGeom.size(); ++i)
	    {
	      ipoi1=pGeom[i].Id();      //counter
	      PSharedE[ipoi1]+=1;       //auxiliar counter: esup2
	    }
	}

      Element::GeometryType& pGeom = rElems.begin()->GetGeometry(); // the first element is taken as reference
      unsigned int Nn= pGeom.size();

      //2.- Reshuffling pass (1)
      unsigned int pn;
      for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); ++in)
	{
	  pn= in->Id();
	  int size= PSharedE[pn]*(Nn-1)+Nn; //it is an estimation of the size like mAverageNodes
	  if(mAverageNodes>size)
	    size=mAverageNodes;

	  (in->GetValue(NEIGHBOUR_ELEMENTS)).resize(PSharedE[pn]);
	  PSurroundN[pn].resize(size);
	}

      //3.- Store the elements in PSurroundE
      for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ++ie)
	{
	  Element::GeometryType& pGeom = ie->GetGeometry();
	  for(unsigned int i = 0; i < pGeom.size(); ++i)
	    {
	      ipoi1=pGeom[i].Id();         //counter
	      PSharedE[ipoi1]-=1;
	      //std::cout<<" NODE "<<ie->Id()<<" "<<PSharedE[ipoi1]<<std::endl;
	      WeakPointerVector<Element >& rNeighElems = rNodes[ipoi1].GetValue(NEIGHBOUR_ELEMENTS);
	      rNeighElems(PSharedE[ipoi1])= Element::WeakPointer( *(ie.base()) );
	    }
	}


      //POINTS SURROUNDING A POINT

      unsigned int nodeID=0;
      unsigned int iel=0;
      unsigned int ipn=0;
      unsigned int rpn=0;

      PSharedN.clear();

      for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); ++in)
	{
	  WeakPointerVector<Element >& rNeighElems = in->GetValue(NEIGHBOUR_ELEMENTS);
	  rpn = in->Id();

	  for (unsigned int sel=0; sel<rNeighElems.size(); ++sel)
	    {

	      iel = rNeighElems[sel].Id();

	      Element::GeometryType& pGeom = rElems[iel].GetGeometry();


	      for (unsigned int nd=0; nd<pGeom.size(); ++nd)
		{
		  ipn = pGeom[nd].Id();

		  if(ipn != rpn )  // Process to find the neighbour points of a point
		    {
		      nodeID=1;

		      if (PSharedN[rpn]!=0)
			{
			  for(unsigned int spn=0; spn<=PSharedN[rpn]; ++spn)
			    {
			      if (PSurroundN[rpn][spn]!=ipn)
				{
				  nodeID=1;
				}
			      else
				{
				  nodeID=0;
				  break;
				}
			    }
			}


		      if (nodeID)
			{
			  PSurroundN[rpn][PSharedN[rpn] ]=ipn;
			  PSharedN[rpn]+=1;
			}
		    }

		}
	    }

	  //Set everything on particles
	  WeakPointerVector<Node<3> >& rN = in->GetValue(NEIGHBOUR_NODES);


	  //std::cout<<" NODE "<<rpn<<" "<<PSharedN[rpn]<<std::endl;


	  for(unsigned int spn=0; spn<PSharedN[rpn]; ++spn)
	    {
	      //std::cout<<" ShNodes "<<PSurroundN[rpn][spn]<<std::endl;
	      Element::NodeType::WeakPointer temp = rNodes(PSurroundN[rpn][spn]);
	      AddUniqueWeakPointer< Node<3> >(rN, temp);

	    }

	}

      return true;
    }


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
    NodalNeighboursSearchProcess& operator=(NodalNeighboursSearchProcess const& rOther);

    /// Copy constructor.
    //NodalNeighboursSearchProcess(NodalNeighboursSearchProcess const& rOther);


    ///@}

  }; // Class NodalNeighboursSearchProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    NodalNeighboursSearchProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const NodalNeighboursSearchProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_NODAL_NEIGHBOURS_SEARCH_PROCESS_H_INCLUDED  defined
