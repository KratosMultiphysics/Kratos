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
class NodalNeighboursSearchProcess
    : public MesherProcess
{
 public:
  ///@name Type Definitions
  ///@{
  /// Pointer definition of NodalNeighboursSearchProcess
  KRATOS_CLASS_POINTER_DEFINITION( NodalNeighboursSearchProcess );

  typedef  ModelPart::NodesContainerType NodesContainerType;
  typedef  ModelPart::ElementsContainerType ElementsContainerType;

  typedef GlobalPointersVector<Node<3> > NodeWeakPtrVectorType;
  typedef GlobalPointersVector<Element> ElementWeakPtrVectorType;
  typedef GlobalPointersVector<Condition> ConditionWeakPtrVectorType;
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

    // double begin_time = OpenMPUtils::GetCurrentTime();

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
    // else
    // {
    //   //print out the mesh generation time
    //   if( mEchoLevel > 1 ){
    //     double end_time = OpenMPUtils::GetCurrentTime();
    //     std::cout<<"  Neighbour Nodes Search time = "<<end_time-begin_time<<std::endl;
    //   }
    //   //PrintNodeNeighbours();
    // }

  };

  void ClearNeighbours()
  {
    for(auto& i_node : mrModelPart.Nodes())
    {
      i_node.GetValue(NEIGHBOUR_ELEMENTS).clear();
      i_node.GetValue(NEIGHBOUR_NODES).clear();
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
  template<class TDataType> void  AddUniquePointer
  (GlobalPointersVector<TDataType>& v, const typename TDataType::WeakPointer candidate)
  {
    typename GlobalPointersVector< TDataType >::iterator i = v.begin();
    typename GlobalPointersVector< TDataType >::iterator endit = v.end();
    while ( i != endit && (i)->Id() != (candidate)->Id())
    {
      i++;
    }
    if( i == endit )
    {
      v.push_back(candidate);
    }

  }

  void CleanNodeNeighbours()
  {
    //*************  Erase old node neighbours  *************//
    for(auto& i_node : mrModelPart.Nodes())
    {
      NodeWeakPtrVectorType& nNodes = i_node.GetValue(NEIGHBOUR_NODES);
      nNodes.clear();
      nNodes.reserve(mAverageNodes);

      ElementWeakPtrVectorType& nElements = i_node.GetValue(NEIGHBOUR_ELEMENTS);
      nElements.clear();
      nElements.reserve(mAverageElements);

      //set fixed nodes as Nodes<3>::STRUCTURE  to not be removed in the meshing
      for(const auto& i_dof : i_node.GetDofs())
      {
        if(i_dof->IsFixed())
        {
          i_node.Set(STRUCTURE);
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
    for(auto& i_node : mrModelPart.Nodes())
    {
      std::cout<<"["<<i_node.Id()<<"]:"<<std::endl;
      std::cout<<"( ";
      ElementWeakPtrVectorType& nElements = i_node.GetValue(NEIGHBOUR_ELEMENTS);
      for(auto& i_nelem : nElements)
      {
        std::cout<< i_nelem.Id()<<", ";
      }
      std::cout<<" )"<<std::endl;
    }

    std::cout<<std::endl;

    std::cout<<" NODES: neighbour nodes: "<<std::endl;

    for(auto& i_node : rNodes)
    {
      std::cout<<"["<<i_node.Id()<<"]:"<<std::endl;
      std::cout<<"( ";
      NodeWeakPtrVectorType& nNodes = i_node.GetValue(NEIGHBOUR_NODES);
      for(auto& i_nnode : nNodes)
      {
        std::cout<< i_nnode.Id()<<", ";
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
    ElementsContainerType& rElements = mrModelPart.Elements();

    //*************  Erase old node neighbours  *************//
    CleanNodeNeighbours();

    //*************  Neighbours of nodes  ************//

    //add the neighbour elements to all the nodes in the mesh
    for(auto i_elem(rElements.begin()); i_elem != rElements.end(); ++i_elem)
    {
      Element::GeometryType& rGeometry = i_elem->GetGeometry();
      for(unsigned int i = 0; i < rGeometry.size(); ++i)
      {
        rGeometry[i].GetValue(NEIGHBOUR_ELEMENTS).push_back(*i_elem.base());
      }
    }

    //adding the neighbouring nodes to all nodes in the mesh
    for(auto& i_node : rNodes)
    {
      ElementWeakPtrVectorType& nElements = i_node.GetValue(NEIGHBOUR_ELEMENTS);
      for(auto& i_nelem : nElements)
      {
        Element::GeometryType& rGeometry = i_nelem.GetGeometry();
        for(unsigned int i = 0; i < rGeometry.size(); ++i)
        {
          if( rGeometry[i].Id() != i_node.Id() )
          {
            NodeWeakPtrVectorType& nNodes = i_node.GetValue(NEIGHBOUR_NODES);
            AddUniquePointer<Node<3> >(nNodes, rGeometry(i));
          }
        }
      }
    }

    return true;
  }


  bool LohnerSearch()
  {
    NodesContainerType&    rNodes = mrModelPart.Nodes();
    ElementsContainerType& rElements = mrModelPart.Elements();

    //*************  Erase old node neighbours  *************//
    CleanNodeNeighbours();

    //*************  Neighbours of nodes  ************//
    unsigned int Ne=rElements.size();
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

    for(auto& i_elem : rElements)
    {
      Element::GeometryType& rGeometry = i_elem.GetGeometry();
      for(unsigned int i = 0; i < rGeometry.size(); ++i)
      {
        ipoi1=rGeometry[i].Id();      //counter
        PSharedE[ipoi1]+=1;       //auxiliar counter: esup2
      }
    }

    Element::GeometryType& rGeometry = rElements.begin()->GetGeometry(); // the first element is taken as reference
    unsigned int Nn= rGeometry.size();

    //2.- Reshuffling pass (1)
    unsigned int pn;
    for(auto& i_node : rNodes)
    {
      pn=i_node.Id();
      int size= PSharedE[pn]*(Nn-1)+Nn; //it is an estimation of the size like mAverageNodes
      if(mAverageNodes>size)
        size=mAverageNodes;

      i_node.GetValue(NEIGHBOUR_ELEMENTS).resize(PSharedE[pn]);
      PSurroundN[pn].resize(size);
    }

    //3.- Store the elements in PSurroundE
    for(auto i_elem(rElements.begin()); i_elem != rElements.end(); ++i_elem)
    {
      Element::GeometryType& rGeometry = i_elem->GetGeometry();
      for(unsigned int i = 0; i < rGeometry.size(); ++i)
      {
        ipoi1=rGeometry[i].Id();         //counter
        PSharedE[ipoi1]-=1;
        //std::cout<<" NODE "<<ie->Id()<<" "<<PSharedE[ipoi1]<<std::endl;
        ElementWeakPtrVectorType& nElements= rNodes[ipoi1].GetValue(NEIGHBOUR_ELEMENTS);
        nElements(PSharedE[ipoi1])= *i_elem.base();
      }
    }

    //POINTS SURROUNDING A POINT
    unsigned int nodeID=0;
    unsigned int ipn=0;
    unsigned int rpn=0;

    PSharedN.clear();

    for(auto& i_node : rNodes)
    {
      ElementWeakPtrVectorType& nElements = i_node.GetValue(NEIGHBOUR_ELEMENTS);

      rpn = i_node.Id();

      for(auto& i_nelem : nElements)
      {

        Element::GeometryType& rGeometry = i_nelem.GetGeometry();

        for (unsigned int nd=0; nd<rGeometry.size(); ++nd)
        {
          ipn = rGeometry[nd].Id();

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
      NodeWeakPtrVectorType& nNodes = i_node.GetValue(NEIGHBOUR_NODES);

      for(unsigned int spn=0; spn<PSharedN[rpn]; ++spn)
      {
        //std::cout<<" ShNodes "<<PSurroundN[rpn][spn]<<std::endl;
        AddUniquePointer<Node<3> >(nNodes, rNodes(PSurroundN[rpn][spn]));
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
