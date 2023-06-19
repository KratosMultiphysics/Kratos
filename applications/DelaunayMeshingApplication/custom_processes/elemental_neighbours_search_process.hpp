//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_ELEMENTAL_NEIGHBOURS_SEARCH_PROCESS_H_INCLUDED )
#define  KRATOS_ELEMENTAL_NEIGHBOURS_SEARCH_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "custom_processes/mesher_process.hpp"
#include "delaunay_meshing_application_variables.h"

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
class ElementalNeighboursSearchProcess
    : public MesherProcess
{
 public:
  ///@name Type Definitions
  ///@{
  /// Pointer definition of ElementalNeighboursSearchProcess
  KRATOS_CLASS_POINTER_DEFINITION( ElementalNeighboursSearchProcess );

  typedef  ModelPart::NodesContainerType NodesContainerType;
  typedef  ModelPart::ElementsContainerType ElementsContainerType;

  typedef Node::WeakPointer NodeWeakPtrType;
  typedef Element::WeakPointer ElementWeakPtrType;
  typedef Condition::WeakPointer ConditionWeakPtrType;

  typedef GlobalPointersVector<Node > NodeWeakPtrVectorType;
  typedef GlobalPointersVector<Element> ElementWeakPtrVectorType;
  typedef GlobalPointersVector<Condition> ConditionWeakPtrVectorType;
  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  /// avg_elems ------ expected number of neighbour elements per node.,
  /// avg_nodes ------ expected number of neighbour Nodes
  /// A better guess for the quantities above -> less memory occupied and faster algorithm

  ElementalNeighboursSearchProcess(ModelPart& rModelPart,
                                   int Dimension,
                                   int EchoLevel = 0,
                                   int AverageElements = 10)
      : mrModelPart(rModelPart)
  {
    mAverageElements = AverageElements;
    mDimension       = Dimension;
    mEchoLevel       = EchoLevel;
  }

  /// Destructor.
  virtual ~ElementalNeighboursSearchProcess()
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

    int method = 0;  //Kratos or Lohner method

    // double begin_time = OpenMPUtils::GetCurrentTime();

    if(method==0)
    {
      //std::cout<<" Kratos Search "<<std::endl;
      success=KratosSearch();
    }
    else
    {
      //std::cout<<" Lohner Search "<<std::endl;
      success=LohnerSearch(); //seems to be worse (needs to be optimized)
    }

    if(!success)
    {
      std::cout<<" ERROR:  Element Neighbours Search FAILED !!! "<<std::endl;
    }
    // else
    // {
    //   //print out the mesh generation time
    //   if( mEchoLevel > 1 ){
    //     double end_time = OpenMPUtils::GetCurrentTime();
    //     std::cout<<"  Neighbour Elements Search time = "<<end_time-begin_time<<std::endl;
    //   }
    //   //PrintElementNeighbours();
    // }


  };


  void ClearNeighbours()
  {
    for(auto& i_node: mrModelPart.Nodes())
    {
      ElementWeakPtrVectorType& nElements = i_node.GetValue(NEIGHBOUR_ELEMENTS);
      nElements.clear();
    }
    for(auto& i_elem : mrModelPart.Elements())
    {
      ElementWeakPtrVectorType& nElements = i_elem.GetValue(NEIGHBOUR_ELEMENTS);
      nElements.clear();
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
    return "ElementalNeighboursSearchProcess";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "ElementalNeighboursSearchProcess";
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
  int mDimension;
  int mEchoLevel;

  ///@}
  ///@name Private Operators
  ///@{
  template<class TDataType> void  AddUniquePointer
  (GlobalPointersVector<TDataType>& v, const typename TDataType::WeakPointer candidate)
  {
    typename GlobalPointersVector< TDataType >::iterator i = v.begin();
    typename GlobalPointersVector< TDataType >::iterator endit = v.end();
    while ( i != endit && (i)->Id() != (candidate.lock())->Id())
    {
      i++;
    }
    if( i == endit )
    {
      v.push_back(candidate);
    }

  }

  ElementWeakPtrType CheckForNeighbourElems1D (unsigned int Id_1, ElementWeakPtrVectorType& nElements, ElementsContainerType::iterator i_elem)
  {
    //look for the faces around node Id_1
    for(auto i_nelem(nElements.begin()); i_nelem != nElements.end(); ++i_nelem)
    {
      //look for the nodes of the neighbour faces
      Geometry<Node >& nGeometry = i_nelem->GetGeometry();
      if(nGeometry.LocalSpaceDimension() == 1){
        for(unsigned int node_i = 0; node_i < nGeometry.size(); ++node_i)
        {
          if(nGeometry[node_i].Id() == Id_1)
          {
            if(i_nelem->Id() != i_elem->Id())
            {
              return Element::WeakPointer(*i_nelem.base());
            }
          }
        }
      }
    }
    return Element::WeakPointer(*i_elem.base());
  }


  ElementWeakPtrType CheckForNeighbourElems2D (unsigned int Id_1, unsigned int Id_2, ElementWeakPtrVectorType& nElements, ElementsContainerType::iterator i_elem)
  {
    //look for the faces around node Id_1
    for(auto i_nelem(nElements.begin()); i_nelem != nElements.end(); ++i_nelem)
    {
      //look for the nodes of the neighbour faces
      Geometry<Node >& nGeometry = i_nelem->GetGeometry();
      if(nGeometry.LocalSpaceDimension() == 2){
        for(unsigned int node_i = 0; node_i < nGeometry.size(); ++node_i)
        {
          if (nGeometry[node_i].Id() == Id_2)
          {
            if(i_nelem->Id() != i_elem->Id())
            {
              return *i_nelem.base();
            }
          }
        }
      }
    }
    return *i_elem.base();
  }

  ElementWeakPtrType CheckForNeighbourElems3D (unsigned int Id_1, unsigned int Id_2, unsigned int Id_3, ElementWeakPtrVectorType& nElements, ElementsContainerType::iterator i_elem)
  {
    //look for the faces around node Id_1
    for(auto i_nelem(nElements.begin()); i_nelem != nElements.end(); ++i_nelem)
    {
      //look for the nodes of the neighbour faces
      Geometry<Node >& nGeometry = i_nelem->GetGeometry();
      if(nGeometry.LocalSpaceDimension() == 3){
        for(unsigned int node_i = 0; node_i < nGeometry.size(); ++node_i)
        {
          if(nGeometry[node_i].Id() == Id_2)
          {
            for(unsigned int node_j = 0; node_j < nGeometry.size(); ++node_j)
            {
              if (nGeometry[node_j].Id() == Id_3)
                if(i_nelem->Id() != i_elem->Id())
                {
                  return *i_nelem.base();
                }
            }
          }
        }
      }
    }
    return *i_elem.base();
  }


  void ResetFlagOptions (Node& rNode)
  {
    rNode.Reset(BOUNDARY);
  }

  void ResetFlagOptions (Element& rElement)
  {
    rElement.Reset(BOUNDARY);
  }


  void CleanElementNeighbours()
  {

    KRATOS_TRY

    //first of all the neighbour nodes and neighbour elements arrays are initialized to the guessed size
    //this cleans the old entries:

    //*************  Erase old node neighbours  *************//
    for(auto& i_node : mrModelPart.Nodes())
    {
      auto& nElements = i_node.GetValue(NEIGHBOUR_ELEMENTS);
      nElements.clear();
      nElements.reserve(mAverageElements);

      ResetFlagOptions(i_node);
    }

    //************* Erase old element neighbours ************//
    for(auto& i_elem : mrModelPart.Elements())
    {
      auto& nElements = i_elem.GetValue(NEIGHBOUR_ELEMENTS);
      nElements.clear();
      nElements.resize(i_elem.GetGeometry().FacesNumber());

      ResetFlagOptions(i_elem);
    }

    KRATOS_CATCH( "" )
  }


  void PrintElementNeighbours()
  {
    KRATOS_TRY

    std::cout<<" NODES: neighbour elems: "<<std::endl;
    for(auto& i_node : mrModelPart.Nodes())
    {
      std::cout<<"["<<i_node.Id()<<"]:"<<std::endl;
      std::cout<<"( ";
      auto& nElements = i_node.GetValue(NEIGHBOUR_ELEMENTS);
      for(const auto& i_nelem : nElements)
      {
        std::cout<< i_nelem.Id()<<", ";
      }
      std::cout<<" )"<<std::endl;
    }

    std::cout<<std::endl;

    std::cout<<" ELEMENTS: neighbour elems: "<<std::endl;

    for(auto& i_elem : mrModelPart.Elements())
    {
      std::cout<<"["<<i_elem.Id()<<"]:"<<std::endl;
      std::cout<<"( ";
      auto& nElements = i_elem.GetValue(NEIGHBOUR_ELEMENTS);
      for(auto& i_nelem : nElements)
      {
        std::cout<< i_nelem.Id()<<", ";
      }
      std::cout<<" )"<<std::endl;
    }


    std::cout<<std::endl;

    KRATOS_CATCH( "" )
  }



  bool KratosSearch()
  {

    KRATOS_TRY

    ElementsContainerType& rElements = mrModelPart.Elements();

    //first of all the neighbour nodes and neighbour elements arrays are initialized to the guessed size
    //this cleans the old entries:

    //*****  Erase old node and element neighbours  *********//
    CleanElementNeighbours();


    //*************  Neigbours of nodes  ************//
    //add the neighbour elements to all the nodes in the mesh
    for(auto i_elem(rElements.begin()); i_elem != rElements.end(); ++i_elem)
    {
      Element::GeometryType& rGeometry = i_elem->GetGeometry();
      for(unsigned int i = 0; i < rGeometry.size(); ++i)
      {
        rGeometry[i].GetValue(NEIGHBOUR_ELEMENTS).push_back(*i_elem.base());
      }
    }

    //*************  Neigbours of elements  *********//
    //add the neighbour elements to all the elements in the mesh

    unsigned int search_performed = false;

    //loop over faces
    if (mDimension==2)
    {
      for(auto i_elem(rElements.begin()); i_elem != rElements.end(); ++i_elem)
      {
        //face nodes
        Geometry<Node >& rGeometry = i_elem->GetGeometry();

        if( rGeometry.FacesNumber() == 3 ){

          auto& nElements = i_elem->GetValue(NEIGHBOUR_ELEMENTS);
          //vector of the 3 faces around the given face
          if(nElements.size() != 3 )
            nElements.resize(3);

          //neighb_face is the vector containing pointers to the three faces around ic:

          // neighbour element over edge 1-2 of element ic;
          nElements(0) = CheckForNeighbourElems2D(rGeometry[1].Id(), rGeometry[2].Id(), rGeometry[1].GetValue(NEIGHBOUR_ELEMENTS), i_elem);
          // neighbour element over edge 2-0 of element ic;
          nElements(1) = CheckForNeighbourElems2D(rGeometry[2].Id(), rGeometry[0].Id(), rGeometry[2].GetValue(NEIGHBOUR_ELEMENTS), i_elem);
          // neighbour element over edge 0-1 of element ic;
          nElements(2) = CheckForNeighbourElems2D(rGeometry[0].Id(), rGeometry[1].Id(), rGeometry[0].GetValue(NEIGHBOUR_ELEMENTS), i_elem);

          unsigned int iface=0;
          for(auto& i_nelem : nElements)
          {
            if (i_nelem.Id() == i_elem->Id())  // If there is no shared element in face nf (the Id coincides)
            {
              i_elem->Set(BOUNDARY);

              DenseMatrix<unsigned int> lpofa; //points that define the faces
              rGeometry.NodesInFaces(lpofa);

              for(unsigned int i = 1; i < rGeometry.FacesNumber(); ++i)
              {
                rGeometry[lpofa(i,iface)].Set(BOUNDARY);  //set boundary particles
              }
            }
            iface++;
          }

        }
        else if( rGeometry.FacesNumber() == 2 ){

          auto& nElements = i_elem->GetValue(NEIGHBOUR_ELEMENTS);

          //vector of the 2 faces around the given face
          if( nElements.size() != 2 )
            nElements.resize(2);

          //neighb_face is the vector containing pointers to the three faces around ic:

          // neighbour element over edge 0 of element ic;
          nElements(0) = CheckForNeighbourElems1D(rGeometry[0].Id(), rGeometry[0].GetValue(NEIGHBOUR_ELEMENTS), i_elem);
          // neighbour element over edge 1 of element ic;
          nElements(1) = CheckForNeighbourElems1D(rGeometry[1].Id(), rGeometry[1].GetValue(NEIGHBOUR_ELEMENTS), i_elem);

          unsigned int iface=0;
          for(auto& i_nelem : nElements)
          {
            if(i_nelem.Id() == i_elem->Id())  // If there is no shared element in face nf (the Id coincides)
            {
              i_elem->Set(BOUNDARY);

              DenseMatrix<unsigned int> lpofa; //points that define the faces
              rGeometry.NodesInFaces(lpofa);

              for(unsigned int i = 1; i < rGeometry.FacesNumber(); ++i)
              {
                rGeometry[lpofa(i,iface)].Set(BOUNDARY);  //set boundary particles
              }
            }
            iface++;
          }
        }
      }

      search_performed = true;
    }

    if (mDimension==3)
    {
      for(auto i_elem(rElements.begin()); i_elem != rElements.end(); ++i_elem)
      {
        //face nodes
        Geometry<Node >& rGeometry = i_elem->GetGeometry();

        if(rGeometry.FacesNumber() == 4){

          //vector of the 4 faces around the given element (3D tetrahedron)
          auto& nElements = i_elem->GetValue(NEIGHBOUR_ELEMENTS);

          if(nElements.size() != 4)
            nElements.resize(4);

          //neighb_face is the vector containing pointers to the three faces around ic:

          // neighbour element over face 1-2-3 of element ic;
          nElements(0) = CheckForNeighbourElems3D(rGeometry[1].Id(), rGeometry[2].Id(), rGeometry[3].Id(), rGeometry[1].GetValue(NEIGHBOUR_ELEMENTS), i_elem);
          // neighbour element over face 2-3-0 of element ic;
          nElements(1) = CheckForNeighbourElems3D(rGeometry[2].Id(), rGeometry[3].Id(), rGeometry[0].Id(), rGeometry[2].GetValue(NEIGHBOUR_ELEMENTS), i_elem);
          // neighbour element over face 3-0-1 of element ic;
          nElements(2) = CheckForNeighbourElems3D(rGeometry[3].Id(), rGeometry[0].Id(), rGeometry[1].Id(), rGeometry[3].GetValue(NEIGHBOUR_ELEMENTS), i_elem);
          // neighbour element over face 0-1-2 of element ic;
          nElements(3) = CheckForNeighbourElems3D(rGeometry[0].Id(), rGeometry[1].Id(), rGeometry[2].Id(), rGeometry[0].GetValue(NEIGHBOUR_ELEMENTS), i_elem);


          unsigned int iface=0;
          for(auto& i_nelem : nElements)
          {
            if(i_nelem.Id() == i_elem->Id())  // If there is no shared element in face nf (the Id coincides)
            {
              i_elem->Set(BOUNDARY);

              DenseMatrix<unsigned int> lpofa; //points that define the faces
              rGeometry.NodesInFaces(lpofa);

              for(unsigned int i = 1; i < rGeometry.FacesNumber(); ++i)
              {
                rGeometry[lpofa(i,iface)].Set(BOUNDARY);  //set boundary particles
                //std::cout<<" SetBoundary ("<<rGeometry[lpofa(i,0)].Id()<<")"<<std::endl;
              }
            }
            iface++;
          }

        }
        else if(rGeometry.FacesNumber() == 3){

          //vector of the 3 faces around the given element (3D triangle)
          auto& nElements = i_elem->GetValue(NEIGHBOUR_ELEMENTS);

          if(nElements.size() != 3)
            nElements.resize(3);

          //neighb_face is the vector containing pointers to the three faces around ic:

          // neighbour element over edge 1-2 of element ic;
          nElements(0) = CheckForNeighbourElems2D(rGeometry[1].Id(), rGeometry[2].Id(), rGeometry[1].GetValue(NEIGHBOUR_ELEMENTS), i_elem);
          // neighbour element over edge 2-0 of element ic;
          nElements(1) = CheckForNeighbourElems2D(rGeometry[2].Id(), rGeometry[0].Id(), rGeometry[2].GetValue(NEIGHBOUR_ELEMENTS), i_elem);
          // neighbour element over edge 0-1 of element ic;
          nElements(2) = CheckForNeighbourElems2D(rGeometry[0].Id(), rGeometry[1].Id(), rGeometry[0].GetValue(NEIGHBOUR_ELEMENTS), i_elem);

          unsigned int iface=0;
          for(auto& i_nelem : nElements)
          {
            if(i_nelem.Id() == i_elem->Id())  // If there is no shared element in face nf (the Id coincides)
            {
              i_elem->Set(BOUNDARY);

              Geometry<Node >& rGeometry = (i_elem)->GetGeometry();

              DenseMatrix<unsigned int> lpofa; //points that define the faces
              rGeometry.NodesInFaces(lpofa);

              for(unsigned int i = 1; i < rGeometry.FacesNumber(); ++i)
              {
                rGeometry[lpofa(i,iface)].Set(BOUNDARY);  //set boundary particles
              }
            }
            iface++;
          }

        }
      }
      search_performed = true;
    }

    if( mrModelPart.NumberOfElements()>0 && search_performed )
      return true;
    else
      return false;

    KRATOS_CATCH( "" )
  }


  bool LohnerSearch()
  {

    KRATOS_TRY

    NodesContainerType& rNodes = mrModelPart.Nodes();
    ElementsContainerType& rElements = mrModelPart.Elements();


    unsigned int Ne=rElements.size();
    unsigned int Np=rNodes.size();

    //*****  Erase old node and element neighbours  *********//
    CleanElementNeighbours();


    //*************  Neigbours of nodes  ************//
    //add the neighbour elements to all the nodes in the mesh
    for(auto i_elem(rElements.begin()); i_elem != rElements.end(); ++i_elem)
    {
      Element::GeometryType& rGeometry = i_elem->GetGeometry();
      if(rGeometry.LocalSpaceDimension() == mrModelPart.GetProcessInfo()[SPACE_DIMENSION]){
        for(unsigned int i = 0; i < rGeometry.size(); ++i)
        {
          rGeometry[i].GetValue(NEIGHBOUR_ELEMENTS).push_back(*i_elem.base());
        }
      }
    }

    //*************  Neigbours of elements  *********//
    //add the neighbour elements to all the elements in the mesh
    //loop over faces

    unsigned int ipoin=0;
    unsigned int nnofa=0;
    unsigned int jelem=0;
    unsigned int icoun=0;
    unsigned int jpoin=0;
    unsigned int nnofj=0;
    unsigned int nface=0;

    DenseVector<unsigned int> lnofa; //number of nodes per face
    DenseMatrix<unsigned int> lpofa; //points that define the faces

    Element::GeometryType& rGeometry = rElements.begin()->GetGeometry(); // the first element is taken as reference
    unsigned int Nf= rGeometry.FacesNumber();     //number of faces

    //lnofa and lpofa defined in Geometry of the element (rGeometry): triangle, quadrilateral, tetrahedron ...
    rGeometry.NumberNodesInFaces(lnofa);
    rGeometry.NodesInFaces(lpofa);

    //Auxiliary vectors
    DenseVector<unsigned int> lhelp (Nf-1); //can be only 2 or 3 nodes per face : Triangles(faces of 2 nodes) Tetrahedra(faces of 3 nodes)
    lhelp.clear();
    DenseVector<unsigned int> lpoin (Np+1);
    lpoin.clear();


    //Elements Surrounding Elements
    int el;
#pragma omp parallel for reduction(+:nface) private(el,ipoin,nnofa,jelem,icoun,jpoin,nnofj) firstprivate(lhelp,lpoin)
    for (el=1; el<(int)Ne+1; ++el) //ELEMENTS START FROM el=1
    {

      for (unsigned int nf=0; nf<Nf; ++nf) //loop over faces
      {
        nnofa=lnofa(nf);

        //Initially assign the same element as a neighbour
        rElements[el].GetValue(NEIGHBOUR_ELEMENTS)(nf) = rElements(el);

        //constant vector, depends on the element
        for (unsigned int t=0; t<nnofa; ++t)
        {
          lhelp(t)=rElements[el].GetGeometry()[lpofa(t,nf)].Id();  //connections of the face
          lpoin(lhelp(t))=1;                                    //mark in lpoin
        }

        ipoin=lhelp(1);   //select a point

        auto& nElements = rNodes[ipoin].GetValue(NEIGHBOUR_ELEMENTS);

        for(auto& i_nelem : nElements)  //loop over elements surronding a point
        {
          jelem=i_nelem.Id();
          unsigned int ielem =rElements[el].Id();

          if(jelem!=ielem)
          {

            for(unsigned int fel=0; fel<Nf; ++fel) //loop over the element faces
            {
              nnofj=lnofa(fel);

              if (nnofj==nnofa)
              {

                icoun=0;
                for (unsigned int jnofa=0; jnofa<nnofa; ++jnofa) //loop to count the number of equal points
                {
                  jpoin= rElements[jelem].GetGeometry()[lpofa(jnofa,fel)].Id();
                  icoun= icoun+lpoin(jpoin);
                }

                if(icoun==nnofa)
                {
                  //store the element
                  rElements[el].GetValue(NEIGHBOUR_ELEMENTS)(nf) = rElements(jelem);
                  //std::cout<<" el "<<el<<" shared "<<jelem<<std::endl;
                }
              }
            }
          }
        }


        if (rElements[el].GetValue(NEIGHBOUR_ELEMENTS)[nf].Id() == rElements[el].Id())  // If there is no shared element in face nf (the Id coincides)
        {

          rElements[el].Set(BOUNDARY);

          //unsigned int nfixed=0;
          for (unsigned int t=0; t<nnofa; ++t) //loop on number of nodes per face
          {
            rNodes[lhelp(t)].Set(BOUNDARY);  //set boundary particles
          }

          nface+=1;

        }
        //loop B is outside to parallelize with omp


        for (unsigned int r=0; r<nnofa; ++r)
        {
          lpoin(lhelp(r))=0;                            //reset lpoin
        }
      }
    }


    //detection of the boundary elements with no face in the boundary and layer elements

    for(auto& i_elem : rElements)
    {
      Element::GeometryType& rGeometry = i_elem.GetGeometry();
      for(unsigned int i = 0; i < rGeometry.size(); ++i)
      {
        if(rGeometry[i].Is(BOUNDARY))
        {
          i_elem.Set(BOUNDARY);
        }
      }

    }

    return true;

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
  ElementalNeighboursSearchProcess& operator=(ElementalNeighboursSearchProcess const& rOther);

  /// Copy constructor.
  //ElementalNeighboursSearchProcess(ElementalNeighboursSearchProcess const& rOther);


  ///@}

}; // Class ElementalNeighboursSearchProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ElementalNeighboursSearchProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ElementalNeighboursSearchProcess& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ELEMENTAL_NEIGHBOURS_SEARCH_PROCESS_H_INCLUDED  defined
