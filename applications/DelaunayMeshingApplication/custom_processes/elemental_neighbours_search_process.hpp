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

///VARIABLES used:
//Data:     NEIGHBOUR_ELEMENTS(set)
//StepData:
//Flags:    (checked) BOUNDARY
//          (set)     BOUNDARY(set)
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
  class ElementalNeighboursSearchProcess
    : public MesherProcess
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ElementalNeighboursSearchProcess
    KRATOS_CLASS_POINTER_DEFINITION( ElementalNeighboursSearchProcess );

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

      double begin_time = OpenMPUtils::GetCurrentTime();

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
      else
        {
	  //print out the mesh generation time
	  if( mEchoLevel > 1 ){
            double end_time = OpenMPUtils::GetCurrentTime();
            std::cout<<"  Neighbour Elements Search time = "<<end_time-begin_time<<std::endl;
          }
	  //PrintElementNeighbours();
        }


    };


    void ClearNeighbours()
    {
      NodesContainerType& rNodes = mrModelPart.Nodes();
      for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); ++in)
        {
	  WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);
	  rE.erase(rE.begin(),rE.end());
        }
      ElementsContainerType& rElems = mrModelPart.Elements();
      for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ++ie)
        {
	  WeakPointerVector<Element >& rE = ie->GetValue(NEIGHBOUR_ELEMENTS);
	  rE.erase(rE.begin(),rE.end());
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


    Element::WeakPointer CheckForNeighbourElems1D (unsigned int Id_1, WeakPointerVector< Element >& neighbour_elem, ElementsContainerType::iterator elem)
    {
      //look for the faces around node Id_1
      for( WeakPointerVector< Element >::iterator i =neighbour_elem.begin(); i != neighbour_elem.end(); ++i)
        {
	  //look for the nodes of the neighbour faces
	  Geometry<Node<3> >& neigh_elem_geometry = (i)->GetGeometry();
          if( neigh_elem_geometry.LocalSpaceDimension() == 1 ){
            for( unsigned int node_i = 0 ; node_i < neigh_elem_geometry.size(); ++node_i)
            {
	      if (neigh_elem_geometry[node_i].Id() == Id_1)
              {
                if(i->Id() != elem->Id())
                {
                  return *(i.base());
                }
              }
            }
          }
        }
      return *(elem.base());
    }


    Element::WeakPointer CheckForNeighbourElems2D (unsigned int Id_1, unsigned int Id_2, WeakPointerVector< Element >& neighbour_elem, ElementsContainerType::iterator elem)
    {
      //look for the faces around node Id_1
      for( WeakPointerVector< Element >::iterator i =neighbour_elem.begin(); i != neighbour_elem.end(); ++i)
        {
	  //look for the nodes of the neighbour faces
	  Geometry<Node<3> >& neigh_elem_geometry = (i)->GetGeometry();
          if( neigh_elem_geometry.LocalSpaceDimension() == 2 ){
            for( unsigned int node_i = 0 ; node_i < neigh_elem_geometry.size(); ++node_i)
            {
	      if (neigh_elem_geometry[node_i].Id() == Id_2)
              {
                if(i->Id() != elem->Id())
                    {
		      return *(i.base());
                    }
              }
            }
          }
        }
      return *(elem.base());
    }

    Element::WeakPointer CheckForNeighbourElems3D (unsigned int Id_1, unsigned int Id_2, unsigned int Id_3, WeakPointerVector< Element >& neighbour_elem, ElementsContainerType::iterator elem)
    {
      //look for the faces around node Id_1
      for( WeakPointerVector< Element >::iterator i = neighbour_elem.begin(); i != neighbour_elem.end(); ++i)
        {
	  //look for the nodes of the neighbour faces
	  Geometry<Node<3> >& neigh_elem_geometry = (i)->GetGeometry();
          if( neigh_elem_geometry.LocalSpaceDimension() == 3 ){
            for( unsigned int node_i = 0 ; node_i < neigh_elem_geometry.size(); ++node_i)
            {
	      if (neigh_elem_geometry[node_i].Id() == Id_2)
              {
                for( unsigned int node_j = 0 ; node_j < neigh_elem_geometry.size(); ++node_j)
                {
                  if (neigh_elem_geometry[node_j].Id() == Id_3)
                    if(i->Id() != elem->Id())
                    {
                      return *(i.base());
                    }
                }
              }
            }
          }
        }
      return *(elem.base());
    }


    void ResetFlagOptions (Node<3>& Node)
    {
      Node.Reset(BOUNDARY);
    }

    void ResetFlagOptions (Element& Elem)
    {
      Elem.Reset(BOUNDARY);
    }


    void CleanElementNeighbours()
    {

      KRATOS_TRY

      NodesContainerType&    rNodes = mrModelPart.Nodes();
      ElementsContainerType& rElems = mrModelPart.Elements();

      //first of all the neighbour nodes and neighbour elements arrays are initialized to the guessed size
      //this cleans the old entries:

      //*************  Erase old node neighbours  *************//
      for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); ++in)
        {
          (in->GetValue(NEIGHBOUR_ELEMENTS)).reserve(mAverageElements);
	  WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);
	  rE.erase(rE.begin(),rE.end());

	  ResetFlagOptions(*in);
        }

      //************* Erase old element neighbours ************//
      for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ++ie)
        {
	  Element::GeometryType& pGeom = ie->GetGeometry();
	  int size= pGeom.FacesNumber();

	  WeakPointerVector<Element >& rE = ie->GetValue(NEIGHBOUR_ELEMENTS);
	  rE.erase(rE.begin(),rE.end() );

	  (ie->GetValue(NEIGHBOUR_ELEMENTS)).resize(size);

	  ResetFlagOptions(*ie);
        }

      KRATOS_CATCH( "" )
    }


    void PrintElementNeighbours()
    {
      KRATOS_TRY

      NodesContainerType& rNodes = mrModelPart.Nodes();
      ElementsContainerType& rElems = mrModelPart.Elements();

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

      std::cout<<" ELEMENTS: neighbour elems: "<<std::endl;

      for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ++ie)
        {
	  std::cout<<"["<<ie->Id()<<"]:"<<std::endl;
	  std::cout<<"( ";
	  WeakPointerVector<Element >& rE = ie->GetValue(NEIGHBOUR_ELEMENTS);
	  for(unsigned int i = 0; i < rE.size(); ++i)
            {
	      std::cout<< rE[i].Id()<<", ";
            }
	  std::cout<<" )"<<std::endl;
        }


      std::cout<<std::endl;

      KRATOS_CATCH( "" )
    }



    bool KratosSearch()
    {

      KRATOS_TRY

      ElementsContainerType& rElems = mrModelPart.Elements();

      //first of all the neighbour nodes and neighbour elements arrays are initialized to the guessed size
      //this cleans the old entries:

      //*****  Erase old node and element neighbours  *********//
      CleanElementNeighbours();


      //*************  Neigbours of nodes  ************//
      //add the neighbour elements to all the nodes in the mesh
      for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ++ie)
        {
	  Element::GeometryType& pGeom = ie->GetGeometry();
          for(unsigned int i = 0; i < pGeom.size(); ++i)
          {
            (pGeom[i].GetValue(NEIGHBOUR_ELEMENTS)).push_back( Element::WeakPointer( *(ie.base()) ) );
          }
        }

      //*************  Neigbours of elements  *********//
      //add the neighbour elements to all the elements in the mesh

      unsigned int search_performed = false;

      //loop over faces
      if (mDimension==2)
        {
	  for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ++ie)
            {
	      //face nodes
	      Geometry<Node<3> >& rGeometry = (ie)->GetGeometry();

	      if( rGeometry.FacesNumber() == 3 ){

		//vector of the 3 faces around the given face
		if( ie->GetValue(NEIGHBOUR_ELEMENTS).size() != 3 )
		  (ie->GetValue(NEIGHBOUR_ELEMENTS)).resize(3);

		WeakPointerVector< Element >& neighb_elems = ie->GetValue(NEIGHBOUR_ELEMENTS);

		//neighb_face is the vector containing pointers to the three faces around ic:

		// neighbour element over edge 1-2 of element ic;
		neighb_elems(0) = CheckForNeighbourElems2D(rGeometry[1].Id(), rGeometry[2].Id(), rGeometry[1].GetValue(NEIGHBOUR_ELEMENTS), ie);
		// neighbour element over edge 2-0 of element ic;
		neighb_elems(1) = CheckForNeighbourElems2D(rGeometry[2].Id(), rGeometry[0].Id(), rGeometry[2].GetValue(NEIGHBOUR_ELEMENTS), ie);
		// neighbour element over edge 0-1 of element ic;
		neighb_elems(2) = CheckForNeighbourElems2D(rGeometry[0].Id(), rGeometry[1].Id(), rGeometry[0].GetValue(NEIGHBOUR_ELEMENTS), ie);

		unsigned int iface=0;
		for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ++ne)
		  {
		    if (ne->Id() == ie->Id())  // If there is no shared element in face nf (the Id coincides)
		      {

			ie->Set(BOUNDARY);

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

		//vector of the 2 faces around the given face
		if( ie->GetValue(NEIGHBOUR_ELEMENTS).size() != 2 )
		  (ie->GetValue(NEIGHBOUR_ELEMENTS)).resize(2);

		WeakPointerVector< Element >& neighb_elems = ie->GetValue(NEIGHBOUR_ELEMENTS);

		//neighb_face is the vector containing pointers to the three faces around ic:

		// neighbour element over edge 0 of element ic;
		neighb_elems(0) = CheckForNeighbourElems1D(rGeometry[0].Id(), rGeometry[0].GetValue(NEIGHBOUR_ELEMENTS), ie);
		// neighbour element over edge 1 of element ic;
		neighb_elems(1) = CheckForNeighbourElems1D(rGeometry[1].Id(), rGeometry[1].GetValue(NEIGHBOUR_ELEMENTS), ie);

		unsigned int iface=0;
		for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ++ne)
		  {
		    if (ne->Id() == ie->Id())  // If there is no shared element in face nf (the Id coincides)
		      {

			ie->Set(BOUNDARY);

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
	  for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ++ie)
            {
	      //face nodes
	      Geometry<Node<3> >& rGeometry = (ie)->GetGeometry();

	      if( rGeometry.FacesNumber() == 4 ){

		//vector of the 4 faces around the given element (3D tetrahedron)
		if( ie->GetValue(NEIGHBOUR_ELEMENTS).size() != 4 )
		  (ie->GetValue(NEIGHBOUR_ELEMENTS)).resize(4);

		WeakPointerVector< Element >& neighb_elems = ie->GetValue(NEIGHBOUR_ELEMENTS);

		//neighb_face is the vector containing pointers to the three faces around ic:

		// neighbour element over face 1-2-3 of element ic;
		neighb_elems(0) = CheckForNeighbourElems3D(rGeometry[1].Id(), rGeometry[2].Id(), rGeometry[3].Id(), rGeometry[1].GetValue(NEIGHBOUR_ELEMENTS), ie);
		// neighbour element over face 2-3-0 of element ic;
		neighb_elems(1) = CheckForNeighbourElems3D(rGeometry[2].Id(), rGeometry[3].Id(), rGeometry[0].Id(), rGeometry[2].GetValue(NEIGHBOUR_ELEMENTS), ie);
		// neighbour element over face 3-0-1 of element ic;
		neighb_elems(2) = CheckForNeighbourElems3D(rGeometry[3].Id(), rGeometry[0].Id(), rGeometry[1].Id(), rGeometry[3].GetValue(NEIGHBOUR_ELEMENTS), ie);
		// neighbour element over face 0-1-2 of element ic;
		neighb_elems(3) = CheckForNeighbourElems3D(rGeometry[0].Id(), rGeometry[1].Id(), rGeometry[2].Id(), rGeometry[0].GetValue(NEIGHBOUR_ELEMENTS), ie);


		unsigned int iface=0;
		for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ++ne)
		  {
		    if (ne->Id() == ie->Id())  // If there is no shared element in face nf (the Id coincides)
                    {

                      ie->Set(BOUNDARY);

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
	      else if( rGeometry.FacesNumber() == 3 ){

		//vector of the 3 faces around the given element (3D triangle)
		if( ie->GetValue(NEIGHBOUR_ELEMENTS).size() != 3 )
		  (ie->GetValue(NEIGHBOUR_ELEMENTS)).resize(3);

		WeakPointerVector< Element >& neighb_elems = ie->GetValue(NEIGHBOUR_ELEMENTS);

		//neighb_face is the vector containing pointers to the three faces around ic:

		// neighbour element over edge 1-2 of element ic;
		neighb_elems(0) = CheckForNeighbourElems2D(rGeometry[1].Id(), rGeometry[2].Id(), rGeometry[1].GetValue(NEIGHBOUR_ELEMENTS), ie);
		// neighbour element over edge 2-0 of element ic;
		neighb_elems(1) = CheckForNeighbourElems2D(rGeometry[2].Id(), rGeometry[0].Id(), rGeometry[2].GetValue(NEIGHBOUR_ELEMENTS), ie);
		// neighbour element over edge 0-1 of element ic;
		neighb_elems(2) = CheckForNeighbourElems2D(rGeometry[0].Id(), rGeometry[1].Id(), rGeometry[0].GetValue(NEIGHBOUR_ELEMENTS), ie);

		unsigned int iface=0;
		for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ++ne)
		  {
		    if (ne->Id() == ie->Id())  // If there is no shared element in face nf (the Id coincides)
		      {

			ie->Set(BOUNDARY);

			Geometry<Node<3> >& rGeometry = (ie)->GetGeometry();

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

      NodesContainerType&    rNodes = mrModelPart.Nodes();
      ElementsContainerType& rElems = mrModelPart.Elements();


      unsigned int Ne=rElems.size();
      unsigned int Np=rNodes.size();

      //*****  Erase old node and element neighbours  *********//
      CleanElementNeighbours();


      //*************  Neigbours of nodes  ************//
      //add the neighbour elements to all the nodes in the mesh
      for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ++ie)
        {
	  Element::GeometryType& pGeom = ie->GetGeometry();
          if(pGeom.LocalSpaceDimension() == mrModelPart.GetProcessInfo()[SPACE_DIMENSION]){
            for(unsigned int i = 0; i < pGeom.size(); ++i)
            {
	      (pGeom[i].GetValue(NEIGHBOUR_ELEMENTS)).push_back( Element::WeakPointer( *(ie.base()) ) );
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

      Element::GeometryType& pGeom = rElems.begin()->GetGeometry(); // the first element is taken as reference
      unsigned int Nf= pGeom.FacesNumber();     //number of faces

      //lnofa and lpofa defined in Geometry of the element (mpGeometry): triangle, quadrilateral, tetrahedron ...
      pGeom.NumberNodesInFaces(lnofa);
      pGeom.NodesInFaces(lpofa);

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
	      rElems[el].GetValue(NEIGHBOUR_ELEMENTS)(nf) =  Element::WeakPointer( rElems(el) );

	      //constant vector, depends on the element
	      for (unsigned int t=0; t<nnofa; ++t)
                {
		  lhelp(t)=rElems[el].GetGeometry()[lpofa(t,nf)].Id();  //connections of the face
		  lpoin(lhelp(t))=1;                                    //mark in lpoin
                }

	      ipoin=lhelp(1);   //select a point

	      WeakPointerVector< Element >& n_elems = rNodes[ipoin].GetValue(NEIGHBOUR_ELEMENTS);

	      for(unsigned int esp=0; esp<n_elems.size(); ++esp)  //loop over elements surronding a point
                {
		  jelem=n_elems[esp].Id();
		  unsigned int iel  =rElems[el].Id();

		  if(jelem!=iel)
                    {

		      for(unsigned int fel=0; fel<Nf; ++fel) //loop over the element faces
                        {
			  nnofj=lnofa(fel);

			  if (nnofj==nnofa)
                            {

			      icoun=0;
			      for (unsigned int jnofa=0; jnofa<nnofa; ++jnofa) //loop to count the number of equal points
                                {
				  jpoin= rElems[jelem].GetGeometry()[lpofa(jnofa,fel)].Id();
				  icoun= icoun+lpoin(jpoin);
                                }

			      if(icoun==nnofa)
                                {
				  //store the element
				  rElems[el].GetValue(NEIGHBOUR_ELEMENTS)(nf) =  Element::WeakPointer( rElems(jelem) );
				  //std::cout<<" el "<<el<<" shared "<<jelem<<std::endl;
                                }
                            }
                        }
                    }
                }


	      if (rElems[el].GetValue(NEIGHBOUR_ELEMENTS)[nf].Id() == rElems[el].Id())  // If there is no shared element in face nf (the Id coincides)
                {

		  rElems[el].Set(BOUNDARY);

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

      for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ++ie)
        {
	  Element::GeometryType& pGeom = ie->GetGeometry();
	  for(unsigned int i = 0; i < pGeom.size(); ++i)
            {
	      if(pGeom[i].Is(BOUNDARY))
                {
		  ie->Set(BOUNDARY);
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
