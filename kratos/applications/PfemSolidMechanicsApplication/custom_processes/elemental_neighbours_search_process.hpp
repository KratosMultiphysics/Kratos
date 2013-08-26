//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_ELEMENTAL_NEIGHBOURS_SEARCH_PROCESS_H_INCLUDED )
#define  KRATOS_ELEMENTAL_NEIGHBOURS_SEARCH_PROCESS_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes
#include <boost/timer.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "pfem_solid_mechanics_application.h"

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
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ElementalNeighboursSearchProcess
    KRATOS_CLASS_POINTER_DEFINITION(ElementalNeighboursSearchProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /// avg_elems ------ expected number of neighbour elements per node.,
    /// avg_nodes ------ expected number of neighbour Nodes
    /// A better guess for the quantities above -> less memory occupied and faster algorithm

    ElementalNeighboursSearchProcess(ModelPart& model_part,
                              unsigned int dim,
                              unsigned int avg_elems = 10,
                              unsigned int method=0)
        : mr_model_part(model_part)
    {
        mavg_elems = avg_elems;
        mdim=dim;
        m_method=method;
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

    virtual void Execute()
    {
        bool success=false;

        boost::timer auxiliary;

        if(m_method==0)
        {
            success=KratosSearch();
        }
        else
        {
            success=LohnerSearch(); //seems to be worse (needs to be optimized)
        }

        if(!success)
        {
            std::cout<<" ERROR:  Element Neighbours Search FAILED !!! "<<std::endl;
        }
        else
        {
            //print out the mesh generation time
            std::cout<<"Neighbour Elements Search time = "<<auxiliary.elapsed()<<std::endl;
            //PrintElementNeighbours();
        }

    };


    void ClearNeighbours()
    {
        NodesContainerType& rNodes = mr_model_part.Nodes();
        for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);
            rE.erase(rE.begin(),rE.end());
        }
        ElementsContainerType& rElems = mr_model_part.Elements();
        for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
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
    virtual std::string Info() const
    {
        return "ElementalNeighboursSearchProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ElementalNeighboursSearchProcess";
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
    ModelPart& mr_model_part;
    unsigned int mavg_elems;
    unsigned int m_method;
    unsigned int mdim;


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
            i++;
        }
        if( i == endit )
        {
            v.push_back(candidate);
        }

    }

    Element::WeakPointer CheckForNeighbourElems2D (unsigned int Id_1, unsigned int Id_2, WeakPointerVector< Element >& neighbour_elem, ElementsContainerType::iterator elem)
    {
        //look for the faces around node Id_1
        for( WeakPointerVector< Element >::iterator i =neighbour_elem.begin(); i != neighbour_elem.end(); i++)
        {
            //look for the nodes of the neighbour faces
            Geometry<Node<3> >& neigh_elem_geometry = (i)->GetGeometry();
            for( unsigned int node_i = 0 ; node_i < neigh_elem_geometry.size(); node_i++)
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
        return *(elem.base());
    }

    Element::WeakPointer CheckForNeighbourElems3D (unsigned int Id_1, unsigned int Id_2, unsigned int Id_3, WeakPointerVector< Element >& neighbour_elem, ElementsContainerType::iterator elem)
    {
        //look for the faces around node Id_1
        for( WeakPointerVector< Element >::iterator i = neighbour_elem.begin(); i != neighbour_elem.end(); i++)
        {
            //look for the nodes of the neighbour faces
            Geometry<Node<3> >& neigh_elem_geometry = (i)->GetGeometry();
            for( unsigned int node_i = 0 ; node_i < neigh_elem_geometry.size(); node_i++)
            {
                if (neigh_elem_geometry[node_i].Id() == Id_2)
                {
                    for( unsigned int node_j = 0 ; node_j < neigh_elem_geometry.size(); node_j++)
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


    void CleanElementNeighbours ()
    {

        NodesContainerType&    rNodes = mr_model_part.Nodes();
        ElementsContainerType& rElems = mr_model_part.Elements();

        //first of all the neighbour nodes and neighbour elements arrays are initialized to the guessed size
        //this cleans the old entries:

        //*************  Erase old node neighbours  *************//
        for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            (in->GetValue(NEIGHBOUR_ELEMENTS)).reserve(mavg_elems);
            WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);
            rE.erase(rE.begin(),rE.end() );

            ResetFlagOptions(*in);
        }

        //************* Erase old element neighbours ************//
        for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
        {
            Element::GeometryType& pGeom = ie->GetGeometry();
            int size= pGeom.EdgesNumber();

            //(ie->GetValue(NEIGHBOUR_ELEMENTS)).reserve(size);

            WeakPointerVector<Element >& rE = ie->GetValue(NEIGHBOUR_ELEMENTS);
            rE.erase(rE.begin(),rE.end() );

            (ie->GetValue(NEIGHBOUR_ELEMENTS)).resize(size);

            ResetFlagOptions(*ie);
        }

    }


    void PrintElementNeighbours ()
    {
        NodesContainerType& rNodes = mr_model_part.Nodes();
        ElementsContainerType& rElems = mr_model_part.Elements();

        std::cout<<" NODES: neighbour elems: "<<std::endl;
        for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            std::cout<<"["<<in->Id()<<"]:"<<std::endl;
            std::cout<<"( ";
            WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);
            for(unsigned int i = 0; i < rE.size(); i++)
            {
                std::cout<< rE[i].Id()<<", ";
            }
            std::cout<<" )"<<std::endl;
        }

        std::cout<<std::endl;

        std::cout<<" ELEMENTS: neighbour elems: "<<std::endl;

        for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
        {
            std::cout<<"["<<ie->Id()<<"]:"<<std::endl;
            std::cout<<"( ";
            WeakPointerVector<Element >& rE = ie->GetValue(NEIGHBOUR_ELEMENTS);
            for(unsigned int i = 0; i < rE.size(); i++)
            {
                std::cout<< rE[i].Id()<<", ";
            }
            std::cout<<" )"<<std::endl;
        }


        std::cout<<std::endl;

    }



    bool KratosSearch()
    {
        //NodesContainerType&    rNodes = mr_model_part.Nodes();
        ElementsContainerType& rElems = mr_model_part.Elements();
        //first of all the neighbour nodes and neighbour elements arrays are initialized to the guessed size
        //this cleans the old entries:

        //*****  Erase old node and element neighbours  *********//
        CleanElementNeighbours();


        //*************  Neigbours of nodes  ************//
        //add the neighbour elements to all the nodes in the mesh
        for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
        {
            Element::GeometryType& pGeom = ie->GetGeometry();
            for(unsigned int i = 0; i < pGeom.size(); i++)
            {
                (pGeom[i].GetValue(NEIGHBOUR_ELEMENTS)).push_back( Element::WeakPointer( *(ie.base()) ) );
            }
        }


        //*************  Neigbours of elements  *********//
        //add the neighbour elements to all the elements in the mesh
        //loop over faces
        if (mdim==2)
        {
            for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
            {
                //face nodes
                Geometry<Node<3> >& geom = (ie)->GetGeometry();
                //vector of the 3 faces around the given face
                (ie->GetValue(NEIGHBOUR_ELEMENTS)).resize(3);
                WeakPointerVector< Element >& neighb_elems = ie->GetValue(NEIGHBOUR_ELEMENTS);

                //neighb_face is the vector containing pointers to the three faces around ic:

                // neighbour element over edge 1-2 of element ic;
                neighb_elems(0) = CheckForNeighbourElems2D(geom[1].Id(), geom[2].Id(), geom[1].GetValue(NEIGHBOUR_ELEMENTS), ie);
                // neighbour element over edge 2-0 of element ic;
                neighb_elems(1) = CheckForNeighbourElems2D(geom[2].Id(), geom[0].Id(), geom[2].GetValue(NEIGHBOUR_ELEMENTS), ie);
                // neighbour element over edge 0-1 of element ic;
                neighb_elems(2) = CheckForNeighbourElems2D(geom[0].Id(), geom[1].Id(), geom[0].GetValue(NEIGHBOUR_ELEMENTS), ie);

		unsigned int counter=0;
                for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ne++)
                {
                    if (ne->Id() == ie->Id())  // If there is no shared element in face nf (the Id coincides)
                    {

                        ie->Set(BOUNDARY);

                        Geometry<Node<3> >& pGeom = (ie)->GetGeometry();
			
			boost::numeric::ublas::matrix<unsigned int> lpofa; //points that define the faces
			pGeom.NodesInFaces(lpofa);
			
                        for(unsigned int i = 0; i < pGeom.size(); i++)
			  {
			    if(i!=counter)
			      pGeom[lpofa(i,0)].Set(BOUNDARY);  //set boundary particles
			  }
			
                    }
		    
		    counter++;
                }
            }
        }
        if (mdim==3)
        {
            for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
            {
                //face nodes
                Geometry<Node<3> >& geom = (ie)->GetGeometry();
                //vector of the 3 faces around the given face
                (ie->GetValue(NEIGHBOUR_ELEMENTS)).resize(4);
                WeakPointerVector< Element >& neighb_elems = ie->GetValue(NEIGHBOUR_ELEMENTS);

                //neighb_face is the vector containing pointers to the three faces around ic:

                // neighbour element over face 1-2-3 of element ic;
                neighb_elems(0) = CheckForNeighbourElems3D(geom[1].Id(), geom[2].Id(), geom[3].Id(), geom[1].GetValue(NEIGHBOUR_ELEMENTS), ie);
                // neighbour element over face 2-3-0 of element ic;
                neighb_elems(1) = CheckForNeighbourElems3D(geom[2].Id(), geom[3].Id(), geom[0].Id(), geom[2].GetValue(NEIGHBOUR_ELEMENTS), ie);
                // neighbour element over face 3-0-1 of element ic;
                neighb_elems(2) = CheckForNeighbourElems3D(geom[3].Id(), geom[0].Id(), geom[1].Id(), geom[3].GetValue(NEIGHBOUR_ELEMENTS), ie);
                // neighbour element over face 0-1-2 of element ic;
                neighb_elems(3) = CheckForNeighbourElems3D(geom[0].Id(), geom[1].Id(), geom[2].Id(), geom[0].GetValue(NEIGHBOUR_ELEMENTS), ie);

		unsigned int counter=0;
                for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ne++)
                {
                    if (ne->Id() == ie->Id())  // If there is no shared element in face nf (the Id coincides)
                    {
                        ie->Set(BOUNDARY);

                        Geometry<Node<3> >& pGeom = (ie)->GetGeometry();

			boost::numeric::ublas::matrix<unsigned int> lpofa; //points that define the faces
			pGeom.NodesInFaces(lpofa);
			
                        for(unsigned int i = 0; i < pGeom.size(); i++)
                        {
			    if(i!=counter)
			      pGeom[lpofa(i,0)].Set(BOUNDARY);  //set boundary particles
			}

                    }
		    counter++;
                }


            }
        }

        return true;
    }


    bool LohnerSearch()
    {

        NodesContainerType&    rNodes = mr_model_part.Nodes();
        ElementsContainerType& rElems = mr_model_part.Elements();


        unsigned int Ne=rElems.size();
        unsigned int Np=rNodes.size();

        //*****  Erase old node and element neighbours  *********//
        CleanElementNeighbours();


        //*************  Neigbours of nodes  ************//
        //add the neighbour elements to all the nodes in the mesh
        for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
        {
            Element::GeometryType& pGeom = ie->GetGeometry();
            for(unsigned int i = 0; i < pGeom.size(); i++)
            {
                (pGeom[i].GetValue(NEIGHBOUR_ELEMENTS)).push_back( Element::WeakPointer( *(ie.base()) ) );
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

        boost::numeric::ublas::vector<unsigned int> lnofa; //number of nodes per face
        boost::numeric::ublas::matrix<unsigned int> lpofa; //points that define the faces

        Element::GeometryType& pGeom = rElems.begin()->GetGeometry(); // the first element is taken as reference
        unsigned int Nf= pGeom.FacesNumber();     //number of faces

        //lnofa and lpofa defined in Geometry of the element (mpGeometry): triangle, quadrilateral, tetrahedron ...
        pGeom.NumberNodesInFaces(lnofa);
        pGeom.NodesInFaces(lpofa);

        //Auxiliary vectors
        boost::numeric::ublas::vector<unsigned int> lhelp (Nf-1); //can be only 2 or 3 nodes per face : Triangles(faces of 2 nodes) Tetrahedra(faces of 3 nodes)
        lhelp.clear();
        boost::numeric::ublas::vector<unsigned int> lpoin (Np+1);
        lpoin.clear();


        //Elements Surrounding Elements
        int el;
        #pragma omp parallel for reduction(+:nface) private(el,ipoin,nnofa,jelem,icoun,jpoin,nnofj) firstprivate(lhelp,lpoin)
        for (el=1; el<int(Ne+1); el++) //ELEMENTS START FROM el=1
        {

            for (unsigned int nf=0; nf<Nf; nf++) //loop over faces
            {
                nnofa=lnofa(nf);

                //Initially assign the same element as a neighbour
                rElems[el].GetValue(NEIGHBOUR_ELEMENTS)(nf) =  Element::WeakPointer( rElems(el) );

                //constant vector, depends on the element
                for (unsigned int t=0; t<nnofa; t++)
                {
                    lhelp(t)=rElems[el].GetGeometry()[lpofa(t,nf)].Id();      //connections of the face
                    lpoin(lhelp(t))=1;                                    //mark in lpoin
                }

                ipoin=lhelp(1);   //select a point

                WeakPointerVector< Element >& n_elems = rNodes[ipoin].GetValue(NEIGHBOUR_ELEMENTS);

                for(unsigned int esp=0; esp<n_elems.size(); esp++)  //loop over elements surronding a point
                {
                    jelem=n_elems[esp].Id();
                    unsigned int iel  =rElems[el].Id();

                    if(jelem!=iel)
                    {

                        for(unsigned int fel=0; fel<Nf; fel++) //loop over the element faces
                        {
                            nnofj=lnofa(fel);

                            if (nnofj==nnofa)
                            {

                                icoun=0;
                                for (unsigned int jnofa=0; jnofa<nnofa; jnofa++) //loop to count the number of equal points
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
                    for (unsigned int t=0; t<nnofa; t++) //loop on number of nodes per face
                    {
                        rNodes[lhelp(t)].Set(BOUNDARY);  //set boundary particles
                    }

                    nface+=1;

                }
                //loop B is outside to parallelize with omp


                for (unsigned int r=0; r<nnofa; r++)
                {
                    lpoin(lhelp(r))=0;                            //reset lpoin
                }
            }
        }


        //detection of the boundary elements with no face in the boundary and layer elements

        for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
        {
            Element::GeometryType& pGeom = ie->GetGeometry();
            for(unsigned int i = 0; i < pGeom.size(); i++)
            {
                if(pGeom[i].Is(BOUNDARY))
                {
                    ie->Set(BOUNDARY);
                }
            }

        }

        return true;

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


