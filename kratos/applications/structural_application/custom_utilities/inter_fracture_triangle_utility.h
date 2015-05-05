/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */
/* *********************************************************
 *
 *   Last Modified by:    $Author: Nelson Lafontaine $
 *   Date:                $Date: 01-27-2010$
 *   Revision:            $Revision: 1.00   $
 *
 * ***********************************************************/

#if !defined(INTER_FRACTURE_TRIANGLE_UTILITY)
#define INTER_FRACTURE_TRIANGLE_UTILITY

#ifdef _OPENMP
#include <omp.h>
#endif

#include "boost/smart_ptr.hpp"
#include <boost/timer.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>



// System includes
#include <string>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <algorithm>



/* Project includes */
#include "structural_application.h"

#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "processes/find_nodal_neighbours_process.h"
#include "processes/find_elements_neighbours_process.h"
#include "processes/find_conditions_neighbours_process.h"
#include "containers/data_value_container.h"
#include "includes/mesh.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/smoothing_utility.h"
#include "geometries/triangle_2d_3.h"
#include "processes/node_erase_process.h"
#include "spatial_containers/spatial_containers.h"
#include "custom_utilities/joint.h"
#include "custom_utilities/disconnect_utility.h"

namespace Kratos
{

class Inter_Fracture_Triangle
{
public:


    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    typedef boost::numeric::ublas::vector<Matrix> Matrix_Order_Tensor;
    typedef boost::numeric::ublas::vector<Vector> Vector_Order_Tensor;
    typedef boost::numeric::ublas::vector<Vector_Order_Tensor> Node_Vector_Order_Tensor;
    typedef Node < 3 > PointType;
    typedef Node < 3 > ::Pointer PointPointerType;
    typedef std::vector<PointType::Pointer> PointVector;
    typedef PointVector::iterator PointIterator;
    typedef Joint<4> Joint2D;

    Inter_Fracture_Triangle(ModelPart& model_part, int domain_size) : mr_model_part(model_part)
    {
        mdomain_size  = domain_size;
        mjoint_create = false;
    }

    ~Inter_Fracture_Triangle()
    {
    }



    ///************************************************************************************************
    ///************************************************************************************************

    bool Detect_And_Split_Elements_Heuristic_Formula(ModelPart& this_model_part)
    {
        //Detect_Node_To_Be_Splitted_Heuristic_Formula(mDTU.Begin(), mDTU.End());
        return false;
    }


    bool Detect_And_Split_Elements(ModelPart& this_model_part)
    {
        KRATOS_TRY

        bool is_split = false;
        array_1d<double, 3 > Failure_Maps;
        WeakPointerVector< Node < 3 > > Nodes_To_Be_Dupplicated;

        unsigned int detect = Detect_Node_To_Be_Splitted(this_model_part, Nodes_To_Be_Dupplicated);
        unsigned  int count = 0;
        if (detect != 0)
        {
            WeakPointerVector< Node < 3 > >::iterator i_begin = Nodes_To_Be_Dupplicated.ptr_begin();
            WeakPointerVector< Node < 3 > >::iterator i_end   = Nodes_To_Be_Dupplicated.ptr_end();

            for(WeakPointerVector< Node < 3 > >::iterator inode = i_begin; inode != i_end; ++inode)
            {
                const Node<3>::Pointer& pNode = (*(inode.base())).lock();
                if(Split_Node(this_model_part, pNode))
                    count++;
            }

            bool split_2             = false;
            unsigned int nodes_added = Finalize(this_model_part, split_2);
            if((count + nodes_added)!=0)
            {
                is_split = true;
                std::cout<< "NUMBER OF NEW NODES CREATED FOR FRACTURING            = "<< count << std::endl;
                std::cout<< "NUMBER OF NEW NODES CREATED FOR SPECIAL SITUATIONS    = "<< nodes_added << std::endl;
                std::cout<< "NUMBER TOTAL OF NEW NODES CREATED                     = "<< count + nodes_added << std::endl;
                //InitializeElementsAndVariables();
                RecomputeNodalMass();
            }
            else
                std::cout<< "NO NEW NODES TO BE CREATED " << std::endl;


        }
        else
        {
            unsigned int zero         = 0;
            unsigned int nodes_added = Finalize(this_model_part, is_split);
            if(is_split==true)
            {
                RecomputeNodalMass();
                std::cout<< "NUMBER OF NEW NODES CREATED FOR FRACTURING            = "<< zero << std::endl;
                std::cout<< "NUMBER OF NEW NODES CREATED FOR SPECIAL SITUATIONS    = "<< nodes_added << std::endl;
                std::cout<< "NUMBER TOTAL OF NEW NODES CREATED                     = "<< zero + nodes_added << std::endl;
            }
            else
                std::cout<< "NO NEW NODES TO BE CREATED " << std::endl;
        }

        CheckAloneElement();
        return is_split;
        KRATOS_CATCH("")
    }


    ///************************************************************************************************
    ///************************************************************************************************



    unsigned int Detect_Node_To_Be_Splitted_Heuristic_Formula(std::vector<Joint2D>::iterator Begin, std::vector<Joint2D>::iterator End)
    {
        KRATOS_TRY
        double dpefa = 0.63;
        double dpefb = 1.8;
        double dpefc = 6.0;
        double dpefm = 0.0;
        double small,sabs,o,s,o1,o2,s1,s2,op,sp,ot,st,z,sigma,tau;
        double e1x,e1y,h,area;
        //int ielem,i0,i1,i2,i3l,m;
        int integ, nfail;
        int nsoft;
        double d1nccx[4];
        double d1nccy[4];

        double dpeft = 3.15e+06;
        double dpepe = 100E+9;
        double dpefs = 3.15e+06;
        double dpegf = 3.0e+2;

        small=1E-6;
        nsoft=0;

        for(std::vector<Joint2D>::iterator Joint = Begin; Joint != End; Joint++)
        {

            array_1d<double,3>& node_rhs_0 =  (*Joint)[0]->FastGetSolutionStepValue(RHS);
            array_1d<double,3>& node_rhs_1 =  (*Joint)[1]->FastGetSolutionStepValue(RHS);
            array_1d<double,3>& node_rhs_2 =  (*Joint)[2]->FastGetSolutionStepValue(RHS);
            array_1d<double,3>& node_rhs_3 =  (*Joint)[3]->FastGetSolutionStepValue(RHS);

            d1nccx[0] = (*Joint)[0]->X();
            d1nccx[1] = (*Joint)[1]->X();
            d1nccx[2] = (*Joint)[2]->X();
            d1nccx[3] = (*Joint)[3]->X();

            d1nccy[0] = (*Joint)[0]->Y();
            d1nccy[1] = (*Joint)[1]->Y();
            d1nccy[2] = (*Joint)[2]->Y();
            d1nccy[3] = (*Joint)[3]->Y();

            e1x=0.50*(d1nccx[1]+d1nccx[2]-d1nccx[0]-d1nccx[3]);
            e1y=0.50*(d1nccy[1]+d1nccy[2]-d1nccy[0]-d1nccy[3]);
            h=std::sqrt(e1x*e1x+e1y*e1y);

            e1x=e1x/(h+small);
            e1y=e1y/(h+small);
            s1=(d1nccy[0]-d1nccy[3])*e1y+(d1nccx[0]-d1nccx[3])*e1x;
            s2=(d1nccy[1]-d1nccy[2])*e1y+(d1nccx[1]-d1nccx[2])*e1x;
            o1=(d1nccy[0]-d1nccy[3])*e1x-(d1nccx[0]-d1nccx[3])*e1y;
            o2=(d1nccy[1]-d1nccy[2])*e1x-(d1nccx[1]-d1nccx[2])*e1y;


            op=2.00*h*dpeft/dpepe;
            sp=2.00*h*dpefs/dpepe;
            ot=std::max((2.00*op),(3.00*dpegf/dpeft));
            st=std::max((2.00*sp),(3.00*dpegf/dpefs));
            nfail=0;

            for(integ=0; integ<3; integ++)
            {
                if(integ==0)
                {
                    o=o1;
                    s=s1;
                }
                else if(integ==2)
                {
                    o=o2;
                    s=s2;
                }
                else
                {
                    o=0.50*(o1+o2);
                    s=0.50*(s1+s2);
                }

                sabs=std::fabs(s);
                if((o>op)&&(sabs>sp))
                {
                    z=std::sqrt(((o-op)/ot)*((o-op)/ot)+((sabs-sp)/st)*((sabs-sp)/st));
                }
                else if(o>op)
                {
                    z=(o-op)/ot;
                }
                else if(sabs>sp)
                {
                    z=(sabs-sp)/st;
                }
                else
                {
                    z=0.00;
                }

                if(z>=1.00)
                {
                    nfail=nfail+1;
                    z=1.00;
                }
                z=(1.00 - ((dpefa+dpefb-1.00)/(dpefa+dpefb))*exp(z*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(1.00-dpefa-dpefb))))*(dpefa*(1.00-z)+dpefb*pow((1.00-z),dpefc));


                if(o<0.00)                 /* normal stress*/
                {
                    sigma=2.00*o*dpeft/op;   /* sigma=R0; */
                }
                else if(o>op)
                {
                    sigma=dpeft*z;
                    nsoft=nsoft+1;
                }
                else
                {
                    sigma=(2.00*o/op-(o/op)*(o/op))*z*dpeft;
                }
                if((sigma>0.00)&&(sabs>sp))           /* shear stress */
                {
                    tau=z*dpefs;
                }
                else if(sigma>0.00)
                {
                    tau=(2.00*(sabs/sp)-(sabs/sp)*(sabs/sp))*z*dpefs;
                }
                else if(sabs>sp)
                {
                    tau=z*dpefs-dpefm*sigma;
                }
                else
                {
                    tau=(2.00*(sabs/sp)-(sabs/sp)*(sabs/sp))*(z*dpefs-dpefm*sigma);
                }
                if(s<0.00)tau=-tau;
                if(integ==0)  /* nodal forces */
                {
                    area=h/6.00; /* area=h/6.0; */
                    node_rhs_0[0] = node_rhs_0[0] - area*(tau*e1x-sigma*e1y);
                    node_rhs_3[0] = node_rhs_3[0] + area*(tau*e1x-sigma*e1y);
                    node_rhs_0[1] = node_rhs_0[0] - area*(tau*e1y+sigma*e1x);
                    node_rhs_3[1] = node_rhs_3[1] + area*(tau*e1y+sigma*e1x);
                }
                else if(integ==1)
                {
                    area=h/3.00;  /* area=h/3.0; */

                    node_rhs_0[0] = node_rhs_0[0] -area*(tau*e1x-sigma*e1y);
                    node_rhs_3[0] = node_rhs_0[3] +area*(tau*e1x-sigma*e1y);
                    node_rhs_0[1] = node_rhs_0[1]-area*(tau*e1y+sigma*e1x);
                    node_rhs_3[1] = node_rhs_3[1]+area*(tau*e1y+sigma*e1x);
                    node_rhs_1[0] = node_rhs_1[0]-area*(tau*e1x-sigma*e1y);
                    node_rhs_2[0] = node_rhs_2[0]+area*(tau*e1x-sigma*e1y);
                    node_rhs_1[1] = node_rhs_1[1]-area*(tau*e1y+sigma*e1x);
                    node_rhs_2[1] = node_rhs_2[1]+area*(tau*e1y+sigma*e1x);
                }
                else
                {
                    area=h/6.00; /* area=h/6.0; */
                    node_rhs_1[0]=node_rhs_1[0]-area*(tau*e1x-sigma*e1y);
                    node_rhs_2[0]=node_rhs_2[0]+area*(tau*e1x-sigma*e1y);
                    node_rhs_1[1]=node_rhs_1[1]-area*(tau*e1y+sigma*e1x);
                    node_rhs_2[1]=node_rhs_2[1]+area*(tau*e1y+sigma*e1x);
                }
            }
        }
        return 0;
        KRATOS_CATCH("")
    }


    ///************************************************************************************************
    ///************************************************************************************************

    unsigned int Detect_Node_To_Be_Splitted(ModelPart& this_model_part, WeakPointerVector< Node < 3 > >& Nodes_To_Be_Dupplicated)
    {
        KRATOS_TRY

        NodesArrayType& pNodes           =  this_model_part.Nodes();
        NodesArrayType::iterator i_begin =  pNodes.ptr_begin();
        NodesArrayType::iterator i_end   =  pNodes.ptr_end();

        for(ModelPart::NodeIterator inode = i_begin; inode != i_end; ++inode)
        {
            double& Condition = inode->GetValue(NODAL_DAMAGE);
            if (Condition > 0.50)
            {
                Nodes_To_Be_Dupplicated.push_back(*(inode.base()));
            }
        }

        return Nodes_To_Be_Dupplicated.size();

        KRATOS_CATCH("")

    }

    ///************************************************************************************************
    ///************************************************************************************************
    ///* Computa el vector, linea, o plano de falla.

    void Calculate_Map_Failure(const Node<3>::Pointer& pNode, array_1d<double, 3 > & Failure_Maps)
    {

        Vector Eigen_Values = ZeroVector(mdomain_size); // Deformaciones o Tensiones principales
        Matrix Eigen_Vector = ZeroMatrix(mdomain_size, mdomain_size); // Direcciones  principales
        Matrix Strain_Tensor = ZeroMatrix(mdomain_size, mdomain_size);

        unsigned int size = 3;
        unsigned int max_iterations = 100;
        double zero_tolerance = 1e-9;
        if (mdomain_size == 3)
        {
            size = 6;
        }


        Vector_Order_Tensor EigenVector;

        EigenVector.resize(3);
        EigenVector[0] = ZeroVector(3);
        EigenVector[1] = ZeroVector(3);
        EigenVector[2] = ZeroVector(3);


        Matrix& Nodal_Values = pNode->GetValue(NODAL_STRAIN);
        Vector temp  = ZeroVector(size);
        for (unsigned int j = 0; j < size; j++)
        {
            temp[j] = Nodal_Values(0, j);
        }

        Strain_Tensor = SD_MathUtils<double>::StrainVectorToTensor(temp);
        SD_MathUtils<double>::EigenVectors(Strain_Tensor, Eigen_Vector, Eigen_Values, zero_tolerance, max_iterations);

        // traccion principal de traccion
        double max = (*std::max_element(Eigen_Values.begin(), Eigen_Values.end()));
        int pos = 0;

        for (unsigned int k = 0; k < Eigen_Values.size(); k++)
        {
            if (max == Eigen_Values(k))
            {
                pos = k;
            }
            for (unsigned int j = 0; j < Eigen_Values.size(); j++)
            {
                EigenVector[k][j] = Eigen_Vector(k, j);
            }
        }

        //Tomo el primer egenvector correspondiente al vector que es perpendicaular al vector de traccion principal mayor
        if (pos == 0)
        {
            pos = 1;
        }
        else
            pos = 0;

        Failure_Maps[0] = EigenVector[pos][0];
        Failure_Maps[1] = EigenVector[pos][1];
        Failure_Maps[2] = EigenVector[pos][2];

    }


    ///************************************************************************************************
    ///************************************************************************************************

    bool Split_Node(ModelPart& this_model_part, const Node<3>::Pointer& pNode)
    {

        KRATOS_TRY

        array_1d<double, 3 > failure_map;
        Calculate_Map_Failure(pNode, failure_map);



        WeakPointerVector< Element > Negative_Elements;
        WeakPointerVector< Element > Positive_Elements;


        Node < 3 > ::Pointer child_node;
        bool node_created = false;

        /// Crea el nuevo nodo y lo inserta al elemento
        node_created = CalculateElements(this_model_part, pNode, child_node, failure_map, Negative_Elements, Positive_Elements);


        if (node_created == true)
        {
            ///UPDATING CONDITIONS
            /// En caso de que el nodo creado modifique la condicion de contorno
            CalculateConditions(this_model_part, pNode, child_node, failure_map);

            /// Updating vecinos
            RecomputeLocalneighbourgs(pNode, Positive_Elements);
            RecomputeLocalneighbourgs(child_node, Negative_Elements);

        }

        return node_created;
        KRATOS_CATCH("")


    }


    ///************************************************************************************************
    ///************************************************************************************************

    void RecomputeLocalneighbourgs(const Node <3>::Pointer& pNode, WeakPointerVector< Element>& Elements)
    {
        KRATOS_TRY

        WeakPointerVector< Element    >& neighb_elems  = pNode->GetValue(NEIGHBOUR_ELEMENTS);
        WeakPointerVector< Node < 3 > >& neighb_nodes  = pNode->GetValue(NEIGHBOUR_NODES);
        //WeakPointerVector< Condition > & neighb_conds  = pNode->GetValue(NEIGHBOUR_CONDITIONS);
        neighb_elems.clear();
        neighb_nodes.clear();
        //neighb_nodes.conds();

        /// nodos y elemntos vecinos de un nodo
        for (WeakPointerVector<Element>::iterator elem = Elements.begin(); elem != Elements.end(); ++elem)
        {
            neighb_elems.push_back(*elem.base());
            Element::GeometryType& geom = elem->GetGeometry();
            for (unsigned int i = 0; i < geom.size(); i++)
            {
                WeakPointerVector< Node < 3 > >::iterator repeated_object = std::find(neighb_nodes.begin(), neighb_nodes.end(), geom[i]);
                if (repeated_object == (neighb_nodes.end()))
                {
                    neighb_nodes.push_back(geom(i));
                }
            }
        }
        //elementos vecinos a vecinos
        for (WeakPointerVector<Element>::iterator elem = Elements.begin(); elem != Elements.end(); ++elem)
        {
            Geometry<Node<3> >& geom = (elem)->GetGeometry();
            (elem->GetValue(NEIGHBOUR_ELEMENTS)).resize(3);
            WeakPointerVector< Element >& neighb_elems = elem->GetValue(NEIGHBOUR_ELEMENTS);
            //neighb_face is the vector containing pointers to the three faces around ic
            //neighb_face[0] = neighbour face over edge 1-2 of element ic;
            //neighb_face[1] = neighbour face over edge 2-0 of element ic;
            //neighb_face[2] = neighbour face over edge 0-1 of element ic;
            neighb_elems(0) = CheckForNeighbourElems(geom[1].Id(), geom[2].Id(), geom[1].GetValue(NEIGHBOUR_ELEMENTS), elem);
            neighb_elems(1) = CheckForNeighbourElems(geom[2].Id(), geom[0].Id(), geom[2].GetValue(NEIGHBOUR_ELEMENTS), elem);
            neighb_elems(2) = CheckForNeighbourElems(geom[0].Id(), geom[1].Id(), geom[0].GetValue(NEIGHBOUR_ELEMENTS), elem);
        }

        KRATOS_CATCH("")

    }

    /// Calculo de los elementos vecinos a un elemento
    Element::WeakPointer CheckForNeighbourElems (unsigned int Id_1, unsigned int Id_2, WeakPointerVector< Element >& neighbour_elem, WeakPointerVector<Element>::iterator elem)
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


    ///************************************************************************************************
    ///************************************************************************************************

    //pNode = the father node
    bool CalculateElements(ModelPart& this_model_part,
                           const Node<3>::Pointer& pNode, // Nodo Padre
                           Node<3>::Pointer& pnode, // Nodo Hijo
                           const array_1d<double, 3 > & failure_map,
                           WeakPointerVector< Element >& Negative_Elements,
                           WeakPointerVector< Element >& Positive_Elements
                          )
    {
        KRATOS_TRY

        double prod = 0.00;
        array_1d<double, 3 > Coord_Point_1;
        array_1d<double, 3 > Coord_Point_2;
        array_1d<double, 3 > normal;
        array_1d<double, 3 > Unit = ZeroVector(3);


        Positive_Elements.reserve(10);
        Negative_Elements.reserve(10);


        WeakPointerVector< Element >& neighb_elems = pNode->GetValue(NEIGHBOUR_ELEMENTS);

        // Normal al plano de fractura
        Calculate_Normal_Faliure_Maps(normal, failure_map);

        noalias(Coord_Point_1) = pNode->Coordinates();
        //Coord_Point_1[1]  = pNode->Y();
        //Coord_Point_1[2]  = pNode->Z();


        for (WeakPointerVector< Element >::iterator neighb_elem = neighb_elems.begin();
                neighb_elem != neighb_elems.end(); neighb_elem++)
        {
            Element::GeometryType& geom = neighb_elem->GetGeometry(); // Nodos del elemento
            noalias(Coord_Point_2) = geom.Center();
            //Find_Coord_Gauss_Points(geom, Coord_Point_2);
            noalias(Unit) = Coord_Point_2 - Coord_Point_1;
            noalias(Unit) = Unit / norm_2(Unit);

            prod = inner_prod(normal, Unit);
            if (prod >= 0.00)
            {
                //std::cout<<"INSERTED_POSITIVE"<<std::endl;
                Positive_Elements.push_back(*(neighb_elem.base()));
            }
            else
            {
                //std::cout<<"INSERTED_NEGATIVE"<std::endl;
                Negative_Elements.push_back(*(neighb_elem.base()));
            }

            Unit = ZeroVector(3);
            prod = 0.00;
        }

        //* Esquinas o superficie extrena donde no hay elementos negativos o positivos
        if (Positive_Elements.size() == 0 || Negative_Elements.size() == 0)
        {
            Positive_Elements.clear();
            Negative_Elements.clear();
            for (WeakPointerVector< Element >::iterator neighb_elem = neighb_elems.begin();
                    neighb_elem != neighb_elems.end(); neighb_elem++)
            {
                Element::GeometryType& geom = neighb_elem->GetGeometry(); // Nodos del elemento
                noalias(Coord_Point_2)      = geom.Center();
                noalias(Unit)               = Coord_Point_2 - Coord_Point_1;
                noalias(Unit)               = Unit / norm_2(Unit);
                if (inner_prod(failure_map, Unit) >= 0.00)
                {
                    Positive_Elements.push_back(*(neighb_elem.base()));
                }
                else
                {
                    Negative_Elements.push_back(*(neighb_elem.base()));
                }
            }
        }


        //* Si el nodo no mas tiene un solo elemtno vecino
        bool& duplicated_pNode = pNode->GetValue(IS_DUPLICATED);
        if ((Positive_Elements.size() == 0 || Negative_Elements.size() == 0) || (Positive_Elements.size() == 1 && Negative_Elements.size() == 0) || (Positive_Elements.size() == 0 && Negative_Elements.size() == 1) || duplicated_pNode == true)
        {
            //std::cout<<"NO INSERTED NODE"<<std::endl;
            return false;
        }
        else
        {


            for (WeakPointerVector< Element >::iterator neighb_elem = neighb_elems.begin();
                    neighb_elem != neighb_elems.end(); neighb_elem++)
            {
                mResetingElements.push_back(*neighb_elem.base());
            }

            unsigned int New_Id = this_model_part.Nodes().size() + 1;
            //Node<3>::Pointer pnode; // the new node
            bool& duplicated_pnode = pNode->GetValue(IS_DUPLICATED);
            duplicated_pnode = true;
            Create_New_Node(this_model_part, failure_map, pnode, New_Id, pNode);

            //* putting de new node to the negative element
            for (WeakPointerVector< Element >::iterator neg_elem = Negative_Elements.begin();
                    neg_elem != Negative_Elements.end(); neg_elem++)
            {
                Element::GeometryType& geom = neg_elem->GetGeometry();
                Insert_New_Node_In_Elements(geom, pnode, pNode->Id());
            }
            return true;
        }

        KRATOS_CATCH("")

    }





    ///************************************************************************************************
    ///************************************************************************************************

    ///Resetea los elementos alrededor de la fratura
    void InitializeElementsAndVariables()
    {
        KRATOS_TRY

        //ElementsArrayType& pElements     = mr_model_part.Elements();
        ProcessInfo& CurrentProcessInfo =  mr_model_part.GetProcessInfo();
        for(WeakPointerVector< Element >::iterator reset_elem = mResetingElements.begin();
                reset_elem != mResetingElements.end(); reset_elem++)
        {
            reset_elem->Initialize();
            reset_elem->InitializeSolutionStep(CurrentProcessInfo);
            reset_elem->FinalizeSolutionStep(CurrentProcessInfo);
        }

        if(mResetingElements.size()!=0)
            mResetingElements.clear();

        KRATOS_CATCH("")
    }



    void CheckAloneElement()
    {

        ProcessInfo& CurrentProcessInfo =  mr_model_part.GetProcessInfo();
        ElementsArrayType& pElements     = mr_model_part.Elements();

#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif



        std::cout<< "Element to be reseted" << std::endl;

        std::vector<double> Variable_Value;
        vector<unsigned int> element_partition;
        CreatePartition(number_of_threads, pElements.size(), element_partition);
        //unsigned int index = 0;
        #pragma omp parallel for private(Variable_Value)
        for(int k=0; k<number_of_threads; k++)
        {
            ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
            ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
            for(ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
            {
                WeakPointerVector<Element>& Neighb_Elem =  it->GetValue(NEIGHBOUR_ELEMENTS);
                if(Neighb_Elem(0).lock()->Id()==it->Id())
                    if(Neighb_Elem(1).lock()->Id()==it->Id())
                        if(Neighb_Elem(2).lock()->Id()==it->Id())
                        {
                            it->GetValueOnIntegrationPoints(DAMAGE, Variable_Value, CurrentProcessInfo);
                            if(Variable_Value[0]>0.00)
                            {
                                KRATOS_WATCH(it->Id() )
                                KRATOS_WATCH(Variable_Value[0])
                                it->ResetConstitutiveLaw();
                                it->Initialize();
                                it->InitializeSolutionStep(CurrentProcessInfo);
                                it->FinalizeSolutionStep(CurrentProcessInfo);
                            }
                        }
            }
        }
    }

    ///************************************************************************************************
    ///************************************************************************************************
    void RecomputeNodalMass()
    {

        KRATOS_TRY
        ProcessInfo& CurrentProcessInfo  = mr_model_part.GetProcessInfo();
        NodesArrayType& pNodes           = mr_model_part.Nodes();
        ElementsArrayType& pElements     = mr_model_part.Elements();


        Matrix MassMatrix;
#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif

        vector<unsigned int> node_partition;
        CreatePartition(number_of_threads, pNodes.size(), node_partition);

        #pragma omp parallel for
        for(int k=0; k<number_of_threads; k++)
        {
            NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
            NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

            for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
            {
                i->FastGetSolutionStepValue(NODAL_MASS) = 0.00;

            }
        }


        vector<unsigned int> element_partition;
        CreatePartition(number_of_threads, pElements.size(), element_partition);

        #pragma omp parallel for private(MassMatrix)
        for(int k=0; k<number_of_threads; k++)
        {
            ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
            ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
            for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
            {
                Element::GeometryType& geom = it->GetGeometry();
                (it)->CalculateMassMatrix(MassMatrix, CurrentProcessInfo);
                unsigned int dim   = geom.WorkingSpaceDimension();
                unsigned int index = 0;
                for (unsigned int i = 0; i <geom.size(); i++)
                {
                    geom[i].SetLock();
                    index = i*dim;
                    geom[i].FastGetSolutionStepValue(NODAL_MASS) += MassMatrix(index,index);
                    geom[i].UnSetLock();
                }
            }
        }
        mr_model_part.GetCommunicator().AssembleCurrentData(NODAL_MASS);

        KRATOS_CATCH("")
    }


    ///************************************************************************************************
    ///************************************************************************************************

    void Create_New_Node(ModelPart& this_model_part, const array_1d<double, 3 > & failure_map, Node < 3 > ::Pointer& pnode, unsigned int& New_Id, const Node<3>::Pointer& pNode)
    {


        //node to get the DOFs from
        int step_data_size   = this_model_part.GetNodalSolutionStepDataSize();
        const double ancho_w = 1E-9;
        array_1d<double, 3 > Aux = -1.00 * (ancho_w * failure_map);

        array_1d<double, 3 > & Coord = pNode->Coordinates();
        noalias(Coord) += Aux;
        Node < 3 > ::DofsContainerType& reference_dofs = (this_model_part.NodesBegin())->GetDofs();

        pnode = this_model_part.CreateNewNode(New_Id, Coord[0], Coord[1], Coord[2]);
        pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize());

        //generating the dofs
        for (Node < 3 > ::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
        {
            Node < 3 > ::DofType& rDof = *iii;
            Node < 3 > ::DofType::Pointer p_new_dof = pnode->pAddDof(rDof);
            if (pNode->IsFixed(iii->GetVariable()) == true)
                (p_new_dof)->FixDof();
            else
            {
                (p_new_dof)->FreeDof();
            }
        }


        unsigned int buffer_size = pnode->GetBufferSize();
        for (unsigned int step = 0; step < buffer_size; step++)
        {
            double* step_data = pnode->SolutionStepData().Data(step);
            double* old_node_data = pNode->SolutionStepData().Data(step);
            //copying this data in the position of the vector we are interested in
            for (signed int j = 0; j < step_data_size; j++)
            {
                step_data[j] = old_node_data[j];
            }
        }


        const array_1d<double, 3 > & disp = pnode->FastGetSolutionStepValue(DISPLACEMENT);
        pnode->X0() = pnode->X() - disp[0];
        pnode->Y0() = pnode->Y() - disp[1];
        pnode->Z0() = pnode->Z() - disp[2];

        const array_1d<double, 3 > & vel_old = pNode->FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3 > & vel_new = pnode->FastGetSolutionStepValue(VELOCITY);
        vel_new = vel_old;

        const array_1d<double, 3 > & accel_old = pNode->FastGetSolutionStepValue(ACCELERATION);
        array_1d<double, 3 > & accel_new = pnode->FastGetSolutionStepValue(ACCELERATION);
        accel_new = accel_old;

    }


    ///************************************************************************************************
    ///************************************************************************************************

    //pNode = the father node
    void CalculateConditions(ModelPart& this_model_part,
                             const Node<3>::Pointer& pNode, // father node
                             Node < 3 > ::Pointer& pnode, // child node
                             const array_1d<double, 3 > & failure_map)
    {
        KRATOS_TRY


        //double prod = 0.00;
        //array_1d<double, 3 > Coord_Point_1;
        //array_1d<double, 3 > Coord_Point_2;
        //array_1d<double, 3 > normal;
        //array_1d<double, 3 > Unit;

        //WeakPointerVector< Condition > Negative_Conditions;
        //Negative_Conditions.reserve(10);
        int j = 0;
        Condition::Pointer rcond;
        WeakPointerVector< Condition >& neighb_conds = pNode->GetValue(NEIGHBOUR_CONDITIONS);
        for(WeakPointerVector< Condition >::iterator neighb_cond = neighb_conds.begin(); neighb_cond != neighb_conds.end(); neighb_cond++)
        {
            if(neighb_cond->GetValue(IS_CONTACT_MASTER)==0)
            {
                rcond = (neighb_conds(j).lock());
                Insert_New_Node_In_Conditions(rcond);
            }
            j++;
        }

        /*
            bool check = CkeckNumberConditionsSurfaces(pNode);
        if(check){
            // Normal al plano de fractura
            Calculate_Normal_Faliure_Maps(normal, failure_map);
            Coord_Point_1 = pNode->Coordinates();

            for (WeakPointerVector< Condition >::iterator neighb_cond = neighb_conds.begin();
                    neighb_cond != neighb_conds.end(); neighb_cond++)
            {

                Condition::GeometryType& geom = neighb_cond->GetGeometry(); // Nodos de las condiciones
                noalias(Coord_Point_2) = geom.Center();
                noalias(Unit) = Coord_Point_2 - Coord_Point_1;
                noalias(Unit) = Unit / norm_2(Unit);
                prod = inner_prod(normal, Unit);
                if (prod < 0.00)
                {
                    Negative_Conditions.push_back(*(neighb_cond.base()));
                }
                Unit = ZeroVector(3);
                prod = 0.00;
            }
            */
        /* Esquinas o superficie extrena donde no hay elementos negativos o positivos */
        /*
            if(Negative_Conditions.size() == 0)
            {
               for (WeakPointerVector< Condition >::iterator neighb_cond = neighb_conds.begin();
                    neighb_cond != neighb_conds.end(); neighb_cond++)
                {

              Condition::GeometryType& geom = neighb_cond->GetGeometry(); // Nodos de las condiciones
              noalias(Coord_Point_2) = geom.Center();
              noalias(Unit) = Coord_Point_2 - Coord_Point_1;
              noalias(Unit) = Unit / norm_2(Unit);
              prod = inner_prod(failure_map, Unit);
              if (prod < 0.00)
              {
        	  Negative_Conditions.push_back(*(neighb_cond.base()));
              }
              Unit = ZeroVector(3);
              prod = 0.00;

                }
              }
        }
            else
        {
          bool exist = false;
          int j =  DetectCondition(pNode, exist);
          KRATOS_WATCH(j)
          if(exist)
          {
        Condition::Pointer rcond = (neighb_conds(j).lock());
        Insert_New_Node_In_Conditions(rcond);

        int a        =  0, b = 1;
        double toler = 1E-8;
        double check = 0.00;
            Condition::GeometryType& geom_cond  = neighb_conds(j).lock()->GetGeometry();
            KRATOS_WATCH(neighb_conds(j).lock()->Id())
        KRATOS_WATCH(geom_cond(0)->Id())
        KRATOS_WATCH(geom_cond(1)->Id())
        KRATOS_WATCH((neighb_conds(j).lock()->GetValue(NEIGHBOUR_ELEMENTS)).size())

            Element::Pointer relem              = (neighb_conds(j).lock()->GetValue(NEIGHBOUR_ELEMENTS))(0).lock();
        KRATOS_WATCH("2222222222")
        Element::GeometryType& geom_elem    = relem->GetGeometry();
        KRATOS_WATCH("33333333")
        array_1d<double,3> cond             = geom_cond.GetPoint(1) - geom_cond.GetPoint(0);
        KRATOS_WATCH("ONEEEEEE")

        array_1d<double,3> elem;
        noalias(cond) = (1.00/norm_2(cond)) * cond;

        for(unsigned int i = 0; i<geom_elem.size(); i++)
        {
          if(i==2) {b = 0; a = 2;}
          elem          = geom_elem.GetPoint(b) - geom_elem.GetPoint(a);
          noalias(elem) = (1.00/norm_2(elem)) * elem;
          check         = std::fabs(inner_prod(elem, cond));
          if( std::fabs(check-1.00)<toler )
             break;
          b++; a++;
        }

         KRATOS_WATCH("TWOOOOOO")
         KRATOS_WATCH(a)
         KRATOS_WATCH(b)

         geom_cond(0) = geom_elem(a);
         geom_cond(1) = geom_elem(b);
         KRATOS_WATCH(pNode->Id())
                 KRATOS_WATCH(j)

         //KRATOS_THROW_ERROR(std::logic_error,  "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" , "");

          }

        }

            if(Negative_Conditions.size()!=0){
            /// putting de new node to the negative conditions
            for (WeakPointerVector< Condition >::iterator neg_cond = Negative_Conditions.begin();
                    neg_cond != Negative_Conditions.end(); neg_cond++)
                       {
        	 Condition::GeometryType& geom = neg_cond->GetGeometry();
                         Insert_New_Node_In_Conditions(geom, pnode, pNode->Id());
                       }
           }

        */
        KRATOS_CATCH("")

    }


    void Insert_New_Node_In_Conditions(Condition::Pointer& rcond)
    {

        //WeakPointerVector<Element>& neighb_elem = rcond->GetValue(NEIGHBOUR_ELEMENTS);
        int a        =  0, b = 1;
        const double toler = 1E-8;
        double check = 0.00;
        a = 0;
        b = 1;
        Condition::GeometryType& geom_cond      = rcond->GetGeometry();
        Element::Pointer relem                  = (rcond->GetValue(NEIGHBOUR_ELEMENTS))(0).lock();
//		int id = rcond->Id();
// 		if(id==1 || id==2 || id==3 || id==4 || id==5 || id==6 || id==7 || id==8)
// 		{
// 		  KRATOS_WATCH(rcond->Id())
// 		  KRATOS_WATCH(relem->Id())
// 		}

        Element::GeometryType& geom_elem        = relem->GetGeometry();
        array_1d<double,3> cond                 = geom_cond.GetPoint(1) - geom_cond.GetPoint(0);
        array_1d<double,3> elem;
        noalias(cond) = (1.00/norm_2(cond)) * cond;
        for(unsigned int i = 0; i<geom_elem.size(); i++)
        {
            if(i==2)
            {
                b = 0;
                a = 2;
            }
            elem          = geom_elem.GetPoint(b) - geom_elem.GetPoint(a);
            noalias(elem) = (1.00/norm_2(elem)) * elem;
            check         = std::fabs(inner_prod(elem, cond));
            if( std::fabs(check-1.00)<toler )
                break;
            b++;
            a++;
        }
        /// Las condiciones se nombran contrario a las manecillas del reloj
        geom_cond(0) = geom_elem(b);
        geom_cond(1) = geom_elem(a);

// 		 KRATOS_WATCH(geom_cond(0)->Id())
// 		 KRATOS_WATCH(geom_cond(1)->Id())
// 		 KRATOS_WATCH("------------------")
    }

    ///************************************************************************************************
    ///************************************************************************************************

    inline void Insert_New_Node_In_Elements(Element::GeometryType& geom, Node < 3 > ::Pointer& pnode, const unsigned int& Node_Id_Old)
    {
        KRATOS_TRY
        for (unsigned int i = 0; i < geom.size(); i++)
            if (geom[i].Id() == Node_Id_Old)
                geom(i) = pnode;
        KRATOS_CATCH("")

    }

    ///************************************************************************************************
    ///************************************************************************************************

    inline void Insert_New_Node_In_Conditions(Condition::GeometryType& geom, Node < 3 > ::Pointer& pnode, const unsigned int& Node_Id_Old)
    {
        KRATOS_TRY
        for (unsigned int i = 0; i < geom.size(); i++)
            if(geom[i].Id() == Node_Id_Old)
                geom(i) = pnode;
        KRATOS_CATCH("")

    }


    ///************************************************************************************************
    ///************************************************************************************************

    inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
    {
        partitions.resize(number_of_threads + 1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for (unsigned int i = 1; i < number_of_threads; i++)
            partitions[i] = partitions[i - 1] + partition_size;
    }

    ///************************************************************************************************
    ///************************************************************************************************

    inline void Find_Coord_Gauss_Points(Element::GeometryType& geom, array_1d<double, 3 > & Coord_Point)
    {
        double x = 0.00;
        double y = 0.00;
        double z = 0.00;
        double fact = 0.33333333333333333333333;

        if (geom.size() == 4)
        {
            fact = 0.25;
        }
        Coord_Point = ZeroVector(3);
        for (unsigned int i = 0; i < geom.size(); i++)
        {

            x = geom[i].X();
            y = geom[i].Y();
            z = geom[i].Z();

            Coord_Point[0] += x;
            Coord_Point[1] += y;
            Coord_Point[2] += z;
        }

        noalias(Coord_Point) = Coord_Point*fact;

    }

    unsigned int Finalize(ModelPart& this_model_part, bool& split)
    {

        KRATOS_TRY

        std::cout<< "FINALIZE BEGIN"<< std::endl;

        unsigned int result              = 0;
        NodesArrayType& pNodes           = this_model_part.Nodes();
        ElementsArrayType& pElements     = this_model_part.Elements();
        //ConditionsArrayType& pConditions = this_model_part.Conditions();

#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif

        vector<unsigned int> element_partition;
        CreatePartition(number_of_threads, pElements.size(), element_partition);
        #pragma omp parallel for
        for(int k=0; k<number_of_threads; k++)
        {
            ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
            ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];
            for (ElementsArrayType::iterator it=it_begin; it!=it_end; ++it)
            {
                it->GetValue(ACTIVATION_LEVEL)=0;
            }
        }

        vector<unsigned int> node_partition;
        number_of_threads  = 1;
        CreatePartition(number_of_threads, pNodes.size(), node_partition);
        //mfail_node.clear();

        unsigned int New_Id      = 0;
        unsigned int level       = 1;
        unsigned int count       = 0;
        array_1d<double, 3> dir  = ZeroVector(3);
        Node < 3 > ::Pointer     pchild_node;
        std::vector< WeakPointerVector<Element> > LevelElements;

        NodesArrayType::iterator i_begin = pNodes.ptr_begin();
        NodesArrayType::iterator i_end   = pNodes.ptr_end();
        for (ModelPart::NodeIterator i = i_begin; i != i_end; ++i)
        {
            i->GetValue(SPLIT_NODAL)       = false;
            if(i->GetValue(IS_BOUNDARY)==1)
            {
                if(CkeckNumberMasterSurfaces(*(i.base())))
                {
                    split = true;
                    WeakPointerVector<Element>& Neighb_Elem =  i->GetValue(NEIGHBOUR_ELEMENTS);
                    level = 0;
                    for(WeakPointerVector<Element>::iterator  neighb = Neighb_Elem.begin(); neighb!= Neighb_Elem.end(); ++neighb)
                    {
                        if(neighb->GetValue(ACTIVATION_LEVEL)==0)
                        {
                            level++;
                            neighb->GetValue(ACTIVATION_LEVEL)=level;
                            WeakPointerVector<Element>& Neighb_Elem_2 =  neighb->GetValue(NEIGHBOUR_ELEMENTS);
                            for(WeakPointerVector<Element>::iterator  neighb_2 = Neighb_Elem_2.begin(); neighb_2!= Neighb_Elem_2.end(); ++neighb_2)
                            {
                                neighb_2->GetValue(ACTIVATION_LEVEL)=level;
                                WeakPointerVector<Element>& Neighb_Elem_3 =  neighb_2->GetValue(NEIGHBOUR_ELEMENTS);
                                for(WeakPointerVector<Element>::iterator  neighb_3 = Neighb_Elem_3.begin(); neighb_3!= Neighb_Elem_3.end(); ++neighb_3)
                                {
                                    if(neighb_3->GetValue(ACTIVATION_LEVEL)==0)
                                        neighb_3->GetValue(ACTIVATION_LEVEL)=level;
                                }
                            }
                        }
                    }

                    if(level!=0)
                    {
// 			KRATOS_WATCH(i->Id())
// 			KRATOS_WATCH(level)
                        count = level;
                        LevelElements.resize(count);
                        for(WeakPointerVector<Element>::iterator  neighb = Neighb_Elem.begin(); neighb!= Neighb_Elem.end(); neighb++)
                        {
                            const int&  the_level = neighb->GetValue(ACTIVATION_LEVEL);
                            LevelElements[the_level-1].push_back(*neighb.base());
                        }

                        Node<3>::Pointer& pfather_node =  *i.base();
// 			KRATOS_WATCH(LevelElements.size())
                        for(unsigned int ii = 0; ii<LevelElements.size()-1; ii++)
                        {
                            /// creo los nodos para cada grupo elementos
                            New_Id  = this_model_part.Nodes().size() + 1;
                            result++;
                            Create_New_Node(this_model_part, dir,  pchild_node, New_Id, pfather_node);
                            WeakPointerVector<Element>& thiswaek = LevelElements[ii];
                            for(WeakPointerVector<Element>::iterator  it = thiswaek.begin(); it != thiswaek.end();  it++)
                            {
// 			    KRATOS_WATCH(it->Id())
// 			    KRATOS_WATCH(it->GetValue(ACTIVATION_LEVEL))
                                Element::GeometryType& geom = it->GetGeometry();
                                for(unsigned int j = 0; j<geom.size(); j++)
                                    if(geom[j].Id()==i->Id())
                                        Insert_New_Node_In_Elements(geom, pchild_node, pfather_node->Id());


                                /// creo los nuevos nodos para las condiciones
                                /*
                                WeakPointerVector<Condition>& Neighb_Conditions   = it->GetValue(NEIGHBOUR_CONDITIONS);
                                KRATOS_WATCH(Neighb_Conditions.size())
                                for(WeakPointerVector<Condition>::iterator  icond = Neighb_Conditions.begin(); icond != Neighb_Conditions.end(); ++icond)
                                  {

                                 KRATOS_WATCH(icond->Id())
                                 Condition::GeometryType& geom_c = icond->GetGeometry();
                                     for(unsigned int j = 0; j<geom_c.size(); j++)
                                       if(geom_c[j].Id()==i->Id())
                                     Insert_New_Node_In_Conditions(geom_c, pchild_node, pfather_node->Id());

                                  }
                                  */
                                //std::cout<< "*************************************"<< std::endl;

                            }

                            int j = 0;
                            Condition::Pointer rcond;
                            WeakPointerVector< Condition >& neighb_conds = pfather_node->GetValue(NEIGHBOUR_CONDITIONS);
                            for(WeakPointerVector< Condition >::iterator neighb_cond = neighb_conds.begin(); neighb_cond != neighb_conds.end(); neighb_cond++)
                            {
                                if(neighb_cond->GetValue(IS_CONTACT_MASTER)==0)
                                {
                                    rcond = (neighb_conds(j).lock());
                                    Insert_New_Node_In_Conditions(rcond);
                                }
                                j++;
                            }
                            /*
                             bool exist = false;
                             WeakPointerVector<Condition>& neighb_conds =  i->GetValue(NEIGHBOUR_CONDITIONS);
                             int j =  DetectCondition(pfather_node, exist);
                             if(exist)
                             {
                            KRATOS_WATCH("AKIIIIIIIIIIIIIIIIIIII")
                            Condition::Pointer rcond = (neighb_conds(j).lock());
                            Insert_New_Node_In_Conditions(rcond);
                                 }
                                 */
                        }

                        for(WeakPointerVector<Element>::iterator  neighb = Neighb_Elem.begin(); neighb!= Neighb_Elem.end(); neighb++)
                        {
                            neighb->GetValue(ACTIVATION_LEVEL)=0;
                            WeakPointerVector<Element>& Neighb_Elem_2 =  neighb->GetValue(NEIGHBOUR_ELEMENTS);
                            for(WeakPointerVector<Element>::iterator  neighb_2 = Neighb_Elem_2.begin(); neighb_2!= Neighb_Elem_2.end(); neighb_2++)
                            {
                                WeakPointerVector<Element>& Neighb_Elem_3 =  neighb_2->GetValue(NEIGHBOUR_ELEMENTS);
                                neighb_2->GetValue(ACTIVATION_LEVEL)=0;
                                for(WeakPointerVector<Element>::iterator  neighb_3 = Neighb_Elem_3.begin(); neighb_3!= Neighb_Elem_3.end(); neighb_3++)
                                {
                                    neighb_3->GetValue(ACTIVATION_LEVEL)=0;
                                }
                            }
                        }

                    }
                }

                LevelElements.clear();
                LevelElements.reserve(10);
            }
        }


        std::cout<< "FINALIZE END"<< std::endl;
        return result;

        KRATOS_CATCH("")
    }



    inline bool CkeckNumberConditionsSurfaces(const Node<3>::Pointer& rnode)
    {
        unsigned int count = 0;
        //unsigned int id = rnode->Id();
        WeakPointerVector<Condition>& Neighb_Conditions = rnode->GetValue(NEIGHBOUR_CONDITIONS);
        for(WeakPointerVector<Condition>::iterator icond = Neighb_Conditions.begin(); icond!=Neighb_Conditions.end(); icond++)
        {
            if( icond->GetValue( IS_CONTACT_MASTER )== 0)
            {
                count++;
            }
        }
        return count>=2;
    }


    inline bool CkeckNumberMasterSurfaces(const Node<3>::Pointer&  rnode)
    {
        unsigned int count = 0;
        WeakPointerVector<Condition>& Neighb_Conditions = rnode->GetValue(NEIGHBOUR_CONDITIONS);
        for(WeakPointerVector<Condition>::iterator icond = Neighb_Conditions.begin(); icond!=Neighb_Conditions.end(); icond++)
            if( icond->GetValue( IS_CONTACT_MASTER )== 1)
                count++;

        return count>2;
    }


    inline int DetectCondition(const Node<3>::Pointer&  rnode, bool& exist )
    {
        unsigned int count = 0;
        exist = false;
        WeakPointerVector<Condition>& Neighb_Conditions = rnode->GetValue(NEIGHBOUR_CONDITIONS);
        for(WeakPointerVector<Condition>::iterator icond = Neighb_Conditions.begin(); icond!=Neighb_Conditions.end(); icond++)
        {
            if( icond->GetValue( IS_CONTACT_MASTER )== 0)
            {
                exist = true;
                break;
            }
            count++;
        }
        return count;
    }


    inline unsigned int ComputeNumNeighbElem(const ElementsArrayType::iterator& relem)
    {
        unsigned int count = 0;
        WeakPointerVector<Element>& Neighb_Elem_Elem = relem->GetValue(NEIGHBOUR_ELEMENTS);
        for(WeakPointerVector<Element>::iterator it  = Neighb_Elem_Elem.begin(); it != Neighb_Elem_Elem.end(); it++)
            if( it->Id()!=relem->Id())
                count++;

        return count;
    }


    inline void Calculate_Normal_Faliure_Maps(array_1d<double, 3 > & normal, const array_1d<double, 3 > & failure_map)
    {
        normal[0] = -failure_map[1];
        normal[1] = failure_map[0];
        normal[2] = 0.00;
        noalias(normal) = normal / norm_2(normal);

    }


    ModelPart& mr_model_part;
    unsigned int mdomain_size;
    bool mjoint_create;
    WeakPointerVector< Element > mResetingElements;



};
}
#endif

