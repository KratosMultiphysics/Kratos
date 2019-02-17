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
*   Last Modified by:    $Author: Nelson $
*   Date:                $Date: 2011-03-29 11:41:31 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/

#if !defined(KRATOS_INTERSECT_TRIANGLES_CASES_INCLUDED)
#define  KRATOS_INTERSECT_TRIANGLES_CASES_INCLUDED
//System includes
//External includes
#include "boost/smart_ptr.hpp"
#include <cmath>

//Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "custom_utilities/sd_math_utils.h"

#include "includes/model_part.h"
#include "includes/mesh.h"
#include "geometries/geometry.h"
#include "includes/element.h"

#include "custom_utilities/segment_2d.h"


namespace Kratos
{
template<class TConfigure>
class IntersectTriangleCases
{
public:

    typedef typename TConfigure::PointType                      PointType;
    typedef typename PointType::CoordinatesArrayType            CoordinatesArrayType;
    typedef typename TConfigure::ContainerType                  ContainerType;
    typedef typename TConfigure::PointerType                    PointerType;
    typedef typename TConfigure::IteratorType                   IteratorType;
    typedef typename TConfigure::ResultContainerType            ResultContainerType;
    typedef typename TConfigure::ResultPointerType              ResultPointerType;
    typedef typename TConfigure::ResultIteratorType             ResultIteratorType;
    typedef typename TConfigure::ContactPairType                ContactPairType;
    typedef typename TConfigure::ContainerContactType           ContainerContactType;
    typedef typename TConfigure::IteratorContactType            IteratorContactType;
    typedef typename TConfigure::PointerContactType             PointerContactType;
    typedef typename TConfigure::PointerTypeIterator            PointerTypeIterator;


    typedef Node<3>                                              NodeType;
    typedef Node<3>::Pointer                                     NodePointerType;
    typedef IntersectionSegment2DToSegment2D                     IntersectionSegments;


    enum IntersectionTriangleCases {Case_0 = 0, /// Default No hay nada
                                    Case_1 = 1,  /// Un nodo dentro   ( Segmento-Punto )
                                    Case_2 = 2,  /// Dos nodos dentro ( Segmento Punto)
                                    Case_3 = 3,  /// Corner dentro    ( Point-Point)
                                    Case_4 = 4,  /// un corner afuera ( Point-Segment distancia larga)
                                    Case_5 = 5,  /// Uno dentro y otro corner afuera
                                    Case_6 = 6
                                   }; /// dos corner afuera ( Point-Segment distancia larga)

    typedef IntersectionTriangleCases Intersect;

    IntersectTriangleCases() {};

    IntersectTriangleCases(ModelPart& model_part) : mr_model_part(model_part) {}
    virtual ~IntersectTriangleCases() {}

    //localiza Casos 1,2,3,5
    Intersect LocateCaseItersection(NodePointerType& NodeOutside,
                                    bool& Change,
                                    const std::vector<Node<3>::Pointer>& InsideNodes,
                                    const PointerType& MasterObject,
                                    const PointerType& SlaveObject)
    {

        Intersect  rCase                                  = Case_0;
        Change                                            = true;
        WeakPointerVector<Condition>& neighb_cond_master  =  MasterObject->GetValue(NEIGHBOUR_CONDITIONS);
        //WeakPointerVector<Condition>& neighb_cond_slave   =  SlaveObject->GetValue(NEIGHBOUR_CONDITIONS);

        const int size    = InsideNodes.size();
        const int size_1  = neighb_cond_master.size();
        //const int size_2  = neighb_cond_slave.size();

// 	    KRATOS_WATCH(size)
// 	    KRATOS_WATCH(size_1)
// 	    KRATOS_WATCH(size_2)
//
//  	    for(unsigned int i = 0; i<size; i++ )
//  	      std::cout <<  "Node_" << i+1 << " = "  << (InsideNodes[i])->Id() << std::endl;

        if(size!=0 && size_1!=0)
        {
            if(size_1==1)
            {
                if(size==1)
                {
                    rCase = Case_1;
                }
                if(size==2)
                {
                    rCase = Case_2;
                }
            }
            else if(size_1 !=1 && (size==2))
            {
                /*if(size_2!=1){
                Intersect  Aux_Case_1 = Case_0;
                Intersect  Aux_Case_2 = Case_0;
                Is_Corner(NodeOutside, Aux_Case_1, InsideNodes[0], SlaveObject, MasterObject);
                //std::cout <<  "Case_a = " << Aux_Case_1 << std::endl;
                //std::cout <<  "Case_b = " << Aux_Case_2 << std::endl;
                Is_Corner(NodeOutside, Aux_Case_2, InsideNodes[1], SlaveObject, MasterObject);
                bool one = (Aux_Case_1==Case_3 && Aux_Case_2==Case_3);
                bool two = (Aux_Case_1==Case_0 && Aux_Case_2==Case_0);
                   if(one==true || two==true)
                       rCase = Case_3;
                }
                else*/
                rCase = Case_2;
            }

            else if(size_1 !=1 && (size==1))  //is corner???
            {
                rCase = Case_1;
                //KRATOS_WATCH("CORNERRRRRRR")
                //Is_Corner(NodeOutside, rCase, InsideNodes[0], SlaveObject, MasterObject);
                //if(rCase == Case_0){rCase = Case_1;}
                //KRATOS_WATCH(rCase)
            }
            else
                rCase = Case_1;
        }
        else if(size!=0 && size_1==0)
        {
            rCase = Case_1;
        }
        else
        {
            rCase = Case_0;
        }

        if(rCase==Case_5 || rCase==Case_3) //uno dentro y otro fuera
            Change = false;


        //std::cout <<  "Vecinos  = " << neighb_cond_master.size() << std::endl;
        //std::cout <<  "Size     = " << size << std::endl;
        //std::cout <<  "Case     = " << rCase << std::endl;
        return rCase;

    }


    // Si el nodo esta dentro de elemento
    void Is_Corner(NodePointerType& NodeOutside,
                   Intersect&  rCase,
                   const NodePointerType& SlaveNode,
                   const PointerType& SlaveObject,
                   const PointerType& MasterObject)
    {

        KRATOS_TRY

        WeakPointerVector<Condition>& neighb_cond_master  = MasterObject->GetValue(NEIGHBOUR_CONDITIONS);
        WeakPointerVector<Condition>& neighb_cond_slave   = SlaveNode->GetValue(NEIGHBOUR_CONDITIONS);
        if(neighb_cond_master.size()!=0 || neighb_cond_slave.size()!=0)
        {
            double distance  = 0.00;
            std::vector<unsigned int>           segment;
            std::vector<double>                 Distances;
            std::vector<unsigned int>::iterator it;
            vector<array_1d<double, 2> >        Points0;
            vector<array_1d<double, 2> >        Points1;
            array_1d<double, 2>                 Point;
            array_1d<double, 2>                 Vect;
            array_1d<double, 2>                 Node_Point;

            Node_Point[0]    = SlaveNode->X();
            Node_Point[1]    = SlaveNode->Y();

            unsigned int I   = 0;
            unsigned int II  = 1;
            unsigned int III = 1;


            Points0.resize(2, false);
            Points1.resize(2, false);


            for(WeakPointerVector<Condition>::iterator cond_slave  = neighb_cond_slave.begin(); cond_slave!= neighb_cond_slave.end(); ++cond_slave)
            {
                Condition::GeometryType& geom = cond_slave->GetGeometry();
                Point[0] = 0.00;
                Point[1] = 0.00;

                Points0(0)[0] = geom[0].X();
                Points0(0)[1] = geom[0].Y();
                Points0(1)[0] = geom[1].X();
                Points0(1)[1] = geom[1].Y();

                I   = 0;
                III = 1;
                for(WeakPointerVector< Condition >::iterator cond  = neighb_cond_master.begin(); cond!= neighb_cond_master.end(); ++cond)
                {
                    Condition::GeometryType& geom_3 = cond->GetGeometry();
                    Points1(0)[0] = geom_3[0].X();
                    Points1(0)[1] = geom_3[0].Y();
                    Points1(1)[0] = geom_3[1].X();
                    Points1(1)[1] = geom_3[1].Y();

                    if(IntersectionSegments::IntersectSegment(Point, Points0, Points1)!=IT_EMPTY)
                    {
                        noalias(Vect) = Point - Node_Point;
                        distance      = std::sqrt(inner_prod(Vect, Vect) );
                        if(segment.size()==0)
                        {
                            Distances.push_back(distance);
                            segment.push_back(I);
                        }
                        else
                        {
                            it = std::find(segment.begin(), segment.end(), I);
                            if(it==segment.end())
                            {
                                Distances.push_back(distance);
                                segment.push_back(I);
                            }
                        }
                    }

                    I++;
                    III++;
                    if(III>neighb_cond_master.size())
                        break;
                }

                II++;
                if(II>neighb_cond_slave.size())
                    break;
            }

            if(segment.size()==0)
            {
                rCase = Case_0;  //Case_1;
                return;
            }
            //is case_1;
            else if(segment.size()==1)
            {
                rCase = Case_0; //Case_1;
                return;
            }
            else
            {
                //busca si el tercer segmento intersecta con el master
                // is case_3 o 5;
                rCase = Is_Case3_Or_Case_5(NodeOutside, SlaveNode, SlaveObject, MasterObject);

                // Verifica si de verdad es corner o si es un nodo dentro. (case_3 or case_1)
                if(rCase==Case_3 && Distances.size()==2)
                {

                    const double toler  = 1E-6;
                    double min   = MinDistancesEdges(SlaveObject, MasterObject);
                    double diff  = std::fabs(Distances[0] - Distances[1]);
                    double ratio = std::fabs(diff/min);
                    //std::cout<< "ratio = " << ratio << std::endl;
                    if(ratio>toler)
                        rCase =Case_1;
                }

// 	       if(rCase==Case_3)
// 		  rCase =Case_1;
            }
        }
        return;
        KRATOS_CATCH("")
    }


    Intersect Is_Case3_Or_Case_5(NodePointerType& NodeOutside,
                                 const NodePointerType& SlaveNode,
                                 const PointerType& SlaveObject,
                                 const PointerType& MasterObject)
    {

        Intersect rCase = Case_3;
        bool is_case5   = false;
        Element::GeometryType& geom_slave  = SlaveObject->GetGeometry();
        //Element::GeometryType& geom_master = MasterObject->GetGeometry();
        vector<array_1d<double, 2> >        Points0;
        vector<array_1d<double, 2> >        Points1;
        array_1d<double, 2>                 Point;
        std::vector<unsigned int>           segment;

        WeakPointerVector<Condition>        rConditions;


        Points0.resize(2, false);
        Points1.resize(2, false);

        unsigned int I   = 0;
        unsigned int III = 1;
        std::vector<NodePointerType> Oposite_Ids(2);
        for(unsigned int i = 0; i<geom_slave.size(); i++)
        {
            if(SlaveNode->Id()!=geom_slave[i].Id())
            {
                Oposite_Ids[I] =  geom_slave(i);
                Points0(I)[0]  =  geom_slave[i].X();
                Points0(I)[1]  =  geom_slave[i].Y();
                I++;
            }
        }

        WeakPointerVector<Condition>& neighb_cond_master  = MasterObject->GetValue(NEIGHBOUR_CONDITIONS);
        for(WeakPointerVector< Condition >::iterator cond  = neighb_cond_master.begin(); cond!= neighb_cond_master.end(); ++cond)
        {
            Condition::GeometryType& geom_3 = cond->GetGeometry();
            Points1(0)[0] = geom_3[0].X();
            Points1(0)[1] = geom_3[0].Y();
            Points1(1)[0] = geom_3[1].X();
            Points1(1)[1] = geom_3[1].Y();

            if(IntersectionSegments::IntersectSegment(Point, Points0, Points1)!=IT_EMPTY)
            {
                rConditions.push_back(*cond.base()); // *neighb_cond_master(I).lock());
            }

            III++;
            if(III>neighb_cond_master.size())
                break;
        }

        if(rConditions.size()==2)
        {
            is_case5 = true;
        }
        if(is_case5 == true)
        {
            III   = 1;
            rCase = Case_5;
            unsigned int j   = 0;
            std::vector<NodePointerType> Ids(4);
            std::vector<unsigned int> segment;

            for( WeakPointerVector< Condition >::iterator rcond = rConditions.begin(); rcond != rConditions.end(); rcond++)
            {
                Condition::GeometryType& geom_3 = (rcond)->GetGeometry();
                Ids[0+j] = geom_3(0);
                Ids[1+j] = geom_3(1);
                j += 2;
                III++;
                if(III>rConditions.size())
                    break;
            }

            for(unsigned int i = 0; i<4; i++)
            {
                NodeOutside=Ids[i];
                for(unsigned int j = 0; j<4 && j!=i; j++)
                {
                    if( NodeOutside->Id() ==Ids[j]->Id() )
                    {
                        break;
                    }
                }
            }

        }
        return rCase;

    }

    // compute the distance triangles edges
    double MinDistancesEdges( const PointerType& SlaveObject,
                              const PointerType& MasterObject)
    {

        Element::GeometryType& geom_slave  = SlaveObject->GetGeometry();
        Element::GeometryType& geom_master = MasterObject->GetGeometry();

        Vector Distances;
        Distances.resize(6, false);
        Distances = ZeroVector(6);

        array_1d<double, 3> Vect;

        noalias(Vect)  = geom_slave[0].Coordinates() - geom_slave[1].Coordinates();
        Distances[0]   = std::sqrt(inner_prod(Vect, Vect ) ) ;

        noalias(Vect)  = geom_slave[1].Coordinates() - geom_slave[2].Coordinates();
        Distances[1]   = std::sqrt(inner_prod(Vect, Vect ) ) ;

        noalias(Vect)  = geom_slave[2].Coordinates() - geom_slave[0].Coordinates();
        Distances[2]   = std::sqrt(inner_prod(Vect, Vect ) ) ;

        noalias(Vect)  = geom_master[0].Coordinates() - geom_master[1].Coordinates();
        Distances[3]   = std::sqrt(inner_prod(Vect, Vect ) ) ;

        noalias(Vect)  = geom_master[1].Coordinates() - geom_master[2].Coordinates();
        Distances[4]   = std::sqrt(inner_prod(Vect, Vect ) ) ;

        noalias(Vect)  = geom_master[2].Coordinates() - geom_master[0].Coordinates();
        Distances[5]   = std::sqrt(inner_prod(Vect, Vect ) ) ;

        double min = (*std::min_element(Distances.begin(), Distances.end() ) );
        return min;

    }


// comparacion con desplazamientos
    bool Test_One(
        const NodePointerType& SlaveNode,
        const PointerType& MasterObject)
    {

        KRATOS_TRY


        vector<array_1d<double, 2> >      Points0;
        vector<array_1d<double, 2> >      Points1;
        array_1d<double, 2>               Point;

        WeakPointerVector<Condition>& neighb_cond_master  = MasterObject->GetValue(NEIGHBOUR_CONDITIONS);
        if(neighb_cond_master.size()!=0)
        {
            array_1d<double,3>& old_pos                       = SlaveNode->FastGetSolutionStepValue(DISPLACEMENT,3);

            Points0.resize(2, false);
            Points1.resize(2, false);

            Points0(0)[0] = SlaveNode->X0() + old_pos[0];
            Points0(0)[1] = SlaveNode->Y0() + old_pos[1];
            Points0(1)[0] = SlaveNode->X();
            Points0(1)[1] = SlaveNode->Y();

            unsigned int JJ = 1;
            for(WeakPointerVector< Condition >::iterator cond  = neighb_cond_master.begin(); cond!= neighb_cond_master.end(); cond++)
            {
                Condition::GeometryType& geom_2 = cond->GetGeometry();

                Points1(0)[0] = geom_2[0].X();
                Points1(0)[1] = geom_2[0].Y();
                Points1(1)[0] = geom_2[1].X();
                Points1(1)[1] = geom_2[1].Y();

                if(IntersectionSegments::IntersectSegment(Point, Points0, Points1)!=IT_EMPTY)
                    return true;

                JJ++;
                if(JJ>neighb_cond_master.size())
                    break;
            }
        }
        return false;
        KRATOS_CATCH("")
    }



private:
    ModelPart mr_model_part;
    static const int IT_POINT   = 0;
    static const int IT_SEGMENT = 1;
    static const int IT_EMPTY   = 2;

};



}//namespace Kratos.

#endif /* KRATOS_INTERSECT_TRIANGLES_CASES_INCLUDED defined */

