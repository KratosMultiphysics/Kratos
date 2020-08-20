//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti.
//

#if !defined(KRATOS_ULF_UTILITIES_INCLUDED )
#define  KRATOS_ULF_UTILITIES_INCLUDED






// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "pfem_2_application_variables.h"
#include "boost/smart_ptr.hpp"
#include "includes/cfd_variables.h"
#include "includes/c2c_variables.h"
#include "includes/deprecated_variables.h"

#include "processes/process.h"
#include "includes/condition.h"
#include "processes/find_nodal_neighbours_process.h"
//#include "incompressible_fluid_application.h"



namespace Kratos
{
class Pfem2Utils
{
public:
    typedef Node<3> NodeType;
    //**********************************************************************************************
    //**********************************************************************************************

    //functions to apply the boundary conditions

    void ApplyBoundaryConditions(ModelPart& ThisModelPart, int laplacian_type)
    {
        KRATOS_TRY;

        if(laplacian_type == 1)
        {
            //identify fluid nodes and mark as boundary the free ones
            for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin(); in!=ThisModelPart.NodesEnd(); in++)
            {
                //marking wet nodes
                if(in->FastGetSolutionStepValue(IS_STRUCTURE) )
                    if( (in->GetValue(NEIGHBOUR_ELEMENTS)).size() != 0)
                        in->FastGetSolutionStepValue(IS_FLUID) = 1.0;
                    else //it is not anymore of fluid
                        in->FastGetSolutionStepValue(IS_FLUID) = 0.0;
                //marking as free surface the lonely nodes
                else if( (in->GetValue(NEIGHBOUR_ELEMENTS)).size() == 0)
                    in->FastGetSolutionStepValue(IS_BOUNDARY) = 1.0;
            }

            //identify the free surface
            for(ModelPart::NodeIterator i = ThisModelPart.NodesBegin() ;
                    i != ThisModelPart.NodesEnd() ; ++i)
            {
                //reset the free surface
                i->FastGetSolutionStepValue(IS_FREE_SURFACE) = 0;
                i->Free(PRESSURE);

                //identify the free surface and fix the pressure accordingly
                if( i->FastGetSolutionStepValue(IS_BOUNDARY) != 0
                        &&
                        i->FastGetSolutionStepValue(IS_STRUCTURE) == 0)
                {
                    i->FastGetSolutionStepValue(IS_FREE_SURFACE) = 1;
                    i->FastGetSolutionStepValue(PRESSURE) = 0.00;
                    i->Fix(PRESSURE);
                }
            }


        }
        else  //case of discrete laplacian
        {
            //identify fluid nodes and mark as boundary the free ones
            for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin(); in!=ThisModelPart.NodesEnd(); in++)
            {
                //marking wet nodes
                if(in->FastGetSolutionStepValue(IS_STRUCTURE) )
                    if( (in->GetValue(NEIGHBOUR_ELEMENTS)).size() != 0)
                        in->FastGetSolutionStepValue(IS_FLUID) = 1.0;
                    else //it is not anymore of fluid
                        in->FastGetSolutionStepValue(IS_FLUID) = 0.0;
                //marking as free surface the lonely nodes
                else if( (in->GetValue(NEIGHBOUR_ELEMENTS)).size() == 0)
                    in->FastGetSolutionStepValue(IS_BOUNDARY) = 1.0;

                if( (in->GetValue(NEIGHBOUR_ELEMENTS)).size() == 0)
                    in->FastGetSolutionStepValue(PRESSURE) = 0.0;
            }

            //identify the free surface
            for(ModelPart::NodeIterator i = ThisModelPart.NodesBegin() ;
                    i != ThisModelPart.NodesEnd() ; ++i)
            {
                //reset the free surface
                i->FastGetSolutionStepValue(IS_FREE_SURFACE) = 0;

                //identify the free surface and fix the pressure accordingly
                if( i->FastGetSolutionStepValue(IS_BOUNDARY) != 0
                        &&
                        i->FastGetSolutionStepValue(IS_STRUCTURE) == 0)
                {
                    i->FastGetSolutionStepValue(IS_FREE_SURFACE) = 1;
                }
            }



        }

        KRATOS_CATCH("");
    }
    //**********************************************************************************************
    //**********************************************************************************************
    void IdentifyFluidNodes(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;

        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                i!=ThisModelPart.ElementsEnd(); i++)
        {
            //calculating shape functions values
            Geometry< Node<3> >& geom = i->GetGeometry();

            for(unsigned int i = 0; i<geom.size(); i++)
                geom[i].FastGetSolutionStepValue(IS_FLUID) = 1;

        }
        KRATOS_CATCH("");
    }

    //**********************************************************************************************
    //**********************************************************************************************
    void ApplyMinimalPressureConditions(std::vector< GlobalPointersVector< Node<3> > >& connected_components)
    {
        KRATOS_TRY;

        //verify that the minimal conditions for pressure are applied
        std::vector<bool> to_be_prescribed(connected_components.size());
        array_1d<double,3> fixed_vel;
        for(unsigned int i = 0; i<connected_components.size(); i++)
        {
            int boundary_nodes = 0;
            int prescribed_vel_nodes = 0;
            GlobalPointersVector< Node<3> >& node_list = connected_components[i];
            for( GlobalPointersVector< Node<3> >::iterator in = node_list.begin();
                    in != node_list.end(); in++)
            {
                //free the pressure
                in->Free(PRESSURE);

                //cound boundary and fixed velocity nodes
                if(in->FastGetSolutionStepValue(IS_BOUNDARY) == 1)
                {
                    //count the nodes added
                    boundary_nodes += 1;

                    //if it is structure add it to the fixed nodes
                    if(in->FastGetSolutionStepValue(IS_STRUCTURE) == 1)
                        prescribed_vel_nodes += 1;

                }
            }

            //if all the boundary nodes have the velocity prescribed => it is needed to
            //prescribe the pressure on 1 node
            if(boundary_nodes == prescribed_vel_nodes)
            {
                bool one_is_prescribed = false;
                for( GlobalPointersVector< Node<3> >::iterator in = node_list.begin();
                        in != node_list.end(); in++)
                {
                    if( one_is_prescribed == false &&
                            in->FastGetSolutionStepValue(IS_BOUNDARY) == 1 )
                    {
                        std::cout << "fixed pressure on node " << in->Id() << std::endl;
                        one_is_prescribed = true;
                        in->Fix(PRESSURE);
                    }
                }
            }
        }
        KRATOS_CATCH("");
    }

    //**********************************************************************************************
    //**********************************************************************************************
    double EstimateDeltaTime(double dt_min, double dt_max, ModelPart& ThisModelPart)
    {
        KRATOS_TRY

        array_1d<double,3> dx, dv;
        double deltatime = dt_max;
        double dvside, lside;

        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                i!=ThisModelPart.ElementsEnd(); i++)
        {
            //calculating shape functions values
            Geometry< Node<3> >& geom = i->GetGeometry();

            for(unsigned int i1 = 0; i1 < geom.size()-1; i1++)
            {
                for(unsigned int i2 = i1 + 1; i2 < geom.size(); i2++)
                {
                    dx[0] = geom[i2].X() - geom[i1].X();
                    dx[1] = geom[i2].Y() - geom[i1].Y();
                    dx[2] = geom[i2].Z() - geom[i1].Z();

                    lside = inner_prod(dx,dx);

                    noalias(dv) = geom[i2].FastGetSolutionStepValue(VELOCITY);
                    noalias(dv) -= geom[i1].FastGetSolutionStepValue(VELOCITY);

                    dvside = inner_prod(dx,dv);

                    double dt;
                    if(dvside < 0.0) //otherwise the side is "getting bigger" so the time step has not to be diminished
                    {
                        dt = fabs( lside/dvside );
                        if(dt < deltatime) deltatime = dt;
                    }

                }
            }
        }

        if(deltatime < dt_min)
        {
            std::cout << "ATTENTION dt_min is being used" << std::endl;
            deltatime = dt_min;
        }


        return deltatime;

        KRATOS_CATCH("")
    }
    //**********************************************************************************************
    //**********************************************************************************************
    void MarkOuterNodes(const array_1d<double,3>& corner1, const array_1d<double,3>& corner2,
                        ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY;

        //add a big number to the id of the nodes to be erased
        int n_erased = 0;
        double xmax, xmin, ymax,ymin, zmax, zmin;



        if(corner1[0] > corner2[0])
        {
            xmax = corner1[0];
            xmin = corner2[0];
        }
        else
        {
            xmax = corner2[0];
            xmin = corner1[0];
        }

        if(corner1[1] > corner2[1])
        {
            ymax = corner1[1];
            ymin = corner2[1];
        }
        else
        {
            ymax = corner2[1];
            ymin = corner1[1];
        }

        if(corner1[2] > corner2[2])
        {
            zmax = corner1[2];
            zmin = corner2[2];
        }
        else
        {
            zmax = corner2[2];
            zmin = corner1[2];
        }



        for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            bool erase = false;
            double& x = in->X();
            double& y = in->Y();
            double& z = in->Z();

            if(x<xmin || x>xmax)		erase = true;
            else if(y<ymin || y>ymax)	erase = true;
            else if(z<zmin || z>zmax)	erase = true;

            if(erase == true)
            {
                n_erased += 1;
                in->Set(TO_ERASE, true);
            }
        }

        KRATOS_CATCH("")
    }

    //**********************************************************************************************
    //**********************************************************************************************
    void MarkExcessivelyCloseNodes(ModelPart::NodesContainerType& rNodes, const double admissible_distance_factor)
    {
        KRATOS_TRY;
        KRATOS_WATCH("ENTERD Mark close nodes")
        double fact2 = admissible_distance_factor*admissible_distance_factor;


        for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            if(in->FastGetSolutionStepValue(IS_STRUCTURE) == 0) //if it is not a wall node i can erase
            {
                double hnode2 = in->FastGetSolutionStepValue(NODAL_H);
                hnode2 *= hnode2; //take the square

                //loop on neighbours and erase if they are too close
                for( GlobalPointersVector< Node<3> >::iterator i = in->GetValue(NEIGHBOUR_NODES).begin();
                        i != in->GetValue(NEIGHBOUR_NODES).end(); i++)
                {
                    if( bool(i->Is(TO_ERASE)) == false) //we can erase the current node only if the neighb is not to be erased
                    {
                        double dx = i->X() - in->X();
                        double dy = i->Y() - in->Y();
                        double dz = i->Z() - in->Z();

                        double dist2 = dx*dx + dy*dy + dz*dz;

                        if(dist2 < fact2 *  hnode2)
                            in->Set(TO_ERASE, true);
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }


    //**********************************************************************************************
    //*		//**********************************************************************************************
    ////////////////////////////////////////////////////////////
    double Length(array_1d<double,3>& Point1, array_1d<double,3>& Point2)
    {
        //KRATOS_WATCH("length calculation")
        return sqrt((Point1[0]-Point2[0])*(Point1[0]-Point2[0]) + (Point1[1]-Point2[1])*(Point1[1]-Point2[1]) +(Point1[2]-Point2[2])*(Point1[2]-Point2[2]));
    }
    double CalculateTriangleArea3D(	array_1d<double,3>& Point1, array_1d<double,3>& Point2, array_1d<double,3>& Point3	)
    {
        //Heron's formula
        double a=Length(Point1, Point2);//sqrt((Point1[0]-Point2[0])*(Point1[0]-Point2[0]) + (Point1[1]-Point2[1])*(Point1[1]-Point2[1]) +(Point1[2]-Point2[2])*(Point1[2]-Point2[2]));
        double b=Length(Point1, Point3);//sqrt((Point3[0]-Point2[0])*(Point3[0]-Point2[0]) + (Point3[1]-Point2[1])*(Point3[1]-Point2[1]) +(Point3[2]-Point2[2])*(Point3[2]-Point2[2]));
        double c=Length(Point2, Point3);//sqrt((Point1[0]-Point3[0])*(Point1[0]-Point3[0]) + (Point1[1]-Point3[1])*(Point1[1]-Point3[1]) +(Point1[2]-Point3[2])*(Point1[2]-Point3[2]));
        double p=0.5*(a+b+c);
        return sqrt(p*(p-a)*(p-b)*(p-c));
    }

    ////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    void CalculateNodalArea(ModelPart& ThisModelPart, int domain_size)
    {

        KRATOS_TRY

        //set to zero the nodal area
        for(ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin();
                in!=ThisModelPart.NodesEnd(); in++)
        {
            in->FastGetSolutionStepValue(NODAL_AREA) = 0.00;
        }

        if(domain_size == 2)
        {
            double area = 0.0;
            for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i!=ThisModelPart.ElementsEnd(); i++)
            {
                //calculating shape functions values
                Geometry< Node<3> >& geom = i->GetGeometry();

                area = GeometryUtils::CalculateVolume2D(geom);
                area *= 0.333333333333333333333333333;


                geom[0].FastGetSolutionStepValue(NODAL_AREA) += area;
                geom[1].FastGetSolutionStepValue(NODAL_AREA) += area;
                geom[2].FastGetSolutionStepValue(NODAL_AREA) += area;

            }
        }
        else if(domain_size == 3)
        {
            for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i!=ThisModelPart.ElementsEnd(); i++)
            {
                double vol;
                //calculating shape functions values
                Geometry< Node<3> >& geom = i->GetGeometry();
                //counting number of structural nodes

                if (geom.size()==4) //not to calculate the nodal area of the membrane which is a 2d element in 3d
                {

                    vol = GeometryUtils::CalculateVolume3D(geom);
                    vol *= 0.25;

                    geom[0].FastGetSolutionStepValue(NODAL_AREA) += vol;
                    geom[1].FastGetSolutionStepValue(NODAL_AREA) += vol;
                    geom[2].FastGetSolutionStepValue(NODAL_AREA) += vol;
                    geom[3].FastGetSolutionStepValue(NODAL_AREA) += vol;




                }
            }
        }
        //r


        KRATOS_CATCH("")
    }
///////////this is not to have nodes close to the free surface
    void MarkNodesCloseToFS(ModelPart& ThisModelPart, int domain_size)
    {
        if (domain_size==2)
        {

            for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i!=ThisModelPart.ElementsEnd(); i++)
            {
                int n_fs=0;

                //counting number on nodes at the wall
                Geometry< Node<3> >& geom = i->GetGeometry();
                n_fs = int(geom[0].FastGetSolutionStepValue(IS_FREE_SURFACE));
                n_fs+= int(geom[1].FastGetSolutionStepValue(IS_FREE_SURFACE));
                n_fs+= int(geom[2].FastGetSolutionStepValue(IS_FREE_SURFACE));
                //if two walls are of freee surface, we check if the third node is close to it or not by passing the alpha-shape
                if (n_fs==2)
                {
                    //if alpha shape tells to preserve
                    if (AlphaShape(1.4, geom)==false)
                    {
                        for (int i=0; i<3; i++)
                        {
                            //if thats not the wall node, remove it
                            if (geom[i].FastGetSolutionStepValue(IS_FREE_SURFACE)==0.0 && geom[i].FastGetSolutionStepValue(IS_FLUID)==1.0  && geom[i].FastGetSolutionStepValue(IS_STRUCTURE)==0.0)
                            {
                                geom[i].Set(TO_ERASE,true);
                                KRATOS_WATCH("NODE CLOSE TO THE FS - WILL BE ERASED!!!!")
                            }
                        }
                    }
                }

            }
        }
        if (domain_size==3)
        {
            KRATOS_WATCH("NOT IMPLEMENTED IN 3D")

        }
    }

    void MarkNodesCloseToWallForBladder(ModelPart& ThisModelPart, const double& crit_distance)

    {
        for(ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin();
                in!=ThisModelPart.NodesEnd(); in++)
        {
            if (in->FastGetSolutionStepValue(DISTANCE)>0.00 && in->FastGetSolutionStepValue(DISTANCE)<crit_distance && in->FastGetSolutionStepValue(IS_STRUCTURE)==0)
            {
                in->Set(TO_ERASE,true);
                KRATOS_WATCH("NODE EXCESSIVELY CLOSE TO WALL!!!!! BLADDER FUNCTION!!!!!!!!!!!!!!!!!!!!!!!!!!111")
            }
        }
    }
    //this is a function originally written by Antonia. It seems to work better then the MarkNodesCloseToWall (its given below)
    void MarkNodesTouchingWall(ModelPart& ThisModelPart, int domain_size, double factor)
    {
        KRATOS_TRY;

        if (domain_size==2)
        {

            for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i!=ThisModelPart.ElementsEnd(); i++)
            {
                double n_str=0;

                //counting number on nodes at the wall
                Geometry< Node<3> >& geom = i->GetGeometry();
                n_str = geom[0].FastGetSolutionStepValue(IS_BOUNDARY);
                n_str+= geom[1].FastGetSolutionStepValue(IS_BOUNDARY);
                n_str+= geom[2].FastGetSolutionStepValue(IS_BOUNDARY);
                //if two walls are at the wall, we check if the third node is close to it or not by passing the alpha-shape
                if (n_str==2.0)
                {
                    BoundedMatrix<double,3,2> sort_coord = ZeroMatrix(3,2);
                    int cnt=1;
                    for (int i=0; i<3; ++i)
                        if(geom[i].FastGetSolutionStepValue(IS_BOUNDARY)==0.0)
                        {
                            sort_coord(0,0) = geom[i].X();
                            sort_coord(0,1) = geom[i].Y();
                        }
                        else
                        {
                            sort_coord(cnt,0) = geom[i].X();
                            sort_coord(cnt,1) = geom[i].Y();
                            cnt++;
                        }

                    array_1d<double,2> vec1 = ZeroVector(2);
                    array_1d<double,2> vec2 = ZeroVector(2);

                    vec1[0] = sort_coord(0,0) - sort_coord(1,0);
                    vec1[1] = sort_coord(0,1) - sort_coord(1,1);

                    vec2[0] = sort_coord(2,0) - sort_coord(1,0);
                    vec2[1] = sort_coord(2,1) - sort_coord(1,1);

                    double outer_prod = 0.0;
                    outer_prod = vec2[1]*vec1[0]-vec1[1]*vec2[0];

                    double length_measure =0.0;
                    length_measure = vec2[0]*vec2[0] + vec2[1]*vec2[1];
                    length_measure = sqrt(length_measure);
                    outer_prod/=length_measure;

                    //KRATOS_WATCH(fabs(outer_prod));
                    //	RATOS_WATCH(factor*length_measure);

                    if (fabs(outer_prod)<factor*length_measure)
                    {
                        for (int i=0; i<3; i++)
                        {
                            //if thats not the wall node, remove it
                            if (geom[i].FastGetSolutionStepValue(IS_BOUNDARY)==0.0)
                            {
                                geom[i].Set(TO_ERASE,true);
                                //KRATOS_WATCH("NODE TOUCHING THE WALL - WILL BE ERASED!!!!")
                            }
                        }
                    }
                }

            }
        }
        else
        {
            for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i!=ThisModelPart.ElementsEnd(); i++)
            {
                double n_str=0;
                //the n_int is just for the bladder example.. otherwise the function works just with is_STR flag
                double n_int=0;

                //counting number on nodes at the wall
                Geometry< Node<3> >& geom = i->GetGeometry();
                if(geom.size() == 4)
                {
                    for(int ii = 0; ii <= domain_size ; ++ii)
                    {
                        n_str += geom[ii].FastGetSolutionStepValue(IS_STRUCTURE);
                        n_int += geom[ii].FastGetSolutionStepValue(IS_INTERFACE);
                    }
                    //if two walls are at the wall, we check if the third node is close to it or not by passing the alpha-shape
                    if (n_str==3.0 && n_int==3.0)
                    {
                        BoundedMatrix<double,4,3> sort_coord = ZeroMatrix(4,3);
                        int cnt=1;
                        for (int i=0; i<4; ++i)
                        {
                            if(geom[i].FastGetSolutionStepValue(IS_STRUCTURE)==0.0)
                            {
                                sort_coord(0,0) = geom[i].X();
                                sort_coord(0,1) = geom[i].Y();
                                sort_coord(0,2) = geom[i].Z();
                            }
                            else
                            {
                                sort_coord(cnt,0) = geom[i].X();
                                sort_coord(cnt,1) = geom[i].Y();
                                sort_coord(cnt,2) = geom[i].Z();
                                cnt++;
                            }
                        }
                        array_1d<double,3> vec1 = ZeroVector(3);
                        array_1d<double,3> vec2 = ZeroVector(3);
                        array_1d<double,3> vec3 = ZeroVector(3);


                        vec1[0] = sort_coord(0,0) - sort_coord(1,0);
                        vec1[1] = sort_coord(0,1) - sort_coord(1,1);
                        vec1[2] = sort_coord(0,2) - sort_coord(1,2);

                        vec2[0] = sort_coord(2,0) - sort_coord(1,0);
                        vec2[1] = sort_coord(2,1) - sort_coord(1,1);
                        vec2[2] = sort_coord(2,2) - sort_coord(1,2);

                        vec3[0] = sort_coord(3,0) - sort_coord(1,0);
                        vec3[1] = sort_coord(3,1) - sort_coord(1,1);
                        vec3[2] = sort_coord(3,2) - sort_coord(1,2);

                        //Control the hight of the thetraedral element
                        //Volume of the tethraedra
                        double vol = (vec2[0]*vec3[1]*vec1[2]-vec2[0]*vec3[2]*vec1[1]+
                                      vec2[1]*vec3[2]*vec1[0]-vec2[1]*vec3[0]*vec1[2]+
                                      vec2[2]*vec3[0]*vec1[1]-vec2[2]*vec3[1]*vec1[0])*0.1666666666667;
                        //Area of the basis
                        array_1d<double,3> outer_prod = ZeroVector(3);
                        outer_prod[0] = vec2[1]*vec3[2]-vec2[2]*vec3[1];
                        outer_prod[1] = vec2[2]*vec3[0]-vec2[0]*vec3[2];
                        outer_prod[2] = vec2[0]*vec3[1]-vec2[1]*vec3[0];
                        double area_base = norm_2(outer_prod);
                        area_base *= 0.5;
                        //height

                        if(area_base >0.0000000001)
                            vol/= area_base;
                        else
                            KRATOS_ERROR<<"error: BAse element has zero area";

                        //vol/=area_base;
                        double length_measure1 = norm_2(vec2);
                        double length_measure = norm_2(vec3);
                        if(length_measure1 < length_measure)
                        {
                            length_measure = length_measure1;
                        }

                        if (fabs(vol)<factor*length_measure)
                        {
                            for (int i=0; i<4; i++)
                            {
                                //if thats not the wall node, remove it
                                //never remove a Lagrangian inlet node
                                if (geom[i].FastGetSolutionStepValue(IS_STRUCTURE)==0.0 && geom[i].FastGetSolutionStepValue(IS_LAGRANGIAN_INLET)==0.0)
                                {
                                    geom[i].Set(TO_ERASE,true);
                                    KRATOS_WATCH("NODE TOUCHING THE WALL - WILL BE ERASED!!!!")
                                }
                            }
                        }
                    }//interface elements
                }//non_shell
            }//all elements loop
        }//domain_size==3
        KRATOS_CATCH("")
    }

    void MarkNodesCloseToWall(ModelPart& ThisModelPart, int domain_size, double alpha_shape)
    {
        if (domain_size==2)
        {

            for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i!=ThisModelPart.ElementsEnd(); i++)
            {
                int n_str=0;

                //counting number on nodes at the wall
                Geometry< Node<3> >& geom = i->GetGeometry();
                n_str = int(geom[0].FastGetSolutionStepValue(IS_STRUCTURE));
                n_str+= int(geom[1].FastGetSolutionStepValue(IS_STRUCTURE));
                n_str+= int(geom[2].FastGetSolutionStepValue(IS_STRUCTURE));
                //if two walls are at the wall, we check if the third node is close to it or not by passing the alpha-shape
                if (n_str==2)
                {
                    //if alpha shape tells to preserve
                    if (AlphaShape(alpha_shape, geom)==false)
                    {
                        for (int i=0; i<3; i++)
                        {
                            //if thats not the wall node, remove it
                            if (geom[i].FastGetSolutionStepValue(IS_STRUCTURE)==0.0)
                            {
                                geom[i].Set(TO_ERASE,true);
                                //KRATOS_WATCH("NODE CLOSE TO THE WALL - WILL BE ERASED!!!!")
                            }
                        }
                    }
                }

            }
        }
        if (domain_size==3)
        {
            KRATOS_WATCH("Inside mark nodes close to wall process")
            for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i!=ThisModelPart.ElementsEnd(); i++)
            {
                int n_str=0;
                int n_fl=0;
                int n_lag=0;
                int n_interf=0;
                //counting number on nodes at the wall
                Geometry< Node<3> >& geom = i->GetGeometry();
                for (unsigned int iii=0; iii<geom.size(); iii++)
                {
                    n_lag += int(geom[iii].FastGetSolutionStepValue(IS_LAGRANGIAN_INLET));
                    n_str += int(geom[iii].FastGetSolutionStepValue(IS_STRUCTURE));
                    n_fl += int(geom[iii].FastGetSolutionStepValue(IS_FLUID));
                    n_interf += int(geom[iii].FastGetSolutionStepValue(IS_INTERFACE));
                }


                //if three nodes are at the wall, we check if the fourth node is close to it or not by passing the alpha-shape
                //if (geom.size()==4.0 && n_str==3.0 && n_fl==4.0)
                if (geom.size()==4 && n_interf==3 && n_lag==0)
                {
                    //if alpha shape tells to preserve
                    if (AlphaShape3D(alpha_shape, geom)==false)
                    {
                        for (unsigned int iii=0; iii<geom.size(); iii++)
                        {
                            //if thats not the wall node, remove it
                            if (geom[iii].FastGetSolutionStepValue(IS_STRUCTURE)==0.0)
                            {
                                geom[iii].Set(TO_ERASE,true);
                                KRATOS_WATCH("NODE CLOSE TO THE WALL - WILL BE ERASED!!!!")
                            }
                        }
                    }
                }


            }
        }



    }



    void MarkLonelyNodesForErasing(ModelPart& ThisModelPart)
    {

//delete lonely nodes //just to try: 19/07/2010
        /*
        for(ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin();
        	in!=ThisModelPart.NodesEnd(); in++)
        {
        in->Set(TO_ERASE,false);

        }
        */

        for(ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin();
                in!=ThisModelPart.NodesEnd(); in++)
        {
            if((in->GetValue(NEIGHBOUR_ELEMENTS)).size() == 0 && in->FastGetSolutionStepValue(IS_STRUCTURE)==0.0 && in->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET)!=1 && in->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET,1)!=1.0)
            {
                in->Set(TO_ERASE,true);
                KRATOS_WATCH("Marking lonelynodes!!!")
            }

        }
        /*
        for(ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin();
        	in!=ThisModelPart.NodesEnd(); in++)
        {
        	if((in->GetValue(NEIGHBOUR_ELEMENTS)).size() == 0 && in->FastGetSolutionStepValue(IS_STRUCTURE)==1.0)
        		{
        		in->FastGetSolutionStepValue(PRESSURE)==0;
        		if ((in)->GetDof(DISPLACEMENT_X).IsFixed())
        			in->FastGetSolutionStepValue(VELOCITY_X)==0;

        		if ((in)->GetDof(DISPLACEMENT_Y).IsFixed())
        			in->FastGetSolutionStepValue(VELOCITY_Y)==0;

        		if ((in)->GetDof(DISPLACEMENT_Z).IsFixed())
        			in->FastGetSolutionStepValue(VELOCITY_Z)==0;


        		}

        }
        */
    }

    bool AlphaShape(double alpha_param, Geometry<Node<3> >& pgeom)
    //bool AlphaShape(double alpha_param, Triangle2D<Node<3> >& pgeom)
    {
        KRATOS_TRY
        BoundedMatrix<double,2,2> J; //local jacobian
        BoundedMatrix<double,2,2> Jinv; //local jacobian
        static array_1d<double,2> c; //center pos
        static array_1d<double,2> rhs; //center pos

        double x0 = pgeom[0].X();
        double x1 = pgeom[1].X();
        double x2 = pgeom[2].X();

        double y0 = pgeom[0].Y();
        double y1 = pgeom[1].Y();
        double y2 = pgeom[2].Y();

        J(0,0)=2.0*(x1-x0);
        J(0,1)=2.0*(y1-y0);
        J(1,0)=2.0*(x2-x0);
        J(1,1)=2.0*(y2-y0);


        double detJ = J(0,0)*J(1,1)-J(0,1)*J(1,0);

        Jinv(0,0) =  J(1,1);
        Jinv(0,1) = -J(0,1);
        Jinv(1,0) = -J(1,0);
        Jinv(1,1) =  J(0,0);

        bounded_matrix<double,2,2> check;


        if(detJ < 1e-12)
        {
            //std::cout << "detJ = " << detJ << std::endl;
            ////mark as boundary
            pgeom[0].GetSolutionStepValue(IS_BOUNDARY) = 1;
            pgeom[1].GetSolutionStepValue(IS_BOUNDARY) = 1;
            pgeom[2].GetSolutionStepValue(IS_BOUNDARY) = 1;
            return false;
        }

        else
        {

            double x0_2 = x0*x0 + y0*y0;
            double x1_2 = x1*x1 + y1*y1;
            double x2_2 = x2*x2 + y2*y2;

            //finalizing the calculation of the inverted matrix
            //std::cout<<"MATR INV"<<MatrixInverse(msJ)<<std::endl;
            Jinv /= detJ;
            //calculating the RHS
            rhs[0] = (x1_2 - x0_2);
            rhs[1] = (x2_2 - x0_2);

            //calculate position of the center
            noalias(c) = prod(Jinv,rhs);

            double radius = sqrt(pow(c[0]-x0,2)+pow(c[1]-y0,2));

            //calculate average h
            double h;
            h =  pgeom[0].FastGetSolutionStepValue(NODAL_H);
            h += pgeom[1].FastGetSolutionStepValue(NODAL_H);
            h += pgeom[2].FastGetSolutionStepValue(NODAL_H);
            h *= 0.333333333;
            if (radius < h*alpha_param)
            {
                return true;
            }
            else
            {
                return false;
            }
        }


        KRATOS_CATCH("")
    }
    bool AlphaShape3D( double alpha_param, Geometry<Node<3> >& geom	)
    {
        KRATOS_TRY

        BoundedMatrix<double,3,3> J; //local jacobian
        BoundedMatrix<double,3,3> Jinv; //local jacobian
        array_1d<double,3> Rhs; //rhs
        array_1d<double,3> xc;
        double radius=0.0;

        const double x0 = geom[0].X();
        const double y0 = geom[0].Y();
        const double z0 = geom[0].Z();
        const double x1 = geom[1].X();
        const double y1 = geom[1].Y();
        const double z1 = geom[1].Z();
        const double x2 = geom[2].X();
        const double y2 = geom[2].Y();
        const double z2 = geom[2].Z();
        const double x3 = geom[3].X();
        const double y3 = geom[3].Y();
        const double z3 = geom[3].Z();

        //calculation of the jacobian
        J(0,0) = x1-x0;
        J(0,1) = y1-y0;
        J(0,2) = z1-z0;
        J(1,0) = x2-x0;
        J(1,1) = y2-y0;
        J(1,2) = z2-z0;
        J(2,0) = x3-x0;
        J(2,1) = y3-y0;
        J(2,2) = z3-z0;
// 			KRATOS_WATCH("33333333333333333333333");

        //inverse of the jacobian
        //first column
        Jinv(0,0) = J(1,1)*J(2,2) - J(1,2)*J(2,1);
        Jinv(1,0) = -J(1,0)*J(2,2) + J(1,2)*J(2,0);
        Jinv(2,0) = J(1,0)*J(2,1) - J(1,1)*J(2,0);
        //second column
        Jinv(0,1) = -J(0,1)*J(2,2) + J(0,2)*J(2,1);
        Jinv(1,1) = J(0,0)*J(2,2) - J(0,2)*J(2,0);
        Jinv(2,1) = -J(0,0)*J(2,1) + J(0,1)*J(2,0);
        //third column
        Jinv(0,2) = J(0,1)*J(1,2) - J(0,2)*J(1,1);
        Jinv(1,2) = -J(0,0)*J(1,2) + J(0,2)*J(1,0);
        Jinv(2,2) = J(0,0)*J(1,1) - J(0,1)*J(1,0);
        //calculation of determinant (of the input matrix)

// 			KRATOS_WATCH("44444444444444444444444444");
        double detJ = J(0,0)*Jinv(0,0)
                      + J(0,1)*Jinv(1,0)
                      + J(0,2)*Jinv(2,0);

        //volume = detJ * 0.16666666667;
// 			KRATOS_WATCH("55555555555555555555555");


        double x0_2 = x0*x0 + y0*y0 + z0*z0;
        double x1_2 = x1*x1 + y1*y1 + z1*z1;
        double x2_2 = x2*x2 + y2*y2 + z2*z2;
        double x3_2 = x3*x3 + y3*y3 + z3*z3;

        //finalizing the calculation of the inverted matrix
        Jinv /= detJ;

        //calculating the RHS
        //calculating the RHS
        Rhs[0] = 0.5*(x1_2 - x0_2);
        Rhs[1] = 0.5*(x2_2 - x0_2);
        Rhs[2] = 0.5*(x3_2 - x0_2);

        //calculate position of the center
        noalias(xc) = prod(Jinv,Rhs);
        //calculate radius
        radius = pow(xc[0] - x0,2);
        radius		  += pow(xc[1] - y0,2);
        radius		  += pow(xc[2] - z0,2);
        radius = sqrt(radius);

        double h;
        h =  geom[0].FastGetSolutionStepValue(NODAL_H);
        h += geom[1].FastGetSolutionStepValue(NODAL_H);
        h += geom[2].FastGetSolutionStepValue(NODAL_H);
        h += geom[3].FastGetSolutionStepValue(NODAL_H);
        h *= 0.250;

        if (radius < h*alpha_param)
        {
            return true;
        }
        else
        {
            return false;
        }

        KRATOS_CATCH("")
    }


    //**********************************************************************************************
    //**********************************************************************************************
    //imposes the velocity that corresponds to a
    void MoveLonelyNodes(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;
        //KRATOS_WATCH("MOVING LONELY NODES")
        double Dt = ThisModelPart.GetProcessInfo()[DELTA_TIME];

        array_1d<double,3> DeltaDisp, acc;

        //const array_1d<double,3> body_force = ThisModelPart.ElementsBegin()->GetProperties()[BODY_FORCE];
        for(ModelPart::NodeIterator i = ThisModelPart.NodesBegin() ;
                i != ThisModelPart.NodesEnd() ; ++i)
        {
            if(
                (i)->FastGetSolutionStepValue(IS_STRUCTURE) == 0 && //if it is not a wall node
                (i)->GetValue(NEIGHBOUR_ELEMENTS).size() == 0 &&//and it is lonely
                ((i)->GetDof(DISPLACEMENT_X).IsFixed() == false || (i)->GetDof(DISPLACEMENT_Y).IsFixed() == false || (i)->GetDof(DISPLACEMENT_Z).IsFixed() == false) //and not the node of the wall, where the 0-displ is prescribed

            )
            {
                //i->Set(TO_ERASE,true);
                //set to zero the pressure
                (i)->FastGetSolutionStepValue(PRESSURE) = 0;

                const array_1d<double,3>& old_vel = (i)->FastGetSolutionStepValue(VELOCITY,1);
                array_1d<double,3>& vel = (i)->FastGetSolutionStepValue(VELOCITY);
                //array_1d<double,3>& acc = (i)->FastGetSolutionStepValue(ACCELERATION);

                noalias(acc) =  (i)->FastGetSolutionStepValue(BODY_FORCE);

                noalias(vel) = old_vel;
                noalias(vel) += Dt * acc ;

                //calculate displacements
                noalias(DeltaDisp) = Dt * vel;

                array_1d<double,3>& disp = i->FastGetSolutionStepValue(DISPLACEMENT);
                noalias(disp) = i->FastGetSolutionStepValue(DISPLACEMENT,1);
                noalias(disp) += DeltaDisp;


            }

        }

        KRATOS_CATCH("")
    }
    //////////////////////
    /*
    void MarkLonelyNodesForErasing(ModelPart& ThisModelPart)
    {
    	KRATOS_TRY;

    	for(ModelPart::NodeIterator i = ThisModelPart.NodesBegin() ;
    		i != ThisModelPart.NodesEnd() ; ++i)
    	{
    		if(
    			(i)->FastGetSolutionStepValue(IS_STRUCTURE) == 0 && //if it is not a wall node
    			(i)->GetValue(NEIGHBOUR_ELEMENTS).size() == 0 &&//and it is lonely
    			((i)->GetDof(DISPLACEMENT_X).IsFixed() == false || (i)->GetDof(DISPLACEMENT_Y).IsFixed() == false || (i)->GetDof(DISPLACEMENT_Z).IsFixed() == false) //and not the node of the wall, where the 0-displ is prescribed

    			)
    		{

    			i->Set(TO_ERASE,true);



    		}

    	}

    	KRATOS_CATCH("")
    }
    */



    //**********************************************************************************************
    //**********************************************************************************************
    //ATTENTION:: returns a negative volume if inverted elements appear
    double CalculateVolume(ModelPart& ThisModelPart, int domain_size)
    {
        KRATOS_TRY;
        KRATOS_WATCH("ENTERED calc vol")
        bool inverted_elements = false;
        //auxiliary vectors
        double Atot = 0.00;
        double Ael;
        //line added on 28th August in order to be able to use membrane elements, i.e. 3 noded elements in 3D
        unsigned int nstruct=0;

        if(domain_size == 2)
        {
            for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i!=ThisModelPart.ElementsEnd(); i++)
            {
                //calculating shape functions values
                Geometry< Node<3> >& geom = i->GetGeometry();
                for (unsigned int kkk=0; kkk<(i->GetGeometry()).size(); kkk++)
                {
                    //n_struct += int( im->GetGeometry()[i].FastGetSolutionStepValue(IS_STRUCTURE) );
                    nstruct+= int(i->GetGeometry()[kkk].FastGetSolutionStepValue(IS_STRUCTURE));
                }

                if (   nstruct!= ( i->GetGeometry()).size() )
                {
                    Ael = GeometryUtils::CalculateVolume2D(geom);
                    Atot += Ael;
                    if(Ael <= 0)
                    {
                        inverted_elements = true;
                    }
                }
                nstruct=0;
            }

        }
        else if (domain_size == 3)
        {
            for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i!=ThisModelPart.ElementsEnd(); i++)
            {
                //calculating shape functions values
                Geometry< Node<3> >& geom = i->GetGeometry();
                for (unsigned int kkk=0; kkk<(i->GetGeometry()).size(); kkk++)
                {
                    //n_struct += int( im->GetGeometry()[i].FastGetSolutionStepValue(IS_STRUCTURE) );
                    nstruct+= int(i->GetGeometry()[kkk].FastGetSolutionStepValue(IS_STRUCTURE));
                }


                //if (   (i->GetGeometry()).size()==4 )
                if (   nstruct!= ( i->GetGeometry()).size() )
                {
                    Ael = GeometryUtils::CalculateVolume3D(geom);
                    Atot += Ael;
                    if(Ael <= 0) inverted_elements = true;
                }
                nstruct=0;
            }

        }

        //set to negative the volume if inverted elements appear
        if( inverted_elements == true)
            Atot = -Atot;
        //return the total area
        KRATOS_WATCH("FINISHED calc vol LAST!!!!!!!! ")
        KRATOS_WATCH(Atot)
        return Atot;

        KRATOS_CATCH("");
    }

    //**********************************************************************************************
    //**********************************************************************************************
    void ReduceTimeStep(ModelPart& ThisModelPart, const double reduction_factor)
    {
        KRATOS_TRY;
        KRATOS_WATCH("INSIDE -REDUCE TIME STEP- process")
        double old_dt = ThisModelPart.GetProcessInfo()[DELTA_TIME];
        double new_dt = reduction_factor * old_dt;
        double time = ThisModelPart.GetProcessInfo()[TIME];

        double new_time =  time - old_dt + new_dt;
        KRATOS_WATCH(time);
        KRATOS_WATCH(new_time);
        (ThisModelPart.GetProcessInfo()).SetCurrentTime(new_time);
        //and now we set the data values to the ones of the previous time

        unsigned int step_data_size = ThisModelPart.GetNodalSolutionStepDataSize();
        KRATOS_WATCH(step_data_size)
        for(ModelPart::NodesContainerType::iterator i = ThisModelPart.NodesBegin();
                i!=ThisModelPart.NodesEnd(); i++)
        {
            //unsigned int buffer_size = i->GetBufferSize();


            //getting the data of the solution step
            double* step_data = i->SolutionStepData().Data(0);
            double* prev_step_data = i->SolutionStepData().Data(1);


            //setting it to the old one (we will intend solve the problem now with the reduced step), for every nodal variable
            for(unsigned int j= 0; j<step_data_size; j++)
            {
                step_data[j] = prev_step_data[j];
            }

        }

        KRATOS_CATCH("");
    }

    //**********************************************************************************************
    //**********************************************************************************************
    void SaveNodalArea(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;

        for(ModelPart::NodesContainerType::iterator i = ThisModelPart.NodesBegin();
                i!=ThisModelPart.NodesEnd(); i++)
        {
            i->FastGetSolutionStepValue(NODAL_AREA,1) = i->FastGetSolutionStepValue(NODAL_AREA);
        }

        KRATOS_CATCH("");
    }

    //**********************************************************************************************
    //**********************************************************************************************

    void SaveReducedPart(ModelPart& full_model_part, ModelPart& reduced_model_part)
    {
    KRATOS_TRY

    /////////////////////////////////////////////////////////////////////////
    //clear reduced_model_part
    reduced_model_part.Conditions().clear();
    reduced_model_part.Elements().clear();
    reduced_model_part.Nodes().clear();
    // unsigned int NumberOfProperties = full_model_part.NumberOfProperties();

    //change the name of the var to another one  - otherwise confusing
    for(ModelPart::NodesContainerType::iterator in = full_model_part.NodesBegin() ; in != full_model_part.NodesEnd() ; ++in)
    {
    in->FastGetSolutionStepValue(ACTIVATION_LEVEL)=false;
    in->FastGetSolutionStepValue(MATERIAL_VARIABLE)=false;
    }

	int n_nodes=0;

        for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; im != full_model_part.ElementsEnd() ; ++im)
        {
	    int n_int=0;
	    int n_nodes=im->GetGeometry().size();

	    for (int i=0; i<n_nodes; i++)
        {

            if(im->GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)==1) n_int+=1.0;//im->GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE);
		}
        if (n_int==n_nodes)
        {

            if(n_nodes==4)
            {
                im->GetGeometry()[0].FastGetSolutionStepValue(MATERIAL_VARIABLE)=true;
                im->GetGeometry()[1].FastGetSolutionStepValue(MATERIAL_VARIABLE)=true;
                im->GetGeometry()[2].FastGetSolutionStepValue(MATERIAL_VARIABLE)=true;
                im->GetGeometry()[3].FastGetSolutionStepValue(MATERIAL_VARIABLE)=true;
            }
            else
            {
                im->GetGeometry()[0].FastGetSolutionStepValue(MATERIAL_VARIABLE)=true;
                im->GetGeometry()[1].FastGetSolutionStepValue(MATERIAL_VARIABLE)=true;
                im->GetGeometry()[2].FastGetSolutionStepValue(MATERIAL_VARIABLE)=true;
                //im->GetGeometry()[3].FastGetSolutionStepValue(MATERIAL_VARIABLE)=true;

            }

	    }

 	if (n_int<n_nodes)
            {
		reduced_model_part.Elements().push_back(*(im.base()));

		for (int i=0;i<n_nodes;i++)
			{
			im->GetGeometry()[i].FastGetSolutionStepValue(ACTIVATION_LEVEL)=true;
			}
            }

            if (n_int>n_nodes)
                KRATOS_ERROR<<"Number of DISABLE flags cant exceed number of the element nodes.... ";

        }



    for(ModelPart::NodesContainerType::iterator in = full_model_part.NodesBegin() ; in != full_model_part.NodesEnd() ; ++in)
        {
            int n_disabled=in->FastGetSolutionStepValue(ACTIVATION_LEVEL);
            if (n_disabled==1)
            {
                reduced_model_part.Nodes().push_back(*(in.base()));
            }

        }

        for(ModelPart::PropertiesContainerType::iterator i_properties = full_model_part.PropertiesBegin() ;
                i_properties != full_model_part.PropertiesEnd() ; ++i_properties)
        {
            reduced_model_part.AddProperties(*(i_properties.base()));

        }





	//find neighbors within reduced model part
	if (n_nodes==3)
		{
		FindNodalNeighboursProcess N_FINDER(reduced_model_part);
		N_FINDER.Execute();
		}
	else if (n_nodes==4)
		{
		FindNodalNeighboursProcess N_FINDER(reduced_model_part);
		N_FINDER.Execute();
		}



        KRATOS_CATCH("")
    }




private:

    //aux vars
    static BoundedMatrix<double,3,3> msJ; //local jacobian
    static BoundedMatrix<double,3,3> msJinv; //inverse jacobian
    static array_1d<double,3> msc; //center pos
    static array_1d<double,3> ms_rhs; //center pos


};

}  // namespace Kratos.

#endif // KRATOS_ULF_UTILITIES_INCLUDED  defined
