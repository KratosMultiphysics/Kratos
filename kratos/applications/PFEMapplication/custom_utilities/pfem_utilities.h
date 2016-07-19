/*
==============================================================================
KratosPFEMApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


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

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2008-11-10 14:23:32 $
//   Revision:            $Revision: 1.9 $
//
//


#if !defined(KRATOS_PFEM_UTILITIES_INCLUDED )
#define  KRATOS_PFEM_UTILITIES_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes
#include "includes/kratos_flags.h"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "PFEM_application.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{

class PfemUtils
{
public:

    //**********************************************************************************************
    //**********************************************************************************************

    //functions to apply the boundary conditions

    void ApplyBoundaryConditions(ModelPart& ThisModelPart, int laplacian_type)
    {
        KRATOS_TRY;

        if (laplacian_type == 1)
        {
            //identify fluid nodes and mark as boundary the free ones
            for (ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin(); in != ThisModelPart.NodesEnd(); in++)
            {
                //marking wet nodes
                if (in->FastGetSolutionStepValue(IS_STRUCTURE))
                {
                    if ((in->GetValue(NEIGHBOUR_ELEMENTS)).size() != 0)
                        in->FastGetSolutionStepValue(IS_FLUID) = 1.0;
                    else //it is not anymore of fluid
                        in->FastGetSolutionStepValue(IS_FLUID) = 0.0;
                }
                //marking as free surface the lonely nodes
                // 					else
                if ((in->GetValue(NEIGHBOUR_ELEMENTS)).size() == 0)
                    in->FastGetSolutionStepValue(IS_BOUNDARY) = 1.0;
            }

            //identify the free surface
            for (ModelPart::NodeIterator i = ThisModelPart.NodesBegin();
                    i != ThisModelPart.NodesEnd(); ++i)
            {
                //reset the free surface
                i->FastGetSolutionStepValue(IS_FREE_SURFACE) = 0;
                i->Free(PRESSURE);

                //identify the free surface and fix the pressure accordingly
                if (i->FastGetSolutionStepValue(IS_BOUNDARY) != 0
                        &&
                        i->FastGetSolutionStepValue(IS_STRUCTURE) == 0)
                {
                    i->FastGetSolutionStepValue(IS_FREE_SURFACE) = 1.0;
                    i->FastGetSolutionStepValue(PRESSURE) = 0.00;
                    i->Fix(PRESSURE);
                }
            }


        }
        else   //case of discrete laplacian
        {
            //identify fluid nodes and mark as boundary the free ones
            for (ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin(); in != ThisModelPart.NodesEnd(); in++)
            {
                in->Free(PRESSURE);
                //marking wet nodes
                if (in->FastGetSolutionStepValue(IS_STRUCTURE))
                {
                    if ((in->GetValue(NEIGHBOUR_ELEMENTS)).size() != 0)
                        in->FastGetSolutionStepValue(IS_FLUID) = 1.0;
                    else //it is not anymore of fluid
                        in->FastGetSolutionStepValue(IS_FLUID) = 0.0;
                }
                //marking as free surface the lonely nodes
                // 					else
                if ((in->GetValue(NEIGHBOUR_ELEMENTS)).size() == 0)
                    in->FastGetSolutionStepValue(IS_BOUNDARY) = 1.0;

                if ((in->GetValue(NEIGHBOUR_ELEMENTS)).size() == 0)
                    in->FastGetSolutionStepValue(PRESSURE) = 0.0;
            }

            //identify the free surface
            for (ModelPart::NodeIterator i = ThisModelPart.NodesBegin();
                    i != ThisModelPart.NodesEnd(); ++i)
            {
                //reset the free surface
                i->FastGetSolutionStepValue(IS_FREE_SURFACE) = 0;

                //identify the free surface and fix the pressure accordingly
                if (i->FastGetSolutionStepValue(IS_BOUNDARY) != 0
                        &&
                        i->FastGetSolutionStepValue(IS_STRUCTURE) == 0)
                {
                    i->FastGetSolutionStepValue(IS_FREE_SURFACE) = 1.0;
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

        for (ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                i != ThisModelPart.ElementsEnd(); i++)
        {
            //calculating shape functions values
            Geometry< Node < 3 > >& geom = i->GetGeometry();

            for (unsigned int i = 0; i < geom.size(); i++)
                geom[i].FastGetSolutionStepValue(IS_FLUID) = 1;

        }
        KRATOS_CATCH("");
    }

    //**********************************************************************************************
    //**********************************************************************************************

    void ApplyMinimalPressureConditions(ModelPart& r_model_part)
    //		void ApplyMinimalPressureConditions(std::vector< WeakPointerVector< Node<3> > >& connected_components)
    {
        KRATOS_TRY;

        /*			//verify that the minimal conditions for pressure are applied
                                std::vector<bool> to_be_prescribed(connected_components.size());
                                array_1d<double,3> fixed_vel;
                                for(unsigned int i = 0; i<connected_components.size(); i++)
                                {
                                        int boundary_nodes = 0;
                                        int prescribed_vel_nodes = 0;
                                        WeakPointerVector< Node<3> >& node_list = connected_components[i];
                                        for( WeakPointerVector< Node<3> >::iterator in = node_list.begin();
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
                                                for( WeakPointerVector< Node<3> >::iterator in = node_list.begin();
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
                                }*/


        typedef Node < 3 > PointType;
        typedef PointerVector<PointType > PointVector;
        //typedef PointVector::iterator PointIterator;

        //generate a container with the layers to be extrapolated
        PointVector work_array;
        work_array.reserve(r_model_part.Nodes().size());

        //reset the counter
        for (ModelPart::NodesContainerType::iterator inode = r_model_part.NodesBegin();
                inode != r_model_part.NodesEnd();
                inode++)
        {
            inode->GetValue(IS_VISITED) = 0.0;
        }

        double connected_component_index = 1;
        // 		PointIterator current_it = work_array.begin();
        unsigned int aux_index = 0;
        unsigned int number_of_connected_components = 0;
        for (ModelPart::NodesContainerType::iterator inode = r_model_part.NodesBegin();
                inode != r_model_part.NodesEnd();
                inode++)
        {

            if (inode->GetValue(IS_VISITED) == 0.0 && (inode->GetValue(NEIGHBOUR_NODES)).size() != 0) //this node is not yet visited - it is in a new
            {
                number_of_connected_components++;
                unsigned int n_free_boundaries = 0;

                //add the node to the work array - this node will be the "seed" of the partition
                inode->GetValue(IS_VISITED) = connected_component_index;
                work_array.push_back(Node < 3 > ::Pointer(*(inode.base())));

                while (work_array.begin() + aux_index < work_array.end())
                {

                    if ((work_array.begin() + aux_index)->FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0)
                    {
                        // KRATOS_WATCH("sono qui");
                        n_free_boundaries++;
                    }

                    //add its neighbours to the work array and mark them
                    WeakPointerVector< Node < 3 > >& neighb_nodes = (work_array.begin() + aux_index)->GetValue(NEIGHBOUR_NODES);
                    // KRATOS_WATCH(neighb_nodes.size() );
                    for (WeakPointerVector< Node < 3 > >::iterator i = neighb_nodes.begin(); i != neighb_nodes.end(); i++)
                    {
                        if (i->GetValue(IS_VISITED) == 0.0)
                        {
                            //add the node to the work array
                            i->GetValue(IS_VISITED) = connected_component_index;
                            work_array.push_back(Node < 3 > ::Pointer(*(i.base())));
                        }

                    }

                    //update the current counter
                    aux_index++;
                }
                // 				KRATOS_WATCH("ccc")

                //fix at least one pressure if needed (the pressure will be fixed on the "seed" node
                if (n_free_boundaries == 0)
                {
                    std::cout << "minimal pressure condition prescribed on node " << inode->Id() << std::endl;
                    inode->Fix(PRESSURE);
                    inode->FastGetSolutionStepValue(PRESSURE) = 0.0;
                }

                //increase the component index counter
                connected_component_index++;
            }


        }
        std::cout << "number of connected components" << number_of_connected_components << std::endl;
        KRATOS_CATCH("");
    }

    //**********************************************************************************************
    //**********************************************************************************************

    double CFLdeltaT(double CFL, double max_dt, ModelPart& ThisModelPart)
    {
        KRATOS_TRY

        array_1d<double, 3 > vel, pt0, pt1;
        double max_dt_cfl = max_dt, a1, b1, a2, b2, X_target, Y_target, top_x, top_y, flag = 0.0;



        for (ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                i != ThisModelPart.ElementsEnd(); i++)
        {
            Geometry< Node < 3 > >& geom = i->GetGeometry();
            double this_node_dt;

            for (int ii = 0; ii < 3; ++ii)
            {
                vel = geom[ii].FastGetSolutionStepValue(VELOCITY);
                double vel_square = inner_prod(vel, vel);
                if (vel_square != 0.0 && geom[ii].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0)
                {
                    a1 = a2 = b1 = b2 = 0.0;
                    top_x = geom[ii].X();
                    top_y = geom[ii].Y();
                    flag = 0.0;
                    //take line points
                    if (ii == 0)
                    {
                        pt1[0] = geom[1].X();
                        pt1[1] = geom[1].Y();
                        pt1[2] = geom[1].Z();


                        pt0[0] = geom[2].X();
                        pt0[1] = geom[2].Y();
                        pt0[2] = geom[2].Z();
                    }
                    else if (ii == 1)
                    {
                        pt1[0] = geom[0].X();
                        pt1[1] = geom[0].Y();
                        pt1[2] = geom[0].Z();


                        pt0[0] = geom[2].X();
                        pt0[1] = geom[2].Y();
                        pt0[2] = geom[2].Z();
                    }
                    else
                    {
                        pt1[0] = geom[1].X();
                        pt1[1] = geom[1].Y();
                        pt1[2] = geom[1].Z();


                        pt0[0] = geom[0].X();
                        pt0[1] = geom[0].Y();
                        pt0[2] = geom[0].Z();
                    }

                    //check possible combination
                    if (vel[0] == 0.0)
                    {
                        if (pt1[0] == pt0[0]) flag = 1.0;
                        else
                        {
                            a2 = (pt1[1] - pt0[1]) / (pt1[0] - pt0[0]);
                            b2 = (pt1[0] * pt0[1] - pt1[1] * pt0[0]) / (pt1[0] - pt0[0]);
                            X_target = top_x;
                            Y_target = a2 * X_target + b2;
                        }
                    }
                    else
                    {
                        a1 = vel[1] / vel[0];
                        b1 = top_y - a1*top_x;

                        if (pt1[0] == pt0[0])
                        {
                            X_target = pt0[0];
                            Y_target = a1 * X_target + b1;
                        }
                        else
                        {
                            a2 = (pt1[1] - pt0[1]) / (pt1[0] - pt0[0]);
                            b2 = (pt1[0] * pt0[1] - pt1[1] * pt0[0]) / (pt1[0] - pt0[0]);

                            if (a1 == a2) flag = 1.0;
                            else
                            {
                                X_target = (b2 - b1) / (a1 - a2);
                                Y_target = a2 * X_target + b2;
                            }
                        }
                    }

                    if (flag == 0.0)
                    {
                        double l_square = (top_x - X_target)*(top_x - X_target) + (top_y - Y_target)*(top_y - Y_target);
                        this_node_dt = l_square / vel_square;
                        this_node_dt = CFL * sqrt(this_node_dt);
                    }
                    else this_node_dt = 1000000;


                    if (this_node_dt < max_dt_cfl)
                    {
                        max_dt_cfl = this_node_dt;

                    }

                }


            }//end of loop over element nodes

        }//end of loop over elements
        return max_dt_cfl;

        KRATOS_CATCH("");
    }
    //**********************************************************************************************
    //**********************************************************************************************

    double ExactDtEstimate(double dt_max, ModelPart& ThisModelPart)
    {
        KRATOS_TRY

        array_1d<double, 3 > dxprim12, dx32, dxprim32, dv12, dv32, dvprim12;
        double deltatime = dt_max;
        double B, A, C;


        for (ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                i != ThisModelPart.ElementsEnd(); i++)
        {
            //calculating velocity and displacement vectors
            Geometry< Node < 3 > >& geom = i->GetGeometry();
            noalias(dv32) = geom[2].FastGetSolutionStepValue(VELOCITY);
            dv32 -= geom[1].FastGetSolutionStepValue(VELOCITY);

            noalias(dv12) = geom[0].FastGetSolutionStepValue(VELOCITY);
            dv12 -= geom[1].FastGetSolutionStepValue(VELOCITY);

            dvprim12[0] = dv12[1];
            dvprim12[1] = -1 * dv12[0];

            dxprim12[0] = geom[0].Y() - geom[1].Y();
            dxprim12[1] = geom[1].X() - geom[0].X();

            dx32[0] = geom[2].X() - geom[1].X();
            dx32[1] = geom[2].Y() - geom[1].Y();

            dxprim32[0] = geom[1].Y() - geom[2].Y();
            dxprim32[1] = geom[2].X() - geom[1].X();

            //Calculate A,B and C to calculate zero´s of dt
            B = inner_prod(dv32, dxprim12);
            B += inner_prod(dv12, dxprim32);

            A = inner_prod(dv32, dvprim12);

            C = inner_prod(dx32, dxprim12);

            double zerodt = dt_max;
            double Zdt1;
            double Zdt2;

            if (A == 0.0)
            {
                if (B != 0.0) zerodt = -1.0 * C / B;

            }
            else
            {
                double delta = 0.0;
                delta = B * B - 4.0 * A*C;
                if (delta > 0.0)
                {
                    Zdt1 = (-1 * B + sqrt(delta)) / (2.0 * A);
                    Zdt2 = (-1 * B - sqrt(delta)) / (2.0 * A);

                    //zerodt = -1.0;
                    if (Zdt1 > 0.0) zerodt = Zdt1;
                    if (zerodt > Zdt2 && Zdt2 > 0.0) zerodt = Zdt2;
                }

            }

            //negative race means evert dt is acceptable
            if (zerodt < 0.0) zerodt = dt_max;

            //global check
            if (zerodt < deltatime)
            {
                deltatime = zerodt;
                KRATOS_WATCH("(((((((((((((((INSIDE ESTIMATE)))))))))))))))");
            }
        }
        return deltatime;

        KRATOS_CATCH("")
    }
    //**********************************************************************************************
    //**********************************************************************************************

    double EstimateDeltaTime(double dt_min, double dt_max, ModelPart& ThisModelPart)
    {
        KRATOS_TRY

        array_1d<double, 3 > dx, dv;
        double deltatime = dt_max;
        double dvside, lside;

        for (ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                i != ThisModelPart.ElementsEnd(); i++)
        {
            //calculating shape functions values
            Geometry< Node < 3 > >& geom = i->GetGeometry();

            for (unsigned int i1 = 0; i1 < geom.size() - 1; i1++)
            {
                for (unsigned int i2 = i1 + 1; i2 < geom.size(); i2++)
                {
                    dx[0] = geom[i2].X() - geom[i1].X();
                    dx[1] = geom[i2].Y() - geom[i1].Y();
                    dx[2] = geom[i2].Z() - geom[i1].Z();

                    lside = inner_prod(dx, dx);

                    noalias(dv) = geom[i2].FastGetSolutionStepValue(VELOCITY);
                    noalias(dv) -= geom[i1].FastGetSolutionStepValue(VELOCITY);

                    dvside = inner_prod(dx, dv);

                    double dt;
                    if (dvside < 0.0) //otherwise the side is "getting bigger" so the time step has not to be diminished
                    {
                        dt = fabs(lside / dvside);
                        if (dt < deltatime) deltatime = dt;
                    }

                }
            }
        }

        if (deltatime < dt_min)
        {
            std::cout << "ATTENTION dt_min is being used" << std::endl;
            deltatime = dt_min;
        }
        //ATTENTION: following code is for testing purposes
        else //check that the time step is nowhere bigger than hnode/vnode
        {
            double dt2 = deltatime*deltatime;
            for (ModelPart::NodeIterator i = ThisModelPart.NodesBegin();
                    i != ThisModelPart.NodesEnd(); ++i)
            {
                const array_1d<double, 3 > & v = i->FastGetSolutionStepValue(VELOCITY);
                const double& hnode = i->FastGetSolutionStepValue(NODAL_H);

                double norm_v2 = inner_prod(v, v);

                double dt_est = dt_max;
                if (norm_v2 > 1e-20)
                    dt_est = hnode * hnode / norm_v2;

                if (dt_est < dt2)
                    dt2 = dt_est;

            }
            double final_estimate = sqrt(dt2);
            if (final_estimate < deltatime)
                deltatime = final_estimate;
            if (deltatime < dt_min)
            {
                std::cout << "ATTENTION dt_min is being used" << std::endl;
                deltatime = dt_min;
            }
        }



        return deltatime;

        KRATOS_CATCH("")
    }
    /*		double EstimateDeltaTime(double dt_min, double dt_max, ModelPart::NodesContainerType& rNodes)
                    {
                            KRATOS_TRY

                                    array_1d<double,3> xc, dist;
                            double deltatime = 1000.0;
                            array_1d<double,3> vaux;

                            for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
                            {
                                    const array_1d<double,3>& vc = in->FastGetSolutionStepValue(VELOCITY);

                                    if((in->GetValue(NEIGHBOUR_NODES)).size() != 0)
                                    {
                                            xc[0] = in->X();
                                            xc[1] = in->Y();
                                            xc[2] = in->Z();

                                            double dt_node = 0.00;
                                            for( WeakPointerVector< Node<3> >::iterator i = in->GetValue(NEIGHBOUR_NODES).begin();
                                                    i != in->GetValue(NEIGHBOUR_NODES).end(); i++)
                                            {
                                                    dist[0] = i->X();
                                                    dist[1] = i->Y();
                                                    dist[2] = i->Z();
                                                    dist -= xc;

                                                    noalias(vaux) =  i->FastGetSolutionStepValue(VELOCITY);
                                                    noalias(vaux) -=  vc;
                                                    double vnorm2 = inner_prod(vaux,vaux);

                                                    if(vnorm2 != 0.0)
                                                    {
                                                            double nom = inner_prod(dist,vaux);
                                                            double dt_side = fabs(nom/vnorm2);
                                                            if(dt_side > dt_node)
                                                                    dt_node = dt_side;
                                                            else if(dt_node < deltatime)
                                                                    deltatime = dt_node;
                                                    }
                                            }
                                    }
                            }

                            if(deltatime > dt_max) deltatime = dt_max;
                            else if (deltatime < dt_min) deltatime = dt_min;

                            return deltatime;

                            KRATOS_CATCH("")
                    }*/

    //**********************************************************************************************
    //**********************************************************************************************

    void MarkOuterNodes(const array_1d<double, 3 > & corner1, const array_1d<double, 3 > & corner2,
                        ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY;

        //add a big number to the id of the nodes to be erased
        int n_erased = 0;
        double xmax, xmin, ymax, ymin, zmax, zmin;



        if (corner1[0] > corner2[0])
        {
            xmax = corner1[0];
            xmin = corner2[0];
        }
        else
        {
            xmax = corner2[0];
            xmin = corner1[0];
        }

        if (corner1[1] > corner2[1])
        {
            ymax = corner1[1];
            ymin = corner2[1];
        }
        else
        {
            ymax = corner2[1];
            ymin = corner1[1];
        }

        if (corner1[2] > corner2[2])
        {
            zmax = corner1[2];
            zmin = corner2[2];
        }
        else
        {
            zmax = corner2[2];
            zmin = corner1[2];
        }



        for (ModelPart::NodesContainerType::iterator in = rNodes.begin(); in != rNodes.end(); in++)
        {
            bool erase = false;
            double& x = in->X();
            double& y = in->Y();
            double& z = in->Z();

            if (x < xmin || x > xmax) erase = true;
            else if (y < ymin || y > ymax) erase = true;
            else if (z < zmin || z > zmax) erase = true;

            if (erase == true)
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
        KRATOS_WATCH("INSSSSSSSSSSIDE MARK EXCESSIVELY");

        const double fact2 = admissible_distance_factor*admissible_distance_factor;



        ModelPart::NodeIterator NodesBegin = rNodes.begin();
        int size = rNodes.size();
        ModelPart::NodeIterator in;
        #pragma omp parallel for firstprivate(size,NodesBegin),private(in)
        for (int k = 0; k < size; k++)
        {
            in = NodesBegin + k;

            if (in->FastGetSolutionStepValue(IS_STRUCTURE) == 0) //if it is not a wall node i can erase
            {
                double hnode2 = in->FastGetSolutionStepValue(NODAL_H);
                hnode2 *= hnode2; //take the square

                //loop on neighbours and erase if they are too close
                for (WeakPointerVector< Node < 3 > >::iterator i = in->GetValue(NEIGHBOUR_NODES).begin();
                        i != in->GetValue(NEIGHBOUR_NODES).end(); i++)
                {
                    if (static_cast<bool> (i->Is(TO_ERASE)) == false) //we can erase the current node only if the neighb is not to be erased
                    {
                        double dx = i->X() - in->X();
                        double dy = i->Y() - in->Y();
                        double dz = i->Z() - in->Z();

                        double dist2 = dx * dx + dy * dy + dz*dz;

                        if (dist2 < fact2 * hnode2)
                            in->Set(TO_ERASE, true);
                    }
                }
            }
        }
        // 	    }

        KRATOS_CATCH("")
    }

    //**********************************************************************************************
    //**********************************************************************************************
    // 		bool AlphaShape3D(double alpha_param, Tetrahedra3D4<Node<3> >& pgeom)
    // 		{
    // 			KRATOS_TRY
    //
    // 				double x0 = pgeom[0].X();
    // 			double x1 = pgeom[1].X();
    // 			double x2 = pgeom[2].X();
    // 			double x3 = pgeom[3].X();
    //
    // 			double y0 = pgeom[0].Y();
    // 			double y1 = pgeom[1].Y();
    // 			double y2 = pgeom[2].Y();
    // 			double y3 = pgeom[3].Y();
    //
    // 			double z0 = pgeom[0].Z();
    // 			double z1 = pgeom[1].Z();
    // 			double z2 = pgeom[2].Z();
    // 			double z3 = pgeom[3].Z();
    //
    // 			//calculation of the jacobian
    // 			msJ(0,0) = x1-x0; msJ(0,1) = y1-y0; msJ(0,2) = z1-z0;
    // 			msJ(1,0) = x2-x0; msJ(1,1) = y2-y0; msJ(1,2) = z2-z0;
    // 			msJ(2,0) = x3-x0; msJ(2,1) = y3-y0; msJ(2,2) = z3-z0;
    //
    // 			//inverse of the jacobian
    // 			//first column
    // 			msJinv(0,0) = msJ(1,1)*msJ(2,2) - msJ(1,2)*msJ(2,1);
    // 			msJinv(1,0) = -msJ(1,0)*msJ(2,2) + msJ(1,2)*msJ(2,0);
    // 			msJinv(2,0) = msJ(1,0)*msJ(2,1) - msJ(1,1)*msJ(2,0);
    // 			//second column
    // 			msJinv(0,1) = -msJ(0,1)*msJ(2,2) + msJ(0,2)*msJ(2,1);
    // 			msJinv(1,1) = msJ(0,0)*msJ(2,2) - msJ(0,2)*msJ(2,0);
    // 			msJinv(2,1) = -msJ(0,0)*msJ(2,1) + msJ(0,1)*msJ(2,0);
    // 			//third column
    // 			msJinv(0,2) = msJ(0,1)*msJ(1,2) - msJ(0,2)*msJ(1,1);
    // 			msJinv(1,2) = -msJ(0,0)*msJ(1,2) + msJ(0,2)*msJ(1,0);
    // 			msJinv(2,2) = msJ(0,0)*msJ(1,1) - msJ(0,1)*msJ(1,0);
    // 			//calculation of determinant (of the input matrix)
    //
    // 			double detJ = msJ(0,0)*msJinv(0,0)
    // 				+ msJ(0,1)*msJinv(1,0)
    // 				+ msJ(0,2)*msJinv(2,0);
    //
    // 			//calculate average h
    // 			double h;
    // 			h =  pgeom[0].FastGetSolutionStepValue(NODAL_H);
    // 			h += pgeom[1].FastGetSolutionStepValue(NODAL_H);
    // 			h += pgeom[2].FastGetSolutionStepValue(NODAL_H);
    // 			h += pgeom[3].FastGetSolutionStepValue(NODAL_H);
    // 			h *= 0.25;
    //
    // 			if(detJ < 1e-6 * h*h*h)  //this is a sliver and we will remove it
    // 			{
    // 				return false;
    // 			}
    // 			else
    // 			{
    //
    // 				double x0_2 = x0*x0 + y0*y0 + z0*z0;
    // 				double x1_2 = x1*x1 + y1*y1 + z1*z1;
    // 				double x2_2 = x2*x2 + y2*y2 + z2*z2;
    // 				double x3_2 = x3*x3 + y3*y3 + z3*z3;
    //
    // 				//finalizing the calculation of the inverted matrix
    // 				msJinv /= detJ;
    //
    // 				//calculating the RHS
    // 				ms_rhs[0] = 0.5*(x1_2 - x0_2);
    // 				ms_rhs[1] = 0.5*(x2_2 - x0_2);
    // 				ms_rhs[2] = 0.5*(x3_2 - x0_2);
    //
    // 				//calculate position of the center
    // 				noalias(msc) = prod(msJinv,ms_rhs);
    //
    // 				//calculate radius
    // 				double radius = pow(msc[0] - x0,2);
    // 				radius		  += pow(msc[1] - y0,2);
    // 				radius		  += pow(msc[2] - z0,2);
    // 				radius = sqrt(radius);
    //
    // 				if (radius < h*alpha_param)
    // 				{
    // 					if(detJ < 1e-4 * h*h*h)  //this is a sliver and we will remove it even if it would be acceptable
    // 					{
    // 						std::cout << "removed sliver = " << detJ << std::endl;
    // 						return false;
    // 					}
    // 					else
    // 						return true;
    // 				}
    // 				else
    // 				{
    // 					return false;
    // 				}
    // 			}
    // 			KRATOS_CATCH("")
    // 		}

    //**********************************************************************************************
    //**********************************************************************************************
    //calculates the ale spatial acceleration, taking in account both the spatial and the time dependent
    //componenents acc = (v-vm)*dv/dx

    template< int TDim > void CalculateSpatialALEAcceleration(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;

        //auxiliary vectors
        boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim> DN_DX;
        array_1d<double, TDim + 1 > N;
        array_1d<double, 3 > vaux3;
        array_1d<double, TDim> vgauss;
        array_1d<double, TDim> aux;
        const array_1d<double, TDim> zero = ZeroVector(TDim);
        const array_1d<double, 3 > zero3 = ZeroVector(3);
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> dv_dx;
        double factor = 1.00 / (TDim + 1);
        const int nnodes = TDim + 1;

        //reset the acceleration on the nodes
        ModelPart::NodesContainerType& rNodes = ThisModelPart.Nodes();
        for (ModelPart::NodesContainerType::iterator in = rNodes.begin(); in != rNodes.end(); in++)
        {
            noalias(in->FastGetSolutionStepValue(ACCELERATION)) = zero3;
        }

        //compute the projection (first step
        for (ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                i != ThisModelPart.ElementsEnd(); i++)
        {
            //calculating shape functions values
            Geometry< Node < 3 > >& geom = i->GetGeometry();
            double Volume;

            GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

            //constructing the gradient
            noalias(dv_dx) = ZeroMatrix(TDim, TDim);
            noalias(vaux3) = zero3;
            for (int I = 0; I < nnodes; I++)
            {
                array_1d<double, 3 > & v = geom[I].FastGetSolutionStepValue(VELOCITY);
                noalias(vaux3) += v;
                noalias(vaux3) -= geom[I].FastGetSolutionStepValue(MESH_VELOCITY);
                for (int i = 0; i < TDim; i++)
                    for (int j = 0; j < TDim; j++)
                        dv_dx(i, j) += DN_DX(I, j) * v[i];
            }
            for (int iii = 0; iii < TDim; iii++)
                vgauss[iii] = factor * vaux3[iii];
            //				vgauss *= factor;

            noalias(aux) = prod(dv_dx, vgauss);

            //project vgauss * dv/dx on the nodes
            for (int I = 0; I < nnodes; I++)
            {
                array_1d<double, 3 > & acc = geom[I].FastGetSolutionStepValue(ACCELERATION);
                for (int iii = 0; iii < TDim; iii++)
                    acc[iii] += N[I] * Volume * aux[iii];
                //noalias(acc) += N[I] * Volume * aux;
            }
        }

        //finalize the calculation of the convective acceleration
        for (ModelPart::NodesContainerType::iterator in = rNodes.begin(); in != rNodes.end(); in++)
        {
            array_1d<double, 3 > & acc = in->FastGetSolutionStepValue(ACCELERATION);
            acc *= in->FastGetSolutionStepValue(DENSITY);
            double mass = in->FastGetSolutionStepValue(NODAL_MASS);
            acc /= mass;
        }

        KRATOS_CATCH("")
    }

    //**********************************************************************************************
    //*		//**********************************************************************************************

    void CalculateNodalMass(ModelPart& ThisModelPart, int domain_size)
    {
        //set to zero the nodal area
        for (ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin();
                in != ThisModelPart.NodesEnd(); in++)
        {
            in->FastGetSolutionStepValue(NODAL_MASS) = 0.00;
        }

        if (domain_size == 2)
        {
            double mass = 0.0;
            for (ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i != ThisModelPart.ElementsEnd(); i++)
            {
                //calculating shape functions values
                Geometry< Node < 3 > >& geom = i->GetGeometry();

                double density = geom[0].FastGetSolutionStepValue(DENSITY)
                                 + geom[1].FastGetSolutionStepValue(DENSITY)
                                 + geom[2].FastGetSolutionStepValue(DENSITY);
                density *= 0.333333333333333333333333333;

                mass = GeometryUtils::CalculateVolume2D(geom); //calculating the elemental area
                mass *= 0.333333333333333333333333333 * density;

                geom[0].FastGetSolutionStepValue(NODAL_MASS) += mass;
                geom[1].FastGetSolutionStepValue(NODAL_MASS) += mass;
                geom[2].FastGetSolutionStepValue(NODAL_MASS) += mass;
            }
        }
        else if (domain_size == 3)
        {
            for (ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i != ThisModelPart.ElementsEnd(); i++)
            {
                double mass;
                //calculating shape functions values
                Geometry< Node < 3 > >& geom = i->GetGeometry();

                double density = geom[0].FastGetSolutionStepValue(DENSITY)
                                 + geom[1].FastGetSolutionStepValue(DENSITY)
                                 + geom[2].FastGetSolutionStepValue(DENSITY)
                                 + geom[3].FastGetSolutionStepValue(DENSITY);
                density *= 0.25;

                mass = GeometryUtils::CalculateVolume3D(geom);
                mass *= 0.25 * density;

                geom[0].FastGetSolutionStepValue(NODAL_MASS) += mass;
                geom[1].FastGetSolutionStepValue(NODAL_MASS) += mass;
                geom[2].FastGetSolutionStepValue(NODAL_MASS) += mass;
                geom[3].FastGetSolutionStepValue(NODAL_MASS) += mass;
            }
        }

        //prevent zero mass
        for (ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin();
                in != ThisModelPart.NodesEnd(); in++)
        {
            if (in->FastGetSolutionStepValue(NODAL_MASS) == 0.00)
            {
                in->FastGetSolutionStepValue(NODAL_MASS) = in->FastGetSolutionStepValue(DENSITY);
            }
        }


    }

    //**********************************************************************************************
    //**********************************************************************************************

    void QuasiLagrangianMove(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;

        double dt = ThisModelPart.GetProcessInfo()[DELTA_TIME];
        double accfactor = 0.5 * dt*dt;
        array_1d<double, 3 > acc;
        for (ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin();
                in != ThisModelPart.NodesEnd(); in++)
        {
            if ((in->GetValue(NEIGHBOUR_ELEMENTS)).size() != 0
                    && in->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
            {
                array_1d<double, 3 > & vel = in->FastGetSolutionStepValue(VELOCITY);
                array_1d<double, 3 > & vel_old = in->FastGetSolutionStepValue(VELOCITY, 1);
                //calculating time derivative of velocity
                noalias(acc) = vel;
                noalias(acc) -= vel_old;
                acc /= dt;

                //adding the spatial acceleration to the acceleration
                noalias(acc) += in->FastGetSolutionStepValue(ACCELERATION, 1);

                //calculate displacements
                array_1d<double, 3 > & disp = in->FastGetSolutionStepValue(DISPLACEMENT);
                noalias(disp) = in->FastGetSolutionStepValue(DISPLACEMENT, 1);
                noalias(disp) += dt * vel;
                noalias(disp) += accfactor * acc;

                //calculate corresponding velocity
                noalias(vel) = vel_old;
                noalias(vel) += dt * acc;



                //calculate mesh velocity
                //array_1d<double,3>& mesh_vel = in->FastGetSolutionStepValue(MESH_VELOCITY);
                //noalias(mesh_vel) = disp;
                //noalias(mesh_vel) -= in->FastGetSolutionStepValue(DISPLACEMENT,1);
                //mesh_vel /= dt;
            }
        }

        KRATOS_CATCH("")
    }


    //**********************************************************************************************
    //**********************************************************************************************

    void Predict(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;

        double dt = ThisModelPart.GetProcessInfo()[DELTA_TIME];
        array_1d<double, 3 > acc;
        for (ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin();
                in != ThisModelPart.NodesEnd(); in++)
        {
            if ((in->GetValue(NEIGHBOUR_ELEMENTS)).size() != 0
                    && in->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
            {
                //calculate displacements
                array_1d<double, 3 > & disp = in->FastGetSolutionStepValue(DISPLACEMENT);
                noalias(disp) = in->FastGetSolutionStepValue(DISPLACEMENT, 1);
                noalias(disp) += dt * in->FastGetSolutionStepValue(VELOCITY, 1);
            }
        }

        KRATOS_CATCH("")
    }

    //**********************************************************************************************
    //**********************************************************************************************
    //imposes the velocity that corresponds to a

    void MoveLonelyNodes(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;

        double Dt = ThisModelPart.GetProcessInfo()[DELTA_TIME];

        array_1d<double, 3 > DeltaDisp, acc;

        KRATOS_WATCH("in MoveLonelyNodes")

        //const array_1d<double,3> body_force = ThisModelPart.ElementsBegin()->GetProperties()[BODY_FORCE];
        for (ModelPart::NodeIterator i = ThisModelPart.NodesBegin();
                i != ThisModelPart.NodesEnd(); ++i)
        {
            if (
                (i)->FastGetSolutionStepValue(IS_STRUCTURE) == 0 && //if it is not a wall node
                (i)->GetValue(NEIGHBOUR_ELEMENTS).size() == 0 //and it is lonely
            )
            {
                //                                     KRATOS_WATCH(i->Id());
                //set to zero the pressure
                (i)->FastGetSolutionStepValue(PRESSURE) = 0;

                const array_1d<double, 3 > & old_vel = (i)->FastGetSolutionStepValue(VELOCITY, 1);
                array_1d<double, 3 > & vel = (i)->FastGetSolutionStepValue(VELOCITY);
                //array_1d<double,3>& acc = (i)->FastGetSolutionStepValue(ACCELERATION);

                noalias(acc) = (i)->FastGetSolutionStepValue(BODY_FORCE);

                noalias(vel) = old_vel;
                noalias(vel) += Dt * acc;

                //calculate displacements
                noalias(DeltaDisp) = Dt * vel;

                array_1d<double, 3 > & disp = i->FastGetSolutionStepValue(DISPLACEMENT);
                noalias(disp) = i->FastGetSolutionStepValue(DISPLACEMENT, 1);
                noalias(disp) += DeltaDisp;

            }

        }

        KRATOS_CATCH("")
    }
    //*********************************************************************************************
    //*********************************************************************************************

    void MarkNodesTouchingWall(ModelPart& ThisModelPart, int domain_size, double factor)
    {
        KRATOS_TRY;

        if (domain_size == 2)
        {

            for (ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i != ThisModelPart.ElementsEnd(); i++)
            {
                double n_str = 0;

                //counting number on nodes at the wall
                Geometry< Node < 3 > >& geom = i->GetGeometry();
                n_str = geom[0].FastGetSolutionStepValue(IS_BOUNDARY);
                n_str += geom[1].FastGetSolutionStepValue(IS_BOUNDARY);
                n_str += geom[2].FastGetSolutionStepValue(IS_BOUNDARY);
                //if two walls are at the wall, we check if the third node is close to it or not by passing the alpha-shape
                if (n_str == 2.0)
                {
                    boost::numeric::ublas::bounded_matrix<double, 3, 2 > sort_coord = ZeroMatrix(3, 2);
                    int cnt = 1;
                    for (int i = 0; i < 3; ++i)
                        if (geom[i].FastGetSolutionStepValue(IS_BOUNDARY) == 0.0)
                        {
                            sort_coord(0, 0) = geom[i].X();
                            sort_coord(0, 1) = geom[i].Y();
                        }
                        else
                        {
                            sort_coord(cnt, 0) = geom[i].X();
                            sort_coord(cnt, 1) = geom[i].Y();
                            cnt++;
                        }

                    array_1d<double, 2 > vec1 = ZeroVector(2);
                    array_1d<double, 2 > vec2 = ZeroVector(2);

                    vec1[0] = sort_coord(0, 0) - sort_coord(1, 0);
                    vec1[1] = sort_coord(0, 1) - sort_coord(1, 1);

                    vec2[0] = sort_coord(2, 0) - sort_coord(1, 0);
                    vec2[1] = sort_coord(2, 1) - sort_coord(1, 1);

                    double outer_prod = 0.0;
                    outer_prod = vec2[1] * vec1[0] - vec1[1] * vec2[0];

                    double length_measure = 0.0;
                    length_measure = vec2[0] * vec2[0] + vec2[1] * vec2[1];
                    length_measure = sqrt(length_measure);
                    outer_prod /= length_measure;

                    //KRATOS_WATCH(fabs(outer_prod));
                    //	RATOS_WATCH(factor*length_measure);

                    if (fabs(outer_prod) < factor * length_measure)
                    {
                        for (int i = 0; i < 3; i++)
                        {
                            //if thats not the wall node, remove it
                            if (geom[i].FastGetSolutionStepValue(IS_BOUNDARY) == 0.0)
                            {
                                geom[i].Set(TO_ERASE, true);
                                //KRATOS_WATCH("NODE TOUCHING THE WALL - WILL BE ERASED!!!!")
                            }
                        }
                    }
                }

            }
        }
        else
        {
            for (ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i != ThisModelPart.ElementsEnd(); i++)
            {
                double n_str = 0;

                //counting number on nodes at the wall
                Geometry< Node < 3 > >& geom = i->GetGeometry();
                if (geom.size() == 4)
                {
                    for (int ii = 0; ii <= domain_size; ++ii)
                        n_str += geom[ii].FastGetSolutionStepValue(IS_STRUCTURE);

                    //if two walls are at the wall, we check if the third node is close to it or not by passing the alpha-shape
                    if (n_str == 3.0)
                    {
                        boost::numeric::ublas::bounded_matrix<double, 4, 3 > sort_coord = ZeroMatrix(4, 3);
                        int cnt = 1;
                        for (int i = 0; i < 4; ++i)
                        {
                            if (geom[i].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0)
                            {
                                sort_coord(0, 0) = geom[i].X();
                                sort_coord(0, 1) = geom[i].Y();
                                sort_coord(0, 2) = geom[i].Z();
                            }
                            else
                            {
                                sort_coord(cnt, 0) = geom[i].X();
                                sort_coord(cnt, 1) = geom[i].Y();
                                sort_coord(cnt, 2) = geom[i].Z();
                                cnt++;
                            }
                        }
                        array_1d<double, 3 > vec1 = ZeroVector(3);
                        array_1d<double, 3 > vec2 = ZeroVector(3);
                        array_1d<double, 3 > vec3 = ZeroVector(3);
                        array_1d<double, 3 > vec4 = ZeroVector(3);


                        vec1[0] = sort_coord(0, 0) - sort_coord(1, 0);
                        vec1[1] = sort_coord(0, 1) - sort_coord(1, 1);
                        vec1[2] = sort_coord(0, 2) - sort_coord(1, 2);

                        vec2[0] = sort_coord(2, 0) - sort_coord(1, 0);
                        vec2[1] = sort_coord(2, 1) - sort_coord(1, 1);
                        vec2[2] = sort_coord(2, 2) - sort_coord(1, 2);

                        vec3[0] = sort_coord(3, 0) - sort_coord(1, 0);
                        vec3[1] = sort_coord(3, 1) - sort_coord(1, 1);
                        vec3[2] = sort_coord(3, 2) - sort_coord(1, 2);

                        //this last vectoris only needed to have all the edges of the base
                        vec4[0] = sort_coord(1, 0) - sort_coord(3, 0);
                        vec4[1] = sort_coord(1, 1) - sort_coord(3, 1);
                        vec4[2] = sort_coord(1, 2) - sort_coord(3, 2);

                        //Control the hight of the thetraedral element
                        //Volume of the tethraedra
                        double vol = (vec2[0] * vec3[1] * vec1[2] - vec2[0] * vec3[2] * vec1[1] +
                                      vec2[1] * vec3[2] * vec1[0] - vec2[1] * vec3[0] * vec1[2] +
                                      vec2[2] * vec3[0] * vec1[1] - vec2[2] * vec3[1] * vec1[0])*0.1666666666667;
                        //Area of the basis
                        array_1d<double, 3 > outer_prod = ZeroVector(3);
                        outer_prod[0] = vec2[1] * vec3[2] - vec2[2] * vec3[1];
                        outer_prod[1] = vec2[2] * vec3[0] - vec2[0] * vec3[2];
                        outer_prod[2] = vec2[0] * vec3[1] - vec2[1] * vec3[0];
                        double area_base = norm_2(outer_prod);
                        area_base *= 0.5;
                        //height
                        if (area_base > 0.0000000000001)
                            vol /= area_base;
                        else
                            KRATOS_WATCH("Something strange: area of a wall triangular element --> zero");

                        double length_measure1 = norm_2(vec2);
                        double length_measure = norm_2(vec3);
                        double length_measure2 = norm_2(vec4);


                        if(length_measure1 < length_measure2)
                        {
                            if(length_measure1 < length_measure)
                                length_measure = length_measure1;
                        }
                        else
                        {
                            if(length_measure2 < length_measure)
                                length_measure = length_measure2;
                        }

                        if (fabs(vol) < factor * length_measure)
                        {
                            for (int i = 0; i < 4; i++)
                            {
                                //if thats not the wall node, remove it
                                if (geom[i].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0)
                                {
                                    geom[i].Set(TO_ERASE, true);
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

    //**********************************************************************************************
    //**********************************************************************************************

    void MarkNodesTouchingInterface(ModelPart& ThisModelPart, int domain_size, double factor)
    {
        KRATOS_TRY;

        if (domain_size == 2)
        {

            for (ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i != ThisModelPart.ElementsEnd(); i++)
            {
                double n_intr = 0;

                //counting number on nodes at the Interface
                Geometry< Node < 3 > >& geom = i->GetGeometry();
                n_intr = geom[0].FastGetSolutionStepValue(IS_INTERFACE);
                n_intr += geom[1].FastGetSolutionStepValue(IS_INTERFACE);
                n_intr += geom[2].FastGetSolutionStepValue(IS_INTERFACE);
                //if two interfaces are at the interface, we check if the third node is close to it or not by passing the alpha-shape
                if (n_intr == 2.0)
                {
                    boost::numeric::ublas::bounded_matrix<double, 3, 2 > sort_coord = ZeroMatrix(3, 2);
                    int cnt = 1;
                    for (int i = 0; i < 3; ++i)
                        if (geom[i].FastGetSolutionStepValue(IS_INTERFACE) == 0.0)
                        {
                            sort_coord(0, 0) = geom[i].X();
                            sort_coord(0, 1) = geom[i].Y();
                        }
                        else
                        {
                            sort_coord(cnt, 0) = geom[i].X();
                            sort_coord(cnt, 1) = geom[i].Y();
                            cnt++;
                        }

                    array_1d<double, 2 > vec1 = ZeroVector(2);
                    array_1d<double, 2 > vec2 = ZeroVector(2);

                    vec1[0] = sort_coord(0, 0) - sort_coord(1, 0);
                    vec1[1] = sort_coord(0, 1) - sort_coord(1, 1);

                    vec2[0] = sort_coord(2, 0) - sort_coord(1, 0);
                    vec2[1] = sort_coord(2, 1) - sort_coord(1, 1);

                    double outer_prod = 0.0;
                    outer_prod = vec2[1] * vec1[0] - vec1[1] * vec2[0];

                    double length_measure = 0.0;
                    length_measure = vec2[0] * vec2[0] + vec2[1] * vec2[1];
                    length_measure = sqrt(length_measure);
                    outer_prod /= length_measure;


                    double hnode2 = geom[0].FastGetSolutionStepValue(NODAL_H);
                    //hnode2*=hnode2;


                    //KRATOS_WATCH(fabs(outer_prod));
                    //	RATOS_WATCH(factor*length_measure);

                    if (fabs(outer_prod) < factor * hnode2)
                    {
                        for (int i = 0; i < 3; i++)
                        {
                            //if thats not the wall node, remove it
                            if (geom[i].FastGetSolutionStepValue(IS_BOUNDARY) == 0.0 && geom[i].FastGetSolutionStepValue(IS_INTERFACE) == 0.0)
                            {
                                geom[i].Set(TO_ERASE, true);
                                //KRATOS_WATCH("NODE TOUCHING THE INTERFACE - WILL BE ERASED!!!!")
                                //KRATOS_WATCH(outer_prod)
                            }
                        }
                    }
                }

            }
            //not to delete interface nodes
            /*for(ModelPart::NodeIterator ind = ThisModelPart.NodesBegin(); ind != ThisModelPart.NodesEnd(); ++ind)
                    {
                    if(ind->FastGetSolutionStepValue(IS_INTERFACE) ==1.0)
                            ind->Set(TO_ERASE, false);

                    }*/



        }
        else if (domain_size == 3)
        {
            int elem_size = ThisModelPart.Elements().size();
            ModelPart::ElementIterator ElemBegin = ThisModelPart.ElementsBegin();
            ModelPart::ElementIterator ele;

            #pragma omp parallel for firstprivate(elem_size,ElemBegin),private(ele)
            for (int kk = 0; kk < elem_size; kk++)

            {
                ele = ElemBegin + kk;


                int n_intr = 0;

                //counting number on nodes at the Interface
                Geometry< Node < 3 > >& geom = ele->GetGeometry();

                if (geom.size() == 4)
                {
                    for (int ii = 0; ii <= domain_size; ++ii)
                        n_intr += int(geom[ii].FastGetSolutionStepValue(IS_INTERFACE));

                    //if three interfaces are at the interface, we check if the third node is close to it or not by passing the alpha-shape
                    if (n_intr == domain_size)
                    {
                        boost::numeric::ublas::bounded_matrix<double, 4, 3 > sort_coord = ZeroMatrix(4, 3);
                        int cnt = 1;
                        int non_interface_id = -1;
                        for (int j = 0; j <= domain_size; ++j)
                        {
                            if (geom[j].FastGetSolutionStepValue(IS_INTERFACE) == 0.0)
                            {
                                sort_coord(0, 0) = geom[j].X();
                                sort_coord(0, 1) = geom[j].Y();
                                sort_coord(0, 2) = geom[j].Z();
                                non_interface_id = j;
                            }
                            else
                            {
                                sort_coord(cnt, 0) = geom[j].X();
                                sort_coord(cnt, 1) = geom[j].Y();
                                sort_coord(cnt, 2) = geom[j].Z();
                                cnt++;
                            }
                        }
                        array_1d<double, 3 > outplane = ZeroVector(3);
                        array_1d<double, 3 > inplane1 = ZeroVector(3);
                        array_1d<double, 3 > inplane2 = ZeroVector(3);

                        outplane[0] = sort_coord(0, 0) - sort_coord(1, 0);
                        outplane[1] = sort_coord(0, 1) - sort_coord(1, 1);
                        outplane[2] = sort_coord(0, 2) - sort_coord(1, 2);

                        inplane1[0] = sort_coord(2, 0) - sort_coord(1, 0);
                        inplane1[1] = sort_coord(2, 1) - sort_coord(1, 1);
                        inplane1[2] = sort_coord(2, 2) - sort_coord(1, 2);

                        inplane2[0] = sort_coord(3, 0) - sort_coord(1, 0);
                        inplane2[1] = sort_coord(3, 1) - sort_coord(1, 1);
                        inplane2[2] = sort_coord(3, 2) - sort_coord(1, 2);

                        //calculate normal to the surface
                        array_1d<double, 3 > outerprod = ZeroVector(3);
                        //boost::numeric::ublas::bounded_matrix<double, 1, 3 > outerprod = ZeroMatrix(1, 3);
                        outerprod[0] = inplane1[1] * inplane2[2] - inplane1[2] * inplane2[1];
                        outerprod[1] = inplane1[2] * inplane2[0] - inplane1[0] * inplane2[2];
                        outerprod[2] = inplane1[0] * inplane2[1] - inplane1[1] * inplane2[0];

                        //outerprod = outer_prod(inplane1,inplane2);


                        double len = inner_prod(outerprod, outerprod);
                        len /= sqrt(len);
                        // double len = outerprod(0,0)*outerprod(0,0) + outerprod(1,0)*outerprod(1,0) + outerprod(2,0)*outerprod(2,0);
                        outerprod /= len;

                        //calc distance to surface
                        double dist_to_surf = 0.0;
                        dist_to_surf = inner_prod(outerprod, outplane);
                        // dist_to_surf = outerprod(0,0)*outplane(0)+ outerprod(1,0)*outplane(1)+outerprod(2,0)*outplane(2);

                        //get H
                        double hnode2 = (ThisModelPart.GetProperties(0))[THICKNESS];

                        //				double hnode2 = geom[non_interface_id].FastGetSolutionStepValue(NODAL_H);
                        //                                  hnode2*=hnode2;

                        //mark node to delete
                        if (fabs(dist_to_surf) < factor * hnode2)
                        {
                            geom[non_interface_id].SetLock();
                            geom[non_interface_id].Set(TO_ERASE, true);

                            geom[non_interface_id].UnSetLock();


                            // KRATOS_WATCH("NODE TOUCHING THE SURFACE - WILL BE ERASED!!!!")
                            // KRATOS_WATCH(dist_to_surf)

                        }

                    }//interface elemenets
                }//non_shell
            }//all elements loop
        }//domain_size==3
        KRATOS_CATCH("")
    }
    //**********************************************************************************************
    //**********************************************************************************************

    void InterfaceDetecting(ModelPart& ThisModelPart, int domain_size, double factor)
    {
        KRATOS_TRY;

        for (ModelPart::NodeIterator ind = ThisModelPart.NodesBegin(); ind != ThisModelPart.NodesEnd(); ++ind)
            ind->FastGetSolutionStepValue(IS_INTERFACE) = 0.0;



        for (ModelPart::ElementsContainerType::iterator elem = ThisModelPart.ElementsBegin();
                elem != ThisModelPart.ElementsEnd(); elem++)
        {

            if (elem->GetValue(IS_WATER_ELEMENT) == 0)
            {

                WeakPointerVector< Element >& neighbor_els = elem->GetValue(NEIGHBOUR_ELEMENTS);
                Geometry< Node < 3 > >& geom = elem->GetGeometry();

                for (int ii = 0; ii < (domain_size + 1); ++ii)
                {

                    if (neighbor_els[ii].GetValue(IS_WATER_ELEMENT) == 1 && neighbor_els[ii].Id() != elem->Id())
                    {

                        if (ii == 0) // 1,2 and/or 3
                        {
                            if (geom[1].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0 && geom[2].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0)
                            {


                                geom[1].FastGetSolutionStepValue(IS_INTERFACE) = 1.0;
                                geom[2].FastGetSolutionStepValue(IS_INTERFACE) = 1.0;
                                if (domain_size == 3)
                                    geom[3].FastGetSolutionStepValue(IS_INTERFACE) = 1.0;

                            }
                        }
                        if (ii == 1) // 0,2 and/or 3
                        {
                            if (geom[0].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0 && geom[2].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0)
                            {


                                geom[0].FastGetSolutionStepValue(IS_INTERFACE) = 1.0;
                                geom[2].FastGetSolutionStepValue(IS_INTERFACE) = 1.0;
                                if (domain_size == 3)
                                    geom[3].FastGetSolutionStepValue(IS_INTERFACE) = 1.0;

                            }
                        }
                        if (ii == 2) // 0,1 and/or 3
                        {
                            if (geom[0].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0 && geom[1].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0)
                            {

                                geom[0].FastGetSolutionStepValue(IS_INTERFACE) = 1.0;
                                geom[1].FastGetSolutionStepValue(IS_INTERFACE) = 1.0;
                                if (domain_size == 3)
                                    geom[3].FastGetSolutionStepValue(IS_INTERFACE) = 1.0;

                            }
                        }
                        if (ii == 3) // 0,1 and/or 2
                        {
                            if (geom[0].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0 && geom[1].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0)
                            {

                                geom[0].FastGetSolutionStepValue(IS_INTERFACE) = 1.0;
                                geom[1].FastGetSolutionStepValue(IS_INTERFACE) = 1.0;
                                if (domain_size == 3)
                                    geom[2].FastGetSolutionStepValue(IS_INTERFACE) = 1.0;

                            }
                        }

                    }
                }

            }
        }
        /* for(ModelPart::NodeIterator ind = ThisModelPart.NodesBegin(); ind != ThisModelPart.NodesEnd(); ++ind)
                ind->FastGetSolutionStepValue(IS_INTERFACE) = 0.0;

         for(ModelPart::NodeIterator ind = ThisModelPart.NodesBegin(); ind != ThisModelPart.NodesEnd(); ++ind)
                {
                if(ind->FastGetSolutionStepValue(IS_WATER) == 0.0 && ind->FastGetSolutionStepValue(IS_STRUCTURE) !=1.0)
                        {
                        //to loop over just one kind of fluid (the one that  has boundary node  at the begining)
                          WeakPointerVector< Node<3> >& neighb = ind->GetValue(NEIGHBOUR_NODES);
                          int other_kind = 0;

                          for( WeakPointerVector< Node<3> >::iterator ngh_ind = neighb.begin(); ngh_ind!=neighb.end(); ngh_ind++)
                             {
                                if(ngh_ind->FastGetSolutionStepValue(IS_WATER) ==1.0 && ngh_ind->FastGetSolutionStepValue(IS_STRUCTURE) !=1.0)
                                                other_kind++;
                             }

                           if(other_kind)
                             ind->FastGetSolutionStepValue(IS_INTERFACE) = 1.0;
                        }

                }*/
        //KRATOS_WATCH("NEw interface nodes are detected");

        //free elements having three interface nodes
        /*for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                        i!=ThisModelPart.ElementsEnd(); i++)
                {
                        Geometry< Node<3> >& geom = i->GetGeometry();
                        double n_intr=0;
                        for(int ii=0;ii<3;++ii)
                                if(geom[ii].FastGetSolutionStepValue(IS_INTERFACE) == 1.0)
                                        n_intr++;

                        if(n_intr == 3.0)
                                for(int ii=0;ii<3;++ii)
                                        geom[ii].FastGetSolutionStepValue(IS_INTERFACE) = 0.0;

                }*/
        KRATOS_CATCH("")
    }

    //**********************************************************************************************
    //**********************************************************************************************

    void ChangeWallWaterFlag(ModelPart& ThisModelPart, int domain_size)
    {
        KRATOS_TRY;
        for (ModelPart::NodeIterator ind = ThisModelPart.NodesBegin(); ind != ThisModelPart.NodesEnd(); ++ind)
        {
            if (ind->FastGetSolutionStepValue(IS_STRUCTURE) == 1.0)
            {
                double current_flag = ind->FastGetSolutionStepValue(IS_WATER);
                double opposit_flag;
                if (current_flag == 1.0)
                    opposit_flag = 0.0;
                else
                    opposit_flag = 1.0;
                int same_flag = 0;

                WeakPointerVector< Node < 3 > >& neighb = ind->GetValue(NEIGHBOUR_NODES);

                for (WeakPointerVector< Node < 3 > >::iterator ngh_ind = neighb.begin(); ngh_ind != neighb.end(); ngh_ind++)
                {
                    if (ngh_ind->FastGetSolutionStepValue(IS_STRUCTURE) != 1.0 && ngh_ind->FastGetSolutionStepValue(IS_INTERFACE) != 1.0)
                    {
                        double ngh_flag = ngh_ind->FastGetSolutionStepValue(IS_WATER);
                        if (current_flag == ngh_flag)
                            same_flag++;
                        else
                            opposit_flag = ngh_flag;
                    }
                }
                if (!same_flag)
                    ind->FastGetSolutionStepValue(IS_WATER) = opposit_flag;
            }
        }
        KRATOS_CATCH("")
    }
    //**********************************************************************************************
    //**********************************************************************************************

    void ChangeInterfaceWaterFlag(ModelPart& ThisModelPart, int domain_size)
    {
        KRATOS_TRY;
        for (ModelPart::NodeIterator ind = ThisModelPart.NodesBegin(); ind != ThisModelPart.NodesEnd(); ++ind)
        {
            if (ind->FastGetSolutionStepValue(IS_INTERFACE) == 1.0)
            {
                double current_flag = ind->FastGetSolutionStepValue(IS_WATER);
                double opposit_flag;
                if (current_flag == 1.0)
                    opposit_flag = 0.0;
                else
                    opposit_flag = 1.0;
                int same_flag = 0;

                WeakPointerVector< Node < 3 > >& neighb = ind->GetValue(NEIGHBOUR_NODES);

                for (WeakPointerVector< Node < 3 > >::iterator ngh_ind = neighb.begin(); ngh_ind != neighb.end(); ngh_ind++)
                {
                    if (ngh_ind->FastGetSolutionStepValue(IS_STRUCTURE) != 1.0 && ngh_ind->FastGetSolutionStepValue(IS_INTERFACE) != 1.0)
                    {
                        double ngh_flag = ngh_ind->FastGetSolutionStepValue(IS_WATER);
                        if (current_flag == ngh_flag)
                            same_flag++;
                        else
                            opposit_flag = ngh_flag;
                    }
                }
                if (!same_flag)
                    ind->FastGetSolutionStepValue(IS_WATER) = opposit_flag;
            }
        }
        KRATOS_CATCH("")
    }
    //**********************************************************************************************
    //**********************************************************************************************

    void AssignNearBoundaryH(ModelPart& ThisModelPart, double factor)
    {
        KRATOS_TRY;
        KRATOS_WATCH("NNNNNNNNNNNNNNearBoundaryHHHHHHHHHHHHHHHH");
        for (ModelPart::NodeIterator ind = ThisModelPart.NodesBegin(); ind != ThisModelPart.NodesEnd(); ++ind)
        {
            const double current_flag = ind->FastGetSolutionStepValue(IS_WATER);

            double same_flag = 1.0;
            double opposit_flag = 0.0;

            WeakPointerVector< Node < 3 > >& neighb = ind->GetValue(NEIGHBOUR_NODES);
            for (WeakPointerVector< Node < 3 > >::iterator ngh_ind = neighb.begin(); ngh_ind != neighb.end(); ngh_ind++)
            {
                const double ngh_flag = ngh_ind->FastGetSolutionStepValue(IS_WATER);
                if (current_flag == ngh_flag)
                    same_flag++;
                else
                    opposit_flag++;
            }

            if (opposit_flag / same_flag > 1.0)
            {
                double current_H = ind->FastGetSolutionStepValue(NODAL_H);
                double new_H = current_H / factor;

                ind->FastGetSolutionStepValue(NODAL_H) = new_H;
                /*  for( WeakPointerVector< Node<3> >::iterator ngh_ind = neighb.begin(); ngh_ind!=neighb.end(); ngh_ind++)
                                 ngh_ind->FastGetSolutionStepValue(NODAL_H)=new_H;*/


            }
        }
        KRATOS_WATCH("NNNNNNNNNNNNNNearBoundaryHHHHHHHHHHHHHHHH~~~~~~~~~~~~~EXIT");
        KRATOS_CATCH("")
    }

    //**********************************************************************************************
    //**********************************************************************************************

    void MoveNodes(ModelPart& ThisModelPart)
    {
        KRATOS_TRY

        for (ModelPart::NodeIterator ind = ThisModelPart.NodesBegin(); ind != ThisModelPart.NodesEnd(); ++ind)
        {
            (ind)->X() = (ind)->X0() + ind->GetSolutionStepValue(DISPLACEMENT_X);
            (ind)->Y() = (ind)->Y0() + ind->GetSolutionStepValue(DISPLACEMENT_Y);
            (ind)->Z() = (ind)->Z0() + ind->GetSolutionStepValue(DISPLACEMENT_Z);
        }
        KRATOS_CATCH("")
    }

    //**********************************************************************************************
    //**********************************************************************************************

    void AssignMeshVelocity(ModelPart& ThisModelPart)
    {
        KRATOS_TRY

        for (ModelPart::NodeIterator ind = ThisModelPart.NodesBegin(); ind != ThisModelPart.NodesEnd(); ++ind)
        {

            (ind)->FastGetSolutionStepValue(MESH_VELOCITY_X) = (ind)->FastGetSolutionStepValue(VELOCITY_X);
            (ind)->FastGetSolutionStepValue(MESH_VELOCITY_Y) = (ind)->FastGetSolutionStepValue(VELOCITY_Y);
            (ind)->FastGetSolutionStepValue(MESH_VELOCITY_Z) = (ind)->FastGetSolutionStepValue(VELOCITY_Z);

        }
        KRATOS_CATCH("")
    }

    //**********************************************************************************************
    //**********************************************************************************************
    //ATTENTION:: returns a negative volume if inverted elements appear

    double CalculateVolume(ModelPart& ThisModelPart, int domain_size)
    {
        KRATOS_TRY;

        bool inverted_elements = false;
        //auxiliary vectors
        double Atot = 0.00;
        if (domain_size == 2)
        {
            for (ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i != ThisModelPart.ElementsEnd(); i++)
            {
                //calculating shape functions values
                Geometry< Node < 3 > >& geom = i->GetGeometry();
                double Ael = 0.0;
                if (geom.size() == 3)
                {
                    Ael = GeometryUtils::CalculateVolume2D(geom);
                    Atot += Ael;
                }
                if (Ael < 0) inverted_elements = true;
            }
        }
        else if (domain_size == 3)
        {
            for (ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i != ThisModelPart.ElementsEnd(); i++)
            {
                //calculating shape functions values
                Geometry< Node < 3 > >& geom = i->GetGeometry();

                double Ael = 0.0;
                if (geom.size() == 4)
                {
                    Ael = GeometryUtils::CalculateVolume3D(geom);
                    Atot += Ael;
                }
                if (Ael < 0) inverted_elements = true;
            }
        }

        //set to negative the volume if inverted elements appear
        if (inverted_elements == true)
            Atot = -Atot;

        //return the total area
        return Atot;
        KRATOS_CATCH("");
    }
    //**********************************************************************************************
    bool CheckInvertedElements(ModelPart& ThisModelPart, int domain_size)
    {
        KRATOS_TRY;

        bool inverted_elements = false;
        //auxiliary vectors
        double Atot = 0.00;
        if (domain_size == 2)
        {
            for (ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i != ThisModelPart.ElementsEnd(); i++)
            {
                //calculating shape functions values
                Geometry< Node < 3 > >& geom = i->GetGeometry();
                double Ael = 0.0;
                if (geom.size() == 3)
                {
                    Ael = GeometryUtils::CalculateVolume2D(geom);
                    Atot += Ael;
                }
                if (Ael < 0) inverted_elements = true;
            }
        }
        else if (domain_size == 3)
        {
            for (ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i != ThisModelPart.ElementsEnd(); i++)
            {
                //calculating shape functions values
                Geometry< Node < 3 > >& geom = i->GetGeometry();

                double Ael = 0.0;
                if (geom.size() == 4)
                {
                    Ael = GeometryUtils::CalculateVolume3D(geom);
                    Atot += Ael;
                }
                if (Ael < 0) inverted_elements = true;
            }
        }

        //return the bool inverted elements
        return inverted_elements;
        KRATOS_CATCH("");
    }

    //**********************************************************************************************

    void ColourAirWaterElement(ModelPart& ThisModelPart, int domain_size)
    {
        KRATOS_TRY;
        KRATOS_WATCH("INSSSSSSSSSIDE COLOUR");
        ModelPart::ElementsContainerType::iterator elem_bg = ThisModelPart.ElementsBegin();
        const unsigned int n_elems = ThisModelPart.Elements().size();

        //#pragma omp parallel for firstprivate(n_elems, elem_bg)
        for (unsigned int jj = 0; jj < n_elems; ++jj)
        {
            ModelPart::ElementsContainerType::iterator elem = elem_bg + jj;

            Geometry< Node < 3 > >& geom = elem->GetGeometry();
            int same_colour = 0;

            //count air flag
            for (int ii = 0; ii <= domain_size; ++ii)
                if (geom[ii].FastGetSolutionStepValue(IS_WATER) == 0.0)
                    same_colour++;

            if (same_colour == (domain_size + 1))
            {
                elem->GetValue(IS_WATER_ELEMENT) = 0.0;

            }
            else
            {
                elem->GetValue(IS_WATER_ELEMENT) = 1.0;


            }

        }
        //#pragma omp barrier
        //detecting interface nodes

        //PfemUtils::InterfaceDetecting(ThisModelPart, domain_size, 5.0);




        /*	//chack for apex elements of all interface
                for(ModelPart::ElementsContainerType::iterator elem = ThisModelPart.ElementsBegin();
                        elem!=ThisModelPart.ElementsEnd(); elem++)
                {
                    if(elem->GetValue(IS_WATER_ELEMENT) == 0.0)
                 {
                        Geometry< Node<3> >& geom = elem->GetGeometry();
                        int is_interface = 0;

                        for(int ii= 0; ii<= domain_size; ++ii)
                                if(geom[ii].FastGetSolutionStepValue(IS_INTERFACE) == 1.0)
                                        is_interface++;

                        if(is_interface == (domain_size + 1))
                          {
        WeakPointerVector< Element >& neighbor_els = elem->GetValue(NEIGHBOUR_ELEMENTS);

                           double ref_elem_flag = 0.0; //interface is always 0
                           double water_neghbs = 0;

                            for(int ii=0; ii<(domain_size+1); ++ii)
                                {

                                if(neighbor_els[ii].GetValue(IS_WATER_ELEMENT) != ref_elem_flag && neighbor_els[ii].Id() != elem->Id())
                                   water_neghbs++;
                                }

                            if(water_neghbs == 1.0) //so change the flag to WATER
                                {
                                        elem->GetValue(IS_WATER_ELEMENT) = 1.0;
                                KRATOS_WATCH("CHANGE    CHANGE   CHANGE   CHANGE   CHANGE   CHANGE");
                                KRATOS_WATCH(elem->Id());
                                }
                           }
                 }
                }*/

        KRATOS_CATCH("");
    }

    //**********************************************************************************************
    //**********************************************************************************************

    void ReduceTimeStep(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;

        double dt = ThisModelPart.GetProcessInfo()[DELTA_TIME];
        double time = ThisModelPart.GetProcessInfo()[TIME];

        ThisModelPart.GetProcessInfo()[TIME] = time - dt * 0.5;
        ThisModelPart.GetProcessInfo()[DELTA_TIME] = dt * 0.5;

        KRATOS_CATCH("");
    }
    //**********************************************************************************************
    //**********************************************************************************************

    double ExplicitDeltaT(double CFL, double min_dt, double max_dt, ModelPart& ThisModelPart)
    {
        KRATOS_TRY

        int NumThreads = OpenMPUtils::GetNumThreads();
        std::vector< double > Threads_dt(NumThreads, 10.0);

        ModelPart::ElementsContainerType::iterator elem_bg = ThisModelPart.ElementsBegin();
        int n_elems = ThisModelPart.Elements().size();

        #pragma omp parallel for firstprivate(n_elems, elem_bg)
        for (int ii = 0; ii < n_elems; ++ii)
        {
            //calculate min_dt
            ModelPart::ElementsContainerType::iterator elem = elem_bg + ii;

            double calc_dt = 1.0;
            elem->Calculate(DELTA_TIME, calc_dt, ThisModelPart.GetProcessInfo());

            int k = OpenMPUtils::ThisThread();
            if (calc_dt < Threads_dt[k])
                Threads_dt[k] = calc_dt;

        }
        #pragma omp barrier

        //KRATOS_WATCH(omp_get_thread_num());
        KRATOS_WATCH(NumThreads);

        double DT = max_dt;
        for (int kk = 0; kk < NumThreads; ++kk)
            if (Threads_dt[kk] < DT)
                DT = Threads_dt[kk];



        if (DT < min_dt) DT = min_dt;

        // double DT = 0.00000001;
        DT *= CFL;
        ThisModelPart.GetProcessInfo()[DELTA_TIME] = DT;
        KRATOS_WATCH("ExplicitDeltaT");
        KRATOS_WATCH(DT);
        return DT;


        /*

                            double& DeltaTime = ThisModelPart.GetProcessInfo()[DELTA_TIME];
                            double DT = 10.0;

                            ModelPart::ElementsContainerType::iterator elem_bg = ThisModelPart.ElementsBegin();
                            unsigned int n_elems = ThisModelPart.Elements().size();

                           #pragma omp parallel for firstprivate(n_elems, elem_bg)
                            for(unsigned int ii=0; ii<n_elems; ++ii)
                                {
                                  //calculate min_dt
                                  ModelPart::ElementsContainerType::iterator elem = elem_bg + ii;

                                  double calc_dt = 1.0;
                                  elem->Calculate(DELTA_TIME, calc_dt, ThisModelPart.GetProcessInfo());

                                  if(calc_dt < DT)
                                      {
                                        #pragma omp critical
                                        DT = CFL*calc_dt;

                                      }

                                }

        KRATOS_WATCH("@@@@@@@@@@@@@@@@@@@@@@@@@@  BEFORE @@@@@@@@@@");
        KRATOS_WATCH(DeltaTime);
                             if(DT> max_dt)
                                DT = max_dt;
                             else if(DT < min_dt)
                                DT = min_dt;
        KRATOS_WATCH("@@@@@@@@@@@@@@@@@@@@@@@@@@  DELTA T EXPLICIT@@@@@@@@@@");
                            DeltaTime = DT;
        KRATOS_WATCH(DeltaTime);


                         return DeltaTime;*/
        KRATOS_CATCH("");
    }


private:

    //aux vars
    static boost::numeric::ublas::bounded_matrix<double, 3, 3 > msJ; //local jacobian
    static boost::numeric::ublas::bounded_matrix<double, 3, 3 > msJinv; //inverse jacobian
    static array_1d<double, 3 > msc; //center pos
    static array_1d<double, 3 > ms_rhs; //center pos



};


} // namespace Kratos.

#endif // KRATOS_PFEM_UTILITIES_INCLUDED  defined 


