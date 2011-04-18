/*
==============================================================================
KratosTestApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2010
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
//   Date:                $Date: 2007-03-06 10:30:31 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_LAGRANGIAN_PARTICLES_UTILITIES_INCLUDED )
#define  KRATOS_LAGRANGIAN_PARTICLES_UTILITIES_INCLUDED

#define PRESSURE_ON_EULERIAN_MESH


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
#include "incompressible_fluid_application.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/timer.h"

namespace Kratos
{

    template< class T, std::size_t dim >
    class DistanceCalculator
    {
    public:

        double operator()(T const& p1, T const& p2)
        {
            double dist = 0.0;
            for (std::size_t i = 0; i < dim; i++)
            {
                double tmp = p1[i] - p2[i];
                dist += tmp*tmp;
            }
            return dist; //square distance because it is easier to work without the square root//
        }
    };

    template<std::size_t TDim> class LagrangianParticleUtils
    {
    public:
        KRATOS_CLASS_POINTER_DEFINITION(LagrangianParticleUtils<TDim>);

        //**********************************************************************************************
        //**********************************************************************************************

        void StreamlineMove(array_1d<double, 3 > & body_force, const double density, const double dt, const double subdivisions, ModelPart& rEulerianModelPart, ModelPart& rLagrangianModelPart)
        {
            KRATOS_TRY

            double density_inverse = 1.0/density;

                    //defintions for spatial search
                    typedef Node < 3 > PointType;
            typedef Node < 3 > ::Pointer PointTypePointer;
            typedef std::vector<PointType::Pointer> PointVector;
            typedef std::vector<PointType::Pointer>::iterator PointIterator;
            typedef std::vector<double> DistanceVector;
            typedef std::vector<double>::iterator DistanceIterator;

            //creating an auxiliary list for the new nodes
            PointVector list_of_nodes;

            //*************
            // Bucket types
            typedef Bucket< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
            // 			   typedef Bins< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > StaticBins;
            // 			   typedef BinsDynamic< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > DynamicBins;
            //*************
            // DynamicBins;

            typedef Tree< KDTreePartition<BucketType> > tree; //Kdtree;
            // 			   typedef Tree< OCTreePartition<BucketType> > tree; 		//Octree;
            // 			   typedef Tree< StaticBins > tree; 		     		//Binstree;
            // 			   typedef Tree< KDTreePartition<StaticBins> > tree; 		//KdtreeBins;
            // 			   typedef typename KdtreeBins::Partitions SubPartitions;
            // 			   typedef Tree< OCTreePartition<StaticBins> > tree; 		//OctreeBins;
            /*
                                       typedef Bins< TDim, PointType, stdPointVector> stdBins;
                                       typedef Tree< Bins<TDim,PointType,stdPointVector> > tree; 	//stdStaticBins;*/

            //starting calculating time of construction of the kdtree
            boost::timer kdtree_construction;

            for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin();
                    node_it != rLagrangianModelPart.NodesEnd(); ++node_it)
            {
                PointTypePointer pnode = *(node_it.base());

                pnode->GetValue(ERASE_FLAG) = true;
                node_it->GetValue(IS_VISITED) = 0;

                //putting the nodes of the destination_model part in an auxiliary list
                list_of_nodes.push_back(pnode);
                
                //reset the position to the position at the end of the step
                array_1d<double, 3 > & old_disp = (node_it)->FastGetSolutionStepValue(DISPLACEMENT,1);
                (node_it)->FastGetSolutionStepValue(DISPLACEMENT) = old_disp;

                (node_it)->X() = (node_it)->X0() + old_disp[0];
                (node_it)->Y() = (node_it)->Y0() + old_disp[1];
                (node_it)->Z() = (node_it)->Z0() + old_disp[2];

            }

            std::cout << "kdt constructin time " << kdtree_construction.elapsed() << std::endl;


            //work arrays
            Node < 3 > work_point(0, 0.0, 0.0, 0.0);
            unsigned int MaximumNumberOfResults = 10000;
            PointVector Results(MaximumNumberOfResults);
            DistanceVector ResultsDistances(MaximumNumberOfResults);
            array_1d<double, TDim + 1 > N; //Shape functions vector//
            array_1d<double, TDim + 1 > pressures; //Shape functions vector//
            boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim> DN_DX;

            array_1d<double, TDim> gradp;
            array_1d<double, 3 > acc_particle;
            array_1d<double, 3 > veulerian;

            //create a spatial database with the list of new nodes
            unsigned int bucket_size = 20;

            double small_dt = dt/subdivisions;

            for(unsigned int substep=0; substep<subdivisions; substep++)
            {
                //compute the tree with the position of the nodes
                tree nodes_tree(list_of_nodes.begin(), list_of_nodes.end(), bucket_size);

                //loop over all of the elements in the eulerian mesh to perform the interpolation
                for (ModelPart::ElementsContainerType::iterator el_it = rEulerianModelPart.ElementsBegin();
                        el_it != rEulerianModelPart.ElementsEnd(); el_it++)
                {
                    Geometry<Node < 3 > >&geom = el_it->GetGeometry();

                    //find the center and "radius" of the element
                    double xc, yc, zc, radius;
                    CalculateCenterAndSearchRadius(geom, xc, yc, zc, radius, N);
                    work_point.X() = xc;
                    work_point.Y() = yc;
                    work_point.Z() = zc;

                    //find all of the new nodes within the radius
                    int number_of_points_in_radius;

                    //look between the new nodes which of them is inside the radius of the circumscribed cyrcle
                    number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(),
                            ResultsDistances.begin(), MaximumNumberOfResults);

                    if (number_of_points_in_radius > 0)
                    {
                        //check if inside
                        for (PointIterator it_found = Results.begin(); it_found != Results.begin() + number_of_points_in_radius; it_found++)
                        {

                            bool is_inside = false;
                            is_inside = CalculatePosition(geom, (*it_found)->X(), (*it_found)->Y(), (*it_found)->Z(), N);

                            if (is_inside == true && (*it_found)->GetValue(IS_VISITED) == 0)
                            {
//                                KRATOS_WATCH("219")
                                (*it_found)->GetValue(IS_VISITED) = 1;

                                //move according to the streamline
                                noalias(veulerian) = N[0]*geom[0].FastGetSolutionStepValue(VELOCITY);
                                for(unsigned int k=1; k<geom.size(); k++)
                                    noalias(veulerian) += N[k]*geom[k].FastGetSolutionStepValue(VELOCITY);

                                array_1d<double, 3 > & disp = (*it_found)->FastGetSolutionStepValue(DISPLACEMENT);

                                noalias(disp) += small_dt*veulerian;

                                (*it_found)->GetValue(ERASE_FLAG) = false;

                                //compute particle velocity
                                noalias(acc_particle) = body_force - N[0]*geom[0].FastGetSolutionStepValue(PRESS_PROJ)*density_inverse;
                                for(unsigned int k=1; k<geom.size(); k++)
                                    noalias(acc_particle) -= N[k]*geom[k].FastGetSolutionStepValue(PRESS_PROJ)*density_inverse;

                                array_1d<double, 3 > & vel_particle = (*it_found)->FastGetSolutionStepValue(VELOCITY);

                                if(substep == 0)
                                    noalias(vel_particle) = veulerian;

                                noalias(vel_particle) += small_dt*acc_particle;
                            }
                        }

                    }
                }

                //position is to be updated only after all of the searches!
                for (ModelPart::NodesContainerType::iterator it = rLagrangianModelPart.NodesBegin();
                        it != rLagrangianModelPart.NodesEnd(); it++)
                {
                    noalias(it->Coordinates()) = it->GetInitialPosition();
                    noalias(it->Coordinates()) += it->FastGetSolutionStepValue(DISPLACEMENT);
                }
            }



            KRATOS_CATCH("")
        }

        //**********************************************************************************************
        //**********************************************************************************************
        //function to seed a list of new nodes

        void Reseed(ModelPart& rEulerianModelPart, ModelPart& rLagrangianModelPart)
        {
            KRATOS_TRY;

            unsigned int id = (rEulerianModelPart.Nodes().end() - 1)->Id() + 1;
            rLagrangianModelPart.Nodes().clear();

            for (ModelPart::NodesContainerType::iterator node_it = rEulerianModelPart.NodesBegin();
                    node_it != rEulerianModelPart.NodesEnd(); node_it++)
            {
                int node_id = id++;
                double x = node_it->X();
                double y = node_it->Y();
                double z = node_it->Z();
                Node < 3 > ::Pointer pnode = rLagrangianModelPart.CreateNewNode(node_id, x, y, z);

                pnode->FastGetSolutionStepValue(VELOCITY) = node_it->FastGetSolutionStepValue(VELOCITY);
            }

            boost::numeric::ublas::bounded_matrix<double, TDim + 2, TDim + 1 > pos;
            boost::numeric::ublas::bounded_matrix<double, TDim + 2, TDim + 1 > N;

            for (ModelPart::ElementsContainerType::iterator el_it = rEulerianModelPart.ElementsBegin();
                    el_it != rEulerianModelPart.ElementsEnd(); el_it++)
            {
                Geometry<Node < 3 > >& geom = el_it->GetGeometry();
                ComputeGaussPointPositions(geom, pos, N);

                for (unsigned int i = 0; i < TDim + 2; i++)
                {
                    int node_id = id++;
                    Node < 3 > ::Pointer pnode = rLagrangianModelPart.CreateNewNode(node_id, pos(i, 0), pos(i, 1), pos(i, 2));

                    array_1d<double, 3 > & vel = pnode->FastGetSolutionStepValue(VELOCITY);
                    noalias(vel) = ZeroVector(3);
                    for (unsigned int j = 0; j < TDim + 1; j++)
                        noalias(vel) += N(i, j) * geom[j].FastGetSolutionStepValue(VELOCITY);
                }
            }

            for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin();
                    node_it != rLagrangianModelPart.NodesEnd(); node_it++)
            {
                node_it->FastGetSolutionStepValue(VELOCITY, 1) = node_it->FastGetSolutionStepValue(VELOCITY);
            }

            KRATOS_CATCH("");
        }

        //**********************************************************************************************
        //**********************************************************************************************

        void VisualizationModelPart(ModelPart& rCompleteModelPart, ModelPart& rEulerianModelPart, ModelPart& rLagrangianModelPart)
        {
            KRATOS_TRY;

            rCompleteModelPart.Elements() = rEulerianModelPart.Elements();
            rCompleteModelPart.Nodes() = rEulerianModelPart.Nodes();

            for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin();
                    node_it != rLagrangianModelPart.NodesEnd(); node_it++)
            {
                rCompleteModelPart.AddNode(*(node_it.base()));
            }

            KRATOS_CATCH("");
        }

        //**********************************************************************************************
        //**********************************************************************************************

        void TransferToEulerianMesh(ModelPart& rEulerianModelPart, ModelPart& rLagrangianModelPart)
        {
            KRATOS_TRY

                    //defintions for spatial search
                    typedef Node < 3 > PointType;
            typedef Node < 3 > ::Pointer PointTypePointer;
            typedef std::vector<PointType::Pointer> PointVector;
            typedef std::vector<PointType::Pointer>::iterator PointIterator;
            typedef std::vector<double> DistanceVector;
            typedef std::vector<double>::iterator DistanceIterator;

            //creating an auxiliary list for the new nodes
            PointVector list_of_nodes;

            //*************
            // Bucket types
            typedef Bucket< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
            // 			   typedef Bins< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > StaticBins;
            // 			   typedef BinsDynamic< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > DynamicBins;
            //*************
            // DynamicBins;

            typedef Tree< KDTreePartition<BucketType> > tree; //Kdtree;
            // 			   typedef Tree< OCTreePartition<BucketType> > tree; 		//Octree;
            // 			   typedef Tree< StaticBins > tree; 		     		//Binstree;
            // 			   typedef Tree< KDTreePartition<StaticBins> > tree; 		//KdtreeBins;
            // 			   typedef typename KdtreeBins::Partitions SubPartitions;
            // 			   typedef Tree< OCTreePartition<StaticBins> > tree; 		//OctreeBins;
            /*
                                       typedef Bins< TDim, PointType, stdPointVector> stdBins;
                                       typedef Tree< Bins<TDim,PointType,stdPointVector> > tree; 	//stdStaticBins;*/

            //starting calculating time of construction of the kdtree
            boost::timer kdtree_construction;

            for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin();
                    node_it != rLagrangianModelPart.NodesEnd(); ++node_it)
            {
                PointTypePointer pnode = *(node_it.base());

                //putting the nodes of the destination_model part in an auxiliary list
                list_of_nodes.push_back(pnode);
            }

            std::cout << "kdt constructin time " << kdtree_construction.elapsed() << std::endl;

            //create a spatial database with the list of new nodes
            unsigned int bucket_size = 20;
            tree nodes_tree(list_of_nodes.begin(), list_of_nodes.end(), bucket_size);

            //work arrays
            Node < 3 > work_point(0, 0.0, 0.0, 0.0);
            unsigned int MaximumNumberOfResults = 10000;
            PointVector Results(MaximumNumberOfResults);
            DistanceVector ResultsDistances(MaximumNumberOfResults);


            if(rEulerianModelPart.NodesBegin()->SolutionStepsDataHas(NODAL_H) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----NODAL_H---- variable!!!!!! ERROR", "");

            for (ModelPart::NodesContainerType::iterator node_it = rEulerianModelPart.NodesBegin();
                    node_it != rEulerianModelPart.NodesEnd(); node_it++)
            {
                work_point.X() = node_it->X();
                work_point.Y() = node_it->Y();
                work_point.Z() = node_it->Z();

                double radius = 0.6 * node_it->FastGetSolutionStepValue(NODAL_H);

                //find all of the new nodes within the radius
                int number_of_points_in_radius;

                //look between the new nodes which of them is inside the radius of the circumscribed cyrcle
                number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(),
                        ResultsDistances.begin(), MaximumNumberOfResults);


                if (number_of_points_in_radius > 0)
                {

                    array_1d<double, 3 > & vel = (node_it)->FastGetSolutionStepValue(VELOCITY);

                    array_1d<double, 3 > original_vel = vel;

                    noalias(vel) = ZeroVector(3);

                    double tot_weight = 0.0;

                    for (int k = 0; k < number_of_points_in_radius; k++)
                    {
                        double weight = 1.0;
                        //                        double weight = 1.0 / (sqrt(ResultsDistances[k]) + 1e-9);
                        tot_weight += weight;

                        PointIterator it_found = Results.begin() + k;

                        const array_1d<double,3> particle_velocity = (*it_found)->FastGetSolutionStepValue(VELOCITY);

                        noalias(vel) += weight * particle_velocity;
                    }

                    vel /= tot_weight;

                    if (node_it->IsFixed(VELOCITY_X))
                        noalias(vel) = original_vel;

                } else
                {
                    if (node_it->IsFixed(VELOCITY_X))
                        node_it->FastGetSolutionStepValue(VELOCITY) = ZeroVector(3);
                }


            }

            KRATOS_CATCH("")
        }

        

        //restarting the step from the beginning

        void RestartStep(ModelPart& rModelPart)
        {
            KRATOS_TRY;

            //setting the variables to their value at the beginning of the time step
            rModelPart.OverwriteSolutionStepData(1, 0);

            //setting the coordinates to their value at the beginning of the step
            for (ModelPart::NodesContainerType::iterator node_it = rModelPart.NodesBegin();
                    node_it != rModelPart.NodesEnd(); node_it++)
            {
                array_1d<double, 3 > & coords = node_it->Coordinates();
                const array_1d<double, 3 > & old_disp = node_it->FastGetSolutionStepValue(DISPLACEMENT, 1);

                coords[0] = node_it->X0() + old_disp[0];
                coords[1] = node_it->Y0() + old_disp[1];
                coords[2] = node_it->Z0() + old_disp[2];
            }


            KRATOS_CATCH("");
        }


    private:

        inline void CalculateCenterAndSearchRadius(Geometry<Node < 3 > >&geom,
                double& xc, double& yc, double& zc, double& R, array_1d<double, 3 > & N
                )
        {
            double x0 = geom[0].X();
            double y0 = geom[0].Y();
            double x1 = geom[1].X();
            double y1 = geom[1].Y();
            double x2 = geom[2].X();
            double y2 = geom[2].Y();


            xc = 0.3333333333333333333 * (x0 + x1 + x2);
            yc = 0.3333333333333333333 * (y0 + y1 + y2);
            zc = 0.0;

            double R1 = (xc - x0)*(xc - x0) + (yc - y0)*(yc - y0);
            double R2 = (xc - x1)*(xc - x1) + (yc - y1)*(yc - y1);
            double R3 = (xc - x2)*(xc - x2) + (yc - y2)*(yc - y2);

            R = R1;
            if (R2 > R) R = R2;
            if (R3 > R) R = R3;

            R = 1.01 * sqrt(R);
        }
        //***************************************
        //***************************************

        inline void CalculateCenterAndSearchRadius(Geometry<Node < 3 > >&geom,
                double& xc, double& yc, double& zc, double& R, array_1d<double, 4 > & N

                )
        {
            double x0 = geom[0].X();
            double y0 = geom[0].Y();
            double z0 = geom[0].Z();
            double x1 = geom[1].X();
            double y1 = geom[1].Y();
            double z1 = geom[1].Z();
            double x2 = geom[2].X();
            double y2 = geom[2].Y();
            double z2 = geom[2].Z();
            double x3 = geom[3].X();
            double y3 = geom[3].Y();
            double z3 = geom[3].Z();


            xc = 0.25 * (x0 + x1 + x2 + x3);
            yc = 0.25 * (y0 + y1 + y2 + y3);
            zc = 0.25 * (z0 + z1 + z2 + z3);

            double R1 = (xc - x0)*(xc - x0) + (yc - y0)*(yc - y0) + (zc - z0)*(zc - z0);
            double R2 = (xc - x1)*(xc - x1) + (yc - y1)*(yc - y1) + (zc - z1)*(zc - z1);
            double R3 = (xc - x2)*(xc - x2) + (yc - y2)*(yc - y2) + (zc - z2)*(zc - z2);
            double R4 = (xc - x3)*(xc - x3) + (yc - y3)*(yc - y3) + (zc - z3)*(zc - z3);

            R = R1;
            if (R2 > R) R = R2;
            if (R3 > R) R = R3;
            if (R4 > R) R = R4;

            R = sqrt(R);
        }

        //***************************************
        //***************************************

        inline bool CalculatePosition(Geometry<Node < 3 > >&geom,
                const double xc, const double yc, const double zc,
                array_1d<double, 3 > & N
                )
        {
            double x0 = geom[0].X();
            double y0 = geom[0].Y();
            double x1 = geom[1].X();
            double y1 = geom[1].Y();
            double x2 = geom[2].X();
            double y2 = geom[2].Y();

            double area = CalculateVol(x0, y0, x1, y1, x2, y2);
            double inv_area = 0.0;
            if (area == 0.0)
            {
                KRATOS_ERROR(std::logic_error, "element with zero area found", "");
            } else
            {
                inv_area = 1.0 / area;
            }


            N[0] = CalculateVol(x1, y1, x2, y2, xc, yc) * inv_area;
            N[1] = CalculateVol(x2, y2, x0, y0, xc, yc) * inv_area;
            N[2] = CalculateVol(x0, y0, x1, y1, xc, yc) * inv_area;


            if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
                return true;

            return false;
        }

        //***************************************
        //***************************************

        inline bool CalculatePosition(Geometry<Node < 3 > >&geom,
                const double xc, const double yc, const double zc,
                array_1d<double, 4 > & N
                )
        {

            double x0 = geom[0].X();
            double y0 = geom[0].Y();
            double z0 = geom[0].Z();
            double x1 = geom[1].X();
            double y1 = geom[1].Y();
            double z1 = geom[1].Z();
            double x2 = geom[2].X();
            double y2 = geom[2].Y();
            double z2 = geom[2].Z();
            double x3 = geom[3].X();
            double y3 = geom[3].Y();
            double z3 = geom[3].Z();

            double vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

            double inv_vol = 0.0;
            if (vol < 0.0000000000001)
            {
                KRATOS_ERROR(std::logic_error, "element with zero vol found", "");
            } else
            {
                inv_vol = 1.0 / vol;
            }

            N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc) * inv_vol;
            N[1] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc) * inv_vol;
            N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc) * inv_vol;
            N[3] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc) * inv_vol;


            if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[3] >= 0.0 &&
                    N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0 && N[3] <= 1.0)
                //if the xc yc zc is inside the tetrahedron return true
                return true;

            return false;
        }

        inline double CalculateVol(const double x0, const double y0,
                const double x1, const double y1,
                const double x2, const double y2
                )
        {
            return 0.5 * ((x1 - x0)*(y2 - y0)- (y1 - y0)*(x2 - x0));
        }
        //***************************************
        //***************************************

        inline double CalculateVol(const double x0, const double y0, const double z0,
                const double x1, const double y1, const double z1,
                const double x2, const double y2, const double z2,
                const double x3, const double y3, const double z3
                )
        {
            double x10 = x1 - x0;
            double y10 = y1 - y0;
            double z10 = z1 - z0;

            double x20 = x2 - x0;
            double y20 = y2 - y0;
            double z20 = z2 - z0;

            double x30 = x3 - x0;
            double y30 = y3 - y0;
            double z30 = z3 - z0;

            double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
            return detJ * 0.1666666666666666666667;
        }

        void ComputeGaussPointPositions(Geometry< Node < 3 > >& geom, boost::numeric::ublas::bounded_matrix<double, 4, 3 > & pos, boost::numeric::ublas::bounded_matrix<double, 4, 3 > & N)
        {
            double one_third = 1.0 / 3.0;
            double one_sixt = 1.0 / 6.0;
            double two_third = 2.0 * one_third;

            N(0, 0) = one_sixt;
            N(0, 1) = one_sixt;
            N(0, 2) = two_third;
            N(1, 0) = two_third;
            N(1, 1) = one_sixt;
            N(1, 2) = one_sixt;
            N(2, 0) = one_sixt;
            N(2, 1) = two_third;
            N(2, 2) = one_sixt;
            N(3, 0) = one_third;
            N(3, 1) = one_third;
            N(3, 2) = one_third;


            //first
            pos(0, 0) = one_sixt * geom[0].X() + one_sixt * geom[1].X() + two_third * geom[2].X();
            pos(0, 1) = one_sixt * geom[0].Y() + one_sixt * geom[1].Y() + two_third * geom[2].Y();
            pos(0, 2) = one_sixt * geom[0].Z() + one_sixt * geom[1].Z() + two_third * geom[2].Z();

            //second
            pos(1, 0) = two_third * geom[0].X() + one_sixt * geom[1].X() + one_sixt * geom[2].X();
            pos(1, 1) = two_third * geom[0].Y() + one_sixt * geom[1].Y() + one_sixt * geom[2].Y();
            pos(1, 2) = two_third * geom[0].Z() + one_sixt * geom[1].Z() + one_sixt * geom[2].Z();

            //third
            pos(2, 0) = one_sixt * geom[0].X() + two_third * geom[1].X() + one_sixt * geom[2].X();
            pos(2, 1) = one_sixt * geom[0].Y() + two_third * geom[1].Y() + one_sixt * geom[2].Y();
            pos(2, 2) = one_sixt * geom[0].Z() + two_third * geom[1].Z() + one_sixt * geom[2].Z();

            //fourth
            pos(3, 0) = one_third * geom[0].X() + one_third * geom[1].X() + one_third * geom[2].X();
            pos(3, 1) = one_third * geom[0].Y() + one_third * geom[1].Y() + one_third * geom[2].Y();
            pos(3, 2) = one_third * geom[0].Z() + one_third * geom[1].Z() + one_third * geom[2].Z();

        }

        void ConsistentMassMatrix(const double A, boost::numeric::ublas::bounded_matrix<double, 3, 3 > & M)
        {
            double c1 = A / 12.0;
            double c2 = 2.0 * c1;
            M(0, 0) = c2;
            M(0, 1) = c1;
            M(0, 2) = c1;
            M(1, 0) = c1;
            M(1, 1) = c2;
            M(1, 2) = c1;
            M(2, 0) = c1;
            M(2, 1) = c1;
            M(2, 2) = c2;
        }

    };

} // namespace Kratos.

#endif // KRATOS_LAGRANGIAN_PARTICLES_UTILITIES_INCLUDED  defined


