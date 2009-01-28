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
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2009-01-14 17:15:24 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/

#if !defined(KRATOS_NODE_SNAPPING_UTILITY_INCLUDED )
#define  KRATOS_NODE_SNAPPING_UTILITY_INCLUDED
//System includes
//External includes
#include "boost/smart_ptr.hpp"

//Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "includes/element.h"
#include "integration/integration_point.h"
#include "spaces/ublas_space.h"

#include "processes/find_nodal_neighbours_process.h"
#include "custom_utilities/circle_3d.h"
#include "custom_utilities/cylinder_3d.h"
#include "custom_utilities/closed_cylinder_3d.h"
#include "custom_utilities/sd_math_utils.h"

#include "utilities/math_utils.h"
#include "containers/vector_map.h"

#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "linear_solvers/cg_solver.h"

#include "custom_utilities/variable_transfer_utility.h"

namespace Kratos
{
    class NodeSnappingUtility
    {
        public:
            typedef Dof<double> TDofType;
            typedef PointerVectorSet<TDofType, IndexedObject> DofsArrayType;
            typedef ModelPart::ElementsContainerType ElementsArrayType;
            typedef double* ContainerType;
            typedef Element::DofsVectorType DofsVectorType;
            typedef Geometry<Node<3> >::IntegrationPointsArrayType IntegrationPointsArrayType;
            typedef Geometry<Node<3> >::GeometryType GeometryType;
            typedef Geometry<Node<3> >::CoordinatesArrayType CoordinatesArrayType;
            typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
            typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
            typedef CGSolver<SpaceType,  LocalSpaceType> CGSolverType;
            typedef Geometry<Node<3> >::GeometriesArrayType GeometriesArrayType;

            /** 
             * Constructor.
             */
            NodeSnappingUtility()
            {
                std::cout << "NodeSnappingUtility created" << std::endl;
            }

            /** 
             * Destructor.
             */
            virtual ~NodeSnappingUtility()
            {}

            /**
			 *
             */
            void AdjustToCylinder( ModelPart& model_part, Point<3> center, Point<3> e1, 
                                         Point<3> e2, Point<3> e3)
            {

            }//AdjustToCylinder

            /** 
             * Adjusts the elements to a closed cylinder and deaktivates the elements inside the cylinder
             * @param model_part current model part
             * @param center starting point of the cylinder, must be midpoint of the circle
             * @param e1 point on the circle above the midpoint
             * @param e2 point on the circle above the midpoint and orthogonal to e2
             * @param e3 point on the direction vector of the cylinder
             */
            void AdjustToClosedCylinder( ModelPart& model_part, Point<3> center, Point<3> e1, Point<3> e2, Point<3> e3){
                std::cout<<"####################### AdjustToClosedCylinder - START"<<std::endl;

//                 std::cout << std::endl;
//                 std::cout << "Controll of Gausspoint transfer at the beginning" << std::endl;
//                 ControllGaussTransfer(model_part);

                std::cout<<"Create cylinder"<<std::endl;
                KRATOS_TRY
                //generating 3D cylinder
                ClosedCylinder3D cylinder( center, e1, e2, e3 );

                std::cout << std::endl;
                std::cout << "Identify elements in the vicinity of the cylinder" << std::endl;
                std::vector<int> vicinityelements;
                std::map<int, int> vicinitynodes;//nodeId, position in map
                IdentifyVicinity(model_part, vicinityelements, vicinitynodes, cylinder);

                std::cout << std::endl;
                std::cout << "Transfer variables from Gausspoints to nodes " << std::endl;
                TransferVariablesToNodes(model_part, vicinityelements, vicinitynodes, INSITU_STRESS);

//                 std::cout << std::endl;
//                 std::cout << "Controll of TransferVariablesToNodes()" << std::endl;
//                 ControllTransferVariablesToNodes(model_part, vicinitynodes);

                std::cout << std::endl;
                std::cout << "calculation of the new position of all contact_slave-nodes" << std::endl;
                MoveCylinderWallNodes(model_part, cylinder);

//                 //deactivation of all elements inside the cylinder
//                 IdentifyInsideElements(model_part, vicinityelements, cylinder);
                
                DefineCapNodes(model_part, vicinityelements, cylinder);

                std::cout << std::endl;
                std::cout << "calculation of the new position of all nodes at the front" << std::endl;
                MoveCylinderCapNodes(model_part, cylinder);

                //deactivation of all elements inside the cylinder
                IdentifyInsideElements(model_part, vicinityelements, cylinder);

                std::cout << std::endl;
                std::cout << "Transfer variables from nodes to Gausspoints" << std::endl;
                TransferVariablesToGaussPoints(model_part, vicinityelements, INSITU_STRESS);

//                 std::cout << std::endl;
//                 std::cout << "Controll of Gausspoint transfer at the end" << std::endl;
//                 ControllGaussTransfer(model_part);
                
                cap_nodes.clear();

                KRATOS_CATCH("")
                std::cout<<"####################### AdjustToClosedCylinder - END"<<std::endl;
            }//AdjustToClosedCylinder

            /** 
             * Indentifies all elements whose center is inside a bigger cylinder
             * @param model_part current model_part
             * @param vicinityelements Ids of inside elements
             * @param vicinitynodes Ids of the nodes belonging to the elements
             * @param cylinder closed cylinder
             */
            void IdentifyVicinity(ModelPart& model_part, std::vector<int>& vicinityelements, std::map<int, int>& vicinitynodes, ClosedCylinder3D& cylinder){
                for ( ElementsArrayType::ptr_iterator it=model_part.Elements().ptr_begin(); it!=model_part.Elements().ptr_end(); ++it){
                    Point<3> center = (*it)->GetGeometry().Center();
                    if (cylinder.IsInBigger(center)) vicinityelements.push_back((*it)->Id());
                }

                //vector with all nodes of the vicinityelements
                for(std::vector<int>::iterator it=vicinityelements.begin(); it != vicinityelements.end(); ++it){
                    for (unsigned int i=0; i<model_part.GetElement(*it).GetGeometry().size(); ++i){
                           vicinitynodes.insert( std::pair<int, int>(model_part.GetElement(*it).GetGeometry().GetPoint(i).Id(), vicinitynodes.size()));
                    }
                }
                std::cout << "size of vicinityelements " << vicinityelements.size() << " size of vicinitynodes " << vicinitynodes.size() << std::endl;
            }

            /** 
             * Returns the intersection of a line and the cap of a closed cylinder
             * @param k0 point inside the cylinder
             * @param k1 point outside the cylinder
             * @param cylinder closed cylinder
             * @return intersection
             */
            Point<3> Intersection(Point<3> k0, Point<3> k1, ClosedCylinder3D cylinder){
                Point<3> b(0,0,0);
                Point<3> m = cylinder.GetCenter();
                Point<3> e1 = cylinder.GetR1();
                Point<3> e2 = cylinder.GetR2();
                Point<3> e3 = cylinder.GetR3();
                Matrix A = ZeroMatrix(3,3);
                Matrix invA = ZeroMatrix(3,3);
                double detA = 0;
                Point<3> solution(0,0,0);
                b = k0 - m - e3;
                for (int i=0; i<3; ++i){
                    A(0,i) = e1(i);
                    A(1,i) = e2(i);
                    A(2,i) = -(k1(i)-k0(i));
                }
                MathUtils<double>::InvertMatrix3(A, invA, detA);
                if( fabs(detA) > 1.0e-6 ) solution = InvAb(invA, b);

                return k0 + solution(2) * (k1 - k0);
            }

            /** 
             * Multiplies the invers matrix AA with the vector b
             * @param invAA invers matrix of AA
			 * @param b vector
			 * @return x
             */            
            Point<3> InvAb(Matrix &invA, Point<3> &b){
                Point<3> x(0,0,0);

                for( unsigned int i=0; i<invA.size1(); ++i ){
                    for( unsigned int k=0; k<invA.size2(); ++k ) {
                        x(k) += invA(i,k)*b(i);
                    }
                }

                return x;
            }
            
            void AdjustToCircle( ModelPart& model_part, Point<3> center, Point<3> e1, Point<3> e2 )
            {
             }
           
            void AdjustNodes( ModelPart& model_part, Point<3>& newPosition )
            {
			}
           
            /**
             * Identifies elements that are located inside a cylinder
			 * @param model_part the current model part
             * @param vicinityelements all elements in the vicinity of the cylinder
             * @param cylinder closed cylinder
             */
            void IdentifyInsideElements( ModelPart& model_part,  std::vector<int>& vicinityelements, ClosedCylinder3D& cylinder )
            {
                std::cout << "testing for elements inside the cylinder..." << std::endl;
                //loop over all elements
                for( std::vector<int>::iterator it=vicinityelements.begin(); it != vicinityelements.end(); ++it) {
                    Point<3> center = model_part.GetElement(*it).GetGeometry().Center();
                    if( cylinder.IsInorOn( center ) ){
//                         bool allnodesinside = true;
                        //loop over all nodes
//                         for (unsigned int i=0; i<model_part.GetElement(*it).GetGeometry().size(); ++i)
//                             if (!cylinder.IsInorOn( model_part.GetElement(*it).GetGeometry()[i] ) ) 
//                                 allnodesinside = false;
//                         if (allnodesinside) model_part.GetElement(*it).GetValue(ACTIVATION_LEVEL) = 1;
                        //TODO: ATTENTION!!! switched off deactivation!!!
//                         model_part.GetElement(*it).GetValue(ACTIVATION_LEVEL) = 1;
                    }
                }//end of ELEMENTS
                std::cout << "done" << std::endl;
			}

            /**
             * Searches the element, the given point lies within and its respective local coordinates
			 * @param elements_set the set of elements in which the point sould be inside
			 * @param point checked point
			 * @param elem_id contains the element ID of the element the point lies inside
			 * @param rResult local coordinates of the the given point in the found element
			 * @retval false if the given point could not be found in the given element set
			 * @retval true if an element could be found which contains the point
             */
            bool FindElement( WeakPointerVector<Element>& elements_set, CoordinatesArrayType& point, 
                              int& elem_id, Point<3>&  rResult ) {
                for( WeakPointerVector<Element>::iterator it = elements_set.begin(); it != elements_set.end(); ++it ) {
                    if( it->GetGeometry().IsInside( point, rResult ) ) {
                        elem_id = it->Id();
                        return true;
                    }
                }
                return false;
            } 

			/**
             * moves all contact_slaves nodes in the vicinty of a closed cylinder to its wall
             * @param model_part the current model part
             * @param cylinder closed cylinder to which the contact_slave nodes should be adapted
             */
            void MoveCylinderWallNodes(ModelPart& model_part, ClosedCylinder3D& cylinder){
                //calculation of the new position of all contact_slave-nodes
                //loop over surfaces
                for (ModelPart::ConditionsContainerType::ptr_iterator it = model_part.Conditions().ptr_begin(); it != model_part.Conditions().ptr_end(); ++it){
                    if( (*it)->GetValue( IS_CONTACT_SLAVE ) ){
                    //the current surface is part of the contact_slave which contains the nodes to be moved on the cylinder wall
                        //loop over all nodes of the current surface
                        for (unsigned int inode=0; inode < (*it)->GetGeometry().size(); ++inode){
                            Point<3> old_position = (*it)->GetGeometry()[inode].GetInitialPosition();
                            Point<3> new_position = cylinder.ClosestPointOnWall(old_position);
                            Point<3> distance = new_position - old_position;
                            double dist = pow(distance[0],2.0)+pow(distance[1],2.0)+pow(distance[2],2.0);
                            if ( dist > 0.000001){
                            //the new node position is different from the old one
                                //neighbour elements of the node
                                WeakPointerVector<Element>& rneigh_el = (*it)->GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS);
                                
                                //move the node
                                if (!MoveNode( model_part, rneigh_el, (*it)->GetGeometry()[inode], new_position)) std::cout << "node " << (*it)->GetGeometry()[inode].Id() << " could not be moved from " << old_position << " to " << new_position << std::endl;
                                else{
                                    //move inner nodes
                                    (*it)->GetGeometry()[inode].FastGetSolutionStepValue( IS_VISITED ) = 3;
                                    MoveInnerNodes(model_part, rneigh_el);
                                }

                            }
                        }
                    }
                }
            }

			/**
             * moves all nodes in the vicinty of the cap of a closed cylinder to its cap
             * @param model_part the current model part
             * @param cylinder closed cylinder to which the nodes should be adapted
             */
            void MoveCylinderCapNodes(ModelPart& model_part, ClosedCylinder3D& cylinder){
                //calculation of the new position of all nodes at the front
                //loop over all nodes
                for (std::vector<int>::iterator it=cap_nodes.begin(); it != cap_nodes.end(); ++it){
                    Point<3> old_position = model_part.GetNode(*it).GetInitialPosition();
                    Point<3> new_position = cylinder.ClosestPointOnCap(old_position);
                    Point<3> distance = new_position - old_position;
                    double dist = pow(distance[0],2.0)+pow(distance[1],2.0)+pow(distance[2],2.0);
                    if ( dist > 0.000001){
                            //the new node position is different from the old one
                        Point<3> step = (new_position - old_position) / 10;

                        //neighbour elements of the node
                        WeakPointerVector<Element>& rneigh_el = model_part.GetNode(*it).GetValue(NEIGHBOUR_ELEMENTS);

                        for (int i=0; i<10; ++i){
                            Point<3> spstep = old_position + step*(i+1);
                            //move the node
                            if (!MoveNode( model_part, rneigh_el, model_part.GetNode(*it), spstep)) std::cout << "node " << model_part.GetNode(*it).Id() << " could not be moved from " << model_part.GetNode(*it).GetInitialPosition() << " to " << spstep << std::endl;
                            else{
                                //move inner nodes
                                model_part.GetNode(*it).FastGetSolutionStepValue( IS_VISITED ) = 3;
                                MoveInnerNodes(model_part, rneigh_el);
                            }
                        }
                    }
                }
            }

			/**
             * should eliminate edges running through the cylinder
             * @param model_part the current model part
             * @param cylinder closed cylinder to which the nodes should be adapted
			 * 
			 * TODO: doesn't work
             */
            void MoveCylinderCircleNodes(ModelPart& model_part, ClosedCylinder3D &cylinder){
                //loop over all elements
                for( ElementsArrayType::iterator ielem = model_part.ElementsBegin(); ielem != model_part.ElementsEnd(); ++ielem ){
                    if (ielem->GetValue(ACTIVATION_LEVEL) != 1){
                        //the current element is still active
                        GeometriesArrayType edges = ielem->GetGeometry().Edges();

                        bool inters = false;
                        int node0 = 0;//node on the cylinderwall
                        int node1 = 1;//node on the cylindercap
                        int helpnode = node1;
                        GeometriesArrayType::iterator edge;

                        //loop over all edges to find if the element is running through cylinder
                        for( GeometriesArrayType::iterator iedge = edges.begin(); iedge != edges.end(); ++iedge ){
                            if (iedge->size() == 3) node1 = 2;//quadratic element
                            if (cylinder.IsOnWall(iedge->pGetPoint(node0)->GetInitialPosition()) && cylinder.IsOnCap(iedge->pGetPoint(node1)->GetInitialPosition())){
                                //edge must run through the cylinder because one node is on the cylinder wall and one on the cap
                                inters = true;
                                edge = iedge;
                            }
                            else if (cylinder.IsOnCap(iedge->pGetPoint(node0)->GetInitialPosition()) && cylinder.IsOnWall(iedge->pGetPoint(node1)->GetInitialPosition())){
                                //edge must run through the cylinder because one node is on the cylinder wall and one on the cap
                                inters = true;
                                edge = iedge;
                                node0 = node1;
                                node1 = helpnode;
                            }
                        }
                        if (inters){
                        //element runs through cylinder
                            std::cout << "element is running through cylinder with nodes " << edge->pGetPoint(node0)->Id() << " and " << edge->pGetPoint(node1)->Id() << std::endl;
                            Point<3> old_position = edge->pGetPoint(node0)->GetInitialPosition();
                            Point<3> new_position = cylinder.ClosestPointOnCap(old_position);

                            //check if there is already a node at the new_position
                            WeakPointerVector< Node<3> >& rneigh_nodes = edge->pGetPoint(node0)->GetValue(NEIGHBOUR_NODES);
                            //loop pver all neighbouring nodes to find a node on the new_position
                            for (unsigned int i=0; i<rneigh_nodes.size(); ++i){
                                if(new_position == rneigh_nodes[i].GetInitialPosition()){
                                    //the current node is on the new_position
                                    Point<3> old_position2 = rneigh_nodes[i].GetInitialPosition();
                                    Point<3> new_position2 = old_position2 + 0.5*(old_position2 - old_position);
                                    WeakPointerVector<Element>& rneigh_el2 = rneigh_nodes[i].GetValue(NEIGHBOUR_ELEMENTS);
                                    //move the node
                                    if (!MoveNode( model_part, rneigh_el2, rneigh_nodes[i], new_position2)) std::cout << "node " << rneigh_nodes[i].Id() << " could not be moved from " << rneigh_nodes[i].GetInitialPosition() << " to " << new_position2 << std::endl;
                                    else{
                                        //move inner nodes
                                        rneigh_nodes[i].FastGetSolutionStepValue( IS_VISITED ) = 1;
                                        MoveInnerNodes(model_part, rneigh_el2);
                                    }
                                }
                            }

                            Point<3> step = (new_position - old_position) / 10;

                            WeakPointerVector<Element>& rneigh_el = edge->pGetPoint(node0)->GetValue(NEIGHBOUR_ELEMENTS);

                            //loop to move the node in 10 steps so that no flipping elements exist
                            for (int i=0; i<10; ++i){
                               //calculate new node position
                                Point<3> spstep = old_position + step*(i+1);

                                //move the node onto the cap
                                if (!MoveNode( model_part, rneigh_el, *(edge->pGetPoint(node0)), spstep)) std::cout << "node " << edge->pGetPoint(node0)->Id() << " could not be moved from " << edge->pGetPoint(node0)->GetInitialPosition() << " to " << spstep << std::endl;
                                else{
                                    //move inner nodes
                                    edge->pGetPoint(node0)->FastGetSolutionStepValue( IS_VISITED ) = 1;
                                    MoveInnerNodes(model_part, rneigh_el);
                                }
                            }
                        }
                    }
                }
            }

			/**
             * moves all inner nodes of the elements rneigh_el if a corner node is moved
             * @param model_part the current model part
             * @param rneigh_el all neighbouring elements of a moved node
             */
            void MoveInnerNodes(ModelPart& model_part, WeakPointerVector<Element>& rneigh_el){
                //loop over all neighbouring elements
                for (unsigned int i=0; i<rneigh_el.size(); ++i){
                    if (rneigh_el[i].GetValue(ACTIVATION_LEVEL) != 1){
                    //the current element is still active
                        GeometriesArrayType edges_local = rneigh_el[i].GetGeometry().Edges();
                        //loop over all edges of an element
                        for (GeometriesArrayType::iterator iedd = edges_local.begin(); iedd != edges_local.end(); ++iedd){
                            if (iedd->size() == 3){
                            //quadratic element with inner node
                                if ( iedd->pGetPoint(1)->FastGetSolutionStepValue( IS_VISITED ) != 3  && 
                                     (iedd->pGetPoint(0)->FastGetSolutionStepValue( IS_VISITED ) == 1 || 
                                      iedd->pGetPoint(2)->FastGetSolutionStepValue( IS_VISITED ) == 1 || 
                                      iedd->pGetPoint(0)->FastGetSolutionStepValue( IS_VISITED ) == 3 ||
                                      iedd->pGetPoint(2)->FastGetSolutionStepValue( IS_VISITED ) == 3 ) ){
                                //one of the edge's nodes is visited
                                    Point<3> new_position = 0.5 * (iedd->pGetPoint(0)->GetInitialPosition() + iedd->pGetPoint(2)->GetInitialPosition());
                                    Point<3> point = new_position - iedd->pGetPoint(1)->GetInitialPosition();
                                    double distance = pow(point[0],2.0)+pow(point[1],2.0)+pow(point[2],2.0);

                                    if ( distance > 0.000001 ){
                                    //the new node position is different from the old one
                                        WeakPointerVector<Element>& rneigh = iedd->pGetPoint(1)->GetValue(NEIGHBOUR_ELEMENTS);

                                        if (!MoveNode( model_part, rneigh, *(iedd->pGetPoint(1)), new_position)){std::cout << "could not move inner node " << iedd->pGetPoint(1)->Id() << " from " << iedd->pGetPoint(1)->GetInitialPosition() << " to " << new_position << std::endl;}
                                        iedd->pGetPoint(1)->FastGetSolutionStepValue( IS_VISITED ) = 2.0;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            /**
             * moves a node to a given position, preserving its deformation state
             * @param model_part the current model part
             * @param rNode the node to be moved
             * @param newPosition the new position of the node
             * @ref Han-Wen Nienhuys and A. Frank van der Stappen: Supporting cuts and
             * finite element deformation in interactive surgery simulation
             * UU-CS-2001-16 (2001)
             */
            bool MoveNode( ModelPart& model_part, WeakPointerVector<Element>& adjacent_elems, Node<3>& rNode, Point<3>& newPosition )
            {
                int elem_id = 0;

                Point<3> newLocalPoint;
                //find element, the new point lies within
                if( FindElement( adjacent_elems, newPosition, elem_id, newLocalPoint ) )
                {
                    Element::Pointer elem = model_part.pGetElement(elem_id);
                    //map all variables to the new position
                    MapNodalValues( rNode, elem, newLocalPoint );
                    //determine current deformation state in newLocalPoint
                    Point<3> undeformed_point( 0.0, 0.0, 0.0 );
                    undeformed_point = newPosition - rNode.GetSolutionStepValue(DISPLACEMENT);
                    //move node to new undeformed_point
                    rNode.SetInitialPosition( undeformed_point );
                    //reset element
                    ResetElements( rNode.GetValue(NEIGHBOUR_ELEMENTS) );
                    return true;
                }
                else
                    return false;
            }

            /**
             * Maps all DOF-values at the given local point to the given node 
			 * @param rNode given node
			 * @param elem element the local point lies inside
			 * @param localPoint point from which all values are maped
             */
            void MapNodalValues( Node<3>& rNode, Element::Pointer& elem, Point<3> localPoint )
            {
                if(rNode.HasDofFor(DISPLACEMENT_X)
                   || rNode.HasDofFor(DISPLACEMENT_Y)
                   || rNode.HasDofFor(DISPLACEMENT_Z))
                {
                    rNode.GetSolutionStepValue(DISPLACEMENT)=
                            MappedValue(elem, localPoint, DISPLACEMENT);
                    rNode.GetSolutionStepValue(DISPLACEMENT_NULL)=
                            MappedValue(elem, localPoint, DISPLACEMENT_NULL);
                    rNode.GetSolutionStepValue(DISPLACEMENT_EINS)=
                            MappedValue(elem, localPoint, DISPLACEMENT_EINS);
                    rNode.GetSolutionStepValue(DISPLACEMENT_NULL_DT)=
                            MappedValue(elem, localPoint, DISPLACEMENT_NULL_DT);
                    rNode.GetSolutionStepValue(ACCELERATION_NULL)=
                            MappedValue(elem, localPoint, ACCELERATION_NULL);
                    rNode.GetSolutionStepValue(DISPLACEMENT_OLD)=
                            MappedValue(elem, localPoint, DISPLACEMENT_OLD);
                }
                if(rNode.HasDofFor(WATER_PRESSURE))
                {
                    rNode.GetSolutionStepValue(WATER_PRESSURE)=
                            MappedValue(elem, localPoint, WATER_PRESSURE);
                    rNode.GetSolutionStepValue(WATER_PRESSURE_NULL)=
                            MappedValue(elem, localPoint, WATER_PRESSURE_NULL);
                    rNode.GetSolutionStepValue(WATER_PRESSURE_EINS)=
                            MappedValue(elem, localPoint, WATER_PRESSURE_EINS);
                    rNode.GetSolutionStepValue(WATER_PRESSURE_NULL_DT)=
                            MappedValue(elem, localPoint, WATER_PRESSURE_NULL_DT);
                    rNode.GetSolutionStepValue(WATER_PRESSURE_NULL_ACCELERATION)=
                            MappedValue(elem, localPoint, WATER_PRESSURE_NULL_ACCELERATION);
                }
                if(rNode.HasDofFor(AIR_PRESSURE))
                {
                    rNode.GetSolutionStepValue(AIR_PRESSURE)=
                            MappedValue(elem, localPoint, AIR_PRESSURE);
                    rNode.GetSolutionStepValue(AIR_PRESSURE_NULL)=
                            MappedValue(elem, localPoint, AIR_PRESSURE_NULL);
                    rNode.GetSolutionStepValue(WATER_PRESSURE_EINS)=
                            MappedValue(elem, localPoint, AIR_PRESSURE_EINS);
                    rNode.GetSolutionStepValue(AIR_PRESSURE_NULL_DT)=
                            MappedValue(elem, localPoint, AIR_PRESSURE_NULL_DT);
                    rNode.GetSolutionStepValue(AIR_PRESSURE_NULL_ACCELERATION)=
                            MappedValue(elem, localPoint, AIR_PRESSURE_NULL_ACCELERATION);
                }
                //for insitu_stress
                if(rNode.SolutionStepsDataHas(INSITU_STRESS) )
                {
                    rNode.GetSolutionStepValue(INSITU_STRESS)=
                            MappedValue(elem, localPoint, INSITU_STRESS);
                }
                

            }

            /**
             * Returns value of the variable at the local point
             * @param elem element the local point lies inside
             * @param localPoint point for which the variable is calculated
             * @param rThisVariable variabletype which should be calculated at the local point
             * @return newValue
             */
            //double MappedValue( Node<3>& node, Element::Pointer& elem, Point<3>& localPoint,
            double MappedValue( Element::Pointer& elem, Point<3>& localPoint,
                        const Variable<double>& rThisVariable )
            {
                double newValue = 0.0;

                for(unsigned int i= 0; i< elem->GetGeometry().size(); i++)
                {
                    newValue+=
                            elem->GetGeometry().ShapeFunctionValue( i, localPoint )
                            *elem->GetGeometry()[i].GetSolutionStepValue(rThisVariable);
                }
                return newValue;
            }

            /**
             * Returns value of the variable at the local point
             * @param elem element the local point lies inside
             * @param localPoint point for which the variable is calculated
             * @param rThisVariable variabletype which should be calculated at the local point
             * @return newValue
             */
            Vector MappedValue( Element::Pointer& elem, Point<3>& localPoint,
                                const Variable<array_1d<double, 3 > >& rThisVariable)
            {
                Vector newValue = ZeroVector(3);

                for(unsigned int i=0; i < elem->GetGeometry().size(); i++)
                {
                    newValue+=
                            elem->GetGeometry().ShapeFunctionValue( i, localPoint )
                            *elem->GetGeometry()[i].GetSolutionStepValue(rThisVariable);
                }
                return newValue;
            }

            /**
             * Returns value of the variable at the local point
             * @param elem element the local point lies inside
             * @param localPoint point for which the variable is calculated
             * @param rThisVariable variabletype which should be calculated at the local point
             * @return newValue
           */
            Vector MappedValue( Element::Pointer& elem, Point<3>& localPoint, const Variable<Kratos::Vector>& rThisVariable){
                Vector newValue = ZeroVector(6);

                for(unsigned int i=0; i < elem->GetGeometry().size(); i++) {
                    newValue+=
                            elem->GetGeometry().ShapeFunctionValue( i, localPoint )
                            *elem->GetGeometry()[i].GetSolutionStepValue(rThisVariable);
                }
                return newValue;
            }

            /**
             * Resets an element after the reference configuration has been changed
             * @param elements_set contains the elements which should be reseted
             */
            void ResetElements( WeakPointerVector<Element>& elements_set )
            {
                for( WeakPointerVector<Element>::iterator it = elements_set.begin(); 
                     it != elements_set.end(); ++it )
                {
                    //reset nodes to reference configuration
                    for( unsigned int i=0; i<(*it).GetGeometry().size(); i++ )
                    {
                        (*it).GetGeometry()[i].X() = (*it).GetGeometry()[i].X0();
                        (*it).GetGeometry()[i].Y() = (*it).GetGeometry()[i].Y0();
                        (*it).GetGeometry()[i].Z() = (*it).GetGeometry()[i].Z0();
                    }
                    //re-initialize element
                    it->Initialize();

                    //restore node positions
                    for( unsigned int i=0; i<(*it).GetGeometry().size(); i++ )
                    {
                        (*it).GetGeometry()[i].X() = (*it).GetGeometry()[i].X0()
                                +(*it).GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X );
                        (*it).GetGeometry()[i].Y() = (*it).GetGeometry()[i].Y0()
                                +(*it).GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y );
                        (*it).GetGeometry()[i].Z() = (*it).GetGeometry()[i].Z0()
                                +(*it).GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z );
                    }
                }
            }

            /**
             * Transfer of rThisVariable defined on integration points to corresponding 
             * nodal values. The transformation is done in a form that ensures a minimization
             * of L_2-norm error (/sum{rThisVariable- f(x)) whereas 
             * f(x)= /sum{shape_func_i*rThisVariable_i}
             * @param model_part current model_part
             * @param vicinityelements contains the Ids of the elements on which the transfer should be done
             * @param vicinitynodes the nodes belonging to vicinityelements
             * @param rThisVariable Matrix-Variable which should be transferred
             * @see TransferVariablesToNodes(ModelPart& model_part, std::vector<int>& vicinityelements, std::map<int, int>& vicinitynodes, Variable<Kratos::Vector>& rThisVariable)
             * @see TransferVariablesToNodes(ModelPart& model_part, std::vector<int>& vicinityelements, std::map<int, int>& vicinitynodes, Variable<double>& rThisVariable)
             * @ref Jiao&Heath: "Common-refinement-based data transfer...", Int.
             * Journal for numer. meth. in eng. 61 (2004) 2402--2427
             * WARNING: this may cause segmentation faults as the respective variables
             * will be created on nodal level while they are originally intended to be 
             * stored on integration points!
            */
            void TransferVariablesToNodes(ModelPart& model_part, std::vector<int>& vicinityelements, std::map<int, int>& vicinitynodes, Variable<Kratos::Matrix>& rThisVariable)
            {
                //loop over all vicinitynodes
                for (std::map<int, int>::iterator it=vicinitynodes.begin(); it != vicinitynodes.end(); ++it)
                    model_part.GetNode(it->first).GetSolutionStepValue(rThisVariable) = ZeroMatrix(3,3);

                //SetUpEquationSystem
                SpaceType::MatrixType M(vicinitynodes.size(), vicinitynodes.size());
                SpaceType::VectorType g(vicinitynodes.size());
                SpaceType::VectorType b(vicinitynodes.size());
                noalias(M)= ZeroMatrix(vicinitynodes.size(), vicinitynodes.size());

                //loop over all vicinityelements
                for( std::vector<int>::iterator it=vicinityelements.begin(); it != vicinityelements.end(); ++it) {
                    const IntegrationPointsArrayType& integration_points =  model_part.GetElement(*it).GetGeometry().IntegrationPoints(model_part.GetElement(*it).GetIntegrationMethod());

                    GeometryType::JacobiansType J(integration_points.size());
                    J = model_part.GetElement(*it).GetGeometry().Jacobian(J, model_part.GetElement(*it).GetIntegrationMethod());

                    const Matrix& Ncontainer = model_part.GetElement(*it).GetGeometry().ShapeFunctionsValues(model_part.GetElement(*it).GetIntegrationMethod());

                    Matrix InvJ(3,3);
                    double DetJ;
                    for(unsigned int point=0; point< integration_points.size(); point++){
                        MathUtils<double>::InvertMatrix(J[point],InvJ,DetJ);

                        double dV= DetJ*integration_points[point].Weight();

                        for(unsigned int prim=0; prim<model_part.GetElement(*it).GetGeometry().size() ; prim++) {
                            for(unsigned int sec=0; sec<model_part.GetElement(*it).GetGeometry().size() ; sec++) {
                                int ii = vicinitynodes.find(model_part.GetElement(*it).GetGeometry()[prim].Id())->second;
                                int jj = vicinitynodes.find(model_part.GetElement(*it).GetGeometry()[sec].Id())->second;
                                M(ii,jj) += Ncontainer(point, prim)*Ncontainer(point, sec)*dV;
                            }
                        }
                    }
                }

                for(unsigned int firstvalue=0; firstvalue<3; firstvalue++){
                    for(unsigned int secondvalue=0; secondvalue<3; secondvalue++){
                        noalias(g)= ZeroVector(vicinitynodes.size());
                        noalias(b)= ZeroVector(vicinitynodes.size());
                        //Transfer of GaussianVariables to Nodal Variablias via L_2-Minimization 
                        // see Jiao + Heath "Common-refinement-based data tranfer ..."
                        // International Journal for numerical methods in engineering 61 (2004) 2402--2427 
                        // for general description of L_2-Minimization
                        for( std::vector<int>::iterator it=vicinityelements.begin(); it != vicinityelements.end(); ++it) {
                            const IntegrationPointsArrayType& integration_points = model_part.GetElement(*it).GetGeometry().IntegrationPoints(model_part.GetElement(*it).GetIntegrationMethod());

                            GeometryType::JacobiansType J(integration_points.size());
                            J = model_part.GetElement(*it).GetGeometry().Jacobian(J, model_part.GetElement(*it).GetIntegrationMethod());
                            std::vector<Matrix> ValuesOnIntPoint(integration_points.size());

                            model_part.GetElement(*it).GetValueOnIntegrationPoints(rThisVariable, ValuesOnIntPoint, model_part.GetProcessInfo());

                            const Matrix& Ncontainer = model_part.GetElement(*it).GetGeometry().ShapeFunctionsValues(model_part.GetElement(*it).GetIntegrationMethod());

                            Matrix InvJ(3,3);
                            double DetJ;
                            for(unsigned int point=0; point< integration_points.size(); point++){
                                MathUtils<double>::InvertMatrix(J[point],InvJ,DetJ);

                                double dV= DetJ*integration_points[point].Weight();

                                for(unsigned int prim=0; prim < model_part.GetElement(*it).GetGeometry().size(); prim++){
                                    int ii = vicinitynodes.find(model_part.GetElement(*it).GetGeometry()[prim].Id())->second;
                                    b(ii) += ValuesOnIntPoint[point](firstvalue, secondvalue)*Ncontainer(point, prim)*dV;
                                }
                            }
                        }
                        SkylineLUFactorizationSolver<SpaceType, SpaceType>().Solve(M, g, b);
//                         CGSolverType(1.0e-8, 15000).Solve(M, g, b);
                        for (std::map<int, int>::iterator it=vicinitynodes.begin(); it != vicinitynodes.end(); ++it)
                            model_part.GetNode(it->first).GetSolutionStepValue(rThisVariable)(firstvalue, secondvalue) = g(it->second);
                    }//END firstvalue
                }//END secondvalue
            }

            /**
             * Transfer of rThisVariable defined on integration points to corresponding 
             * nodal values. The transformation is done in a form that ensures a minimization
             * of L_2-norm error (/sum{rThisVariable- f(x)) whereas 
             * f(x)= /sum{shape_func_i*rThisVariable_i}
             * @param model_part current model_part
             * @param vicinityelements contains the Ids of the elements on which the transfer should be done
             * @param vicinitynodes the nodes belonging to vicinityelements
             * @param rThisVariable Vector-Variable which should be transferred
               * @see TransferVariablesToNodes(ModelPart& model_part, std::vector<int>& vicinityelements, std::map<int, int>& vicinitynodes, Variable<Kratos::Matrix>& rThisVariable)
               * @see TransferVariablesToNodes(ModelPart& model_part, std::vector<int>& vicinityelements, std::map<int, int>& vicinitynodes, Variable<double>& rThisVariable)
             * @ref Jiao&Heath: "Common-refinement-based data transfer...", Int.
             * Journal for numer. meth. in eng. 61 (2004) 2402--2427
             * WARNING: this may cause segmentation faults as the respective variables
             * will be created on nodal level while they are originally intended to be 
             * stored on integration points!
             */
            void TransferVariablesToNodes(ModelPart& model_part, std::vector<int>& vicinityelements, std::map<int,int>& vicinitynodes, Variable<Kratos::Vector>& rThisVariable)
            {
                //loop over all vicinitynodes
                for (std::map<int, int>::iterator it=vicinitynodes.begin(); it != vicinitynodes.end(); ++it)
                    model_part.GetNode(it->first).GetSolutionStepValue(rThisVariable) = ZeroVector(6);

                //SetUpEquationSystem
                SpaceType::MatrixType M(vicinitynodes.size(), vicinitynodes.size());
                SpaceType::VectorType g(vicinitynodes.size());
                SpaceType::VectorType b(vicinitynodes.size());
                noalias(M)= ZeroMatrix(vicinitynodes.size(), vicinitynodes.size());

                //loop over all vicinityelements
                boost::progress_display show_progress( vicinityelements.size() );
                for( std::vector<int>::iterator it=vicinityelements.begin(); it != vicinityelements.end(); ++it) {
                    const IntegrationPointsArrayType& integration_points = model_part.GetElement(*it).GetGeometry().IntegrationPoints(model_part.GetElement(*it).GetIntegrationMethod());

                    GeometryType::JacobiansType J(integration_points.size());
                    J = model_part.GetElement(*it).GetGeometry().Jacobian(J, model_part.GetElement(*it).GetIntegrationMethod());

                    const Matrix& Ncontainer = model_part.GetElement(*it).GetGeometry().ShapeFunctionsValues(model_part.GetElement(*it).GetIntegrationMethod());

                    Matrix InvJ(3,3);
                    double DetJ;
                    for(unsigned int point=0; point< integration_points.size(); point++){
                        MathUtils<double>::InvertMatrix(J[point],InvJ,DetJ);


                        double dV= DetJ*integration_points[point].Weight();

                        for(unsigned int prim=0; prim < model_part.GetElement(*it).GetGeometry().size() ; prim++)
                            for(unsigned int sec=0; sec < model_part.GetElement(*it).GetGeometry().size(); sec++){
                                int ii = vicinitynodes.find(model_part.GetElement(*it).GetGeometry()[prim].Id())->second;
                                int jj = vicinitynodes.find(model_part.GetElement(*it).GetGeometry()[sec].Id())->second;
                                M(ii,jj) += Ncontainer(point, prim)*Ncontainer(point, sec)*dV;
                            }
                    }
                    ++show_progress;
                }

                for(unsigned int firstvalue=0; firstvalue<6; firstvalue++)
                {
                    noalias(g)= ZeroVector(vicinitynodes.size());
                    noalias(b)= ZeroVector(vicinitynodes.size());
                    //Transfer of GaussianVariables to Nodal Variablias via L_2-Minimization 
                    // see Jiao + Heath "Common-refinement-based data tranfer ..."
                    // International Journal for numerical methods in engineering 61 (2004) 2402--2427 
                    // for general description of L_2-Minimization
                    for( std::vector<int>::iterator it=vicinityelements.begin(); it != vicinityelements.end(); ++it) {
                        const IntegrationPointsArrayType& integration_points = model_part.GetElement(*it).GetGeometry().IntegrationPoints(model_part.GetElement(*it).GetIntegrationMethod());

                        GeometryType::JacobiansType J(integration_points.size());
                        J = model_part.GetElement(*it).GetGeometry().Jacobian(J, model_part.GetElement(*it).GetIntegrationMethod());

                        std::vector<Vector> ValuesOnIntPoint(integration_points.size());

                        model_part.GetElement(*it).GetValueOnIntegrationPoints(rThisVariable, ValuesOnIntPoint, model_part.GetProcessInfo());

                        const Matrix& Ncontainer = model_part.GetElement(*it).GetGeometry().ShapeFunctionsValues(model_part.GetElement(*it).GetIntegrationMethod());

                        Matrix InvJ(3,3);
                        double DetJ;
                        for(unsigned int point=0; point< integration_points.size(); point++){
                            MathUtils<double>::InvertMatrix(J[point],InvJ,DetJ);

                            double dV= DetJ*integration_points[point].Weight();

                            for(unsigned int prim=0; prim < model_part.GetElement(*it).GetGeometry().size(); prim++){
                                int ii = vicinitynodes.find(model_part.GetElement(*it).GetGeometry()[prim].Id())->second;
                                b(ii) += ValuesOnIntPoint[point](firstvalue)*Ncontainer(point, prim)*dV;
                            }
                        }
                    }
                    
                    SkylineLUFactorizationSolver<SpaceType, SpaceType>().Solve(M, g, b);
//                     CGSolverType(1.0e-8, 15000).Solve(M, g, b);
                    for (std::map<int, int>::iterator it=vicinitynodes.begin(); it != vicinitynodes.end(); ++it)
                        model_part.GetNode(it->first).GetSolutionStepValue(rThisVariable)(firstvalue)= g(it->second);
                }//END firstvalue
            }

            /**
             * Transfer of rThisVariable defined on integration points to corresponding 
             * nodal values. The transformation is done in a form that ensures a minimization
             * of L_2-norm error (/sum{rThisVariable- f(x)) whereas 
             * f(x)= /sum{shape_func_i*rThisVariable_i}
             * @param model_part current model_part
             * @param vicinityelements contains the Ids of the elements on which the transfer should be done
             * @param vicinitynodes the nodes belonging to vicinityelements
             * @param rThisVariable double-Variable which should be transferred
             * @see TransferVariablesToNodes(ModelPart& model_part, std::vector<int>& vicinityelements, std::map<int, int>& vicinitynodes, Variable<Kratos::Matrix>& rThisVariable)
             * @see TransferVariablesToNodes(ModelPart& model_part, std::vector<int>& vicinityelements, std::map<int, int>& vicinitynodes, Variable<Kratos::Vector>& rThisVariable)
             * @ref Jiao&Heath: "Common-refinement-based data transfer...", Int.
             * Journal for numer. meth. in eng. 61 (2004) 2402--2427
             * WARNING: this may cause segmentation faults as the respective variables
             * will be created on nodal level while they are originally intended to be 
             * stored on integration points!
             */
            void TransferVariablesToNodes(ModelPart& model_part, std::vector<int>& vicinityelements, std::map<int, int>& vicinitynodes, Variable<double>& rThisVariable)
            {
                //loop over all vicinitynodes
                for (std::map<int, int>::iterator it=vicinitynodes.begin(); it != vicinitynodes.end(); ++it)
                    model_part.GetNode(it->first).GetSolutionStepValue(rThisVariable) = 0.0;

                //SetUpEquationSystem
                SpaceType::MatrixType M(vicinitynodes.size(), vicinitynodes.size());
                SpaceType::VectorType g(vicinitynodes.size());
                SpaceType::VectorType b(vicinitynodes.size());
                noalias(M)= ZeroMatrix(vicinitynodes.size(), vicinitynodes.size());
                noalias(g)= ZeroVector(vicinitynodes.size());
                noalias(b)= ZeroVector(vicinitynodes.size());

                //Transfer of GaussianVariables to Nodal Variablias via L_2-Minimization 
                // see Jiao + Heath "Common-refinement-based data tranfer ..."
                // International Journal for numerical methods in engineering 61 (2004) 2402--2427 
                // for general description of L_2-Minimization
                for( std::vector<int>::iterator it=vicinityelements.begin(); it != vicinityelements.end(); ++it) {
                    const IntegrationPointsArrayType& integration_points = model_part.GetElement(*it).GetGeometry().IntegrationPoints(model_part.GetElement(*it).GetIntegrationMethod());

                    GeometryType::JacobiansType J(integration_points.size());
                    J = model_part.GetElement(*it).GetGeometry().Jacobian(J, model_part.GetElement(*it).GetIntegrationMethod());

                    std::vector<double> ValuesOnIntPoint(integration_points.size());

                    model_part.GetElement(*it).GetValueOnIntegrationPoints(rThisVariable, ValuesOnIntPoint, model_part.GetProcessInfo());

                    const Matrix& Ncontainer = model_part.GetElement(*it).GetGeometry().ShapeFunctionsValues(model_part.GetElement(*it).GetIntegrationMethod());

                    Matrix InvJ(3,3);
                    double DetJ;
                    for(unsigned int point=0; point< integration_points.size(); point++){
                        MathUtils<double>::InvertMatrix(J[point],InvJ,DetJ);

                        double dV= DetJ*integration_points[point].Weight();

                        for(unsigned int prim=0; prim < model_part.GetElement(*it).GetGeometry().size(); prim++){ 
                            int ii = vicinitynodes.find(model_part.GetElement(*it).GetGeometry()[prim].Id())->second;
                            b(ii) += ValuesOnIntPoint[point]*Ncontainer(point, prim)*dV;
                            for(unsigned int sec=0 ; sec < model_part.GetElement(*it).GetGeometry().size(); sec++){
                                int jj = vicinitynodes.find(model_part.GetElement(*it).GetGeometry()[sec].Id())->second;
                                M(ii,jj) += Ncontainer(point, prim)*Ncontainer(point, sec)*dV;
                            }
                        }
                    }
                }
                SkylineLUFactorizationSolver<SpaceType, SpaceType>().Solve(M, g, b);
//                 CGSolverType(1.0e-8, 15000).Solve(M, g, b);
                for (std::map<int, int>::iterator it=vicinitynodes.begin(); it != vicinitynodes.end(); ++it)
                    model_part.GetNode(it->first).GetSolutionStepValue(rThisVariable)= g(it->second);
            }

            /**
             * Transfer of rThisVariable stored on nodes to integration point via
             * approximation by shape functions
             * @param model_part current model part
             * @param vicinityelements contains the Ids of the elements on which the transfer should be done
             * @param rThisVariable Matrix-Variable which should be transferred
             * @see TransferVariablesToGaussPoints(ModelPart& model_part, Variable<Kratos::Vector>& rThisVariable)
             * @see TransferVariablesToGaussPoints(ModelPart& model_part, Variable<double>& rThisVariable)
             */
            void TransferVariablesToGaussPoints(ModelPart& model_part, std::vector<int>& vicinityelements, Variable<Kratos::Matrix>& rThisVariable){
                for( std::vector<int>::iterator it=vicinityelements.begin(); it != vicinityelements.end(); ++it) {
                    const IntegrationPointsArrayType& integration_points =
                            model_part.GetElement(*it).GetGeometry().IntegrationPoints(model_part.GetElement(*it).GetIntegrationMethod());

                    std::vector<Matrix> ValuesOnIntPoint(integration_points.size());

                    const Matrix& Ncontainer = model_part.GetElement(*it).GetGeometry().ShapeFunctionsValues(model_part.GetElement(*it).GetIntegrationMethod());

                    for( unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++){
                        ValuesOnIntPoint[PointNumber].resize(3,3,false);

                        noalias(ValuesOnIntPoint[PointNumber])= ZeroMatrix(3,3);

                        for(unsigned int node= 0; node < model_part.GetElement(*it).GetGeometry().size(); node++){
                            ValuesOnIntPoint[PointNumber]                                   +=Ncontainer(PointNumber, node) * model_part.GetElement(*it).GetGeometry()[node].GetSolutionStepValue(rThisVariable);
                        }
                    }

                    model_part.GetElement(*it).SetValueOnIntegrationPoints( rThisVariable, ValuesOnIntPoint, model_part.GetProcessInfo());
                }
            }
            /**
             * Transfer of rThisVariable stored on nodes to integration point via
             * approximation by shape functions
             * @param model_part current model part
             * @param vicinityelements contains the Ids of the elements on which the transfer should be done
             * @param rThisVariable Matrix-Variable which should be transferred
             * @see TransferVariablesToGaussPoints(ModelPart& model_part, Variable<Kratos::Matrix>& rThisVariable)
             * @see TransferVariablesToGaussPoints(ModelPart& model_part, Variable<double>& rThisVariable)
             */
            void TransferVariablesToGaussPoints(ModelPart& model_part, std::vector<int>& vicinityelements, Variable<Kratos::Vector>& rThisVariable){
                for( std::vector<int>::iterator it=vicinityelements.begin(); it != vicinityelements.end(); ++it) {
                    const IntegrationPointsArrayType& integration_points =
                            model_part.GetElement(*it).GetGeometry().IntegrationPoints(model_part.GetElement(*it).GetIntegrationMethod());

                    std::vector<Vector> ValuesOnIntPoint(integration_points.size());

                    const Matrix& Ncontainer = model_part.GetElement(*it).GetGeometry().ShapeFunctionsValues(model_part.GetElement(*it).GetIntegrationMethod());

                    for( unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++){
                        ValuesOnIntPoint[PointNumber].resize(6,false);

                        noalias(ValuesOnIntPoint[PointNumber])= ZeroVector(6);

                        for(unsigned int node= 0; node < model_part.GetElement(*it).GetGeometry().size(); node++){
                            ValuesOnIntPoint[PointNumber]                                   +=Ncontainer(PointNumber, node) * model_part.GetElement(*it).GetGeometry()[node].GetSolutionStepValue(rThisVariable);
                        }
                    }

                    model_part.GetElement(*it).SetValueOnIntegrationPoints( rThisVariable, ValuesOnIntPoint, model_part.GetProcessInfo());
                }
            }
            /**
             * Transfer of rThisVariable stored on nodes to integration point via
             * approximation by shape functions
             * @param model_part model_part on which the transfer should be done
             * @param vicinityelements contains the Ids of the elements on which the transfer should be done
             * @param rThisVariable double-Variable which should be transferred
             * @see TransferVariablesToGaussPoints(ModelPart& model_part, Variable<Kratos::Matrix>& rThisVariable)
             * @see TransferVariablesToGaussPoints(ModelPart& model_part, Variable<Kratos::Vector>& rThisVariable)
             */
            void TransferVariablesToGaussPoints(ModelPart& model_part, std::vector<int>& vicinityelements, Variable<double>& rThisVariable){
                for( std::vector<int>::iterator it=vicinityelements.begin(); it != vicinityelements.end(); ++it) {
                    const IntegrationPointsArrayType& integration_points =
                            model_part.GetElement(*it).GetGeometry().IntegrationPoints(model_part.GetElement(*it).GetIntegrationMethod());

                    std::vector<double> ValuesOnIntPoint(integration_points.size());

                    const Matrix& Ncontainer = model_part.GetElement(*it).GetGeometry().ShapeFunctionsValues(model_part.GetElement(*it).GetIntegrationMethod());

                    for( unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++){
                        ValuesOnIntPoint[PointNumber]= 0.0;

                        for(unsigned int node= 0; node < model_part.GetElement(*it).GetGeometry().size(); node++){
                            ValuesOnIntPoint[PointNumber]                                   +=Ncontainer(PointNumber, node) * model_part.GetElement(*it).GetGeometry()[node].GetSolutionStepValue(rThisVariable);
                        }
                    }

                    model_part.GetElement(*it).SetValueOnIntegrationPoints( rThisVariable, ValuesOnIntPoint, model_part.GetProcessInfo());
                }
            }

            void SetInsituStress(ModelPart& model_part){
                //predicting insitu stress
//                 std::cout << " set insitu stresses " << std::endl;
//                 for ( ElementsArrayType::ptr_iterator it=model_part.Elements().ptr_begin(); it!=model_part.Elements().ptr_end(); ++it){
//                     //get integration points
//                     GeometryType::IntegrationPointsArrayType integration_points = (*it)->GetGeometry().IntegrationPoints( (*it)->GetIntegrationMethod() );
//                     Point<3> coords;
//                     std::vector<Vector> InsituStresses;
//                     Vector InsituStress = ZeroVector(6);
//                     for( unsigned int i=0; i< integration_points.size(); i++ ){
//                         (*it)->GetGeometry().GlobalCoordinates(coords, integration_points[i]);
//                         for (int j=0; j<3; ++j)
//                             InsituStress[j] = coords[j];
//                         InsituStresses.push_back( InsituStress );
//                     }
//                     (*it)->SetValueOnIntegrationPoints(INSITU_STRESS,InsituStresses,model_part.GetProcessInfo() );
//                 }

                std::cout << " allocate memory for insitu stresses at nodes " << std::endl;
                Vector vec = ZeroVector(6);
                for (ModelPart::NodeIterator i=model_part.NodesBegin(); i<model_part.NodesEnd(); ++i)
                    i->GetSolutionStepValue(INSITU_STRESS)= vec;
            }
            
            void ControllGaussTransfer(ModelPart& model_part){
                double error;
                std::cout << "controll of Gauss transfer " << std::endl;
                for ( ElementsArrayType::ptr_iterator it=model_part.Elements().ptr_begin(); it!=model_part.Elements().ptr_end(); ++it){
                    //get integration points
                    GeometryType::IntegrationPointsArrayType integration_points = (*it)->GetGeometry().IntegrationPoints( (*it)->GetIntegrationMethod() );
                    std::vector<Vector> ValuesOnIntPoint(integration_points.size());
                    (*it)->GetValueOnIntegrationPoints(INSITU_STRESS, ValuesOnIntPoint, model_part.GetProcessInfo());
                    Point<3> coords;
                    for( unsigned int i=0; i<integration_points.size(); i++ ){
                        (*it)->GetGeometry().GlobalCoordinates(coords, integration_points[i]);
                        for (int j=0; j<3; ++j){
                            error = ValuesOnIntPoint[i](j) - coords[j];
                            if (error >= 0.000001 || error <= -0.000001)
                                std::cout << "error " << error << std::endl;
                        }
                    }
                }
            }
            
            void ControllTransferVariablesToNodes(ModelPart& model_part, std::map<int, int>& vicinitynodes){
                double error;

                for (std::map<int, int>::iterator it=vicinitynodes.begin(); it != vicinitynodes.end(); ++it){
                    Point<3> point = model_part.GetNode(it->first).GetInitialPosition();
                    for (int i=0; i<3; ++i){
                        error = model_part.GetNode(it->first).GetSolutionStepValue(INSITU_STRESS)[i] -  point[i];
                        if (error >= 0.000001 || error <= -0.000001)
                            std::cout << "error  " << error << std::endl;
                    }
                }
            }



            void ExtractCapNodes( boost::python::list surface_nodes, int len_surface_nodes )
            {
                for( int it=0; it < len_surface_nodes; it++ ){
                    boost::python::extract<int> x( surface_nodes[it] );
                    if( x.check() ) cap_nodes.push_back( x() );
                    else break;
                }
            }
            
            void DefineCapNodes(ModelPart& model_part, std::vector<int>& vicinityelements, ClosedCylinder3D& cylinder){
                if (cap_nodes.size() != 0) return;
                for( std::vector<int>::iterator it=vicinityelements.begin(); it != vicinityelements.end(); ++it){
                    if (model_part.GetElement(*it).GetValue(ACTIVATION_LEVEL) != 1){
                    //current element is still active
                        GeometriesArrayType edges = model_part.GetElement(*it).GetGeometry().Edges();

                        //corner nodes of edge
                        int node0 = 0;
                        int node1 = 1;

                        //loop over all edges to find if the element is running through cylinder
                        for( GeometriesArrayType::iterator iedge = edges.begin(); iedge != edges.end(); ++iedge ){
                            if (iedge->size() == 3) node1 = 2;//quadratic element
                            if (cylinder.IsIn(iedge->pGetPoint(node0)->GetInitialPosition())) cap_nodes.push_back(iedge->pGetPoint(node0)->Id());
                            if (cylinder.IsIn(iedge->pGetPoint(node1)->GetInitialPosition())) cap_nodes.push_back(iedge->pGetPoint(node1)->Id());
                        }
                    }
                }
                std::sort(cap_nodes.begin(), cap_nodes.end());
                std::vector<int>::iterator newend = std::unique(cap_nodes.begin(), cap_nodes.end());
                cap_nodes.erase(newend, cap_nodes.end());
            }
            
            void TestElements(ModelPart& model_part){
                std::cout << "in TestElements()" << std::endl; 
                for ( ElementsArrayType::ptr_iterator it=model_part.Elements().ptr_begin(); it!=model_part.Elements().ptr_end(); ++it){
                    std::cout << "elem " << (*it)->Id() << std::endl;
                    GeometriesArrayType edges = (*it)->GetGeometry().Edges();
                    for( GeometriesArrayType::iterator iedge = edges.begin(); iedge != edges.end(); ++iedge ) 
                        std::cout << " edge with nodes " << iedge->pGetPoint(0)->Id() << " " << iedge->pGetPoint(1)->Id() << " " << iedge->pGetPoint(2)->Id() << std::endl;
                }
            }

        private:
            
            std::vector<int> cap_nodes;
    };//Class NodeSnappingUtility
}//namespace Kratos.

#endif /* KRATOS_NODE_SNAPPING_UTILITY  defined */
