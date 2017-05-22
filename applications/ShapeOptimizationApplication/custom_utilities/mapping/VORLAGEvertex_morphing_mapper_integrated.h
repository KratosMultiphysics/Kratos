// ==============================================================================
/*
 KratosShapeOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

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
*/
//==============================================================================
//
//   Project Name:        KratosShape                            $
//   Created by:          $Author:    daniel.baumgaertner@tum.de $
//   Date:                $Date:                   December 2016 $
//   Revision:            $Revision:                         0.0 $
//
// ==============================================================================

#ifndef VERTEX_MORPHING_MAPPER_INTEGRATED_H
#define VERTEX_MORPHING_MAPPER_INTEGRATED_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <iomanip>      // for std::setprecision

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/python.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "../../kratos/includes/define.h"
#include "../../kratos/processes/process.h"
#include "../../kratos/includes/node.h"
#include "../../kratos/includes/element.h"
#include "../../kratos/includes/model_part.h"
#include "../../kratos/includes/kratos_flags.h"
#include "../../kratos/spatial_containers/spatial_containers.h"
#include "../../kratos/utilities/timer.h"
#include "../../kratos/processes/node_erase_process.h"
#include "../../kratos/utilities/binbased_fast_point_locator.h"
#include "../../kratos/utilities/normal_calculation_utils.h"
#include "../../kratos/spaces/ublas_space.h"
#include "../../kratos/processes/find_conditions_neighbours_process.h"
#include "../../kratos/utilities/math_utils.h"
#include "shape_optimization_application.h"
#include "filter_function.h"

// ==============================================================================

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

class VertexMorphingMapperIntegration
{
public:
    ///@name Type Definitions
    ///@{

    // ==========================================================================
    // Type definitions for better reading later
    // ==========================================================================
    typedef array_1d<double,3> array_3d;
    typedef Node < 3 > PointType;
    typedef Node < 3 > ::Pointer PointTypePointer;
    typedef std::vector<PointType::Pointer> PointVector;
    typedef std::vector<PointType::Pointer>::iterator PointIterator;
    typedef std::vector<double> DistanceVector;
    typedef std::vector<double>::iterator DistanceIterator;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    // ==========================================================================
    // Type definitions for linear algebra including sparse systems
    // ==========================================================================
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef typename SparseSpaceType::MatrixType SparseMatrixType;
    typedef typename SparseSpaceType::VectorType VectorType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    /// Pointer definition of VertexMorphingMapperIntegration
    KRATOS_CLASS_POINTER_DEFINITION(VertexMorphingMapperIntegration);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VertexMorphingMapperIntegration( ModelPart& model_part,
                             std::string filter_function_type,
                             bool use_mesh_preserving_mapping_matrix,
                             double filter_size,
                             const int max_nodes_affected )
        : mr_opt_model_part(model_part),
          m_filter_function_type(filter_function_type),
          m_use_mesh_preserving_mapping_matrix(use_mesh_preserving_mapping_matrix),
          m_filter_size(filter_size),
          m_max_nodes_affected(max_nodes_affected),
          m_number_of_design_variables(model_part.Nodes().size())
    {
        // Set precision for output
        std::cout.precision(12);

        // Create map to obtain local mapping matrix Id from global node Id
        unsigned int i = 0;
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            // Store local mapping matrix Id on the node
            node_i->SetValue(MAPPING_MATRIX_ID,i);

            // iterator design variable iterator i
            i++;
        }

        // store neighbouring information
        std::cout << "############# integration Mode" << std::endl;
        FindConditionsNeighboursProcess find_conditions_neighbours_process(mr_opt_model_part,
                                                        mr_opt_model_part.GetProcessInfo()[DOMAIN_SIZE]);
        find_conditions_neighbours_process.Execute();
        std::cout << "############# integration Mode" << std::endl;
    }

    /// Destructor.
    virtual ~VertexMorphingMapperIntegration()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    void compute_mapping_matrix()
    {
        KRATOS_TRY;

        m_mapping_matrix.clear();

        // Measure time
        boost::timer mapping_time;

        Element::IntegrationMethod m_integration_method = GeometryData::GI_GAUSS_2;

        std::cout << "\n> Starting computation of mapping matrix..." << std::endl;

        // Initialize filter matrix
        m_mapping_matrix.resize(m_number_of_design_variables,m_number_of_design_variables);

        // Initialize mapping matrix (zero possible entries)
        m_mapping_matrix.clear();

        // Creating an auxiliary list for the nodes to be searched on
        PointVector list_of_nodes;

        // Start constructing and computing the kdtree
        typedef Bucket< 3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
        typedef Tree< KDTreePartition<BucketType> > tree;

/*
        // create three with center coordinates of Conditions
        for (ModelPart::ConditionsContainerType::iterator condition_it = mr_opt_model_part.ConditionsBegin(); condition_it != mr_opt_model_part.ConditionsEnd(); ++condition_it)
        {
            Geometry< Node<3> >& rGeometry = condition_it->GetGeometry();
            list_of_nodes.push_back(*(rGeometry.Center()));
        }


        // loop nodes i

            // search for neighbour Conditions

            //loop neighbour Conditions
            for(ConditionsContainerType::iterator ic = rConds.begin(); ic!=rConds.end(); ic++)
            {

                // get geometry

                // filter function has to be scaled according r so that A = 1

                // do gauss integration of geometry times filter function

                // distribute that value to the nodes
            }
*/

        // starting calculating time of construction of the kdtree
        for (ModelPart::NodesContainerType::iterator node_it = mr_opt_model_part.NodesBegin(); node_it != mr_opt_model_part.NodesEnd(); ++node_it)
        {
            PointTypePointer pnode = *(node_it.base());

            // Putting the nodes of interest in an auxiliary list
            list_of_nodes.push_back(pnode);
        }

        // Arrays needed for spatial search
        PointVector nodes_affected(m_max_nodes_affected);
        DistanceVector resulting_squared_distances(m_max_nodes_affected);

        // Compute tree with the node positions
        tree nodes_tree(list_of_nodes.begin(), list_of_nodes.end(), m_max_nodes_affected);

        // Prepare Weighting function to be used in the mapping
        FilterFunction filter_function( m_filter_function_type, m_filter_size );

        // Loop over all design variables
        for (ModelPart::NodeIterator node_itr = mr_opt_model_part.NodesBegin(); node_itr != mr_opt_model_part.NodesEnd(); ++node_itr)
        {
            // Get node information
            int i_ID = node_itr->Id();
            ModelPart::NodeType& node_i = mr_opt_model_part.Nodes()[i_ID];
            array_3d i_coord = node_i.Coordinates();

            // Get tID of the node in the mapping matrix
            int i = node_i.GetValue(MAPPING_MATRIX_ID);

            // Perform spatial search for current node
            unsigned int number_of_nodes_affected;
            number_of_nodes_affected = nodes_tree.SearchInRadius(node_i, m_filter_size, nodes_affected.begin(),resulting_squared_distances.begin(), m_max_nodes_affected);

            // loop neighbour nodes and get the conditions
            // first the geometries have to be stored at the nodes once
            std::vector<int> condition_ids;
            for(unsigned int j_itr = 0 ; j_itr<number_of_nodes_affected ; j_itr++)
            {
                // Get node information
                int j_ID = nodes_affected[j_itr]->Id();
                ModelPart::NodeType& node_j = mr_opt_model_part.Nodes()[j_ID];
                const WeakPointerVector<Condition>& rC = node_j.GetValue(NEIGHBOUR_CONDITIONS);
                //KRATOS_WATCH(rC.size());
                for (unsigned int i=0; i< rC.size(); i++)
                {
                    condition_ids.push_back(rC[i].Id());
                }
            }

            // remove duplicated entries from vector
            sort(condition_ids.begin(),condition_ids.end());
            condition_ids.erase( unique( condition_ids.begin(), condition_ids.end() ), condition_ids.end() );

            double sum_weights = 0.0;

            // loop conditions

            for(std::vector<int>::iterator c_itr=condition_ids.begin() ; c_itr!=condition_ids.end() ; c_itr++)
            {
                // Get geometry information of current condition
                Condition::GeometryType& geom_i = mr_opt_model_part.Conditions()[*c_itr].GetGeometry();
                unsigned int n_nodes = geom_i.size();

                //KRATOS_WATCH(n_nodes);

                // Evaluate shape functions of design surface according specified integration methdod
                const Condition::GeometryType::IntegrationPointsArrayType& integration_points = geom_i.IntegrationPoints(m_integration_method);
                const unsigned int number_of_integration_points = integration_points.size();
                const Matrix& N_container = geom_i.ShapeFunctionsValues(m_integration_method);

                // TODO determinant is not defined for elements im able to use
                //Kratos::Vector det_J;
                //geom_i.DeterminantOfJacobian(det_J, m_integration_method);
                Condition::GeometryType::JacobiansType J;
                geom_i.Jacobian(J,m_integration_method);

                //KRATOS_WATCH(number_of_integration_points);

                for ( unsigned int point_number = 0; point_number < number_of_integration_points; point_number++ )
                {
                    //KRATOS_WATCH(number_of_integration_points);
                    // Get weight for integration
                    double integration_weight = integration_points[point_number].Weight();
                    //KRATOS_WATCH(integration_weight);

                    // Get FEM-shape-function-value for current integration point
                    Vector N_FEM_GPi = row( N_container, point_number);
                    //KRATOS_WATCH(N_FEM_GPi);

                    // Get node information
                    PointType::CoordinatesArrayType gp_i_coord;
                    geom_i.GlobalCoordinates(gp_i_coord, integration_points[point_number].Coordinates());
                    //KRATOS_WATCH(gp_i_coord);

                    // TODO use gramian to determine det_J
                    Kratos::Matrix J_Ti = trans(J[point_number]);
                    Kratos::Matrix JJ_T = prod(J_Ti,J[point_number] );
                    //KRATOS_WATCH(J[point_number]);
                    //KRATOS_WATCH(J_Ti);
                    //KRATOS_WATCH(JJ_T);
                    double detGram;
                    if (mr_opt_model_part.GetProcessInfo()[DOMAIN_SIZE] == 2)
                        detGram = JJ_T(0,0);
                    else
                        detGram = MathUtils<double>::Det(JJ_T);
                    //KRATOS_WATCH(detGram);
                    double det_J = sqrt(detGram);
                    //KRATOS_WATCH(det_J);

                    for (unsigned int node_ctr=0; node_ctr<n_nodes; node_ctr++)
                    {
                        // Get the id of the node in the mapping matrix
                        int j = geom_i[node_ctr].GetValue(MAPPING_MATRIX_ID);

                        // evaluate shape function at gauss point
                        double N_j_gpi = geom_i.ShapeFunctionValue(point_number,node_ctr,m_integration_method);

                        // Computation of weight according specified weighting function
                        // Note that we did not compute the square root of the distances to save this expensive computation (it is not needed here)
                        double Aij = filter_function.compute_weight(gp_i_coord,node_i.Coordinates());
                        Aij *= N_j_gpi;
                        //KRATOS_WATCH(Aij);
                        Aij *= integration_weight;
                        if (mr_opt_model_part.GetProcessInfo()[DOMAIN_SIZE] == 2)
                            Aij /= m_filter_size;

                        // consider jacobian
                        Aij *= det_J;

                        // In this way we implicitly preserve the in-plane mesh quality (pure Heuristic)
                        m_mapping_matrix(i+0,j) += Aij;

                        // Computed for integration of weighting function later using post-scaling
                        sum_weights += Aij;
                    }
                }
            }

            // Here we perform the post scaling by the sum of all weights

            for(unsigned int j_itr = 0 ; j_itr<mr_opt_model_part.Nodes().size() ; j_itr++)
            {
                // Get the id of the node in the mapping matrix
                int j = j_itr;

                // Scaling of the weights
                //m_mapping_matrix(3*i+0,j) /= sum_weights;
                //m_mapping_matrix(3*i+1,j) /= sum_weights;
                if (m_mapping_matrix(i,j)>0.0)
                    m_mapping_matrix(i,j) /= sum_weights;
            }

        }
        if (false)
        {
            for (unsigned int i=0; i<m_mapping_matrix.size1()/3;i++)
            {
                for (unsigned int j = 0; j<m_mapping_matrix.size2(); j++)
                    std::cout << m_mapping_matrix(i*3+1,j) <<",";
                std::cout << std::endl;
            }
        }

        // Console output for information
        std::cout << "> Time needed for computation of mapping matrix: " << mapping_time.elapsed() << " s" << std::endl;

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void map_sensitivities_to_design_space( bool constraint_given )
    {
        KRATOS_TRY;

        // Measure time
        boost::timer mapping_time;
        std::cout << "\n> Starting mapping of sensitivities to design space..." << std::endl;

        // use scaling with M matrix
        // TODO
        //SparseMatrixType M;
        //SparseMatrixType M_inv;
        //double detM;
        //M.resize(m_number_of_design_variables,m_number_of_design_variables);
        //M_inv.resize(m_number_of_design_variables,m_number_of_design_variables);
        // ask Daniel for multiplication and solution/inverse

        //SparseSpaceType::TransposeMult(m_mapping_matrix,m_mapping_matrix,M_inv);
        //MathUtils<double>::InvertMatrix( M, M_inv, detM);

        // Map objective sensitivities

        // Assign nodal sensitivities to vector used for mapping
        VectorType dJdX;
        dJdX.resize(m_number_of_design_variables*3);
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            int i = node_i->GetValue(MAPPING_MATRIX_ID);
            array_3d node_sens = node_i->FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY);
            dJdX[i*3+0] = node_sens[0];
            dJdX[i*3+1] = node_sens[1];
            dJdX[i*3+2] = node_sens[2];
        }
        VectorType dCdX;
        if(constraint_given)
        {
            std::cout << "NOT IMPLEMENTED ARMIN " << std::endl;
            exit(-1);
            dCdX.resize(m_number_of_design_variables*3);
            for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
            {
                int i = node_i->GetValue(MAPPING_MATRIX_ID);
                array_3d node_sens = node_i->FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY);
                dCdX[i*3+0] = node_sens[0];
                dCdX[i*3+1] = node_sens[1];
                dCdX[i*3+2] = node_sens[2];
            }
        }

        // Perform mapping according to the specified mapping matrix
        if(m_use_mesh_preserving_mapping_matrix)
        {
            VectorType dJds;
            VectorType dJdn;
            dJds.resize(m_number_of_design_variables);
            dJdn.resize(m_number_of_design_variables);

            // project values on normals
            for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
            {
                int i = node_i->GetValue(MAPPING_MATRIX_ID);
                array_3d node_normal = node_i->FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);
                dJdn[i] = dJdX[i*3+0] * node_normal[0]+
                          dJdX[i*3+1] * node_normal[1]+
                          dJdX[i*3+2] * node_normal[2];
            }

            SparseSpaceType::Mult(m_mapping_matrix,dJdn,dJds);
            // SparseSpaceType::TransposeMult(m_mapping_matrix,dJdn,dJds);

            // Assign results to nodal variables
            for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
            {
                int i = node_i->GetValue(MAPPING_MATRIX_ID);
                VectorType dJds_i = ZeroVector(3);
                dJds_i(0) = dJds[i];
                node_i->FastGetSolutionStepValue(MAPPED_OBJECTIVE_SENSITIVITY) = dJds_i;
            }

            // Repeat mapping to map constraint sensitivities
            if(constraint_given)
            {
                std::cout << "NOT IMPLEMENTED ARMIN " << std::endl;
                exit(-1);
                // Perform mapping
                VectorType dCds;
                dCds.resize(m_number_of_design_variables);
                SparseSpaceType::TransposeMult(m_mapping_matrix,dCdX,dCds);

                // Assign results to nodal variables
                for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
                {
                    int i = node_i->GetValue(MAPPING_MATRIX_ID);
                    VectorType dCds_i = ZeroVector(3);
                    dCds_i(0) = dCds[i];
                    node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY) = dCds_i;
                }
            }
        }
        else // if complete matrix is used for mapping
        {
            VectorType dJds(m_number_of_design_variables*3);
            VectorType dJdX_dir(m_number_of_design_variables);
            VectorType dJds_dir(m_number_of_design_variables);
            for (int dir=0; dir<3; dir++){
                // project values on normals
                for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
                {
                    int i = node_i->GetValue(MAPPING_MATRIX_ID);
                    dJdX_dir[i] = dJdX[i*3+dir];
                }
                SparseSpaceType::Mult(m_mapping_matrix,dJdX_dir,dJds_dir);
                //SparseSpaceType::TransposeMult(m_mapping_matrix,dJdX_dir,dJds_dir);
                for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
                {
                    int i = node_i->GetValue(MAPPING_MATRIX_ID);
                    dJds[i*3+dir] = dJds_dir[i];
                }
            }

            // Assign results to nodal variables
            for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
            {
                int i = node_i->GetValue(MAPPING_MATRIX_ID);
                VectorType dJds_i = ZeroVector(3);
                dJds_i(0) = dJds[i*3+0];
                dJds_i(1) = dJds[i*3+1];
                dJds_i(2) = dJds[i*3+2];
                node_i->FastGetSolutionStepValue(MAPPED_OBJECTIVE_SENSITIVITY) = dJds_i;
            }

            // Repeat mapping to map constraint sensitivities
            if(constraint_given)
            {
                std::cout << "NOT IMPLEMENTED ARMIN " << std::endl;
                exit(-1);
                // Perform mapping
                VectorType dCds;
                dCds.resize(m_number_of_design_variables*3);
                SparseSpaceType::TransposeMult(m_mapping_matrix,dCdX,dCds);

                // Assign results to nodal variables
                for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
                {
                    int i = node_i->GetValue(MAPPING_MATRIX_ID);
                    VectorType dCds_i = ZeroVector(3);
                    dCds_i(0) = dCds[i*3+0];
                    dCds_i(1) = dCds[i*3+1];
                    dCds_i(2) = dCds[i*3+2];
                    node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY) = dCds_i;
                }
            }
        }

        // Console output for information
        std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void map_design_update_to_geometry_space()
    {
        KRATOS_TRY;

        // Measure time of mapping
        boost::timer mapping_time;
        std::cout << "\n> Starting mapping of design update to geometry space..." << std::endl;

        VectorType ds;
        VectorType dx;
        dx.resize(m_number_of_design_variables*3);

        // temporary vectors
        VectorType dx_dir(m_number_of_design_variables);
        VectorType ds_dir(m_number_of_design_variables);

        // Assign design update to vector which shall be mapped (depending specified mapping matrix )
        if(m_use_mesh_preserving_mapping_matrix)
        {
            for (int dir=0; dir<3; dir++)
        	{
                for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
                {
                    int i = node_i->GetValue(MAPPING_MATRIX_ID);
                    array_3d node_normal = node_i->FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);
                    ds_dir[i] = node_i->FastGetSolutionStepValue(DESIGN_UPDATE_X)*node_normal[dir];
                }

                SparseSpaceType::Mult(m_mapping_matrix,ds_dir,dx_dir);

                for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
                {
                    int i = node_i->GetValue(MAPPING_MATRIX_ID);
                    dx[i*3+dir] = dx_dir[i];
                }
            }
        }
        else
        {
            for (int dir=0; dir<3; dir++)
        	{
                ds.resize(m_number_of_design_variables*3);
                for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
                {
                    int i = node_i->GetValue(MAPPING_MATRIX_ID);
                    // TODO ARMIN
                    ds_dir[i] = node_i->FastGetSolutionStepValue(DESIGN_UPDATE)[dir];
                    //ds[3*i+0] = 0.;//node_i->FastGetSolutionStepValue(DESIGN_UPDATE_X);
                    //ds[3*i+1] = 1.;//node_i->FastGetSolutionStepValue(DESIGN_UPDATE_Y);
                    //ds[3*i+2] = 0.;//node_i->FastGetSolutionStepValue(DESIGN_UPDATE_Z);
                }
                // Perform mapping to compute shape update
                SparseSpaceType::Mult(m_mapping_matrix,ds_dir,dx_dir);

                for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
                {
                    int i = node_i->GetValue(MAPPING_MATRIX_ID);
                    dx[i*3+dir] = dx_dir[i];
                }
            }
        }


        // Assign dx as nodal shape updates and update coordinates accordingly
        for (ModelPart::NodeIterator node_i = mr_opt_model_part.NodesBegin(); node_i != mr_opt_model_part.NodesEnd(); ++node_i)
        {
            // If shape update deactivated, set it to zero
            if(node_i->FastGetSolutionStepValue(SHAPE_UPDATES_DEACTIVATED))
            {
                array_3d zero_array(3,0.0);
                noalias(node_i->FastGetSolutionStepValue(SHAPE_UPDATE)) = zero_array;
            }
            // In case it is not deactivated, it is checked if it is on a specified boundary beyond which no update is wanted
            else
            {
                // Read shape update from solution vector
                int i = node_i->GetValue(MAPPING_MATRIX_ID);
                array_3d shape_update;
                shape_update[0] = dx[3*i+0];
                shape_update[1] = dx[3*i+1];
                shape_update[2] = dx[3*i+2];

                // If node is on specified boundary, project shape update on specified boundary plane
                // I.e. remove component that is normal to the boundary plane
                if(node_i->FastGetSolutionStepValue(IS_ON_BOUNDARY))
                {
                    array_3d boundary_plane = node_i->FastGetSolutionStepValue(BOUNDARY_PLANE);
                    shape_update = shape_update - (inner_prod(shape_update,boundary_plane))*boundary_plane/norm_2(boundary_plane);
                }

                // Assign shape update to nodal variable
                noalias(node_i->FastGetSolutionStepValue(SHAPE_UPDATE)) = shape_update;

                // Update coordinates
                node_i->X() += shape_update[0];
                node_i->Y() += shape_update[1];
                node_i->Z() += shape_update[2];

                // Add final shape update to previous updates
                noalias(node_i->FastGetSolutionStepValue(SHAPE_CHANGE_ABSOLUTE)) += shape_update;
            }
        }

        // Console output for information
        std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;

        KRATOS_CATCH("");
    }

    // ==============================================================================

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
        return "VertexMorphingMapperIntegration";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "VertexMorphingMapperIntegration";
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

    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    ModelPart& mr_opt_model_part;
    std::string m_filter_function_type;
    bool m_use_mesh_preserving_mapping_matrix;
    const double m_filter_size;
    const int m_max_nodes_affected;
    const unsigned int m_number_of_design_variables;
    // ==============================================================================
    // General working arrays
    // ==============================================================================
    SparseMatrixType m_mapping_matrix;

    ///@}
    ///@name Private Operators
    ///@{


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
//      VertexMorphingMapperIntegration& operator=(VertexMorphingMapperIntegration const& rOther);

    /// Copy constructor.
//      VertexMorphingMapperIntegration(VertexMorphingMapperIntegration const& rOther);


    ///@}

}; // Class VertexMorphingMapperIntegration

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // VERTEX_MORPHING_MAPPER_INTEGRATED_H
