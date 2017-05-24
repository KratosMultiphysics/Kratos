// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef MAPPER_VERTEX_MORPHING_INTEGRATION_H
#define MAPPER_VERTEX_MORPHING_INTEGRATION_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

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

class MapperVertexMorphingIntegration
{
public:
    ///@name Type Definitions
    ///@{

    // Type definitions for better reading later
    typedef array_1d<double,3> array_3d;
    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef std::vector<double>::iterator DoubleVectorIterator;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    // Type definitions for linear algebra including sparse systems
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef SparseSpaceType::MatrixType SparseMatrixType;
    typedef SparseSpaceType::VectorType VectorType;

    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;

    /// Pointer definition of MapperVertexMorphingIntegration
    KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphingIntegration);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperVertexMorphingIntegration( ModelPart& designSurface, Parameters& optimizationSettings )
        : mrDesignSurface( designSurface ),
          mNumberOfDesignVariables(designSurface.Nodes().size()),
          mFilterType( optimizationSettings["design_variables"]["filter"]["filter_function_type"].GetString() ),
          mFilterRadius( optimizationSettings["design_variables"]["filter"]["filter_radius"].GetDouble() ),
          mMaxNumberOfNeighbors( optimizationSettings["design_variables"]["filter"]["max_nodes_in_filter_radius"].GetInt() )
    {
        CreateListOfNodesOfDesignSurface();
        CreateFilterFunction();
        InitializeMappingVariables();
        AssignMappingIds();
        SetIntegrationMethod(optimizationSettings);
        FindNeighbourConditions();

        // optional flags
        try {
            mConsistentBackwardMapping = optimizationSettings["design_variables"]["consistent_backward_mapping"].GetBool();
        }
        catch (...)
        {
            mConsistentBackwardMapping = false;
        }

        ComputeMappingMatrix(); // TODO this should not happen in constructor
    }

    /// Destructor.
    virtual ~MapperVertexMorphingIntegration()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    void CreateListOfNodesOfDesignSurface()
    {
        for (ModelPart::NodesContainerType::iterator node_it = mrDesignSurface.NodesBegin(); node_it != mrDesignSurface.NodesEnd(); ++node_it)
        {
            NodeTypePointer pnode = *(node_it.base());
            mListOfNodesOfDesignSurface.push_back(pnode);
        }
    }

    // --------------------------------------------------------------------------
    void CreateFilterFunction()
    {
        mpFilterFunction = boost::shared_ptr<FilterFunction>(new FilterFunction(mFilterType, mFilterRadius));
    }

    // --------------------------------------------------------------------------
    void SetIntegrationMethod( Parameters& optimizationSettings )
    {
        int integrationMethod = optimizationSettings["design_variables"]["integration_method"].GetInt();
        mAreaWeightedNodeSum = false;
        if (integrationMethod == 0)
            mAreaWeightedNodeSum = true;
        else if (integrationMethod == 1)
            mIntegrationMethod = GeometryData::GI_GAUSS_1;
        else if (integrationMethod == 2)
            mIntegrationMethod = GeometryData::GI_GAUSS_2;
        else if (integrationMethod == 3)
            mIntegrationMethod = GeometryData::GI_GAUSS_3;
        else if (integrationMethod == 4)
            mIntegrationMethod = GeometryData::GI_GAUSS_4;
        else if (integrationMethod == 5)
            mIntegrationMethod = GeometryData::GI_GAUSS_5;
        else
        {
            std::cout << "\n> Integration method " << integrationMethod << " not valid! USING DEFAULT: 2 " << std::endl;
            mIntegrationMethod = GeometryData::GI_GAUSS_2;
        }
    }

    // --------------------------------------------------------------------------
    void FindNeighbourConditions()
    {

            // store neighbouring information
        std::cout << "> Computing neighbour conditions ..." << std::endl;
        FindConditionsNeighboursProcess find_conditions_neighbours_process(mrDesignSurface,
                                                        mrDesignSurface.GetProcessInfo()[DOMAIN_SIZE]);
        find_conditions_neighbours_process.Execute();
    }

    // --------------------------------------------------------------------------
    void InitializeMappingVariables()
    {
        mMappingMatrix.resize(mNumberOfDesignVariables,mNumberOfDesignVariables);
        mMappingMatrix.clear();

        x_variables_in_design_space.resize(mNumberOfDesignVariables,0.0);
        y_variables_in_design_space.resize(mNumberOfDesignVariables,0.0);
        z_variables_in_design_space.resize(mNumberOfDesignVariables,0.0);

        x_variables_in_geometry_space.resize(mNumberOfDesignVariables,0.0);
        y_variables_in_geometry_space.resize(mNumberOfDesignVariables,0.0);
        z_variables_in_geometry_space.resize(mNumberOfDesignVariables,0.0);
    }

    // --------------------------------------------------------------------------
    void AssignMappingIds()
    {
        unsigned int i = 0;
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
            node_i->SetValue(MAPPING_ID,i++);
    }

    // --------------------------------------------------------------------------
        void ComputeMappingMatrix()
   {
        boost::timer timer;
        std::cout << "> Computing mapping matrix to perform mapping..." << std::endl;

        CreateSearchTreeWithAllNodesOnDesignSurface();
        ComputeEntriesOfMappingMatrix();

        std::cout << "> Mapping matrix computed in: " << timer.elapsed() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    void CreateSearchTreeWithAllNodesOnDesignSurface()
    {
        mpSearchTree = boost::shared_ptr<KDTree>(new KDTree(mListOfNodesOfDesignSurface.begin(), mListOfNodesOfDesignSurface.end(), mBucketSize));
    }

    // --------------------------------------------------------------------------
    void ComputeEntriesOfMappingMatrix()
    {
        // Arrays needed for spatial search
        NodeVector neighbor_nodes(mMaxNumberOfNeighbors);
        std::vector<double> resulting_squared_distances(mMaxNumberOfNeighbors);

        // Loop over all design variables
        for (ModelPart::NodeIterator node_itr = mrDesignSurface.NodesBegin(); node_itr != mrDesignSurface.NodesEnd(); ++node_itr)
        {
            ModelPart::NodeType& node_i = *node_itr; //TODO renaming to geometryNode -->see design_node in orig version

            unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( node_i,
                                                                    mFilterRadius,
                                                                    neighbor_nodes.begin(),
                                                                    resulting_squared_distances.begin(),
                                                                    mMaxNumberOfNeighbors );

            std::vector<double> list_of_weights( number_of_neighbors, 0.0 );
            double sum_of_weights = 0.0;

            ThrowWarningIfMaxNodeNeighborsReached( node_i, number_of_neighbors );
            ComputeWeightForAllNeighbors( node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );
            FillMappingMatrixWithWeights( node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );
        }
    }

    // --------------------------------------------------------------------------
    void ThrowWarningIfMaxNodeNeighborsReached( ModelPart::NodeType& given_node, unsigned int number_of_neighbors )
    {
        if(number_of_neighbors >= mMaxNumberOfNeighbors)
            std::cout << "\n> WARNING!!!!! For node " << given_node.Id() << " and specified filter radius, maximum number of neighbor nodes (=" << mMaxNumberOfNeighbors << " nodes) reached!" << std::endl;
    }

    // --------------------------------------------------------------------------
    void ComputeWeightForAllNeighbors(  ModelPart::NodeType& node_i,
                                        NodeVector& neighbor_nodes,
                                        unsigned int number_of_neighbors,
                                        std::vector<double>& list_of_weights,
                                        double& sum_of_weights )
    {
        for(unsigned int j_itr = 0 ; j_itr<number_of_neighbors ; j_itr++)
        {
            // Get node information
            ModelPart::NodeType& node_j = *neighbor_nodes[j_itr];

            // Get all neighbour conditions
            const WeakPointerVector<Condition>& rConditions = node_j.GetValue(NEIGHBOUR_CONDITIONS);

            // loop conditions
            for(unsigned int c_itr=0; c_itr<rConditions.size(); c_itr++)
            {
                // Get geometry of current condition
                Condition rCondition = rConditions[c_itr];
                Condition::GeometryType& geom_i = rCondition.GetGeometry();

                if (mAreaWeightedNodeSum){

                    // Computation of weight according specified weighting function
                    // Note that we did not compute the square root of the distances to save this expensive computation (it is not needed here)
                    double Aij = mpFilterFunction->compute_weight(node_j.Coordinates(),node_i.Coordinates());
                    Aij *= geom_i.DomainSize();

                    // Add values to list
                    list_of_weights[j_itr] += Aij;

                    // Computed for integration of weighting function later using post-scaling
                    sum_of_weights += Aij;
                }
                else
                {
                    // Get geometry information of current condition
                    unsigned int n_nodes = geom_i.size();
                    int localNodeIndex = -1;
                    for (unsigned int node_ctr=0; node_ctr<n_nodes; node_ctr++)
                    {
                        if (geom_i[node_ctr].Id() == node_j.Id())
                            localNodeIndex = node_ctr;
                    }

                    // Evaluate shape functions of design surface according specified integration method
                    const Condition::GeometryType::IntegrationPointsArrayType& integrationPoints = geom_i.IntegrationPoints(mIntegrationMethod);
                    const unsigned int numberOfIntegrationPoints = integrationPoints.size();
                    const Matrix& N_container = geom_i.ShapeFunctionsValues(mIntegrationMethod);

                    for ( unsigned int pointNumber = 0; pointNumber < numberOfIntegrationPoints; pointNumber++ )
                    {

                        // Get FEM-shape-function-value for current integration point
                        Vector N_FEM_GPi = row( N_container, pointNumber);

                        // Get gp coordinates
                        NodeType::CoordinatesArrayType gp_i_coord;
                        geom_i.GlobalCoordinates(gp_i_coord, integrationPoints[pointNumber].Coordinates());

                        // Computation of weight according specified weighting function
                        // Note that we did not compute the square root of the distances to save this expensive computation (it is not needed here)
                        double Aij = mpFilterFunction->compute_weight(gp_i_coord,node_i.Coordinates());

                        // multiply with evaluation of shape function at gauss point
                        Aij *= geom_i.ShapeFunctionValue(pointNumber,localNodeIndex,mIntegrationMethod);;

                        // Get weight for integration
                        Aij *= integrationPoints[pointNumber].Weight();

                        // consider jacobian
                        Aij *= geom_i.DeterminantOfJacobian(pointNumber,mIntegrationMethod);

                        // Add values to list
                        list_of_weights[j_itr] += Aij;

                        // Computed for integration of weighting function later using post-scaling
                        sum_of_weights += Aij;
                    }
                }
            }
        }
    }

    // --------------------------------------------------------------------------
    void FillMappingMatrixWithWeights(  ModelPart::NodeType& design_node,
                                        NodeVector& neighbor_nodes,
                                        unsigned int number_of_neighbors,
                                        std::vector<double>& list_of_weights,
                                        double& sum_of_weights )
    {
        unsigned int row_id = design_node.GetValue(MAPPING_ID);
        for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
        {
            ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            int collumn_id = neighbor_node.GetValue(MAPPING_ID);

            double weight = list_of_weights[neighbor_itr] / sum_of_weights;
            mMappingMatrix.push_back(row_id,collumn_id,weight);
        }
    }

    // --------------------------------------------------------------------------
    void MapToDesignSpace( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rNodalVariableInDesignSpace )
    {
        boost::timer mapping_time;
        std::cout << "\n> Starting to map " << rNodalVariable.Name() << " to design space..." << std::endl;

        RecomputeMappingMatrixIfGeometryHasChanged();
        PrepareVectorsForMappingToDesignSpace( rNodalVariable );
        if (mConsistentBackwardMapping)
            MultiplyVectorsWithConsistentBackwardMappingMatrix();
        else
            MultiplyVectorsWithTransposeMappingMatrix();
        AssignResultingDesignVectorsToNodalVariable( rNodalVariableInDesignSpace );

        std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    void MapToGeometrySpace( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rNodalVariableInGeometrySpace )
    {
        boost::timer mapping_time;
        std::cout << "\n> Starting to map " << rNodalVariable.Name() << " to geometry space..." << std::endl;

        RecomputeMappingMatrixIfGeometryHasChanged();
        PrepareVectorsForMappingToGeometrySpace( rNodalVariable );
        MultiplyVectorsWithMappingMatrix();
        AssignResultingGeometryVectorsToNodalVariable( rNodalVariableInGeometrySpace );

        std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    void RecomputeMappingMatrixIfGeometryHasChanged()
    {
        if(HasGeometryChanged())
        {
            InitializeComputationOfMappingMatrix();
            ComputeMappingMatrix();
        }
    }

    // --------------------------------------------------------------------------
    void PrepareVectorsForMappingToDesignSpace( const Variable<array_3d> &rNodalVariable )
    {
        x_variables_in_design_space.clear();
        y_variables_in_design_space.clear();
        z_variables_in_design_space.clear();
        x_variables_in_geometry_space.clear();
        y_variables_in_geometry_space.clear();
        z_variables_in_geometry_space.clear();

        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
            int i = node_i->GetValue(MAPPING_ID);
            array_3d& nodal_variable = node_i->FastGetSolutionStepValue(rNodalVariable);
            x_variables_in_geometry_space[i] = nodal_variable[0];
            y_variables_in_geometry_space[i] = nodal_variable[1];
            z_variables_in_geometry_space[i] = nodal_variable[2];
        }
    }
    // --------------------------------------------------------------------------
    void PrepareVectorsForMappingToGeometrySpace( const Variable<array_3d> &rNodalVariable )
    {
        x_variables_in_design_space.clear();
        y_variables_in_design_space.clear();
        z_variables_in_design_space.clear();
        x_variables_in_geometry_space.clear();
        y_variables_in_geometry_space.clear();
        z_variables_in_geometry_space.clear();

        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
            int i = node_i->GetValue(MAPPING_ID);
            array_3d& nodal_variable = node_i->FastGetSolutionStepValue(rNodalVariable);
            x_variables_in_design_space[i] = nodal_variable[0];
            y_variables_in_design_space[i] = nodal_variable[1];
            z_variables_in_design_space[i] = nodal_variable[2];
        }
    }

    // --------------------------------------------------------------------------
    void MultiplyVectorsWithTransposeMappingMatrix()
    {
        SparseSpaceType::TransposeMult(mMappingMatrix,x_variables_in_geometry_space,x_variables_in_design_space);
        SparseSpaceType::TransposeMult(mMappingMatrix,y_variables_in_geometry_space,y_variables_in_design_space);
        SparseSpaceType::TransposeMult(mMappingMatrix,z_variables_in_geometry_space,z_variables_in_design_space);
    }

    // --------------------------------------------------------------------------
    void MultiplyVectorsWithConsistentBackwardMappingMatrix()
    {
        noalias(x_variables_in_design_space) = prod(mMappingMatrix,x_variables_in_geometry_space);
        noalias(y_variables_in_design_space) = prod(mMappingMatrix,y_variables_in_geometry_space);
        noalias(z_variables_in_design_space) = prod(mMappingMatrix,z_variables_in_geometry_space);
    }

    // --------------------------------------------------------------------------
    void MultiplyVectorsWithMappingMatrix()
    {
        noalias(x_variables_in_geometry_space) = prod(mMappingMatrix,x_variables_in_design_space);
        noalias(y_variables_in_geometry_space) = prod(mMappingMatrix,y_variables_in_design_space);
        noalias(z_variables_in_geometry_space) = prod(mMappingMatrix,z_variables_in_design_space);
    }

    // --------------------------------------------------------------------------
    void AssignResultingDesignVectorsToNodalVariable( const Variable<array_3d> &rNodalVariable )
    {
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
            int i = node_i->GetValue(MAPPING_ID);

            Vector node_vector = ZeroVector(3);
            node_vector(0) = x_variables_in_design_space[i];
            node_vector(1) = y_variables_in_design_space[i];
            node_vector(2) = z_variables_in_design_space[i];
            node_i->FastGetSolutionStepValue(rNodalVariable) = node_vector;
        }
    }

    // --------------------------------------------------------------------------
    void AssignResultingGeometryVectorsToNodalVariable( const Variable<array_3d> &rNodalVariable )
    {
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
            int i = node_i->GetValue(MAPPING_ID);

            Vector node_vector = ZeroVector(3);
            node_vector(0) = x_variables_in_geometry_space[i];
            node_vector(1) = y_variables_in_geometry_space[i];
            node_vector(2) = z_variables_in_geometry_space[i];
            node_i->FastGetSolutionStepValue(rNodalVariable) = node_vector;
        }
    }

    // --------------------------------------------------------------------------
    bool HasGeometryChanged()
    {
        double sumOfAllCoordinates = 0.0;
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
            array_3d& coord = node_i->Coordinates();
            sumOfAllCoordinates += coord[0] + coord[1] + coord[2];
        }

        if(IsFirstMappingOperation())
        {
            mControlSum = sumOfAllCoordinates;
            return false;
        }
        else if (mControlSum == sumOfAllCoordinates)
            return false;
        else
        {
            mControlSum = sumOfAllCoordinates;
            return true;
        }
    }

    // --------------------------------------------------------------------------
    void InitializeComputationOfMappingMatrix()
    {
        mpSearchTree.reset();
        mMappingMatrix.clear();
    }


    // --------------------------------------------------------------------------
    bool IsFirstMappingOperation()
    {
        if(mControlSum == 0.0)
            return true;
        else
             return false;
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
        return "MapperVertexMorphingIntegration";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MapperVertexMorphingIntegration";
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
    ModelPart& mrDesignSurface;
    const unsigned int mNumberOfDesignVariables;
    std::string mFilterType;
    double mFilterRadius;
    unsigned int mMaxNumberOfNeighbors;
    FilterFunction::Pointer mpFilterFunction;
    bool mConsistentBackwardMapping;
    Element::IntegrationMethod mIntegrationMethod;
    bool mAreaWeightedNodeSum;

    // ==============================================================================
    // Variables for spatial search
    // ==============================================================================
    unsigned int mBucketSize = 100;
    NodeVector mListOfNodesOfDesignSurface;
    KDTree::Pointer mpSearchTree;

    // ==============================================================================
    // Variables for mapping
    // ==============================================================================
    SparseMatrixType mMappingMatrix;
    Vector x_variables_in_design_space, y_variables_in_design_space, z_variables_in_design_space;
    Vector x_variables_in_geometry_space, y_variables_in_geometry_space, z_variables_in_geometry_space;
    double mControlSum = 0.0;

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
//      MapperVertexMorphingIntegration& operator=(MapperVertexMorphingIntegration const& rOther);

    /// Copy constructor.
//      MapperVertexMorphingIntegration(MapperVertexMorphingIntegration const& rOther);


    ///@}

}; // Class MapperVertexMorphingIntegration

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MAPPER_VERTEX_MORPHING_INTEGRATION_H
