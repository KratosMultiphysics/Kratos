// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Geiser Armin, https://github.com/armingeiser
//
// ==============================================================================

#ifndef MAPPER_VERTEX_MORPHING_MESH_INDEPENDENT_NORMAL_H
#define MAPPER_VERTEX_MORPHING_MESH_INDEPENDENT_NORMAL_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "mapper_vertex_morphing_normal.h"
#include "processes/find_conditions_neighbours_process.h"
#include "utilities/math_utils.h"

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

class MapperVertexMorphingMeshIndependentNormal : public MapperVertexMorphingNormal
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperVertexMorphingMeshIndependentNormal
    KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphingMeshIndependentNormal);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperVertexMorphingMeshIndependentNormal( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters MapperSettings )
        : MapperVertexMorphingNormal(rOriginModelPart, rDestinationModelPart, MapperSettings)
    {
    }

    /// Destructor.
    virtual ~MapperVertexMorphingMeshIndependentNormal()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        if (mIsMappingInitialized == false)
        {
            FindNeighbourConditions(mrOriginModelPart);
            FindNeighbourConditions(mrDestinationModelPart);
        }

        MapperVertexMorphingNormal::Initialize();
    }

    void ComputeMappingMatrix() override
    {
        mDiagonalMassMatrix.resize(mrOriginModelPart.Nodes().size(), false);
        noalias(mDiagonalMassMatrix) = ZeroVector(mDiagonalMassMatrix.size());

        MapperVertexMorphingNormal::ComputeMappingMatrix();

        mSqrtOfInverseDiagonalMassMatrix.resize(mDiagonalMassMatrix.size(), false);

        #pragma omp parallel for
        for (int i=0; i<mDiagonalMassMatrix.size(); ++i)
        {
            mSqrtOfInverseDiagonalMassMatrix[i] = std::sqrt(1.0/mDiagonalMassMatrix[i]);
        }

        // mSqrtOfInverseDiagonalMassMatrix *= 1.0/norm_inf(mSqrtOfInverseDiagonalMassMatrix);
    }

    void Map( const Variable<array_3d> &rOriginVariable, const Variable<array_3d> &rDestinationVariable) override
    {
        KRATOS_ERROR << "Not implemented!" << std::endl;
    }

    void Map( const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable) override
    {
        KRATOS_ERROR << "Not implemented!" << std::endl;
    }

    void Map( const Variable<double> &rOriginVariable, const Variable<array_3d> &rDestinationVariable) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        ScaleOriginValues(rOriginVariable);
        MapperVertexMorphingNormal::Map(rOriginVariable, rDestinationVariable);
    }

    void InverseMap( const Variable<array_3d> &rDestinationVariable, const Variable<array_3d> &rOriginVariable) override
    {
        KRATOS_ERROR << "Not implemented!" << std::endl;
    }

    void InverseMap(const Variable<double> &rDestinationVariable, const Variable<double> &rOriginVariable) override
    {
        KRATOS_ERROR << "Not implemented!" << std::endl;
    }

    void InverseMap(const Variable<array_3d> &rDestinationVariable, const Variable<double> &rOriginVariable) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        MapperVertexMorphingNormal::InverseMap(rDestinationVariable, rOriginVariable);
        ScaleOriginValues(rOriginVariable);
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
        return "MapperVertexMorphingMeshIndependentNormal";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MapperVertexMorphingMeshIndependentNormal";
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
            double Aij = mpFilterFunction->compute_weight(node_j.Coordinates(),node_i.Coordinates());
            Aij *= mOriginNodalAreas[node_j.GetValue(MAPPING_ID)];

            // Add values to list
            list_of_weights[j_itr] += Aij;

            // Computed for integration of weighting function later using post-scaling
            sum_of_weights += Aij;
        }
    }

    void InitializeComputationOfMappingMatrix() override
    {
        // from base class
        MapperVertexMorphingNormal::InitializeComputationOfMappingMatrix();

        ComputeNodalAreas(mrOriginModelPart, mOriginNodalAreas);
        ComputeNodalAreas(mrDestinationModelPart, mDestinationNodalAreas);
    }

    // --------------------------------------------------------------------------
    void FillMappingMatrixWithWeights(  ModelPart::NodeType& destination_node,
                                        NodeVector& neighbor_nodes,
                                        unsigned int number_of_neighbors,
                                        std::vector<double>& list_of_weights,
                                        double& sum_of_weights )
    {
        unsigned int row_id = destination_node.GetValue(MAPPING_ID);
        const double area = mDestinationNodalAreas[row_id];
        for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
        {
            ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            int collumn_id = neighbor_node.GetValue(MAPPING_ID);
            const array_3d& normal = neighbor_node.FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);

            array_3d weight = normal * list_of_weights[neighbor_itr] / sum_of_weights;
            mMappingMatrix.insert_element(row_id*3+0,collumn_id,weight[0]);
            mMappingMatrix.insert_element(row_id*3+1,collumn_id,weight[1]);
            mMappingMatrix.insert_element(row_id*3+2,collumn_id,weight[2]);

            mDiagonalMassMatrix[collumn_id] += list_of_weights[neighbor_itr] / sum_of_weights * area;
        }
    }

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

    Vector mOriginNodalAreas;
    Vector mDestinationNodalAreas;
    Vector mSqrtOfInverseDiagonalMassMatrix;
    Vector mDiagonalMassMatrix;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void FindNeighbourConditions(ModelPart& rModelPart)
    {
        FindConditionsNeighboursProcess find_conditions_neighbours_process(rModelPart, rModelPart.GetProcessInfo()[DOMAIN_SIZE]);
        find_conditions_neighbours_process.Execute();
    }

    void ComputeNodalAreas(ModelPart& rModelPart, Vector& rNodalAreas){
        rNodalAreas.resize(rModelPart.Nodes().size(), false);
        noalias(rNodalAreas) = ZeroVector(rNodalAreas.size());
        for(auto& node_i : rModelPart.Nodes())
        {
            const int& i = node_i.GetValue(MAPPING_ID);

            // Get all neighbour conditions
            const GlobalPointersVector<Condition>& rConditions = node_i.GetValue(NEIGHBOUR_CONDITIONS);

            // loop conditions
            for(unsigned int c_itr=0; c_itr<rConditions.size(); c_itr++)
            {
                // Get geometry of current condition
                Condition rCondition = rConditions[c_itr];
                Condition::GeometryType& geom_i = rCondition.GetGeometry();
                rNodalAreas[i] += geom_i.DomainSize() / geom_i.size();
            }
        }
    }

    void ScaleOriginValues(const Variable<double>& rOriginVariable) {
        for (auto& node_i : mrOriginModelPart.Nodes())
        {
            const int mapping_id = node_i.GetValue(MAPPING_ID);
            node_i.FastGetSolutionStepValue(rOriginVariable) *= mSqrtOfInverseDiagonalMassMatrix[mapping_id];
        }
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class MapperVertexMorphingMeshIndependentNormal

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MAPPER_VERTEX_MORPHING_MESH_INDEPENDENT_NORMAL_H
