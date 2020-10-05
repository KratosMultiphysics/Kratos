// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Geiser Armin, https://github.com/armingeiser
//
// ==============================================================================

#ifndef MAPPER_VERTEX_MORPHING_MESH_INDEPENDENT_H
#define MAPPER_VERTEX_MORPHING_MESH_INDEPENDENT_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "mapper_vertex_morphing.h"
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

class MapperVertexMorphingMeshIndependent : public MapperVertexMorphing
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperVertexMorphingMeshIndependent
    KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphingMeshIndependent);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperVertexMorphingMeshIndependent( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters MapperSettings )
        : MapperVertexMorphing(rOriginModelPart, rDestinationModelPart, MapperSettings)
    {
    }

    /// Destructor.
    virtual ~MapperVertexMorphingMeshIndependent()
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

        MapperVertexMorphing::Initialize();
    }

    void ComputeMappingMatrix() override
    {
        MapperVertexMorphing::ComputeMappingMatrix();

        KRATOS_WATCH("diagonal")

        // compute the scaling for the mesh indpendency
        Vector diagonal_mass_matrix(mrOriginModelPart.Nodes().size());
        SparseSpaceType::TransposeMult(mMappingMatrix, mDestinationNodalAreas, diagonal_mass_matrix);

        mSqrtOfInverseDiagonalMassMatrix.resize(diagonal_mass_matrix.size(), false);

        // diagonal_mass_matrix *= 1.0/norm_inf(diagonal_mass_matrix);

        #pragma omp parallel for
        for (int i=0; i<diagonal_mass_matrix.size(); ++i)
        {
            mSqrtOfInverseDiagonalMassMatrix[i] = std::sqrt(1.0/diagonal_mass_matrix[i]);
        }
    }

    void Map( const Variable<array_3d> &rOriginVariable, const Variable<array_3d> &rDestinationVariable) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        ScaleOriginValues(rOriginVariable);
        MapperVertexMorphing::Map(rOriginVariable, rDestinationVariable);
    }

    void Map( const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        ScaleOriginValues(rOriginVariable);
        MapperVertexMorphing::Map(rOriginVariable, rDestinationVariable);
    }

    void InverseMap( const Variable<array_3d> &rDestinationVariable, const Variable<array_3d> &rOriginVariable) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        MapperVertexMorphing::InverseMap(rDestinationVariable, rOriginVariable);
        ScaleOriginValues(rOriginVariable);
    }

    void InverseMap(const Variable<double> &rDestinationVariable, const Variable<double> &rOriginVariable) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        MapperVertexMorphing::InverseMap(rDestinationVariable, rOriginVariable);
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
        return "MapperVertexMorphingMeshIndependent";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MapperVertexMorphingMeshIndependent";
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
                                        double& sum_of_weights ) override
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
        MapperVertexMorphing::InitializeComputationOfMappingMatrix();

        ComputeNodalAreas(mrOriginModelPart, mOriginNodalAreas);
        ComputeNodalAreas(mrDestinationModelPart, mDestinationNodalAreas);
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

    template <typename T>
    void ScaleOriginValues(const T& rOriginVariable) {
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

}; // Class MapperVertexMorphingMeshIndependent

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MAPPER_VERTEX_MORPHING_MESH_INDEPENDENT_H
