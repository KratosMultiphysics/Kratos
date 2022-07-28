// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Geiser Armin, https://github.com/armingeiser
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "mapper_vertex_morphing_improved_integration.h"
#include "mapper_vertex_morphing.h"
#include "processes/find_conditions_neighbours_process.h"
#include "utilities/math_utils.h"

// ==============================================================================
namespace Kratos
{

void MapperVertexMorphingImprovedIntegration::Initialize()
{
    SetIntegrationMethod();
    FindNeighbourConditions();

    MapperVertexMorphing::Initialize();
}

void MapperVertexMorphingImprovedIntegration::Update()
{
    FindNeighbourConditions();

    MapperVertexMorphing::Update();
}

void MapperVertexMorphingImprovedIntegration::SetIntegrationMethod()
{
    std::string integration_method = mMapperSettings["integration_method"].GetString();
    int number_of_gauss_points = mMapperSettings["number_of_gauss_points"].GetInt();

    mAreaWeightedNodeSum = true;

}

void MapperVertexMorphingImprovedIntegration::FindNeighbourConditions()
{
    KRATOS_INFO("ShapeOpt") << "Computing neighbour conditions ..." << std::endl;
    FindConditionsNeighboursProcess find_conditions_neighbours_process(mrOriginModelPart, mrOriginModelPart.GetProcessInfo()[DOMAIN_SIZE]);
    find_conditions_neighbours_process.Execute();
}

void MapperVertexMorphingImprovedIntegration::ComputeWeightForAllNeighbors(  ModelPart::NodeType& node_i,
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
        if (mAreaWeightedNodeSum){
            // Computation of weight according specified weighting function
            // Note that we did not compute the square root of the distances to save this expensive computation (it is not needed here)
            double Aij = mpFilterFunction->ComputeWeight(node_j.Coordinates(),node_i.Coordinates());
            // array_3d& node_j_damp = node_j.GetValue(DAMPING_FACTOR);
            // array_3d& node_i_damp = node_i.GetValue(DAMPING_FACTOR);

            list_of_weights[j_itr] += (Aij * nodalAreas[node_j.GetValue(MAPPING_ID)]);
            sum_of_weights += (Aij * nodalAreas[node_j.GetValue(MAPPING_ID)]);            

            // if(node_i_damp[0]<0.00000000001){
            //     list_of_weights[j_itr] = 0.0;
            //     sum_of_weights = 1.0;
            // }
            // else
            // {
            //     list_of_weights[j_itr] += Aij  * pow(node_j_damp[0],4);
            //     sum_of_weights += (Aij * nodalAreas[node_j.GetValue(MAPPING_ID)]);
            // }            
        }
    }
}

void MapperVertexMorphingImprovedIntegration::Map( const Variable<array_3d> &rOriginVariable, const Variable<array_3d> &rDestinationVariable)
{

    for(auto& node_i : mrOriginModelPart.Nodes())
    {
        // int i = node_i.GetValue(MAPPING_ID);
        array_3d& r_nodal_variable = node_i.FastGetSolutionStepValue(rOriginVariable);
        r_nodal_variable[0] *= nodalAreas[node_i.GetValue(MAPPING_ID)];
        r_nodal_variable[1] *= nodalAreas[node_i.GetValue(MAPPING_ID)];
        r_nodal_variable[2] *= nodalAreas[node_i.GetValue(MAPPING_ID)];
    }
    MapperVertexMorphing::Map(rOriginVariable,rDestinationVariable);
}

void MapperVertexMorphingImprovedIntegration::InverseMap( const Variable<array_3d> &rDestinationVariable, const Variable<array_3d> &rOriginVariable)
{
    MapperVertexMorphing::InverseMap(rDestinationVariable,rOriginVariable);
}

void MapperVertexMorphingImprovedIntegration::InverseMap( const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable)
{
}

void MapperVertexMorphingImprovedIntegration::Map( const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable)
{
}

void MapperVertexMorphingImprovedIntegration::InitializeComputationOfMappingMatrix()
{
    // from base class
    MapperVertexMorphing::InitializeComputationOfMappingMatrix();

    // necessary for this class
    if (mAreaWeightedNodeSum)
    {
        nodalAreas.resize(mrOriginModelPart.Nodes().size(),0.0);
        for(auto& node_i : mrOriginModelPart.Nodes())
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
                nodalAreas[i] += geom_i.DomainSize() / geom_i.size();
            }
        }
    }
}

}  // namespace Kratos.
