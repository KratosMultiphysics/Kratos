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

#ifndef MAPPER_VERTEX_MORPHING_IMPROVED_INTEGRATION_H
#define MAPPER_VERTEX_MORPHING_IMPROVED_INTEGRATION_H

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

class MapperVertexMorphingImprovedIntegration : public MapperVertexMorphing
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperVertexMorphingImprovedIntegration
    KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphingImprovedIntegration);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperVertexMorphingImprovedIntegration( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters MapperSettings )
        : MapperVertexMorphing(rOriginModelPart, rDestinationModelPart, MapperSettings)
    {
    }

    /// Destructor.
    virtual ~MapperVertexMorphingImprovedIntegration()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // --------------------------------------------------------------------------
    void Initialize() override
    {
        if (mIsMappingInitialized == false)
        {
            SetIntegrationMethod();
            FindNeighbourConditions();
        }

        MapperVertexMorphing::Initialize();
    }
    // --------------------------------------------------------------------------

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
        return "MapperVertexMorphingImprovedIntegration";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MapperVertexMorphingImprovedIntegration";
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

    Element::IntegrationMethod mIntegrationMethod;
    bool mAreaWeightedNodeSum;
    std::vector<double> nodalAreas;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // --------------------------------------------------------------------------
    void SetIntegrationMethod()
    {
        std::string integration_method = mMapperSettings["integration_method"].GetString();
        int number_of_gauss_points = mMapperSettings["number_of_gauss_points"].GetInt();

        if (integration_method.compare("area_weighted_sum") == 0)
            mAreaWeightedNodeSum = true;
        else if (integration_method.compare("gauss_integration") == 0)
        {
            mAreaWeightedNodeSum = false;
            if (number_of_gauss_points == 1)
                mIntegrationMethod = GeometryData::GI_GAUSS_1;
            else if (number_of_gauss_points == 2)
                mIntegrationMethod = GeometryData::GI_GAUSS_2;
            else if (number_of_gauss_points == 3)
                mIntegrationMethod = GeometryData::GI_GAUSS_3;
            else if (number_of_gauss_points == 4)
                mIntegrationMethod = GeometryData::GI_GAUSS_4;
            else if (number_of_gauss_points == 5)
                mIntegrationMethod = GeometryData::GI_GAUSS_5;
            else
            {
                KRATOS_WARNING("ShapeOpt::MapperVertexMorphingImprovedIntegration") << "\n> Number_of_gauss_points: " << number_of_gauss_points << " not valid! Using default: 2 " << std::endl;
                mIntegrationMethod = GeometryData::GI_GAUSS_2;
            }
        }
        else{
            KRATOS_ERROR << "\n> Integration method " << integration_method << " unknown!" << std::endl;
        }
    }

    // --------------------------------------------------------------------------
    void FindNeighbourConditions()
    {
        KRATOS_INFO("ShapeOpt") << "Computing neighbour conditions ..." << std::endl;
        FindConditionsNeighboursProcess find_conditions_neighbours_process(mrOriginModelPart, mrOriginModelPart.GetProcessInfo()[DOMAIN_SIZE]);
        find_conditions_neighbours_process.Execute();
    }

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
            if (mAreaWeightedNodeSum){
                // Computation of weight according specified weighting function
                // Note that we did not compute the square root of the distances to save this expensive computation (it is not needed here)
                double Aij = mpFilterFunction->compute_weight(node_j.Coordinates(),node_i.Coordinates());
                Aij *= nodalAreas[node_j.GetValue(MAPPING_ID)];

                // Add values to list
                list_of_weights[j_itr] += Aij;

                // Computed for integration of weighting function later using post-scaling
                sum_of_weights += Aij;
            }
            else
            {
                const GlobalPointersVector<Condition>& rConditions = node_j.GetValue(NEIGHBOUR_CONDITIONS);

                // loop conditions
                for(unsigned int c_itr=0; c_itr<rConditions.size(); c_itr++)
                {
                    // Get geometry of current condition
                    Condition rCondition = rConditions[c_itr];
                    Condition::GeometryType& geom_i = rCondition.GetGeometry();

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
    void InitializeComputationOfMappingMatrix() override
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
    // --------------------------------------------------------------------------

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
//      MapperVertexMorphingImprovedIntegration& operator=(MapperVertexMorphingImprovedIntegration const& rOther);

    /// Copy constructor.
//      MapperVertexMorphingImprovedIntegration(MapperVertexMorphingImprovedIntegration const& rOther);


    ///@}

}; // Class MapperVertexMorphingImprovedIntegration

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MAPPER_VERTEX_MORPHING_IMPROVED_INTEGRATION_H
