//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes

#include "includes/checks.h"
#include "includes/model_part.h"
#include "utilities/divide_triangle_2d_3.h"
#include "utilities/divide_tetrahedra_3d_4.h"

// Application includes
#include "embedded_skin_visualization_process.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/* Public functions *******************************************************/
EmbeddedSkinVisualizationProcess::EmbeddedSkinVisualizationProcess(
    const ModelPart& rModelPart,
    ModelPart& rVisualizationModelPart)
    : Process(), mrModelPart(rModelPart), mrVisualizationModelPart(rVisualizationModelPart) {
}

EmbeddedSkinVisualizationProcess::EmbeddedSkinVisualizationProcess(
    const ModelPart& rModelPart,
    ModelPart& rVisualizationModelPart,
    Parameters& rParameters)
    : Process(), mrModelPart(rModelPart), mrVisualizationModelPart(rVisualizationModelPart) {

    Parameters default_parameters( R"(
    {
        "model_part_name"               : "default_model_part_name"
        "visualization_model_part_name" : "visualization_default_model_part_name"
    })");

    rParameters.ValidateAndAssignDefaults(default_parameters);

}

void EmbeddedSkinVisualizationProcess::ExecuteInitialize() {

    KRATOS_TRY;

    // Required variables check
    const auto& r_node = *mrModelPart.NodesBegin();
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);

    KRATOS_CATCH("");
}

void EmbeddedSkinVisualizationProcess::ExecuteAfterOutputStep() {

    #pragma omp parallel for
    for (int i_elem = 0; i_elem < static_cast<int>(mrModelPart.NumberOfElements()); ++i_elem)
    {

        auto it_elem = mrModelPart.ElementsBegin() + i_elem;

        // Get element geometry
        auto r_geometry = it_elem->GetGeometry();
        const unsigned int n_nodes = r_geometry.PointsNumber();

        // Check if the element is split
        double n_pos (0.0), n_neg(0.0);
        Vector nodal_distances(n_nodes);
        for (unsigned int i_node = 0; i_node < n_nodes; ++i_node){
            const double dist = r_geometry[i_node].FastGetSolutionStepValue(DISTANCE);
            nodal_distances[i_node] = dist;
            if (dist > 0.0)
                n_pos++;
            else
                n_neg++;
        }
        bool is_split = (n_pos > 0 && n_neg > 0) ? true : false;

        if (is_split){

            // Set the split utility and compute the splitting pattern
            DivideGeometry::Pointer p_split_utility = this->GetGeometrySplitUtility(r_geometry, nodal_distances);
            p_split_utility->GenerateDivision();
            p_split_utility->GenerateIntersectionsSkin();

            // Save the geometries from the splitting pattern in the visualization model part

        }



    }
}

/* Protected functions ****************************************************/

/* Private functions ******************************************************/
DivideGeometry::Pointer EmbeddedSkinVisualizationProcess::GetGeometrySplitUtility(
    const Geometry<Node<3>> &rGeometry,
    const Vector& rNodalDistances){

    // Get the geometry type
    const GeometryData::KratosGeometryType geometry_type = rGeometry.GetGeometryType();
    
    // Return the split utility according to the geometry tyoe
    switch (geometry_type){
        case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
            return Kratos::make_shared<DivideTriangle2D3>(rGeometry, rNodalDistances);
        case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
            return Kratos::make_shared<DivideTetrahedra3D4>(rGeometry, rNodalDistances);
        default:
            KRATOS_ERROR << "Asking for non-implemented geometry splitting utility.";
    }   
}

};  // namespace Kratos.
