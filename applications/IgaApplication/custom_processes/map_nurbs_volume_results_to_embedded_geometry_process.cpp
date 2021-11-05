//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Manuel Messmer

// Project includes
#include "map_nurbs_volume_results_to_embedded_geometry_process.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{

    MapNurbsVolumeResultsToEmbeddedGeometryProcess::MapNurbsVolumeResultsToEmbeddedGeometryProcess(
        Model& rModel, Parameters ThisParameters) : mrModel(rModel), mThisParameters(ThisParameters)
    {
        mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
        KRATOS_ERROR_IF_NOT( rModel.HasModelPart( mThisParameters["main_model_part_name"].GetString()) )
            << "MapNurbsVolumeResultsToEmbeddedGeometryProcess: Model Part '" <<  mThisParameters["main_model_part_name"].GetString() << "' does not exist." << std::endl;

        KRATOS_ERROR_IF_NOT( rModel.HasModelPart( mThisParameters["embedded_model_part_name"].GetString()) )
            << "MapNurbsVolumeResultsToEmbeddedGeometryProcess: Model Part '" <<  mThisParameters["embedded_model_part_name"].GetString() << "' does not exist." << std::endl;

        ModelPart& main_model_part = mrModel.GetModelPart(mThisParameters["main_model_part_name"].GetString());
        KRATOS_ERROR_IF_NOT( main_model_part.HasGeometry(mThisParameters["nurbs_volume_name"].GetString()) )
            << "MapNurbsVolumeResultsToEmbeddedGeometryProcess: Model Part '" <<  mThisParameters["main_model_part_name"].GetString() << "' does not have Geometry: '"
                << mThisParameters["nurbs_volume_name"].GetString() << "'. " << std::endl;
        ModelPart::GeometryType::Pointer p_geometry = main_model_part.pGetGeometry(mThisParameters["nurbs_volume_name"].GetString());

        KRATOS_ERROR_IF( p_geometry->GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Nurbs_Volume)
            << "MapNurbsVolumeResultsToEmbeddedGeometryProcess: Geometry: '" <<  mThisParameters["nurbs_volume_name"].GetInt() << "' is no 'Kratos_Nurbs_Volume-Geometry'." << std::endl;
   }

    void MapNurbsVolumeResultsToEmbeddedGeometryProcess::MapNodalValues(const Variable<array_1d<double,3>>& rVariable){

        // Get Model Parts
        ModelPart& main_model_part = mrModel.GetModelPart(mThisParameters["main_model_part_name"].GetString());
        ModelPart& embedded_model_part = mrModel.GetModelPart(mThisParameters["embedded_model_part_name"].GetString());

        // Get Nurbs Volume Geometry
        GeometryPointerType p_geometry = main_model_part.pGetGeometry(mThisParameters["nurbs_volume_name"].GetString());

        const SizeType number_nodes_embedded = embedded_model_part.NumberOfNodes();
        IntegrationPointsArrayType integration_points(number_nodes_embedded);

        const CoordinatesArrayType lower_point = p_geometry->begin()->GetInitialPosition();
        const CoordinatesArrayType upper_point = (p_geometry->end()-1)->GetInitialPosition();

        const auto node_itr_begin = embedded_model_part.NodesBegin();
        IndexPartition<IndexType>(embedded_model_part.NumberOfNodes()).for_each([&](IndexType i) {
            auto node_itr = node_itr_begin + i;
            // Map point into parameter space
            CoordinatesArrayType local_point;
            local_point[0] = (node_itr->X() - lower_point[0]) / std::abs( lower_point[0] - upper_point[0]);
            local_point[1] = (node_itr->Y() - lower_point[1]) / std::abs( lower_point[1] - upper_point[1]);
            local_point[2] = (node_itr->Z() - lower_point[2]) / std::abs( lower_point[2] - upper_point[2]);
            integration_points[i] = IntegrationPoint<3>(local_point, 0.0);
        });

        IntegrationInfo integration_info = p_geometry->GetDefaultIntegrationInfo();
        GeometriesArrayType geometry_list;
        p_geometry->CreateQuadraturePointGeometries(geometry_list, 1, integration_points, integration_info);

        IndexPartition<IndexType>(embedded_model_part.NumberOfNodes()).for_each([&](IndexType i) {
            auto node_itr = node_itr_begin + i;
            auto& quadrature_point = geometry_list[i];
            Vector N = row(quadrature_point.ShapeFunctionsValues(), 0);
            array_1d<double, 3> value = ZeroVector(3);
            for( IndexType j = 0; j < quadrature_point.size(); ++j){
                value += quadrature_point[j].FastGetSolutionStepValue(rVariable,0) * N(j);
            }
            // Write value onto embedded geometry
            node_itr->FastGetSolutionStepValue(rVariable, 0) = value;
        });
    }
} // End namespace Kratos
