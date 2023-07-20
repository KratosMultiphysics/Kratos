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

        // We get the list of nodal variables.
        for (const std::string& r_variable_name : mThisParameters["nodal_results"].GetStringArray()){
            if (KratosComponents<Variable<double>>::Has(r_variable_name)){
                mDoubleVariableNode.push_back(&(KratosComponents< Variable<double>>::Get(r_variable_name)));
            } else if (KratosComponents<Variable<array_1d<double, 3> >>::Has(r_variable_name)) {
                mArrayVariableNode.push_back(&(KratosComponents<Variable<array_1d<double, 3>>>::Get(r_variable_name)));
            } else {
                KRATOS_ERROR << "Only double, and array variables are allowed in the variables list:"
                    << " nodal_results.\n";
            }
        }

        // We get the list of gauss point variables.
        for (const std::string& r_variable_name : mThisParameters["gauss_point_results"].GetStringArray()){
            if (KratosComponents<Variable<double>>::Has(r_variable_name)){
                mDoubleVariableGauss.push_back(&(KratosComponents< Variable<double>>::Get(r_variable_name)));
            } else if (KratosComponents<Variable<array_1d<double, 3> >>::Has(r_variable_name)) {
                mArrayVariableGauss.push_back(&(KratosComponents<Variable<array_1d<double, 3>>>::Get(r_variable_name)));
            } else if (KratosComponents<Variable<Vector>>::Has(r_variable_name)) {
                mVectorVariableGauss.push_back(&(KratosComponents<Variable<Vector>>::Get(r_variable_name)));
            } else if (KratosComponents<Variable<Matrix>>::Has(r_variable_name)) {
                mMatrixVariableGauss.push_back(&(KratosComponents<Variable<Matrix>>::Get(r_variable_name)));
            } else {
                KRATOS_ERROR << "Only double, array, vector and matrix variables are allowed in the variables list:"
                    << " gauss_point_results.\n";
            }
        }

   }

    void MapNurbsVolumeResultsToEmbeddedGeometryProcess::MapVariables(){

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
            integration_points[i] = IntegrationPoint<3>(local_point, 1.0);
        });

        IntegrationInfo integration_info = p_geometry->GetDefaultIntegrationInfo();
        GeometriesArrayType geometry_list;
        // Make sure one qudrature point geometry is created for each integration point.
        integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::EXTENDED_GAUSS);
        integration_info.SetQuadratureMethod(1, IntegrationInfo::QuadratureMethod::EXTENDED_GAUSS);
        integration_info.SetQuadratureMethod(2, IntegrationInfo::QuadratureMethod::EXTENDED_GAUSS);
        p_geometry->CreateQuadraturePointGeometries(geometry_list, 2, integration_points, integration_info);

        const auto p_properties = main_model_part.pGetProperties(1);
        const auto p_background_element = main_model_part.ElementsBegin();
        const auto& r_process_info = main_model_part.GetProcessInfo();

        IndexPartition<IndexType>(embedded_model_part.NumberOfNodes()).for_each([&](IndexType i) {
            auto node_itr = node_itr_begin + i;
            auto& r_quadrature_point = geometry_list[i];

            //// Map nodal results
            // Scalar
            const Vector N = row(r_quadrature_point.ShapeFunctionsValues(), 0);
            for ( const auto p_var : mDoubleVariableNode) {
                double value = 0.0;
                for( IndexType j = 0; j < r_quadrature_point.size(); ++j){
                    value += r_quadrature_point[j].FastGetSolutionStepValue(*p_var) * N(j);
                }
                // Write value onto embedded geometry
                node_itr->FastGetSolutionStepValue(*p_var) = value;
            }

            // Array
            for ( const auto p_var : mArrayVariableNode) {
                array_1d<double, 3> value = ZeroVector(3);
                for( IndexType j = 0; j < r_quadrature_point.size(); ++j){
                    value += r_quadrature_point[j].FastGetSolutionStepValue(*p_var) * N(j);
                }
                // Write value onto embedded geometry
                node_itr->FastGetSolutionStepValue(*p_var) = value;
            }


            //// Map gauss point results
            auto p_quadrature_point = geometry_list(i);
            Element::Pointer p_element = p_background_element->Create(1, p_quadrature_point, p_properties);
            p_element->Initialize(r_process_info);

            // Scalar
            for ( const auto p_var : mDoubleVariableGauss) {
                std::vector<double> aux_result(1);
                p_element->CalculateOnIntegrationPoints(*p_var, aux_result, r_process_info);
                node_itr->SetValue(*p_var, aux_result[0]);
            }
            // Array
            for ( const auto p_var : mArrayVariableGauss) {
                std::vector<array_1d<double, 3>> aux_result(1);
                p_element->CalculateOnIntegrationPoints(*p_var, aux_result, r_process_info);
                node_itr->SetValue(*p_var, aux_result[0]);
            }
            // Vector
            for ( const auto p_var : mVectorVariableGauss) {
                std::vector<Vector> aux_result(1);
                p_element->CalculateOnIntegrationPoints(*p_var, aux_result, r_process_info);
                node_itr->SetValue(*p_var, aux_result[0]);
            }
            // Matrix
            for ( const auto p_var : mMatrixVariableGauss) {
                std::vector<Matrix> aux_result(1);
                p_element->CalculateOnIntegrationPoints(*p_var, aux_result, r_process_info);
                node_itr->SetValue(*p_var, aux_result[0]);
            }
        });
    }
} // End namespace Kratos
