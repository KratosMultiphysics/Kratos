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
#include "assign_integration_points_to_background_elements_process.h"
// #include "utilities/parallel_utilities.h"

namespace Kratos
{

    AssignIntegrationPointsToBackgroundElementsProcess::AssignIntegrationPointsToBackgroundElementsProcess(
        Model& rModel, Parameters ThisParameters) : mrModel(rModel), mThisParameters(ThisParameters)
    {
        mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
        KRATOS_ERROR_IF_NOT( rModel.HasModelPart( mThisParameters["main_model_part_name"].GetString()) )
            << "AssignIntegrationPointsToBackgroundElementsProcess: Model Part '" <<  mThisParameters["main_model_part_name"].GetString() << "' does not exist." << std::endl;

        KRATOS_ERROR_IF_NOT( rModel.HasModelPart( mThisParameters["embedded_model_part_name"].GetString()) )
            << "AssignIntegrationPointsToBackgroundElementsProcess: Model Part '" <<  mThisParameters["embedded_model_part_name"].GetString() << "' does not exist." << std::endl;

        ModelPart& main_model_part = mrModel.GetModelPart(mThisParameters["main_model_part_name"].GetString());
        KRATOS_ERROR_IF_NOT( main_model_part.HasGeometry(mThisParameters["nurbs_volume_name"].GetString()) )
            << "AssignIntegrationPointsToBackgroundElementsProcess: Model Part '" <<  mThisParameters["main_model_part_name"].GetString() << "' does not have Geometry: '"
                << mThisParameters["nurbs_volume_name"].GetString() << "'. " << std::endl;
        ModelPart::GeometryType::Pointer p_geometry = main_model_part.pGetGeometry(mThisParameters["nurbs_volume_name"].GetString());

        KRATOS_ERROR_IF( p_geometry->GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Nurbs_Volume)
            << "AssignIntegrationPointsToBackgroundElementsProcess: Geometry: '" <<  mThisParameters["nurbs_volume_name"].GetInt() << "' is no 'Kratos_Nurbs_Volume-Geometry'." << std::endl;
   }

    void AssignIntegrationPointsToBackgroundElementsProcess::AssignIntegrationPoints(){

        // Get Model Parts
        ModelPart& main_model_part = mrModel.GetModelPart(mThisParameters["main_model_part_name"].GetString());
        ModelPart& embedded_model_part = mrModel.GetModelPart(mThisParameters["embedded_model_part_name"].GetString());

        // Get Nurbs Volume Geometry
        GeometryPointerType p_geometry = main_model_part.pGetGeometry(mThisParameters["nurbs_volume_name"].GetString());

        const CoordinatesArrayType lower_point = p_geometry->begin()->GetInitialPosition();
        const CoordinatesArrayType upper_point = (p_geometry->end()-1)->GetInitialPosition();

        const auto element_itr_begin = embedded_model_part.ElementsBegin();
        std::cout << "Number of Elements: " << embedded_model_part.NumberOfElements() << std::endl;
        // Get element type
        IntegrationPointsArrayType integration_points(embedded_model_part.NumberOfElements());

        for( IndexType i = 0; i < embedded_model_part.NumberOfElements(); ++i){
            auto element_itr = element_itr_begin + i;
            auto& tmp_integration_points = element_itr->GetGeometry().IntegrationPoints();
            // Check if size is one!
            auto ip_physical_space = element_itr->GetGeometry().Center();

            // Map point into parameter space
            CoordinatesArrayType local_point;
            local_point[0] = (ip_physical_space[0] - lower_point[0]) / std::abs( lower_point[0] - upper_point[0]);
            local_point[1] = (ip_physical_space[1] - lower_point[1]) / std::abs( lower_point[1] - upper_point[1]);
            local_point[2] = (ip_physical_space[2] - lower_point[2]) / std::abs( lower_point[2] - upper_point[2]);

            integration_points[i] = IntegrationPoint<3>(local_point, tmp_integration_points[0].Weight());
        }

        IntegrationInfo integration_info = p_geometry->GetDefaultIntegrationInfo();
        GeometriesArrayType geometry_list;
        p_geometry->CreateQuadraturePointGeometries(geometry_list, 2, integration_points, integration_info);


        embedded_model_part.AddProperties(main_model_part.pGetProperties(0));

        auto& sub_model_part = embedded_model_part.GetSubModelPart("FormFinding");
        sub_model_part.AddProperties(main_model_part.pGetProperties(0));

        auto p_properties = main_model_part.pGetProperties(0);
        std::cout << "main_model_part.pGetProperties(0) "  << *p_properties << std::endl;


        const auto element_ptr_begin = embedded_model_part.Elements().ptr_begin();
        const auto& r_proces_info = main_model_part.GetProcessInfo();

        for( IndexType i = 0; i < embedded_model_part.NumberOfElements(); ++i){
            auto element_ptr = element_ptr_begin + i;
            IndexType current_id = (*element_ptr)->Id();

            auto& p_quadrature_point = geometry_list(i);
            //create the new element
            std::string element_name = "UpdatedLagrangianElement3D8N";
            ElementType const& r_clone_element = KratosComponents<ElementType>::Get(element_name);
            Element::Pointer p_element = r_clone_element.Create(current_id, p_quadrature_point, p_properties);
            p_element->Initialize(r_proces_info);
            (*element_ptr) = std::move(p_element);
        }
        auto points = embedded_model_part.ElementsBegin()->GetGeometry().IntegrationPoints();
        int id = points.size();
        std::cout << "embedded_model_part.: " << embedded_model_part.NumberOfElements() << std::endl;
        // const SizeType number_nodes_embedded = embedded_model_part.NumberOfNodes();
        // IntegrationPointsArrayType integration_points(number_nodes_embedded);

        // const CoordinatesArrayType lower_point = p_geometry->begin()->GetInitialPosition();
        // const CoordinatesArrayType upper_point = (p_geometry->end()-1)->GetInitialPosition();

        // const auto node_itr_begin = embedded_model_part.NodesBegin();
        // IndexPartition<IndexType>(embedded_model_part.NumberOfNodes()).for_each([&](IndexType i) {
        //     auto node_itr = node_itr_begin + i;
        //     // Map point into parameter space
        //     CoordinatesArrayType local_point;
        //     local_point[0] = (node_itr->X() - lower_point[0]) / std::abs( lower_point[0] - upper_point[0]);
        //     local_point[1] = (node_itr->Y() - lower_point[1]) / std::abs( lower_point[1] - upper_point[1]);
        //     local_point[2] = (node_itr->Z() - lower_point[2]) / std::abs( lower_point[2] - upper_point[2]);
        //     integration_points[i] = IntegrationPoint<3>(local_point, 0.0);
        // });


        // IndexPartition<IndexType>(embedded_model_part.NumberOfNodes()).for_each([&](IndexType i) {
        //     auto node_itr = node_itr_begin + i;
        //     auto& quadrature_point = geometry_list[i];
        //     Vector N = row(quadrature_point.ShapeFunctionsValues(), 0);
        //     array_1d<double, 3> value = ZeroVector(3);
        //     for( IndexType j = 0; j < quadrature_point.size(); ++j){
        //         value += quadrature_point[j].FastGetSolutionStepValue(rVariable,0) * N(j);
        //     }
        //     // Write value onto embedded geometry
        //     node_itr->FastGetSolutionStepValue(rVariable, 0) = value;
        // });
    }
} // End namespace Kratos
