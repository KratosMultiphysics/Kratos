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
#include "geometries/nurbs_volume_geometry.h"

namespace Kratos
{

    AssignIntegrationPointsToBackgroundElementsProcess::AssignIntegrationPointsToBackgroundElementsProcess(
        Model& rModel, Parameters ThisParameters) : mrModel(rModel), mThisParameters(ThisParameters)
    {
        mIsAssignedFlag = false;

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

    void AssignIntegrationPointsToBackgroundElementsProcess::ExecuteBeforeOutputStep(){

        if( !mIsAssignedFlag ){
            mIsAssignedFlag = true;
            // Get Model Parts
            ModelPart& main_model_part = mrModel.GetModelPart(mThisParameters["main_model_part_name"].GetString());
            ModelPart& embedded_model_part = mrModel.GetModelPart(mThisParameters["embedded_model_part_name"].GetString());

            // Get Nurbs Volume Geometry
            GeometryPointerType p_geometry = main_model_part.pGetGeometry(mThisParameters["nurbs_volume_name"].GetString());

            // Create quadrature points in parameter space of nurbs volume
            IntegrationPointsArrayType integration_points(embedded_model_part.NumberOfElements());
            const auto element_itr_begin = embedded_model_part.ElementsBegin();
            for( IndexType i = 0; i < embedded_model_part.NumberOfElements(); ++i){
                auto element_itr = element_itr_begin + i;
                auto& tmp_integration_points = element_itr->GetGeometry().IntegrationPoints();
                KRATOS_ERROR_IF_NOT( tmp_integration_points.size() == 1 )
                    << "AssignIntegrationPointsToBackgroundElementsProcess: " << "Number of integration points does not equal 1." << std::endl;
                // Get location of integration points in physical space
                auto ip_physical_space = element_itr->GetGeometry().Center();

                // Map point into parameter space
                CoordinatesArrayType local_point;
                p_geometry->ProjectionPointGlobalToLocalSpace(ip_physical_space, local_point);

                integration_points[i] = IntegrationPoint<3>(local_point, tmp_integration_points[0].Weight());
            }
            IntegrationInfo integration_info = p_geometry->GetDefaultIntegrationInfo();
            GeometriesArrayType geometry_list;
            // Make sure one quadrature point geometry is created for each integration point.
            integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::CUSTOM);
            integration_info.SetQuadratureMethod(1, IntegrationInfo::QuadratureMethod::CUSTOM);
            integration_info.SetQuadratureMethod(2, IntegrationInfo::QuadratureMethod::CUSTOM);
            p_geometry->CreateQuadraturePointGeometries(geometry_list, 2, integration_points, integration_info);

            // Get properties from main model part
            auto p_properties = main_model_part.pGetProperties(1);

            const auto p_background_element = main_model_part.ElementsBegin();
            const auto& r_proces_info = main_model_part.GetProcessInfo();
            const auto element_ptr_begin = embedded_model_part.Elements().ptr_begin();
            for( IndexType i = 0; i < embedded_model_part.NumberOfElements(); ++i){
                auto element_ptr = element_ptr_begin + i;
                // Get id and geometry
                IndexType current_id = (*element_ptr)->Id();
                auto p_quadrature_point = geometry_list(i);
                // Create and initialize new element - same element type as background element
                Element::Pointer p_element = p_background_element->Create(current_id, p_quadrature_point, p_properties);
                p_element->Initialize(r_proces_info);
                // Point to newly created element
                (*element_ptr) = std::move(p_element);
            }
        }
    }
} // End namespace Kratos
