//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


// Project includes
#include "mapping_geometries_modeler.h"


namespace Kratos
{
    ///@name Stages
    ///@{

    void MappingGeometriesModeler::SetupGeometryModel()
    {
        KRATOS_ERROR_IF_NOT(mParameters.Has("origin_model_part_name"))
            << "Missing \"origin_model_part_name\" in MappingGeometriesModeler Parameters." << std::endl;
        const std::string origin_model_part_name = mParameters["origin_model_part_name"].GetString();
        ModelPart& origin_model_part = mpModel->HasModelPart(origin_model_part_name)
            ? mpModel->GetModelPart(origin_model_part_name)
            : mpModel->CreateModelPart(origin_model_part_name);

        KRATOS_ERROR_IF_NOT(mParameters.Has("destination_model_part_name"))
            << "Missing \"destination_model_part_name\" in MappingGeometriesModeler Parameters." << std::endl;
        const std::string destination_model_part_name = mParameters["destination_model_part_name"].GetString();
        ModelPart& destination_model_part = mpModel->HasModelPart(destination_model_part_name)
            ? mpModel->GetModelPart(destination_model_part_name)
            : mpModel->CreateModelPart(destination_model_part_name);

        ModelPart& coupling_model_part = mpModel->HasModelPart("coupling")
            ? mpModel->GetModelPart("coupling")
            : mpModel->CreateModelPart("coupling");

        ModelPart& model_part_coupling_quadrature_points = mpModel->HasModelPart("coupling_quadrature_points")
            ? mpModel->GetModelPart("coupling_quadrature_points")
            : mpModel->CreateModelPart("coupling_quadrature_points");


        MappingIntersectionUtilities::FindIntersection1DGeometries2D(
            origin_model_part, destination_model_part, coupling_model_part, 1e-6);
        MappingIntersectionUtilities::CreateQuadraturePointsCoupling1DGeometries2D(
            coupling_model_part, model_part_coupling_quadrature_points, 1e-6);
    }

    void MappingGeometriesModeler::PrepareGeometryModel()
    {
        // TODO add other functions

    }

    ///@}
}
