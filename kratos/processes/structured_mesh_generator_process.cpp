//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// System includes
#include <vector>
#include <numeric>

// External includes

// Project includes
#include "processes/structured_mesh_generator_process.h"
#include "geometries/geometry.h"
#include "geometries/point.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "includes/checks.h"
#include "processes/skin_detection_process.h"

namespace Kratos
{
StructuredMeshGeneratorProcess::StructuredMeshGeneratorProcess(const GeometryType& rGeometry, ModelPart& rOutputModelPart, Parameters TheParameters)
    : Process()
    , mrGeometry(rGeometry)
    , mrOutputModelPart(rOutputModelPart) {

    TheParameters["element_name"]; // Should be given by caller! if not thorws an error

    ValidateTheDefaultParameters(TheParameters);

    mStartNodeId = TheParameters["start_node_id"].GetInt();
    mStartElementId = TheParameters["start_element_id"].GetInt();
    mStartConditionId = TheParameters["start_condition_id"].GetInt();

    mNumberOfDivisions[0] = TheParameters["number_of_divisions_X"].GetInt();
    mNumberOfDivisions[1] = TheParameters["number_of_divisions_Y"].GetInt();
    mNumberOfDivisions[2] = TheParameters["number_of_divisions_Z"].GetInt();

    mElementPropertiesId = TheParameters["elements_properties_id"].GetInt();
    mConditiongPropertiesId = TheParameters["conditions_properties_id"].GetInt();
    mElementName = TheParameters["element_name"].GetString();
    mConditionName = TheParameters["condition_name"].GetString();
    mCreateBodySubModelPart = TheParameters["create_body_sub_model_part"].GetBool();
    mCreateSkinSubModelPart = TheParameters["create_skin_sub_model_part"].GetBool();
    if (mCreateBodySubModelPart) {
        mBodySubModelPartName = TheParameters["body_sub_model_part_name"].GetString();
    }
    if (mCreateSkinSubModelPart) {
        mSkinSubModelPartName = TheParameters["skin_sub_model_part_name"].GetString();
    }

    Check();
}

StructuredMeshGeneratorProcess::~StructuredMeshGeneratorProcess()
{

}

void StructuredMeshGeneratorProcess::Execute()
{
    // Generate 3D mesh
    if (mrGeometry.LocalSpaceDimension() == 2) {
        Generate2DMesh();
    } else if (mrGeometry.LocalSpaceDimension() == 3) {
        Generate3DMesh();
    } else {
        KRATOS_ERROR << "Not supported geometry is given" << std::endl;
    }

    // Generate body model part if required
    if (mCreateBodySubModelPart) {
        // Create the body model part
        auto& r_body_sub_model_part = mrOutputModelPart.CreateSubModelPart(mBodySubModelPartName);
        // Set the nodal ids array
        std::vector<ModelPart::IndexType> ids_nodes(mrOutputModelPart.NumberOfNodes());
        std::iota(ids_nodes.begin(), ids_nodes.end(), mStartNodeId);
        r_body_sub_model_part.AddNodes(ids_nodes);
        // Set the element ids array
        std::vector<ModelPart::IndexType> ids_elems(mrOutputModelPart.NumberOfElements());
        std::iota(ids_elems.begin(), ids_elems.end(), mStartElementId);
        r_body_sub_model_part.AddElements(ids_elems);
    }

    // Generate skin if required
    if (mCreateSkinSubModelPart) {
        const Parameters skin_parameters = Parameters(R"(
        {
            "name_auxiliar_model_part" : "",
            "name_auxiliar_condition"  : "Condition"
        })" );
        skin_parameters["name_auxiliar_model_part"].SetString(mSkinSubModelPartName);

        if (mConditionName != "PLEASE SPECIFY IT") {
            skin_parameters["name_auxiliar_condition"].SetString(mConditionName);
        } else {
            KRATOS_WARNING("StructuredMeshGeneratorProcess") << "Condition name not specified. Default geometrical conditions will be generated" << std::endl;
        }
        if (mrGeometry.LocalSpaceDimension() == 2) {
            SkinDetectionProcess<2>(mrOutputModelPart, skin_parameters).Execute();
        } else {
            SkinDetectionProcess<3>(mrOutputModelPart, skin_parameters).Execute();
        }
    }

}

void StructuredMeshGeneratorProcess::ValidateTheDefaultParameters(Parameters TheParameters)
{
    if(TheParameters.Has("number_of_divisions")){
        if(TheParameters["number_of_divisions"].IsInt() && (!TheParameters.Has("number_of_divisions_X") && !TheParameters.Has("number_of_divisions_Y") && !TheParameters.Has("number_of_divisions_Z"))) 
        {
            int ndivisions = TheParameters["number_of_divisions"].GetInt();
            TheParameters.RemoveValue("number_of_divisions");
            TheParameters.AddInt("number_of_divisions_X", ndivisions);
            TheParameters.AddInt("number_of_divisions_Y", ndivisions);
            TheParameters.AddInt("number_of_divisions_Z", ndivisions);
        } else {
            KRATOS_THROW_ERROR(std::invalid_argument, "Please specify number_of_divisions as an int or the component of each direction","")
        }
    }
    TheParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

const Parameters StructuredMeshGeneratorProcess::GetDefaultParameters() const
{
    const Parameters default_parameters(R"(
    {
        "create_skin_sub_model_part" : true,
        "create_body_sub_model_part" : false,
        "skin_sub_model_part_name"   : "Skin",
        "body_sub_model_part_name"   : "Body",
        "start_node_id"              : 1,
        "start_element_id"           : 1,
        "start_condition_id"         : 1,
        "number_of_divisions_X"      : 1,
        "number_of_divisions_Y"      : 1,
        "number_of_divisions_Z"      : 1,
        "elements_properties_id"     : 0,
        "conditions_properties_id"   : 0,
        "element_name"               : "PLEASE SPECIFY IT",
        "condition_name"             : "PLEASE SPECIFY IT"
    }  )");
    return default_parameters;
}

std::string StructuredMeshGeneratorProcess::Info() const
{
    return "StructuredMeshGeneratorProcess";
}

void StructuredMeshGeneratorProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

void StructuredMeshGeneratorProcess::PrintData(std::ostream& rOStream) const
{

}

void StructuredMeshGeneratorProcess::Generate2DMesh()
{
    Point min_point(1.00, 1.00, 1.00);
    Point max_point(-1.00, -1.00, -1.00);
    GetLocalCoordinatesRange(min_point, max_point);

    GenerateNodes2D(min_point, max_point);

    GenerateTriangularElements();
}

void StructuredMeshGeneratorProcess::Generate3DMesh()
{
    Point min_point(1.00, 1.00, 1.00);
    Point max_point(-1.00, -1.00, -1.00);
    GetLocalCoordinatesRange(min_point, max_point);

    GenerateNodes3D(min_point, max_point);

    GenerateTetrahedraElements();
}

void StructuredMeshGeneratorProcess::GenerateNodes2D(Point const& rMinPoint, Point const& rMaxPoint)
{
    GeometryType::CoordinatesArrayType local_element_size = rMaxPoint - rMinPoint;
    local_element_size[0] /= mNumberOfDivisions[0];
    local_element_size[1] /= mNumberOfDivisions[1];
    //const std::size_t local_space_dimension = mrGeometry.LocalSpaceDimension();
    Point local_coordinates = rMinPoint;
    auto global_coordinates = Point{ZeroVector(3)};
    std::size_t node_id = mStartNodeId;

    for (std::size_t j = 0; j <= mNumberOfDivisions[1]; j++) {
        for (std::size_t i = 0; i <= mNumberOfDivisions[0]; i++) {
            local_coordinates[0] = rMinPoint[0] + (i * local_element_size[0]);
            local_coordinates[1] = rMinPoint[1] + (j * local_element_size[1]);
            mrGeometry.GlobalCoordinates(global_coordinates, local_coordinates);
            mrOutputModelPart.CreateNewNode(node_id++, global_coordinates[0],
                                                            global_coordinates[1],
                                                            global_coordinates[2]);
        }
    }
}

void StructuredMeshGeneratorProcess::GenerateNodes3D(Point const& rMinPoint, Point const& rMaxPoint)
{
    GeometryType::CoordinatesArrayType local_element_size = rMaxPoint - rMinPoint;
    local_element_size[0] /= mNumberOfDivisions[0];
    local_element_size[1] /= mNumberOfDivisions[1];
    local_element_size[2] /= mNumberOfDivisions[2];
    Point local_coordinates = rMinPoint;
    auto global_coordinates = Point{ZeroVector(3)};
    std::size_t node_id = mStartNodeId;

    for (std::size_t k = 0; k <= mNumberOfDivisions[2]; k++) {
        for (std::size_t j = 0; j <= mNumberOfDivisions[1]; j++) {
            for (std::size_t i = 0; i <= mNumberOfDivisions[0]; i++) {
                local_coordinates[0] = rMinPoint[0] + (i * local_element_size[0]);
                local_coordinates[1] = rMinPoint[1] + (j * local_element_size[1]);
                local_coordinates[2] = rMinPoint[2] + (k * local_element_size[2]);
                mrGeometry.GlobalCoordinates(global_coordinates, local_coordinates);
                mrOutputModelPart.CreateNewNode(node_id++, global_coordinates[0],
                                                                global_coordinates[1],
                                                                global_coordinates[2]);
            }
        }
    }
}

void StructuredMeshGeneratorProcess::GenerateTriangularElements()
{
    std::size_t element_id = mStartElementId;

    Properties::Pointer p_properties = mrOutputModelPart.CreateNewProperties(mElementPropertiesId);
    std::vector<ModelPart::IndexType> element_connectivity(3);

    for (std::size_t j = 0; j < mNumberOfDivisions[1]; j++) {
        for (std::size_t i = 0; i < mNumberOfDivisions[0]; i++) {
            element_connectivity = { GetNodeId(i,j,0), GetNodeId(i + 1,j + 1,0), GetNodeId(i + 1,j,0) };
            mrOutputModelPart.CreateNewElement(mElementName, element_id++, element_connectivity, p_properties);

            element_connectivity = { GetNodeId(i,j,0), GetNodeId(i,j + 1,0), GetNodeId(i + 1,j + 1,0) };
            mrOutputModelPart.CreateNewElement(mElementName, element_id++, element_connectivity, p_properties);
        }
    }
}

void StructuredMeshGeneratorProcess::GenerateTetrahedraElements()
{
    Properties::Pointer p_properties = mrOutputModelPart.CreateNewProperties(mElementPropertiesId);

    for (std::size_t k = 0; k < mNumberOfDivisions[2]; k++) {
        for (std::size_t j = 0; j < mNumberOfDivisions[1]; j++) {
            for (std::size_t i = 0; i < mNumberOfDivisions[0]; i++) {
                CreateCellTetrahedra(i, j, k, p_properties);
            }
        }
    }
}

void  StructuredMeshGeneratorProcess::CreateCellTetrahedra(std::size_t I, std::size_t J, std::size_t K, Properties::Pointer pProperties)
{
    using point_in_cell_position_type = std::array<std::size_t, 3>;
    using tetrahedra_connectivity_in_cell_type = std::array<std::size_t, 4>;
    constexpr std::size_t number_of_cases = 6;
    constexpr point_in_cell_position_type cell_points[8] = { {{ 0,0,0 }},{{ 1,0,0 }},{{ 1,1,0 }},{{ 0,1,0 }},
                                                              {{ 0,0,1 }},{{ 1,0,1 }},{{ 1,1,1 }},{{ 0,1,1 }}  };

    constexpr tetrahedra_connectivity_in_cell_type connectivity_cases[number_of_cases] = { {{ 0,3,6,2 }},{{ 3,6,7,0 }},{{ 4,7,6,0 }},
                                                                                            {{ 0,4,5,6 }},{{ 0,1,2,6 }},{{ 1,5,6,0 }} };
    std::vector<ModelPart::IndexType> element_connectivity(4);

    for (std::size_t i_case = 0; i_case < number_of_cases; i_case++) {
        auto connectivity = connectivity_cases[i_case];
        for (std::size_t i_position = 0; i_position < 4; i_position++)
        {
            auto& cell_point = cell_points[connectivity[i_position]];
            element_connectivity[i_position] = GetNodeId(I + cell_point[0], J + cell_point[1], K + cell_point[2]);
        }
        mrOutputModelPart.CreateNewElement(mElementName, mStartElementId++, element_connectivity, pProperties);
    }
}

std::size_t StructuredMeshGeneratorProcess::GetNodeId(std::size_t I, std::size_t J, std::size_t K)
{
    return mStartNodeId + (K * (mNumberOfDivisions[1] + 1) * (mNumberOfDivisions[0] + 1)) + (J * (mNumberOfDivisions[0] + 1)) + I;
}

void StructuredMeshGeneratorProcess::GetLocalCoordinatesRange(Point& rMinPoint, Point& rMaxPoint)
{
    const std::size_t local_space_dimension = mrGeometry.LocalSpaceDimension();
    Matrix geometry_points_local_coordinates;
    mrGeometry.PointsLocalCoordinates(geometry_points_local_coordinates);

    const std::size_t number_of_points = mrGeometry.size();

    for (std::size_t i_point = 0; i_point < number_of_points; i_point++) {
        for (std::size_t i_dimension = 0; i_dimension < local_space_dimension; i_dimension++) {
            rMinPoint[i_dimension] = std::min(rMinPoint[i_dimension], geometry_points_local_coordinates(i_point, i_dimension));
            rMaxPoint[i_dimension] = std::max(rMaxPoint[i_dimension], geometry_points_local_coordinates(i_point, i_dimension));
        }
    }
}

int StructuredMeshGeneratorProcess::Check()
{
    KRATOS_TRY

    KRATOS_CHECK(CheckDomainGeometry());
    KRATOS_CHECK(KratosComponents<Element>::Has(mElementName));

    if ((mrGeometry.GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4) &&
        (mrGeometry.GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Hexahedra3D8))
        KRATOS_ERROR << "An unsupported geometry was given. Only Quadrilateral2D4 and Hexahedra3D8 are supported and given geometry is : " << mrGeometry << std::endl;

    KRATOS_CHECK_NOT_EQUAL(mNumberOfDivisions[0], 0);
    KRATOS_CHECK_NOT_EQUAL(mNumberOfDivisions[1], 0);
    KRATOS_CHECK_NOT_EQUAL(mNumberOfDivisions[2], 0);

    return 0;

    KRATOS_CATCH("")
}

bool StructuredMeshGeneratorProcess::CheckDomainGeometry()
{
    if (mrGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4) {
        return CheckDomainGeometryConnectivityForQuadrilateral2D4();
    } else if (mrGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8) {
        return CheckDomainGeometryConnectivityForHexahedra3D8();
    }

    return true;
}

bool StructuredMeshGeneratorProcess::CheckDomainGeometryConnectivityForQuadrilateral2D4()
{
    using triangle_connectivity_in_cell_type = std::array<std::size_t, 3>;
    constexpr std::size_t number_of_cases = 2;

    std::vector<std::array<double, 3> > cell_points;
    for (std::size_t i = 0; i < 4; ++i){
        std::array<double, 3> coordinates{{mrGeometry[i][0], mrGeometry[i][1], mrGeometry[i][2]}};
        cell_points.push_back(coordinates);
    }

    constexpr triangle_connectivity_in_cell_type connectivity_cases[number_of_cases] = { {{ 0,2,1 }},{{ 0,3,2 }} };

    std::vector<Point::Pointer> my_points(3);
    double min_area = 1.0;

    for (std::size_t i_case = 0; i_case < number_of_cases; i_case++) {
        auto connectivity = connectivity_cases[i_case];
        for (std::size_t i_position = 0; i_position < 3; i_position++)
        {
            auto& cell_point = cell_points[connectivity[i_position]];
            Point::Pointer pPi(new Point(cell_point[0], cell_point[1], cell_point[2]));
            my_points[i_position] = pPi;
        }

        Triangle2D3<Point > trial_triangle(my_points[0], my_points[1], my_points[2]);
        min_area = std::min(min_area, trial_triangle.DomainSize());
    }

    bool all_triangles_have_positive_area = min_area > 0.0;
    return all_triangles_have_positive_area;
}


bool StructuredMeshGeneratorProcess::CheckDomainGeometryConnectivityForHexahedra3D8()
{
    using tetrahedra_connectivity_in_cell_type = std::array<std::size_t, 4>;
    constexpr std::size_t number_of_cases = 6;

    std::vector<std::array<double, 3> > cell_points;
    for (std::size_t i = 0; i < 8; ++i){
        std::array<double, 3> coordinates{{mrGeometry[i][0], mrGeometry[i][1], mrGeometry[i][2]}};
        cell_points.push_back(coordinates);
    }

    constexpr tetrahedra_connectivity_in_cell_type connectivity_cases[number_of_cases] = { {{ 0,3,6,2 }},{{ 3,6,7,0 }},{{ 4,7,6,0 }},
                                                                                            {{ 0,4,5,6 }},{{ 0,1,2,6 }},{{ 1,5,6,0 }} };
    std::vector<Point::Pointer> my_points(4);
    double min_volume = 1.0;

    for (std::size_t i_case = 0; i_case < number_of_cases; i_case++) {
        auto connectivity = connectivity_cases[i_case];
        for (std::size_t i_position = 0; i_position < 4; i_position++)
        {
            auto& cell_point = cell_points[connectivity[i_position]];
            Point::Pointer pPi(new Point(cell_point[0], cell_point[1], cell_point[2]));
            my_points[i_position] = pPi;
        }

        Tetrahedra3D4<Point > trial_tetra(my_points[0], my_points[1], my_points[2], my_points[3]);
        min_volume = std::min(min_volume, trial_tetra.DomainSize());
    }

    bool all_tetrahedra_have_positive_volume = min_volume > 0.0;
    return all_tetrahedra_have_positive_volume;
}

}  // namespace Kratos.
