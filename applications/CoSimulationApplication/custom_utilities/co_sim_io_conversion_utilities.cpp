//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes

// Project includes
#include "includes/parallel_environment.h"
#include "co_sim_io_conversion_utilities.h"

namespace Kratos {
namespace {

//  TODO refactor with switch?
const std::map<GeometryData::KratosGeometryType, CoSimIO::ElementType> elem_type_map {
    {GeometryData::KratosGeometryType::Kratos_Hexahedra3D20,    CoSimIO::ElementType::Hexahedra3D20},
    {GeometryData::KratosGeometryType::Kratos_Hexahedra3D27,    CoSimIO::ElementType::Hexahedra3D27},
    {GeometryData::KratosGeometryType::Kratos_Hexahedra3D8,     CoSimIO::ElementType::Hexahedra3D8},
    {GeometryData::KratosGeometryType::Kratos_Prism3D15,        CoSimIO::ElementType::Prism3D15},
    {GeometryData::KratosGeometryType::Kratos_Prism3D6,         CoSimIO::ElementType::Pyramid3D13},
    {GeometryData::KratosGeometryType::Kratos_Pyramid3D13,      CoSimIO::ElementType::Pyramid3D5},
    {GeometryData::KratosGeometryType::Kratos_Pyramid3D5,       CoSimIO::ElementType::Prism3D6},
    {GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4, CoSimIO::ElementType::Quadrilateral2D4},
    {GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8, CoSimIO::ElementType::Quadrilateral2D8},
    {GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9, CoSimIO::ElementType::Quadrilateral2D9},
    {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4, CoSimIO::ElementType::Quadrilateral3D4},
    {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8, CoSimIO::ElementType::Quadrilateral3D8},
    {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9, CoSimIO::ElementType::Quadrilateral3D9},
    {GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10,   CoSimIO::ElementType::Tetrahedra3D10},
    {GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4,    CoSimIO::ElementType::Tetrahedra3D4},
    {GeometryData::KratosGeometryType::Kratos_Triangle2D3,      CoSimIO::ElementType::Triangle2D3},
    {GeometryData::KratosGeometryType::Kratos_Triangle2D6,      CoSimIO::ElementType::Triangle2D6},
    {GeometryData::KratosGeometryType::Kratos_Triangle3D3,      CoSimIO::ElementType::Triangle3D3},
    {GeometryData::KratosGeometryType::Kratos_Triangle3D6,      CoSimIO::ElementType::Triangle3D6},
    {GeometryData::KratosGeometryType::Kratos_Line2D2,          CoSimIO::ElementType::Line2D2},
    {GeometryData::KratosGeometryType::Kratos_Line2D3,          CoSimIO::ElementType::Line2D3},
    {GeometryData::KratosGeometryType::Kratos_Line3D2,          CoSimIO::ElementType::Line3D2},
    {GeometryData::KratosGeometryType::Kratos_Line3D3,          CoSimIO::ElementType::Line3D3},
    {GeometryData::KratosGeometryType::Kratos_Point2D,          CoSimIO::ElementType::Point2D},
    {GeometryData::KratosGeometryType::Kratos_Point3D,          CoSimIO::ElementType::Point3D}
};

//  TODO refactor with switch?
const std::map<CoSimIO::ElementType, std::string> elem_name_map {
    {CoSimIO::ElementType::Hexahedra3D20, "Element3D20N"},
    {CoSimIO::ElementType::Hexahedra3D27, "Element3D27N"},
    {CoSimIO::ElementType::Hexahedra3D8, "Element3D8N"},
    {CoSimIO::ElementType::Prism3D15, "Element3D15N"},
    {CoSimIO::ElementType::Prism3D6, "Element3D6N"},
    {CoSimIO::ElementType::Quadrilateral2D4, "Element2D4N"},
    {CoSimIO::ElementType::Quadrilateral2D8, "Element2D8N"},
    {CoSimIO::ElementType::Quadrilateral2D9, "Element2D9N"},
    {CoSimIO::ElementType::Quadrilateral3D8, "Element3D8N"},
    {CoSimIO::ElementType::Tetrahedra3D10, "Element3D10N"},
    {CoSimIO::ElementType::Tetrahedra3D4, "Element3D4N"},
    {CoSimIO::ElementType::Triangle2D3, "Element2D3N"},
    {CoSimIO::ElementType::Triangle2D6, "Element2D6N"},
    {CoSimIO::ElementType::Triangle3D3, "Element3D3N"},
    {CoSimIO::ElementType::Line2D2, "Element2D2N"},
    {CoSimIO::ElementType::Line3D2, "Element3D2N"},
    {CoSimIO::ElementType::Point2D, "Element2D1N"},
    {CoSimIO::ElementType::Point3D, "Element3D1N"}
};

} // anonymous namespace

void CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(
    const CoSimIO::ModelPart& rCoSimIOModelPart,
    Kratos::ModelPart& rKratosModelPart,
    const DataCommunicator& rDataComm)
{
    KRATOS_TRY

    // fill ModelPart from received Mesh
    KRATOS_ERROR_IF(rKratosModelPart.NumberOfNodes() > 0) << "ModelPart is not empty, it has nodes!" << std::endl;
    KRATOS_ERROR_IF(rKratosModelPart.NumberOfProperties() > 0) << "ModelPart is not empty, it has properties!" << std::endl;
    KRATOS_ERROR_IF(rKratosModelPart.IsDistributed()) << "ModelPart cannot be distributed!" << std::endl;

    const bool is_distributed = rDataComm.IsDistributed();
    const int my_rank = rDataComm.Rank();

    if (is_distributed) {
        rKratosModelPart.AddNodalSolutionStepVariable(PARTITION_INDEX); // to be on the safe side
    }

    // check if nodes are ordered consecutively
    // this is unfortunately necessary for several resons, e.g. AddElements and the ParallelFillCommunicator
    std::size_t max_node_id = 0;
    for (const auto& r_node : rCoSimIOModelPart.Nodes()) {
        KRATOS_ERROR_IF(max_node_id >= static_cast<std::size_t>(r_node.Id())) << "The nodes must be consecutively ordered!" << std::endl;
        max_node_id = r_node.Id();
    }

    for (const auto& r_node : rCoSimIOModelPart.LocalNodes()) {
        auto p_node = rKratosModelPart.CreateNewNode(
            r_node.Id(),
            r_node.X(),
            r_node.Y(),
            r_node.Z()
        );

        if (is_distributed) {
            // alternatively use VariableUtils
            p_node->FastGetSolutionStepValue(PARTITION_INDEX) = my_rank;
        }
    };

    KRATOS_ERROR_IF(!is_distributed && rCoSimIOModelPart.GetPartitionModelParts().size()>0) << "Ghost entities exist in CoSimIO ModelPart in serial simulation!" << std::endl;

    for (const auto& r_partition_pair : rCoSimIOModelPart.GetPartitionModelParts()) {
        const int partition_index = r_partition_pair.first;
        const std::unique_ptr<CoSimIO::ModelPart>& rp_partition_model_part = r_partition_pair.second;

        max_node_id = 0;
        for (const auto& r_node : rp_partition_model_part->Nodes()) {
            KRATOS_ERROR_IF(max_node_id >= static_cast<std::size_t>(r_node.Id())) << "The nodes must be consecutively ordered!" << std::endl;
            max_node_id = r_node.Id();

            auto p_node = rKratosModelPart.CreateNewNode(
                r_node.Id(),
                r_node.X(),
                r_node.Y(),
                r_node.Z()
            );

            p_node->FastGetSolutionStepValue(PARTITION_INDEX) = partition_index;
        };
    }

    Properties::Pointer p_props;
    if (rCoSimIOModelPart.NumberOfElements() > 0) {
        p_props = rKratosModelPart.CreateNewProperties(0);
        rKratosModelPart.Elements().reserve(rCoSimIOModelPart.NumberOfElements());
    }

    std::vector<IndexType> conn;
    std::size_t max_elem_id = 0;
    for (auto elem_it=rCoSimIOModelPart.ElementsBegin(); elem_it!=rCoSimIOModelPart.ElementsEnd(); ++elem_it) {
        const auto& r_elem = **elem_it;

        KRATOS_ERROR_IF(max_elem_id >= static_cast<std::size_t>(r_elem.Id())) << "The elements must be consecutively ordered!" << std::endl;
        max_elem_id = r_elem.Id();

        if (conn.size() != r_elem.NumberOfNodes()) {
            conn.resize(r_elem.NumberOfNodes());
        };

        const auto nodes_begin = r_elem.NodesBegin();
        for (std::size_t i=0; i<r_elem.NumberOfNodes(); ++i) {
            conn[i] = (*(nodes_begin+i))->Id();
        };

        auto elem_name_it = elem_name_map.find(r_elem.Type());
        if (elem_name_it == elem_name_map.end()) {
            std::stringstream err;
            err << "No Kratos element found for this element type (" << static_cast<int>(r_elem.Type()) << ")!\nOnly the following types are available:";
            for (const auto& r_type_name_pair : elem_name_map) {
                err << "\n\t" << r_type_name_pair.second;
            }
            KRATOS_ERROR << err.str();
        }

        rKratosModelPart.CreateNewElement(
            elem_name_it->second,
            r_elem.Id(),
            conn,
            p_props
        );
    };

    if (rDataComm.IsDistributed()) {
        // this calls the ParallelFillCommunicator
        ParallelEnvironment::CreateFillCommunicatorFromGlobalParallelism(rKratosModelPart, rDataComm)->Execute();
    }

    KRATOS_CATCH("")
}

void CoSimIOConversionUtilities::KratosModelPartToCoSimIOModelPart(
    const Kratos::ModelPart& rKratosModelPart,
    CoSimIO::ModelPart& rCoSimIOModelPart)
{
    KRATOS_TRY

    const int my_rank = rKratosModelPart.GetCommunicator().MyPID();
    const bool is_distributed = rKratosModelPart.IsDistributed();

    const auto is_local_node = [my_rank, is_distributed](const Kratos::Node<3>& rNode) -> bool {
        if (is_distributed) {
            const int node_rank = rNode.FastGetSolutionStepValue(PARTITION_INDEX);
            return node_rank == my_rank;
        } else {
            return true;
        }
    };

    for (const auto& r_node : rKratosModelPart.Nodes()) {
        // must be done in one loop to preserve order
        if (is_local_node(r_node)) {
            rCoSimIOModelPart.CreateNewNode(
                r_node.Id(),
                // TODO: use initial or current coordinates?
                r_node.X0(),
                r_node.Y0(),
                r_node.Z0()
            );
        } else {
            rCoSimIOModelPart.CreateNewGhostNode(
                r_node.Id(),
                // TODO: use initial or current coordinates?
                r_node.X0(),
                r_node.Y0(),
                r_node.Z0(),
                r_node.FastGetSolutionStepValue(PARTITION_INDEX)
            );
        }
    }

    CoSimIO::ConnectivitiesType conn;
    for (const auto& r_elem : rKratosModelPart.Elements()) {
        const auto& r_geom = r_elem.GetGeometry();
        if (conn.size() != r_geom.PointsNumber()) {
            conn.resize(r_geom.PointsNumber());
        };

        for (std::size_t i=0; i<r_geom.PointsNumber(); ++i) {
            conn[i] = r_geom[i].Id();
        }

        auto elem_type_it = elem_type_map.find(r_geom.GetGeometryType());
        KRATOS_ERROR_IF(elem_type_it == elem_type_map.end()) << "No CoSimIO element type found for this Kratos element type (" << static_cast<int>(r_geom.GetGeometryType()) << ")!" << std::endl;

        rCoSimIOModelPart.CreateNewElement(
            r_elem.Id(),
            elem_type_it->second,
            conn
        );
    };

    KRATOS_CATCH("")
}

CoSimIO::Info CoSimIOConversionUtilities::InfoFromParameters(const Parameters rSettings)
{
    KRATOS_TRY

    CoSimIO::Info info;

    for (auto it = rSettings.begin(); it != rSettings.end(); ++it) {
        if      (it->IsString())       info.Set<std::string>(it.name(),   it->GetString());
        else if (it->IsInt())          info.Set<int>(it.name(),           it->GetInt());
        else if (it->IsBool())         info.Set<bool>(it.name(),          it->GetBool());
        else if (it->IsDouble())       info.Set<double>(it.name(),        it->GetDouble());
        else if (it->IsSubParameter()) info.Set<CoSimIO::Info>(it.name(), InfoFromParameters(*it));
        else KRATOS_WARNING("Kratos-CoSimIO") << "Setting with name \"" << it.name() << "\" cannot be converted to CoSimIO::Info and is ignored!" << std::endl;
    }

    return info;

    KRATOS_CATCH("")
}

}  // namespace Kratos.
