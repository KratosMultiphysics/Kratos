//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Andrea Gorgi

// Project includes
#include "iga_contact_process_gap_sbm.h"
#include "includes/define.h"
#include "iga_application_variables.h"
#include "includes/global_pointer_variables.h"
#include "spatial_containers/bins_dynamic.h"
#include "containers/pointer_vector.h"

#include <limits>

namespace Kratos
{

    IgaContactProcessGapSbm::IgaContactProcessGapSbm(
        Model& rModel, Parameters ThisParameters) :
        Process(),
        mpModel(&rModel),
        mParameters(ThisParameters)
    {
        mEchoLevel = mParameters.Has("echo_level") ? mParameters["echo_level"].GetInt() : 0;

        KRATOS_ERROR_IF_NOT(ThisParameters.Has("analysis_model_part_name")) << "::[IgaContactProcessGapSbm]::"
                            << " Missing \"analysis_model_part_name\" parameter. "<< std::endl;

        KRATOS_ERROR_IF_NOT(ThisParameters.Has("contact_sub_model_part_name")) << "::[IgaContactProcessGapSbm]::"
                            << " Missing \"contact_sub_model_part_name\" parameter. "<< std::endl;

        KRATOS_ERROR_IF_NOT(ThisParameters.Has("contact_parameters")) << "::[IgaContactProcessGapSbm]::"
                            << " Missing \"contact_parameters\" parameter. "<< std::endl;

        KRATOS_ERROR_IF_NOT(ThisParameters["contact_parameters"].Has("slave_model_part")) << "::[IgaContactProcessGapSbm]::"
                            << " Missing \"contact_parameters/slave_model_part\" parameter. "<< std::endl;

        KRATOS_ERROR_IF_NOT(ThisParameters["contact_parameters"].Has("master_model_part")) << "::[IgaContactProcessGapSbm]::"
                            << " Missing \"contact_parameters/master_model_part\" parameter. "<< std::endl;

        //-----------------------------------------------------------------------------------------------
        // Obtain SLAVE interface b_reps
        std::string slave_model_part_name = mParameters["contact_parameters"]["slave_model_part"]["sub_model_part_name"].GetString();
        const std::string slave_layer_name = mParameters["contact_parameters"]["slave_model_part"]["layer_name"].GetString();
        slave_model_part_name += "." + slave_layer_name;

        KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(slave_model_part_name)) << "ERROR: SLAVE MODEL PART "
                                                << slave_model_part_name << " NOT CREATED" << std::endl;

        mrSlaveModelPart = &(mpModel->GetModelPart(slave_model_part_name));

        const IndexType slave_property_id = mParameters["contact_parameters"]["slave_model_part"]["property_id"].GetInt();
        mpPropSlave = mrSlaveModelPart->pGetProperties(slave_property_id);

        //--------------------------------------------------------------------------------------------------
        // Obtain MASTER interface b_reps
        std::string master_model_part_name = mParameters["contact_parameters"]["master_model_part"]["sub_model_part_name"].GetString();
        const std::string master_layer_name = mParameters["contact_parameters"]["master_model_part"]["layer_name"].GetString();
        master_model_part_name += "." + master_layer_name;

        KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(master_model_part_name)) << "ERROR: MASTER MODEL PART "
                                                << master_model_part_name << " NOT CREATED" << std::endl;

        mrMasterModelPart = &(mpModel->GetModelPart(master_model_part_name));

        const IndexType master_property_id = mParameters["contact_parameters"]["master_model_part"]["property_id"].GetInt();
        mpPropMaster = mrMasterModelPart->pGetProperties(master_property_id);

        //------------------------------------------------------------------------
        std::string analysis_model_part_name = mParameters["analysis_model_part_name"].GetString();
        std::string contact_sub_model_part_name = mParameters["contact_sub_model_part_name"].GetString();

        std::string contact_model_part_name = analysis_model_part_name + ".ContactInterface." + contact_sub_model_part_name;

        mrContactModelPart = &(mpModel->CreateModelPart(contact_model_part_name));
    }

    void IgaContactProcessGapSbm::Execute()
    {
        // TODO: create or retrieve "contact" submodelpart of the contact model part.
        // NOTE: this is the target where all coupling quadrature geometries will be stored.
        ModelPart& r_contact_sub_model_part = mrContactModelPart->HasSubModelPart("contact")
            ? mrContactModelPart->GetSubModelPart("contact")
            : mrContactModelPart->CreateSubModelPart("contact");

        if (mrMasterModelPart->NumberOfConditions() > 0) {
            std::vector<IndexType> condition_ids;
            condition_ids.reserve(mrMasterModelPart->NumberOfConditions());
            for (auto& r_condition : mrMasterModelPart->Conditions()) {
                condition_ids.push_back(r_condition.Id());
            }
            r_contact_sub_model_part.AddConditions(condition_ids);
        }

        if (!mrSlaveModelPart->HasSubModelPart("contact")) {
            KRATOS_WARNING("IgaContactProcessGapSbm")
                << "Slave model part has no 'contact' submodelpart. Nothing to do." << std::endl;
            return;
        }

        using PointType = NodeType;
        using PointTypePointer = NodeType::Pointer;
        using PointVector = std::vector<PointType::Pointer>;
        using PointIterator = std::vector<PointType::Pointer>::iterator;
        using DistanceVector = std::vector<double>;
        using DistanceIterator = std::vector<double>::iterator;
        using DynamicBins = BinsDynamic<3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using PointerType = DynamicBins::PointerType;

        ModelPart& r_slave_contact = mrSlaveModelPart->GetSubModelPart("contact");

        // Build bins for slave contact geometries (based on centers).
        // TODO: use deformed coordinates for the centers once available.
        PointVector slave_center_points;
        slave_center_points.reserve(r_slave_contact.NumberOfGeometries());
        for (auto& r_slave_geometry : r_slave_contact.Geometries()) {
            const auto& r_center = r_slave_geometry.Center();
            slave_center_points.push_back(
                Kratos::make_intrusive<PointType>(r_slave_geometry.Id(), r_center.X(), r_center.Y(), r_center.Z()));
        }

        if (slave_center_points.empty()) {
            return;
        }

        DynamicBins slave_bins(slave_center_points.begin(), slave_center_points.end());

        IndexType next_projection_node_id = 1;
        if (mrContactModelPart->GetRootModelPart().NumberOfNodes() > 0) {
            next_projection_node_id = (mrContactModelPart->GetRootModelPart().NodesEnd() - 1)->Id() + 1;
        }

        for (auto& r_master_condition : mrMasterModelPart->Conditions()) {
            r_master_condition.SetValue(IDENTIFIER, "MASTER");
            const auto& r_master_center = r_master_condition.GetGeometry().Center();
            // TODO: use deformed coordinates for the master quadrature point once available.
            PointType master_query_point(0, r_master_center.X(), r_master_center.Y(), r_master_center.Z());

            PointerType p_nearest = nullptr;
            double nearest_distance = std::numeric_limits<double>::max();
            slave_bins.SearchNearestPoint(master_query_point, p_nearest, nearest_distance);
            if (!p_nearest) {
                continue;
            }

            auto p_slave_geometry = r_slave_contact.pGetGeometry(p_nearest->Id());
            if (!p_slave_geometry) {
                continue;
            }

            auto p_slave_projection_node = Kratos::make_intrusive<NodeType>(
                next_projection_node_id++, p_nearest->X(), p_nearest->Y(), p_nearest->Z());

            r_master_condition.GetGeometry().SetValue(PROJECTION_NODE, p_slave_projection_node);

            std::vector<Geometry<Node>::Pointer> neighbour_geometries;
            if (p_slave_geometry->Has(NEIGHBOUR_GEOMETRIES)) {
                neighbour_geometries = p_slave_geometry->GetValue(NEIGHBOUR_GEOMETRIES);
            }
            neighbour_geometries.push_back(p_slave_geometry);
            p_slave_projection_node->SetValue(NEIGHBOUR_GEOMETRIES, neighbour_geometries);
        }

        for (auto& r_slave_geometry : r_slave_contact.Geometries()) {
            r_slave_geometry.SetValue(IDENTIFIER, "SLAVE");
        }

        // NOTE: no extra work is required on the slave model part beyond the projection.

        EntitiesUtilities::InitializeEntities<Condition>(mrContactModelPart->GetParentModelPart());
    }

} // namespace Kratos
