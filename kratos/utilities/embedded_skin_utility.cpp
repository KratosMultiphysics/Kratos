//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "variable_utils.h"
#include "embedded_skin_utility.h"
#include "modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "utilities/divide_triangle_2d_3.h"
#include "utilities/divide_tetrahedra_3d_4.h"

namespace Kratos
{
    template<std::size_t TDim>
    void EmbeddedSkinUtility<TDim>::GenerateSkin()
    {
        // Erase all the geometrical entities of embedded skin model part
        this->Clear();

        // Initialize the ids. for the new entries
        // For the temporal ids. there is no necessity of synchronizing between processors
        unsigned int temp_node_id = (mrModelPart.NumberOfNodes() > 0) ? ((mrModelPart.NodesEnd()-1)->Id()) + 1 : 1;
        unsigned int temp_cond_id = (mrModelPart.NumberOfConditions() > 0) ? ((mrModelPart.ConditionsEnd()-1)->Id()) + 1 : 1;

        // Auxiliar vectors to store pointers to the new entities
        // These vectors will be use when renumbering the entities ids.
        ModelPart::NodesContainerType new_nodes_vect;
        ModelPart::ConditionsContainerType new_conds_vect;

        // Set the new condition properties
        Properties::Pointer p_cond_prop = this->SetSkinEntitiesProperties();

        int n_elems = mrModelPart.NumberOfElements();
        for (int i_elem = 0; i_elem < n_elems; ++i_elem){
            const auto it_elem = mrModelPart.ElementsBegin() + i_elem;

            // Get element geometry
            const auto &r_geometry = it_elem->GetGeometry();
            const Vector nodal_distances = this->SetDistancesVector(*it_elem);

            // If split, create the intersection geometry
            if (this->ElementIsSplit(r_geometry, nodal_distances)) {
                const auto p_elem = Kratos::make_intrusive<Element>(*it_elem);
                this->ComputeElementSkin(p_elem, nodal_distances, temp_node_id, temp_cond_id, p_cond_prop, new_nodes_vect, new_conds_vect);
            }
        }

        this->RenumberAndAddSkinEntities(new_nodes_vect, new_conds_vect);
    }

    template<std::size_t TDim>
    void EmbeddedSkinUtility<TDim>::InterpolateMeshVariableToSkin(
        const Variable<double> &rMeshVariable,
		const Variable<double> &rSkinVariable)
    {
        this->InterpolateMeshVariableToSkinSpecialization<double>(rMeshVariable, rSkinVariable);
    }

    template<std::size_t TDim>
    void EmbeddedSkinUtility<TDim>::InterpolateMeshVariableToSkin(
        const Variable<array_1d<double,3>> &rMeshVariable,
		const Variable<array_1d<double,3>> &rSkinVariable)
    {
        this->InterpolateMeshVariableToSkinSpecialization<array_1d<double,3>>(rMeshVariable, rSkinVariable);
    }

    template<std::size_t TDim>
    void EmbeddedSkinUtility<TDim>::InterpolateDiscontinuousMeshVariableToSkin(
        const Variable<double> &rMeshVariable,
		const Variable<double> &rSkinVariable,
        const std::string &rInterfaceSide)
    {
        this->InterpolateMeshVariableToSkinSpecialization<double>(rMeshVariable, rSkinVariable, rInterfaceSide);
    }

    template<std::size_t TDim>
    void EmbeddedSkinUtility<TDim>::InterpolateDiscontinuousMeshVariableToSkin(
        const Variable<array_1d<double,3>> &rMeshVariable,
		const Variable<array_1d<double,3>> &rSkinVariable,
        const std::string &rInterfaceSide)
    {
        this->InterpolateMeshVariableToSkinSpecialization<array_1d<double,3>>(rMeshVariable, rSkinVariable, rInterfaceSide);
    }

    template<std::size_t TDim>
    void EmbeddedSkinUtility<TDim>::Clear()
    {
        // Remove all the entries of the edge nodes map
        mEdgeNodesMap.clear();

        // Flag all the geometrical entities with the TO_ERASE flag
        VariableUtils().SetFlag<ModelPart::NodesContainerType>(TO_ERASE, true, mrSkinModelPart.Nodes());
        VariableUtils().SetFlag<ModelPart::ElementsContainerType>(TO_ERASE, true, mrSkinModelPart.Elements());
        VariableUtils().SetFlag<ModelPart::ConditionsContainerType>(TO_ERASE, true, mrSkinModelPart.Conditions());

        // Remove all the skin model part geometrical entities
        mrSkinModelPart.RemoveNodesFromAllLevels(TO_ERASE);
        mrSkinModelPart.RemoveElementsFromAllLevels(TO_ERASE);
        mrSkinModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
    }

    template<std::size_t TDim>
    void EmbeddedSkinUtility<TDim>::ComputeElementSkin(
        const Element::Pointer pElement,
        const Vector &rNodalDistances,
        unsigned int &rTempNodeId,
        unsigned int &rTempCondId,
        Properties::Pointer pCondProp,
        ModelPart::NodesContainerType &rNewNodesVect,
        ModelPart::ConditionsContainerType &rNewCondsVect)
    {
        // Set the split utility and compute the splitting pattern
        const auto &r_geom = pElement->GetGeometry();
        DivideGeometry::Pointer p_split_utility = this->SetDivideGeometryUtility(r_geom, rNodalDistances);
        p_split_utility->GenerateDivision();
        p_split_utility->GenerateIntersectionsSkin();

        // Get the interface geometries from the splitting pattern (consider only the positive side)
        const unsigned int n_interface_geom = (p_split_utility->mPositiveInterfaces).size();

        std::vector<DivideGeometry::IndexedPointGeometryPointerType> split_interface_geometries;
        split_interface_geometries.reserve(n_interface_geom);
        split_interface_geometries.insert(split_interface_geometries.end(), (p_split_utility->mPositiveInterfaces).begin(), (p_split_utility->mPositiveInterfaces).end());

        // Create the split interface geometries in the skin model part
        for (unsigned int i_int_geom = 0; i_int_geom < n_interface_geom; ++i_int_geom){
            const auto p_int_sub_geom = split_interface_geometries[i_int_geom];
            const auto p_int_sub_geom_type = p_int_sub_geom->GetGeometryType();
            const unsigned int sub_int_geom_n_nodes = p_int_sub_geom->PointsNumber();

            // Fill the new condition nodes array
            Condition::NodesArrayType sub_int_geom_nodes_array;
            for (unsigned int i_node = 0; i_node < sub_int_geom_n_nodes; ++i_node){
                // Create the new node
                const auto &r_sub_int_geom_node = (*p_int_sub_geom)[i_node];
                auto p_new_node = mrSkinModelPart.CreateNewNode(
                    rTempNodeId,
                    r_sub_int_geom_node.X(),
                    r_sub_int_geom_node.Y(),
                    r_sub_int_geom_node.Z());
                rTempNodeId++;
                rNewNodesVect.push_back(p_new_node);
                sub_int_geom_nodes_array.push_back(p_new_node);

                // Add the new node data to the map
                // Note that we take advantage of the Id numbering of the intersection nodes to identify the intersected edge
                const unsigned int intersected_edge_id = r_sub_int_geom_node.Id() - r_geom.PointsNumber();
                const auto new_info = std::make_tuple(pElement, intersected_edge_id);
                mEdgeNodesMap.insert(EdgeNodesMapType::value_type(p_new_node, new_info));
            }

            // Set the new condition pointer
            auto p_new_cond = this->pCreateNewCondition(
                p_int_sub_geom_type,
                sub_int_geom_nodes_array,
                rTempCondId,
                pCondProp);

            // Create the new condition
            rNewCondsVect.push_back(p_new_cond);

            // Update the new elements id. counter
            rTempCondId++;
        }
    }

    template<std::size_t TDim>
    void EmbeddedSkinUtility<TDim>::RenumberAndAddSkinEntities(
        const ModelPart::NodesContainerType &rNewNodesVect,
        const ModelPart::ConditionsContainerType &rNewCondsVect)
    {
        // Once all the entities have been created, renumber the ids.
        // Created entities local number partial reduction
        int n_nodes_local = rNewNodesVect.size();
        int n_conds_local = rNewCondsVect.size();
        int n_nodes_local_scansum, n_conds_local_scansum;
        const auto& r_comm = mrModelPart.GetCommunicator().GetDataCommunicator();
        n_nodes_local_scansum = r_comm.ScanSum(n_nodes_local);
        n_conds_local_scansum = r_comm.ScanSum(n_conds_local);

        // Origin model part number of entities
        int n_nodes_orig = mrModelPart.NumberOfNodes();
        int n_conds_orig = mrModelPart.NumberOfConditions();
        n_nodes_orig = r_comm.SumAll(n_nodes_orig);
        n_conds_orig = r_comm.SumAll(n_conds_orig);

        // Initialize the new ids. values
        std::size_t new_node_id(n_nodes_orig + n_nodes_local_scansum - n_nodes_local + 1);
        std::size_t new_cond_id(n_conds_orig + n_conds_local_scansum - n_conds_local + 1);

        // Set the new entities ids. values
        auto new_nodes_begin = rNewNodesVect.begin();
        auto new_conds_begin = rNewCondsVect.begin();

        for (int i_node = 0; i_node < static_cast<int>(rNewNodesVect.size()); ++i_node){
            auto it_node = new_nodes_begin + i_node;
            const unsigned int new_id = new_node_id + i_node;
            it_node->SetId(new_id);
        }

        for (int i_cond = 0; i_cond < static_cast<int>(rNewCondsVect.size()); ++i_cond){
            auto it_cond = new_conds_begin + i_cond;
            const unsigned int new_id = new_cond_id + i_cond;
            it_cond->SetId(new_id);
        }

        // Add the created conditions to the skin model part
        mrSkinModelPart.AddConditions(rNewCondsVect.begin(), rNewCondsVect.end());

        // Wait for all nodes to renumber its nodes
        r_comm.Barrier();
    }

    template<std::size_t TDim>
    Geometry< Node<3> >::Pointer EmbeddedSkinUtility<TDim>::pCreateNewConditionGeometry(
        const GeometryData::KratosGeometryType &rOriginGeometryType,
        const Condition::NodesArrayType &rNewNodesArray)
    {
        switch(rOriginGeometryType){
            case GeometryData::KratosGeometryType::Kratos_Line2D2:
                return Kratos::make_shared<Line2D2< Node<3> > >(rNewNodesArray);
            case GeometryData::KratosGeometryType::Kratos_Triangle3D3:
                return Kratos::make_shared<Triangle3D3< Node<3> > >(rNewNodesArray);
            default:
                KRATOS_ERROR << "Implement the skin generation for the intersection geometry type: " << rOriginGeometryType;
        }
    }

    template<std::size_t TDim>
    Condition::Pointer EmbeddedSkinUtility<TDim>::pCreateNewCondition(
        const GeometryData::KratosGeometryType &rOriginGeometryType,
        const Condition::NodesArrayType &rNewNodesArray,
        const unsigned int &rConditionId,
        const Properties::Pointer pConditionProperties)
    {
        auto p_new_geom = this->pCreateNewConditionGeometry(rOriginGeometryType, rNewNodesArray);
        return mrConditionPrototype.Create(rConditionId, p_new_geom, pConditionProperties);
    }

    template<>
    const std::string EmbeddedSkinUtility<2>::GetConditionType()
    {
        return "LineCondition2D2N";
    }

    template<>
    const std::string EmbeddedSkinUtility<3>::GetConditionType()
    {
        return "SurfaceCondition3D3N";
    }

    template<std::size_t TDim>
    Properties::Pointer EmbeddedSkinUtility<TDim>::SetSkinEntitiesProperties()
    {
        // Set the properties for the new elements depending if the
        // element is in the positive or negative side of the cut.
        // In this way, two layers will appear in GiD.
        unsigned int max_prop_id = 0;
        for (auto it_prop = mrModelPart.GetRootModelPart().PropertiesBegin(); it_prop < mrModelPart.GetRootModelPart().PropertiesEnd(); ++it_prop){
            if (max_prop_id < it_prop->Id()){
                max_prop_id = it_prop->Id();
            }
        }

        Properties::Pointer p_new_prop = Kratos::make_shared<Properties>(max_prop_id + 1);
        mrModelPart.AddProperties(p_new_prop);

        return p_new_prop;
    }

    template<std::size_t TDim>
    bool inline EmbeddedSkinUtility<TDim>::ElementIsSplit(
        const Geometry<Node<3>> &rGeometry,
        const Vector &rNodalDistances)
    {
        unsigned int n_pos (0), n_neg(0);
        for (unsigned int i_node = 0; i_node < rGeometry.PointsNumber(); ++i_node) {
            if (rNodalDistances[i_node] > 0.0) {
                n_pos++;
            } else {
                n_neg++;
            }
        }

        return (n_pos > 0 && n_neg > 0) ? true : false;
    }

    template<std::size_t TDim>
    const Vector EmbeddedSkinUtility<TDim>::SetDistancesVector(const Element &rElement)
    {
        const auto &r_geom = rElement.GetGeometry();
        Vector nodal_distances(r_geom.PointsNumber());

        if (mLevelSetType == LevelSetTypeEnum::Continuous) {
            // Continuous nodal distance function case
            for (unsigned int i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
                nodal_distances[i_node] = r_geom[i_node].FastGetSolutionStepValue(DISTANCE);
            }
        } else if (mLevelSetType == LevelSetTypeEnum::Discontinuous) {
            // Discontinuous elemental distance function case
            nodal_distances = rElement.GetValue(ELEMENTAL_DISTANCES);
        } else {
            KRATOS_ERROR << "Level set type must be either \'continuous\' or \'discontinuous\'";
        }

        return nodal_distances;
    }

    template<std::size_t TDim>
    DivideGeometry::Pointer EmbeddedSkinUtility<TDim>::SetDivideGeometryUtility(
        const Geometry<Node<3>> &rGeometry,
        const Vector& rNodalDistances)
    {
        // Get the geometry type
        const GeometryData::KratosGeometryType geometry_type = rGeometry.GetGeometryType();

        // Return the divide geometry utility
        switch (geometry_type){
            case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
                return Kratos::make_shared<DivideTriangle2D3>(rGeometry, rNodalDistances);
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
                return Kratos::make_shared<DivideTetrahedra3D4>(rGeometry, rNodalDistances);
            default:
                KRATOS_ERROR << "Asking for a non-implemented divide geometry utility.";
        }
    }

    template<std::size_t TDim>
    ModifiedShapeFunctions::UniquePointer EmbeddedSkinUtility<TDim>::pCreateModifiedShapeFunctions(
        const Geometry<Node<3>>::Pointer pGeometry,
        const Vector& rNodalDistances)
    {
        // Get the geometry type
        const GeometryData::KratosGeometryType geometry_type = pGeometry->GetGeometryType();

        // Return the modified shape functions utility
        if (mLevelSetType == LevelSetTypeEnum::Continuous) {
            switch (geometry_type) {
                case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
                    return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rNodalDistances);
                case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
                    return Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(pGeometry, rNodalDistances);
                default:
                    KRATOS_ERROR << "Asking for a non-implemented modified shape functions geometry.";
            }
        } else if (mLevelSetType == LevelSetTypeEnum::Discontinuous) {
            switch (geometry_type) {
                case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
                    return Kratos::make_unique<Triangle2D3AusasModifiedShapeFunctions>(pGeometry, rNodalDistances);
                case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
                    return Kratos::make_unique<Tetrahedra3D4AusasModifiedShapeFunctions>(pGeometry, rNodalDistances);
                default:
                    KRATOS_ERROR << "Asking for a non-implemented Ausas modified shape functions geometry.";
            }
        } else {
            KRATOS_ERROR << "Level set type must be either \'continuous\' or \'discontinuous\'.";
        }
    }

    template<std::size_t TDim>
    Matrix EmbeddedSkinUtility<TDim>::GetModifiedShapeFunctionsValues(
        const ModifiedShapeFunctions::UniquePointer &rpModifiedShapeFunctions,
        const std::string &rInterfaceSide) const
    {
        Vector w_int;
        Matrix int_sh_func;
        ModifiedShapeFunctions::ShapeFunctionsGradientsType int_grads;

        if (rInterfaceSide == "positive") {
            rpModifiedShapeFunctions->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
                int_sh_func,
                int_grads,
                w_int,
                GeometryData::GI_GAUSS_2);
        } else if (rInterfaceSide == "negative") {
            rpModifiedShapeFunctions->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
                int_sh_func,
                int_grads,
                w_int,
                GeometryData::GI_GAUSS_2);
        } else {
            KRATOS_ERROR << "Interface side must be either 'positive' or 'negative'. Got " << rInterfaceSide;
        }

        return int_sh_func;
    }

    template<std::size_t TDim>
    Matrix EmbeddedSkinUtility<TDim>::GetModifiedShapeFunctionsValuesOnEdge(
        const ModifiedShapeFunctions::UniquePointer &rpModifiedShapeFunctions,
        const std::string &rInterfaceSide) const
    {
        Matrix edge_N_values;
        if (rInterfaceSide == "positive") {
            rpModifiedShapeFunctions->ComputeShapeFunctionsOnPositiveEdgeIntersections(edge_N_values);
        } else if (rInterfaceSide == "negative") {
            rpModifiedShapeFunctions->ComputeShapeFunctionsOnNegativeEdgeIntersections(edge_N_values);
        } else {
            KRATOS_ERROR << "Interface side must be either 'positive' or 'negative'. Got " << rInterfaceSide;
        }

        return edge_N_values;
    }

    template class Kratos::EmbeddedSkinUtility<2>;
	template class Kratos::EmbeddedSkinUtility<3>;
}
