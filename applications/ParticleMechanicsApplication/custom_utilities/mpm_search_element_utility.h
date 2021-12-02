//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#ifndef KRATOS_MPM_SEARCH_ELEMENT_UTILITY
#define KRATOS_MPM_SEARCH_ELEMENT_UTILITY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/quadrature_points_utility.h"

#include "particle_mechanics_application_variables.h"

#include "geometries/geometry.h"
#include "includes/model_part.h"

#include "pqmpm_partition_utilities.h"

namespace Kratos
{
namespace MPMSearchElementUtility
{
    // Standard types
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef Node<3> NodeType;
    typedef typename ModelPart::GeometryType GeometryType;

    inline double CrossProductDet2D(array_1d<double, 3> VectorA, array_1d<double, 3> VectorB)
    {
        return (VectorA[0] * VectorB[1] - VectorB[0] * VectorA[1]);
    }

    inline bool CheckIsInside(const GeometryType& rGeom, array_1d<double, 3>& LocalCoords, const array_1d<double, 3>& Coords, const double Tolerance, const bool IsCalcLocalCoords = true)
    {

        // TODO some optimisation for simple 2D shapes.

        return rGeom.IsInside(Coords, LocalCoords, Tolerance);
    }

    inline void ConstructNeighbourRelations(GeometryType& rGeom, const ModelPart& rBackgroundGridModelPart)
    {
        std::vector<typename Geometry<Node<3>>::Pointer> geometry_neighbours;
        for (IndexType j = 0; j < rBackgroundGridModelPart.NumberOfElements(); j++)
        {
            auto p_geometry_neighbour = (rBackgroundGridModelPart.ElementsBegin() + j)->pGetGeometry();
            if (p_geometry_neighbour->Id() != rGeom.Id()) // dont add the parent as its own neighbour
            {
                for (IndexType n = 0; n < p_geometry_neighbour->size(); n++)
                {
                    for (IndexType k = 0; k < rGeom.size(); k++)
                    {
                        if (rGeom[k].Id() == (*p_geometry_neighbour)[n].Id()) {
                            // Prevent duplicate additions
                            bool add_entry = true;
                            for (size_t i = 0; i < geometry_neighbours.size(); i++)
                            {
                                if (geometry_neighbours[i]->Id() == p_geometry_neighbour->Id())
                                {
                                    add_entry = false;
                                    break;
                                }
                            }
                            if (add_entry)
                            {
                                geometry_neighbours.push_back(p_geometry_neighbour);
                            }
                            break;
                        }
                    }
                }
            }

        }
        #pragma omp critical
        rGeom.SetValue(GEOMETRY_NEIGHBOURS, geometry_neighbours);
    }


    inline bool IsExplicitAndNeedsCorrection(GeometryType::Pointer pQuadraturePoint, const ProcessInfo& rProcessInfo)
    {
        if (rProcessInfo.Has(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)) {
            if (rProcessInfo.GetValue(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)) {
                if (pQuadraturePoint->IntegrationPointsNumber() == 1)
                {
                    for (size_t i = 0; i < pQuadraturePoint->ShapeFunctionsValues().size2(); ++i)
                    {
                        if (pQuadraturePoint->ShapeFunctionsValues()(0, i) < std::numeric_limits<double>::epsilon()) return true;
                    }
                }
            }
        }

        return false;
    }

    inline GeometryType& FindGridGeom(GeometryType& rParentGeom,
        const ModelPart& rBackgroundGridModelPart,
        const double Tolerance,
        const array_1d<double, 3>& xg,
        array_1d<double, 3>& rLocalCoords,
        const ProcessInfo& rProcessInfo,
        bool& IsFound)
    {
        IsFound = false;

        if (CheckIsInside(rParentGeom, rLocalCoords, xg, Tolerance)) {
            IsFound = true;
            return rParentGeom;
        }
        else
        {
            if (!rParentGeom.Has(GEOMETRY_NEIGHBOURS))
                ConstructNeighbourRelations(rParentGeom, rBackgroundGridModelPart);

            auto& geometry_neighbours = rParentGeom.GetValue(GEOMETRY_NEIGHBOURS);
            for (IndexType k = 0; k < geometry_neighbours.size(); ++k) {
                if (CheckIsInside(*geometry_neighbours[k], rLocalCoords, xg, Tolerance)) {
                    IsFound = true;
                    return *(geometry_neighbours[k].get());
                }
            }
        }

        return rParentGeom;
    }


    inline void UpdatePartitionedQuadraturePoint(const ModelPart& rBackgroundGridModelPart,
        const array_1d<double, 3>& rCoordinates,
        Element& rMasterMaterialPoint,
        typename GeometryType::Pointer pQuadraturePointGeometry,
        const double Tolerance)
    {
        KRATOS_TRY;

        array_1d<double, 3> local_coords;
        pQuadraturePointGeometry->IsInside(rCoordinates, local_coords, Tolerance);
        PQMPMPartitionUtilities::PartitionMasterMaterialPointsIntoSubPoints(rBackgroundGridModelPart, rCoordinates,
            local_coords, rMasterMaterialPoint, pQuadraturePointGeometry, Tolerance);

        KRATOS_CATCH("");
    }


    inline void NeighbourSearchElements(const ModelPart& rMPMModelPart,
        const ModelPart& rBackgroundGridModelPart,
        std::vector<typename Element::Pointer>& rMissingElements,
        const double Tolerance)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rMPMModelPart.Elements().size()); ++i) {
            auto element_itr = (rMPMModelPart.ElementsBegin() + i);
            array_1d<double, 3> local_coordinates;
            bool is_found = false;
            std::vector<array_1d<double, 3>> xg;
            element_itr->CalculateOnIntegrationPoints(MP_COORD, xg, rBackgroundGridModelPart.GetProcessInfo());

            GeometryType& r_found_geom = FindGridGeom(element_itr->GetGeometry().GetGeometryParent(0),
                rBackgroundGridModelPart, Tolerance, xg[0], local_coordinates,
                rMPMModelPart.GetProcessInfo(), is_found);

            if (is_found)
            {
                const bool is_pqmpm = (rBackgroundGridModelPart.GetProcessInfo().Has(IS_PQMPM))
                    ? rBackgroundGridModelPart.GetProcessInfo().GetValue(IS_PQMPM) : false;
                if (is_pqmpm)
                {
                    // Updates the quadrature point geometry.
                    (*element_itr).GetGeometry().SetGeometryParent(&r_found_geom);
                    PQMPMPartitionUtilities::PartitionMasterMaterialPointsIntoSubPoints(rBackgroundGridModelPart, xg[0],
                        local_coordinates, *element_itr, element_itr->pGetGeometry(), Tolerance);
                }
                else
                {
                    CreateQuadraturePointsUtility<Node<3>>::UpdateFromLocalCoordinates(
                        element_itr->pGetGeometry(), local_coordinates,
                        element_itr->GetGeometry().IntegrationPoints()[0].Weight(), r_found_geom);
                }

                if (IsExplicitAndNeedsCorrection(element_itr->pGetGeometry(), rBackgroundGridModelPart.GetProcessInfo()))
                    is_found = false;
                else {
                    for (IndexType j = 0; j < r_found_geom.PointsNumber(); ++j)
                        r_found_geom.Points()[j].Set(ACTIVE);
                }
            }
            if(!is_found)
            {
                #pragma omp critical
                rMissingElements.push_back(&*element_itr);
            }
        }
    }


    //


    inline void NeighbourSearchConditions(const ModelPart& rMPMModelPart,
        const ModelPart& rBackgroundGridModelPart,
        std::vector<typename Condition::Pointer>& rMissingConditions,
        const double Tolerance)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rMPMModelPart.Conditions().size()); ++i) {
            auto condition_itr = rMPMModelPart.Conditions().begin() + i;

            std::vector<array_1d<double, 3>> xg;
            condition_itr->CalculateOnIntegrationPoints(MPC_COORD, xg, rMPMModelPart.GetProcessInfo());

            if (xg.size() > 0 && condition_itr->Is(BOUNDARY))
            {
                array_1d<double, 3> local_coordinates;
                bool is_found = false;

                GeometryType& r_found_geom = FindGridGeom(condition_itr->GetGeometry().GetGeometryParent(0),
                    rBackgroundGridModelPart, Tolerance, xg[0], local_coordinates,
                    rMPMModelPart.GetProcessInfo(), is_found);

                if (is_found)
                {
                    CreateQuadraturePointsUtility<Node<3>>::UpdateFromLocalCoordinates(
                        condition_itr->pGetGeometry(), local_coordinates,
                        condition_itr->GetGeometry().IntegrationPoints()[0].Weight(), r_found_geom);
                        
                    for (IndexType j = 0; j < r_found_geom.PointsNumber(); ++j)
                        r_found_geom[j].Set(ACTIVE);
                }
                else
                {
                    #pragma omp critical
                    rMissingConditions.push_back(&*condition_itr);
                }
            }
        }
    }


    inline bool IsFixExplicitAndOnElementEdge(const Vector& N, const ProcessInfo& rProcessInfo)
    {
        if (rProcessInfo.Has(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)) {
            if (rProcessInfo.GetValue(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)) {
                // check if MP is exactly on the edge of the element, this gives spurious strains in explicit
                for (SizeType i = 0; i < N.size(); ++i) {
                    if (std::abs(N[i]) < std::numeric_limits<double>::epsilon()) {
                        return true;
                    }
                }
            }
        }

        return false;
    }


    template <std::size_t TDimension>
    void BinBasedSearchElementsAndConditions(ModelPart& rMPMModelPart,
        ModelPart& rBackgroundGridModelPart,
        std::vector<typename Element::Pointer>& rMissingElements,
        std::vector<typename Condition::Pointer>& rMissingConditions,
        const std::size_t MaxNumberOfResults, const double Tolerance)
    {
        const ProcessInfo& r_process_info = rBackgroundGridModelPart.GetProcessInfo();
        bool is_pqmpm = (r_process_info.Has(IS_PQMPM))
            ? r_process_info.GetValue(IS_PQMPM) : false;

        // Search background grid and make element active
        Vector N;
        const int max_result = 1000;

        #pragma omp parallel
        {
            BinBasedFastPointLocator<TDimension> SearchStructure(rBackgroundGridModelPart);
            SearchStructure.UpdateSearchDatabase();
            typename BinBasedFastPointLocator<TDimension>::ResultContainerType results(max_result);

            // Element search and assign background grid
            #pragma omp for
            for (int i = 0; i < static_cast<int>(rMissingElements.size()); ++i) {
                auto element_itr = *(rMissingElements.begin() + i);
                std::vector<array_1d<double, 3>> xg;
                element_itr->CalculateOnIntegrationPoints(MP_COORD, xg, rMPMModelPart.GetProcessInfo());
                typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();
                Element::Pointer pelem;

                // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                bool is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);

                if (is_found == true) {
                    if (IsFixExplicitAndOnElementEdge(N, r_process_info) && !is_pqmpm) {
                        // MP is exactly on the edge. Now we give it a little 'nudge'
                        array_1d<double, 3> xg_nudged = array_1d<double, 3>(xg[0]);
                        std::vector<array_1d<double, 3>> mp_vel;
                        element_itr->CalculateOnIntegrationPoints(MP_VELOCITY, mp_vel, rMPMModelPart.GetProcessInfo());
                        xg_nudged += r_process_info[DELTA_TIME] / 1000.0 * mp_vel[0];
                        if (SearchStructure.FindPointOnMesh(xg_nudged, N, pelem, result_begin, MaxNumberOfResults, Tolerance)) {
                            element_itr->SetValuesOnIntegrationPoints(MP_COORD, { xg_nudged }, rMPMModelPart.GetProcessInfo());
                            KRATOS_INFO("MPMSearchElementUtility") << "WARNING: To prevent spurious explicit stresses, Material Point "
                                << element_itr->Id() << " was nudged." << std::endl;
                        } else {
                            is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);
                            KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Material Point " << element_itr->Id()
                                << " lies exactly on an element edge and may give spurious results." << std::endl;
                        }
                    }
                    pelem->Set(ACTIVE);

                    const bool is_pqmpm = (rBackgroundGridModelPart.GetProcessInfo().Has(IS_PQMPM))
                        ? rBackgroundGridModelPart.GetProcessInfo().GetValue(IS_PQMPM) : false;
                    if (is_pqmpm)
                    {
                        // Updates the quadrature point geometry.
                        (*element_itr).GetGeometry().SetGeometryParent((pelem->pGetGeometry().get()));
                        UpdatePartitionedQuadraturePoint(rBackgroundGridModelPart, xg[0],
                            *element_itr, pelem->pGetGeometry(), Tolerance);
                    }
                    else
                    {
                        auto p_quadrature_point_geometry = element_itr->pGetGeometry();
                        array_1d<double, 3> local_coordinates;
                        p_quadrature_point_geometry->PointLocalCoordinates(local_coordinates, xg[0]);
                        CreateQuadraturePointsUtility<Node<3>>::UpdateFromLocalCoordinates(
                            p_quadrature_point_geometry, local_coordinates,
                            p_quadrature_point_geometry->IntegrationPoints()[0].Weight(), pelem->GetGeometry());
                    }

                    auto& r_geometry = element_itr->GetGeometry();
                    for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j)
                        r_geometry[j].Set(ACTIVE);
                }
                else {
                    KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Search Element for Material Point: " << element_itr->Id()
                        << " is failed. Geometry is cleared." << std::endl;

                    element_itr->GetGeometry().clear();
                    element_itr->Reset(ACTIVE);
                    element_itr->Set(TO_ERASE);
                }
            }

            // Condition search and assign background grid
            #pragma omp for
            for (int i = 0; i < static_cast<int>(rMissingConditions.size()); ++i) {
                auto condition_itr = *(rMissingConditions.begin() + i);
                std::vector<array_1d<double, 3>> xg;
                condition_itr->CalculateOnIntegrationPoints(MPC_COORD, xg, rMPMModelPart.GetProcessInfo());

                if (xg.size() > 0) {
                    // Only search for particle based BCs!
                    // Grid BCs are still applied on MP_model_part but we don't want to search for them.
                    typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();
                    Element::Pointer pelem;

                    // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                    bool is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);

                    if (is_found == true) {
                        auto p_quadrature_point_geometry = condition_itr->pGetGeometry();
                        array_1d<double, 3> local_coordinates;
                        p_quadrature_point_geometry->PointLocalCoordinates(local_coordinates, xg[0]);
                        CreateQuadraturePointsUtility<Node<3>>::UpdateFromLocalCoordinates(
                            p_quadrature_point_geometry, local_coordinates,
                            p_quadrature_point_geometry->IntegrationPoints()[0].Weight(), pelem->GetGeometry());
                        
                        auto& r_geometry = condition_itr->GetGeometry();

                        for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j)
                            r_geometry[j].Set(ACTIVE);
                    } else {
                        KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Search Element for Material Point Condition: " << condition_itr->Id()
                            << " is failed. Geometry is cleared." << std::endl;

                        condition_itr->GetGeometry().clear();
                        condition_itr->Reset(ACTIVE);
                        condition_itr->Set(TO_ERASE);
                    }
                }
            }
        }
    }

    inline void ResetElementsAndNodes(ModelPart& rBackgroundGridModelPart)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rBackgroundGridModelPart.Elements().size()); ++i) {
            auto element_itr = rBackgroundGridModelPart.Elements().begin() + i;
            auto& r_geometry = element_itr->GetGeometry();
            element_itr->Reset(ACTIVE);

            for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j)
                r_geometry[j].Reset(ACTIVE);

        }
    }

    /**
     * @brief Search element connectivity for each particle
     * @details A search is performed to know in which grid element the material point falls.
     * If one or more material points fall in the grid element, the grid element is
     * set to be active and its connectivity is associated to the material point
     * element.
     * STEPS:
     * 1) All the elements are set to be INACTIVE
     * 2) A searching is performed and the grid elements which contain at least a MP are set to be ACTIVE
     *
     */
    template<std::size_t TDimension>
    void SearchElement(ModelPart& rBackgroundGridModelPart, ModelPart& rMPMModelPart, const std::size_t MaxNumberOfResults,
        const double Tolerance)
    {
        ResetElementsAndNodes(rBackgroundGridModelPart);

        std::vector<typename Element::Pointer> missing_elements;
        std::vector<typename Condition::Pointer> missing_conditions;

        NeighbourSearchElements(rMPMModelPart, rBackgroundGridModelPart, missing_elements, Tolerance);
        NeighbourSearchConditions(rMPMModelPart, rBackgroundGridModelPart, missing_conditions, Tolerance);

        if (missing_conditions.size() > 0 || missing_elements.size() > 0)
            BinBasedSearchElementsAndConditions<TDimension>(rMPMModelPart,
                rBackgroundGridModelPart, missing_elements, missing_conditions,
                MaxNumberOfResults, Tolerance);

        

        for (auto& submodelpart : rMPMModelPart.SubModelParts())
        {
            // Create additional Point load conditions if there is a line load
            if (submodelpart.HasSubModelPart("lagrange_condition"))
            {
                auto& lagrange_model_part =    rMPMModelPart.HasSubModelPart("lagrange_condition")
                                        ? rMPMModelPart.GetSubModelPart("lagrange_condition")
                                        : rMPMModelPart.CreateSubModelPart("lagrange_condition");
                

                // Delete old conditions and the corresponding nodes
                std::vector<Kratos::IndexType> toRemove;
                const auto it_condition_begin = lagrange_model_part.Conditions().begin();
                for (int c = 0; c < static_cast<int>(lagrange_model_part.Conditions().size()); ++c) {
                    auto it_condition = it_condition_begin + c;
                    std::vector<int> node_id(1);
                    it_condition->CalculateOnIntegrationPoints(MPC_CORRESPONDING_NODE_ID, node_id, rMPMModelPart.GetProcessInfo() );
                    rBackgroundGridModelPart.RemoveNodeFromAllLevels(node_id[0]);
                    
                    toRemove.push_back(it_condition->Id());
                }
                for (unsigned int c : toRemove) lagrange_model_part.RemoveConditionFromAllLevels(c);

                Vector N;
                const int max_result = 1000;

                BinBasedFastPointLocator<TDimension> SearchStructure(rBackgroundGridModelPart);
                SearchStructure.UpdateSearchDatabase();
                typename BinBasedFastPointLocator<TDimension>::ResultContainerType results(max_result);

                // Create Lagrange Condition
                const Condition& new_condition = KratosComponents<Condition>::Get("MPMParticleLagrangeDirichletCondition");
                IndexType condition_id = submodelpart.GetRootModelPart().Conditions().back().Id() + 1;

                for(auto mpc=submodelpart.ConditionsBegin(); mpc!=submodelpart.ConditionsEnd(); ++mpc){
                    mpc->Set(VISITED,false);
                }

                for(auto mpc=submodelpart.ConditionsBegin(); mpc!=submodelpart.ConditionsEnd(); ++mpc){
                    if (mpc->IsNot(VISITED)){
                        int id_parent = mpc->GetGeometry().GetGeometryParent(0).Id();
                        Properties::Pointer properties = mpc->pGetProperties();
                        std::vector<double> mpc_area(1);
                        mpc_area[0]=0.0;
                        std::vector<array_1d<double, 3>> xg = {ZeroVector(3)};
                        std::vector<array_1d<double, 3>> normal = {ZeroVector(3)};
                        std::vector<array_1d<double, 3>> imposed_disp = {ZeroVector(3)};
                        std::vector<array_1d<double, 3>> mpc_velocity = { ZeroVector(3) };
                        std::vector<array_1d<double, 3>> mpc_imposed_velocity = { ZeroVector(3) };
                        std::vector<array_1d<double, 3>> mpc_acceleration = { ZeroVector(3) };
                        std::vector<array_1d<double, 3>> mpc_imposed_acceleration = { ZeroVector(3) };
                        std::vector<array_1d<double, 3>> mpc_contact_force = { ZeroVector(3) };
                        std::vector<int> mpc_counter(1);
                        mpc_counter[0]=0;
                        int corresponding_condition_id = int(condition_id);
                        int lagrange_node_id = rBackgroundGridModelPart.Nodes().size() + 1;
                        for(auto mpc2=submodelpart.ConditionsBegin(); mpc2!=submodelpart.ConditionsEnd(); ++mpc2){
                            int id_parent2 = mpc2->GetGeometry().GetGeometryParent(0).Id();                          
                            if (id_parent == id_parent2){
                                mpc2->Set(VISITED,true);
                                std::vector<array_1d<double, 3>> xg_tmp;
                                std::vector<double> mpc_area_tmp(1);
                                std::vector<array_1d<double, 3>> normal_tmp;
                                std::vector<array_1d<double, 3>> imposed_disp_tmp;
                                mpc2->CalculateOnIntegrationPoints(MPC_AREA, mpc_area_tmp, rMPMModelPart.GetProcessInfo() );
                                mpc2->CalculateOnIntegrationPoints(MPC_COORD, xg_tmp, rMPMModelPart.GetProcessInfo() );
                                mpc2->CalculateOnIntegrationPoints(MPC_NORMAL, normal_tmp, rMPMModelPart.GetProcessInfo() );
                                mpc2->CalculateOnIntegrationPoints(MPC_IMPOSED_DISPLACEMENT, imposed_disp_tmp, rMPMModelPart.GetProcessInfo() );
                                mpc2->SetValuesOnIntegrationPoints(MPC_CORRESPONDING_NODE_ID, {lagrange_node_id}, rMPMModelPart.GetProcessInfo() );
                                mpc2->SetValuesOnIntegrationPoints(MPC_CORRESPONDING_CONDITION_ID, {corresponding_condition_id}, rMPMModelPart.GetProcessInfo() );
                                mpc_area[0] += mpc_area_tmp[0];
                                xg[0] += xg_tmp[0];
                                normal[0] += normal_tmp[0];
                                imposed_disp[0] += imposed_disp_tmp[0];
                                mpc_counter[0] += 1;
                            }
                        }

                        mpc_area[0] = mpc_area[0]/ mpc_counter[0];
                        for ( IndexType i = 0; i < TDimension; i++ ) {
                            xg[0][i] = xg[0][i]/ mpc_counter[0];
                            
                            imposed_disp[0][i] = imposed_disp[0][i]/ mpc_counter[0];
                            normal[0][i] = normal[0][i]/ mpc_counter[0];
                        }
                        const bool is_slip = mpc->Is(SLIP);
                        const bool is_contact = mpc->Is(CONTACT);
                        const bool is_interface = mpc->Is(INTERFACE);

                        typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();
                        Element::Pointer pelem;
                        bool is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);

                        auto p_new_geometry = CreateQuadraturePointsUtility<Node<3>>::CreateFromCoordinates(
                                        pelem->pGetGeometry(), xg[0],
                                        mpc_area[0]);


                        Condition::Pointer p_condition = new_condition.Create(condition_id, p_new_geometry, properties);

                        p_condition->SetValuesOnIntegrationPoints(MPC_COORD, xg , rMPMModelPart.GetProcessInfo());
                        p_condition->SetValuesOnIntegrationPoints(MPC_AREA,  mpc_area  , rMPMModelPart.GetProcessInfo());
                        p_condition->SetValuesOnIntegrationPoints(MPC_NORMAL, normal, rMPMModelPart.GetProcessInfo());
                        p_condition->SetValuesOnIntegrationPoints(MPC_COUNTER, mpc_counter, rMPMModelPart.GetProcessInfo());

                        // p_condition->SetValuesOnIntegrationPoints(MPC_DISPLACEMENT, { mpc_displacement }, process_info);
                        p_condition->SetValuesOnIntegrationPoints(MPC_IMPOSED_DISPLACEMENT, imposed_disp, rMPMModelPart.GetProcessInfo());
                        p_condition->SetValuesOnIntegrationPoints(MPC_VELOCITY, { mpc_velocity }, rMPMModelPart.GetProcessInfo());
                        p_condition->SetValuesOnIntegrationPoints(MPC_IMPOSED_VELOCITY, { mpc_imposed_velocity }, rMPMModelPart.GetProcessInfo());
                        p_condition->SetValuesOnIntegrationPoints(MPC_ACCELERATION, { mpc_acceleration }, rMPMModelPart.GetProcessInfo());
                        p_condition->SetValuesOnIntegrationPoints(MPC_IMPOSED_ACCELERATION, { mpc_imposed_acceleration }, rMPMModelPart.GetProcessInfo());
                        p_condition->SetValuesOnIntegrationPoints(MPC_CONTACT_FORCE,  mpc_contact_force , rMPMModelPart.GetProcessInfo());

                        auto p_new_node = rBackgroundGridModelPart.CreateNewNode(rBackgroundGridModelPart.Nodes().size() + 1, xg[0][0], xg[0][1], xg[0][2]);
                        p_new_node->AddDof(VECTOR_LAGRANGE_MULTIPLIER_X,WEIGHTED_VECTOR_RESIDUAL_X);
                        p_new_node->AddDof(VECTOR_LAGRANGE_MULTIPLIER_Y,WEIGHTED_VECTOR_RESIDUAL_Y);
                        p_new_node->AddDof(VECTOR_LAGRANGE_MULTIPLIER_Z,WEIGHTED_VECTOR_RESIDUAL_Z);
                        p_new_node->AddDof(DISPLACEMENT_X,REACTION_X);
                        p_new_node->AddDof(DISPLACEMENT_Y,REACTION_Y);
                        p_new_node->AddDof(DISPLACEMENT_Z,REACTION_Z);

                        p_condition->SetValue(MPC_LAGRANGE_NODE, p_new_node);

                        if (is_slip)
                            p_condition->Set(SLIP);
                        if (is_contact)
                            p_condition->Set(CONTACT);
                        if (is_interface)
                            p_condition->Set(INTERFACE);

                        // Add the MP Condition to the model part
                        lagrange_model_part.AddCondition(p_condition);

                        condition_id+=1;
                    }
                }

            }
        }

    }
} // end namespace MPMSearchElementUtility

} // end namespace Kratos

#endif // KRATOS_MPM_SEARCH_ELEMENT_UTILITY

