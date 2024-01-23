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


#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/quadrature_points_utility.h"
#include "mpm_application_variables.h"
#include "geometries/geometry.h"
#include "includes/model_part.h"
#include "pqmpm_partition_utilities.h"

namespace Kratos::MPMSearchElementUtility
{
    // Standard types
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef Node NodeType;
    typedef typename ModelPart::GeometryType GeometryType;


    inline double CrossProductDet2D(
        array_1d<double, 3> VectorA,
        array_1d<double, 3> VectorB
        )
    {
        return (VectorA[0] * VectorB[1] - VectorB[0] * VectorA[1]);
    }


    inline bool CheckIsInside(
        const GeometryType& rGeom,
        array_1d<double, 3>& LocalCoords,
        const array_1d<double, 3>& Coords,
        const double Tolerance
        )
    {
        // TODO some optimisation for simple 2D shapes.
        return rGeom.IsInside(Coords, LocalCoords, Tolerance);
    }


    inline void ConstructNeighbourRelations(
        GeometryType& rGeom,
        const ModelPart& rBackgroundGridModelPart
        )
    {
        std::vector<typename Geometry<Node>::Pointer> geometry_neighbours;
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
                            if (add_entry) {
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


    inline bool IsExplicitAndNeedsCorrection(
        GeometryType::Pointer pQuadraturePoint,
        const ProcessInfo& rProcessInfo
        )
    {
        if (rProcessInfo.Has(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)) {
            if (rProcessInfo.GetValue(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)) {
                if (pQuadraturePoint->IntegrationPointsNumber() == 1) {
                    for (size_t i = 0; i < pQuadraturePoint->ShapeFunctionsValues().size2(); ++i) {
                        if (pQuadraturePoint->ShapeFunctionsValues()(0, i) < std::numeric_limits<double>::epsilon()) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    inline GeometryType& FindGridGeom(
        GeometryType& rParentGeom,
        const ModelPart& rBackgroundGridModelPart,
        const double Tolerance,
        const array_1d<double, 3>& xg,
        array_1d<double, 3>& rLocalCoords,
        const ProcessInfo& rProcessInfo,
        bool& IsFound
        )
    {
        IsFound = false;
        if (CheckIsInside(rParentGeom, rLocalCoords, xg, Tolerance)) {
            IsFound = true;
            return rParentGeom;
        } else {
            if (!rParentGeom.Has(GEOMETRY_NEIGHBOURS)) {
                ConstructNeighbourRelations(rParentGeom, rBackgroundGridModelPart);
            }
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


    inline void UpdatePartitionedQuadraturePoint(
        const ModelPart& rBackgroundGridModelPart,
        const array_1d<double, 3>& rCoordinates,
        Element& rMasterMaterialPoint,
        typename GeometryType::Pointer pQuadraturePointGeometry,
        const double Tolerance
        )
    {
        KRATOS_TRY;

        array_1d<double, 3> local_coords;
        pQuadraturePointGeometry->IsInside(rCoordinates, local_coords, Tolerance);
        PQMPMPartitionUtilities::PartitionMasterMaterialPointsIntoSubPoints(rBackgroundGridModelPart, rCoordinates,
            local_coords, rMasterMaterialPoint, pQuadraturePointGeometry, Tolerance);

        KRATOS_CATCH("");
    }


    inline void NeighbourSearchElements(
        const ModelPart& rMPMModelPart,
        const ModelPart& rBackgroundGridModelPart,
        std::vector<typename Element::Pointer>& rMissingElements,
        const double Tolerance
        )
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

            if (is_found) {
                const bool is_pqmpm = (rBackgroundGridModelPart.GetProcessInfo().Has(IS_PQMPM))
                    ? rBackgroundGridModelPart.GetProcessInfo().GetValue(IS_PQMPM) : false;
                if (is_pqmpm) {
                    // Updates the quadrature point geometry.
                    (*element_itr).GetGeometry().SetGeometryParent(&r_found_geom);
                    PQMPMPartitionUtilities::PartitionMasterMaterialPointsIntoSubPoints(rBackgroundGridModelPart, xg[0],
                        local_coordinates, *element_itr, element_itr->pGetGeometry(), Tolerance);
                } else {
                    CreateQuadraturePointsUtility<Node>::UpdateFromLocalCoordinates(
                        element_itr->pGetGeometry(), local_coordinates,
                        element_itr->GetGeometry().IntegrationPoints()[0].Weight(), r_found_geom);
                }
                if (IsExplicitAndNeedsCorrection(element_itr->pGetGeometry(), rBackgroundGridModelPart.GetProcessInfo())) {
                    is_found = false;
                } else {
                    for (IndexType j = 0; j < r_found_geom.PointsNumber(); ++j) {
                        r_found_geom.Points()[j].Set(ACTIVE);
                    }
                }
            } else {
                #pragma omp critical
                rMissingElements.push_back(&*element_itr);
            }
        }
    }


    inline void NeighbourSearchConditions(
        const ModelPart& rMPMModelPart,
        const ModelPart& rBackgroundGridModelPart,
        std::vector<typename Condition::Pointer>& rMissingConditions,
        const double Tolerance
        )
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rMPMModelPart.Conditions().size()); ++i) {
            auto condition_itr = rMPMModelPart.Conditions().begin() + i;
            std::vector<array_1d<double, 3>> xg;
            condition_itr->CalculateOnIntegrationPoints(MPC_COORD, xg, rMPMModelPart.GetProcessInfo());

            if (xg.size() > 0 && condition_itr->Is(BOUNDARY)) {
                array_1d<double, 3> local_coordinates;
                bool is_found = false;

                GeometryType& r_found_geom = FindGridGeom(condition_itr->GetGeometry().GetGeometryParent(0),
                    rBackgroundGridModelPart, Tolerance, xg[0], local_coordinates,
                    rMPMModelPart.GetProcessInfo(), is_found);

                if (is_found) {
                    CreateQuadraturePointsUtility<Node>::UpdateFromLocalCoordinates(
                        condition_itr->pGetGeometry(), local_coordinates,
                        condition_itr->GetGeometry().IntegrationPoints()[0].Weight(), r_found_geom);

                    for (IndexType j = 0; j < r_found_geom.PointsNumber(); ++j) {
                        r_found_geom[j].Set(ACTIVE);
                    }
                } else {
                    #pragma omp critical
                    rMissingConditions.push_back(&*condition_itr);
                }
            }
        }
    }


    inline bool IsFixExplicitAndOnElementEdge(
        const Vector& N,
        const ProcessInfo& rProcessInfo
        )
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
    void BinBasedSearchElementsAndConditions(
        ModelPart& rMPMModelPart,
        ModelPart& rBackgroundGridModelPart,
        std::vector<typename Element::Pointer>& rMissingElements,
        std::vector<typename Condition::Pointer>& rMissingConditions,
        const std::size_t MaxNumberOfResults,
        const double Tolerance
        )
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

                if (is_found) {
                    pelem->Set(ACTIVE);

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
                    if (is_pqmpm) {
                        // Updates the quadrature point geometry.
                        (*element_itr).GetGeometry().SetGeometryParent((pelem->pGetGeometry().get()));
                        UpdatePartitionedQuadraturePoint(rBackgroundGridModelPart, xg[0],
                            *element_itr, pelem->pGetGeometry(), Tolerance);
                    } else {
                        auto p_quadrature_point_geometry = element_itr->pGetGeometry();
                        array_1d<double, 3> local_coordinates;
                        pelem->pGetGeometry()->PointLocalCoordinates(local_coordinates, xg[0]);
                        CreateQuadraturePointsUtility<Node>::UpdateFromLocalCoordinates(
                            p_quadrature_point_geometry, local_coordinates,
                            p_quadrature_point_geometry->IntegrationPoints()[0].Weight(), pelem->GetGeometry());
                    }
                    auto& r_geometry = element_itr->GetGeometry();
                    for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j) {
                        r_geometry[j].Set(ACTIVE);
                    }
                } else {
                    KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Search Element for Material Point: "
                        << element_itr->Id() << " is failed. Geometry is cleared." << std::endl;
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

                    if (is_found) {
                        auto p_quadrature_point_geometry = condition_itr->pGetGeometry();
                        array_1d<double, 3> local_coordinates;
                        pelem->pGetGeometry()->PointLocalCoordinates(local_coordinates, xg[0]);
                        CreateQuadraturePointsUtility<Node>::UpdateFromLocalCoordinates(
                            p_quadrature_point_geometry, local_coordinates,
                            p_quadrature_point_geometry->IntegrationPoints()[0].Weight(), pelem->GetGeometry());

                        auto& r_geometry = condition_itr->GetGeometry();

                        for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j) {
                            r_geometry[j].Set(ACTIVE);
                        }
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
            element_itr->Reset(ACTIVE);
            auto& r_geometry = element_itr->GetGeometry();
            for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j) {
                r_geometry[j].Reset(ACTIVE);
            }
        }
    }


    /**
     * @brief Search element connectivity for each material point
     * @details A search is performed to know in which grid element the material point falls.
     * If one or more material points fall in the grid element, the grid element is
     * set to be active and its connectivity is associated to the material point
     * element.
     * STEPS:
     * 1) All the elements are set to be INACTIVE
     * 2) A searching is performed and the grid elements which contain at least a MP are set to be ACTIVE
     */
    template<std::size_t TDimension>
    void SearchElement(
        ModelPart& rBackgroundGridModelPart,
        ModelPart& rMPMModelPart,
        const std::size_t MaxNumberOfResults,
        const double Tolerance
        )
    {
        ResetElementsAndNodes(rBackgroundGridModelPart);

        std::vector<typename Element::Pointer> missing_elements;
        std::vector<typename Condition::Pointer> missing_conditions;

        if (!rMPMModelPart.GetProcessInfo()[IS_RESTARTED]) {
            NeighbourSearchElements(rMPMModelPart, rBackgroundGridModelPart, missing_elements, Tolerance);
            NeighbourSearchConditions(rMPMModelPart, rBackgroundGridModelPart, missing_conditions, Tolerance);
        } else {
            missing_elements.resize(rMPMModelPart.Elements().size());
            IndexPartition(rMPMModelPart.Elements().size()).for_each([&](std::size_t i){
                missing_elements[i] = &*(rMPMModelPart.ElementsBegin() + i); // maybe add/remove *&?
            });

            missing_conditions.resize(rMPMModelPart.Conditions().size());
            IndexPartition(rMPMModelPart.Conditions().size()).for_each([&](std::size_t i){
                missing_conditions[i] = &*(rMPMModelPart.Conditions().begin() + i); // maybe add/remove *&?
            });
        }

        if (missing_conditions.size() > 0 || missing_elements.size() > 0) {
            BinBasedSearchElementsAndConditions<TDimension>(rMPMModelPart,
                rBackgroundGridModelPart, missing_elements, missing_conditions,
                MaxNumberOfResults, Tolerance);
        }
    }
} // end namespace Kratos::MPMSearchElementUtility
