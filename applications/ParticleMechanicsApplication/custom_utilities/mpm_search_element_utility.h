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

namespace Kratos
{
namespace MPMSearchElementUtility
{

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef typename ModelPart::GeometryType GeometryType;

    inline double CrossProductDet2D(array_1d<double, 3> VectorA, array_1d<double, 3> VectorB)
    {
        return (VectorA[0] * VectorB[1] - VectorB[0] * VectorA[1]);
    }

    bool CheckIsInside(const GeometryType& rGeom, array_1d<double, 3>& LocalCoords, const array_1d<double, 3>& Coords, const double Tolerance)
    {
        bool is_inside = true;
        if (rGeom.Dimension() == 2)
        {
            is_inside = true;
            // Do walk around method
            Vector cross_products(rGeom.PointsNumber());
            for (size_t i = 0; i < rGeom.PointsNumber(); ++i)
            {
                if (rGeom.Points()[i].Coordinates()[2] != 0.0) {
                    return rGeom.IsInside(Coords, LocalCoords, Tolerance);
                    break;
                }
                cross_products[i] = CrossProductDet2D(Coords - rGeom.Points()[i].Coordinates(),
                    rGeom.Points()[(i+1)% rGeom.PointsNumber()].Coordinates()- rGeom.Points()[i].Coordinates());
            }
            for (size_t i = 1; i < cross_products.size(); ++i)
            {
                if (cross_products[i] * cross_products[0] < 0.0)
                {
                    is_inside = false;
                    break;
                }
            }

        }

        if (is_inside) return rGeom.IsInside(Coords, LocalCoords, Tolerance);

        return false;
    }

    void ConstructNeighbourRelations(ModelPart& rBackgroundGridModelPart)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rBackgroundGridModelPart.NumberOfElements()); ++i)
        {
            auto& r_geometry = (rBackgroundGridModelPart.ElementsBegin() + i)->GetGeometry();
            std::vector<typename Geometry<Node<3>>::Pointer> geometry_neighbours;
            for (IndexType j = 0; j < rBackgroundGridModelPart.NumberOfElements(); j++)
            {
                auto p_geometry_neighbour = (rBackgroundGridModelPart.ElementsBegin() + j)->pGetGeometry();
                for (IndexType n = 0; n < p_geometry_neighbour->size(); n++)
                {
                    for (IndexType k = 0; k < r_geometry.size(); k++)
                    {
                        if (r_geometry[k].Id() == (*p_geometry_neighbour)[n].Id()) {
                            geometry_neighbours.push_back(p_geometry_neighbour);
                            break;
                        }
                    }
                }
            }
            r_geometry.SetValue(GEOMETRY_NEIGHBOURS, geometry_neighbours);
        }
    }

    bool IsExplicitAndNeedsCorrection(GeometryType::Pointer pQuadraturePoint, const ProcessInfo& rProcessInfo)
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

    void NeighbourSearchElements(ModelPart& rMPMModelPart, std::vector<typename Element::Pointer>& rMissingElements,
        const double Tolerance, const ProcessInfo& rProcessInfo)
    {
        #pragma omp for
        for (int i = 0; i < static_cast<int>(rMPMModelPart.Elements().size()); ++i) {
            auto element_itr = rMPMModelPart.Elements().begin() + i;

            std::vector<array_1d<double, 3>> xg;
            element_itr->CalculateOnIntegrationPoints(MP_COORD, xg, rMPMModelPart.GetProcessInfo());

            GeometryType& r_parent_geometry = element_itr->GetGeometry().GetGeometryParent(0);

            array_1d<double, 3> local_coordinates;
            bool is_found = CheckIsInside(r_parent_geometry,local_coordinates,xg[0],Tolerance);
            if (!is_found)
            {
                auto& geometry_neighbours = r_parent_geometry.GetValue(GEOMETRY_NEIGHBOURS);
                for (IndexType k = 0; k < geometry_neighbours.size(); k++)
                {
                    if (CheckIsInside(*geometry_neighbours[k], local_coordinates, xg[0], Tolerance))
                    {
                        auto p_new_geometry = CreateQuadraturePointsUtility<Node<3>>::CreateFromLocalCoordinates(
                            *(geometry_neighbours[k].get()), local_coordinates,
                            element_itr->GetGeometry().IntegrationPoints()[0].Weight());

                        if (!IsExplicitAndNeedsCorrection(p_new_geometry, rProcessInfo)) {
                            is_found = true;
                            // Update geometry of particle element
                            element_itr->SetGeometry(p_new_geometry);

                            for (IndexType j = 0; j < geometry_neighbours[k]->PointsNumber(); ++j)
                                geometry_neighbours[k]->Points()[j].Set(ACTIVE);
                            break;
                        }
                    }
                }
            }
            else {
                //pelem->Set(ACTIVE);

                auto p_new_geometry = CreateQuadraturePointsUtility<Node<3>>::CreateFromLocalCoordinates(
                    r_parent_geometry, local_coordinates,
                    element_itr->GetGeometry().IntegrationPoints()[0].Weight());

                if (IsExplicitAndNeedsCorrection(p_new_geometry, rProcessInfo)) is_found = false;
                else
                {
                    // Update geometry of particle element
                    element_itr->SetGeometry(p_new_geometry);

                    for (IndexType j = 0; j < r_parent_geometry.PointsNumber(); ++j)
                        r_parent_geometry[j].Set(ACTIVE);
                }
            }
            if (!is_found) {
                #pragma omp critical
                rMissingElements.push_back(&*element_itr);
            }
        }
    }

    void NeighbourSearchConditions(ModelPart& rMPMModelPart, std::vector<typename Condition::Pointer>& rMissingConditions,
        const double Tolerance, const ProcessInfo& rProcessInfo)
    {
        #pragma omp for
        for (int i = 0; i < static_cast<int>(rMPMModelPart.Conditions().size()); ++i) {
            auto condition_itr = rMPMModelPart.Conditions().begin() + i;

            std::vector<array_1d<double, 3>> xg;
            condition_itr->CalculateOnIntegrationPoints(MPC_COORD, xg, rMPMModelPart.GetProcessInfo());
            if (xg.size() > 0)
            {
                GeometryType& r_parent_geometry = condition_itr->GetGeometry();

                array_1d<double, 3> local_coordinates;
                //bool is_found = r_parent_geometry.IsInside(xg[0], local_coordinates, Tolerance);
                bool is_found = CheckIsInside(r_parent_geometry, local_coordinates, xg[0], Tolerance);
                if (!is_found)
                {
                    auto& geometry_neighbours = r_parent_geometry.GetValue(GEOMETRY_NEIGHBOURS);
                    for (IndexType k = 0; k < geometry_neighbours.size(); k++)
                    {
                        //if (geometry_neighbours[k]->IsInside(xg[0], local_coordinates, Tolerance))
                        if (CheckIsInside(*geometry_neighbours[k], local_coordinates, xg[0], Tolerance))
                        {
                            is_found = true;

                            condition_itr->GetGeometry() = *(geometry_neighbours[k].get());

                            for (IndexType j = 0; j < geometry_neighbours[k]->PointsNumber(); ++j)
                                geometry_neighbours[k]->Points()[j].Set(ACTIVE);
                            break;
                        }
                    }
                }
                else {
                    //pelem->Set(ACTIVE);

                    condition_itr->GetGeometry() = r_parent_geometry;

                    for (IndexType j = 0; j < r_parent_geometry.PointsNumber(); ++j)
                        r_parent_geometry[j].Set(ACTIVE);
                }
                if (!is_found) {
                    #pragma omp critical
                    rMissingConditions.push_back(&*condition_itr);
                }
            }
        }
    }

    template <std::size_t TDimension>
    void BinBasedSearchElementsAndConditions(ModelPart& rMPMModelPart,
        ModelPart& rBackgroundGridModelPart,
        std::vector<typename Element::Pointer>& rMissingElements,
        std::vector<typename Condition::Pointer>& rMissingConditions,
        const std::size_t MaxNumberOfResults, const double Tolerance)
    {
        const ProcessInfo& r_process_info = rBackgroundGridModelPart.GetProcessInfo();
        const bool is_fix_explicit_mp_on_grid_edge = (r_process_info.Has(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE))
            ? r_process_info.GetValue(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)
            : false;

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


                if (is_found && is_fix_explicit_mp_on_grid_edge) {
                    // check if MP is exactly on the edge of the element, this gives spurious strains in explicit
                    bool isOnEdge = false;
                    for (SizeType i = 0; i < N.size(); ++i) {
                        if (std::abs(N[i]) < std::numeric_limits<double>::epsilon()) {
                            isOnEdge = true;
                            break;
                        }
                    }
                    if (isOnEdge) {
                        // MP is exactly on the edge. Now we give it a little 'nudge'
                        array_1d<double, 3> xg_nudged = array_1d<double, 3>(xg[0]);
                        const double& delta_time = r_process_info[DELTA_TIME];
                        std::vector<array_1d<double, 3>> mp_vel;
                        element_itr->CalculateOnIntegrationPoints(MP_VELOCITY, mp_vel, rMPMModelPart.GetProcessInfo());
                        array_1d<double, 3> nudge_displacement = delta_time / 1000.0 * mp_vel[0];
                        xg_nudged += nudge_displacement;
                        is_found = SearchStructure.FindPointOnMesh(xg_nudged, N, pelem, result_begin, MaxNumberOfResults, Tolerance);
                        // check if the nudged point is found...
                        if (is_found) {
                            // store the nudged MP position
                            element_itr->SetValuesOnIntegrationPoints(MP_COORD, { xg_nudged }, rMPMModelPart.GetProcessInfo());
                            KRATOS_INFO("MPMSearchElementUtility") << "WARNING: To prevent spurious explicit stresses, Material Point " << element_itr->Id()
                                << " was nudged by " << nudge_displacement << std::endl;
                        }
                        else {
                            // find the un-nudged MP again
                            is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);
                            KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Material Point " << element_itr->Id()
                                << " lies exactly on an element edge and may give spurious results." << std::endl;
                        }
                    }
                }


                if (is_found == true) {
                    pelem->Set(ACTIVE);

                    auto p_new_geometry = CreateQuadraturePointsUtility<Node<3>>::CreateFromCoordinates(
                        pelem->pGetGeometry(), xg[0],
                        element_itr->GetGeometry().IntegrationPoints()[0].Weight());

                    // Update geometry of particle element
                    element_itr->SetGeometry(p_new_geometry);

                    for (IndexType j = 0; j < p_new_geometry->PointsNumber(); ++j)
                        (*p_new_geometry)[j].Set(ACTIVE);
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
                        //pelem->Set(ACTIVE);
                        condition_itr->GetGeometry() = pelem->GetGeometry();
                        auto& r_geometry = condition_itr->GetGeometry();

                        for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j)
                            r_geometry[j].Set(ACTIVE);
                    }
                    else {
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


    void ResetElementsAndNodes(ModelPart& rBackgroundGridModelPart)
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
        const ProcessInfo& r_process_info = rBackgroundGridModelPart.GetProcessInfo();

        ResetElementsAndNodes(rBackgroundGridModelPart);

        const bool use_neighbour_search = true; //TODO delete - for testing only

        std::vector<typename Element::Pointer> missing_elements;
        std::vector<typename Condition::Pointer> missing_conditions;

        if (use_neighbour_search)
        {

            if (!rBackgroundGridModelPart.ElementsBegin()->GetGeometry().Has(GEOMETRY_NEIGHBOURS))
                ConstructNeighbourRelations(rBackgroundGridModelPart);

            NeighbourSearchElements(rMPMModelPart, missing_elements, Tolerance, r_process_info);
            NeighbourSearchConditions(rMPMModelPart, missing_conditions, Tolerance, r_process_info);
        }
        else
        {
            missing_elements.resize(rMPMModelPart.NumberOfElements());
            for (size_t i = 0; i < rMPMModelPart.NumberOfElements(); ++i)
            {
                missing_elements[i] = &*(rMPMModelPart.Elements().begin() + i);
            }

            missing_conditions.resize(rMPMModelPart.NumberOfConditions());
            for (size_t i = 0; i < rMPMModelPart.NumberOfConditions(); ++i)
            {
                missing_conditions[i] = &*(rMPMModelPart.Conditions().begin() + i);
            }
        }

        if (missing_conditions.size() > 0 || missing_elements.size() > 0)
            BinBasedSearchElementsAndConditions<TDimension>(rMPMModelPart,
                rBackgroundGridModelPart, missing_elements, missing_conditions,
                MaxNumberOfResults, Tolerance);

    }
} // end namespace MPMSearchElementUtility

} // end namespace Kratos

#endif // KRATOS_MPM_SEARCH_ELEMENT_UTILITY

