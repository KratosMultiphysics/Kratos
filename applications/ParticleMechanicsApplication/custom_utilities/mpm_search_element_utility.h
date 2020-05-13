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
#include "geometries/geometry_shape_function_container.h"
#include "custom_geometries/quadrature_point_partitioned_geometry.h"

#include "geometries/quadrilateral_2d_4.h"
#include "geometries/bounding_box.h"

namespace Kratos
{
namespace MPMSearchElementUtility
{
    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef GeometricalObject::GeometryType GeometryType;

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
        const bool is_explicit = (r_process_info.Has(IS_EXPLICIT))
            ? r_process_info.GetValue(IS_EXPLICIT)
            : false;
        const bool is_pqmpm = (r_process_info.Has(IS_PQMPM))
            ? r_process_info.GetValue(IS_PQMPM)
            : false;

        // Reset elements to inactive
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rBackgroundGridModelPart.Elements().size()); ++i) {
            auto element_itr = rBackgroundGridModelPart.Elements().begin() + i;
            auto& r_geometry = element_itr->GetGeometry();
            element_itr->Reset(ACTIVE);

            for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j)
                r_geometry[j].Reset(ACTIVE);

        }
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
            for (int i = 0; i < static_cast<int>(rMPMModelPart.Elements().size()); ++i) {
                auto element_itr = rMPMModelPart.Elements().begin() + i;

                std::vector<array_1d<double, 3>> xg;
                element_itr->CalculateOnIntegrationPoints(MP_COORD, xg, rMPMModelPart.GetProcessInfo());
                typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();

                Element::Pointer pelem;

                // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                bool is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);


                if (is_found && is_explicit) {
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
                        if (is_found){
                            // store the nudged MP position
                            element_itr->SetValuesOnIntegrationPoints(MP_COORD, { xg_nudged }, rMPMModelPart.GetProcessInfo());
                            KRATOS_INFO("MPMSearchElementUtility") << "WARNING: To prevent spurious explicit stresses, Material Point " << element_itr->Id()
                                << " was nudged by " << nudge_displacement << std::endl;
                        }
                        else {
                            // find the un-nudged MP again
                            is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);
                            KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Material Point " << element_itr->Id()
                                << " lies exactly on an element edge and may give spurious results."<< std::endl;
                        }
                    }
                }


                if (is_found == true) {
                    pelem->Set(ACTIVE);

                    // location: xg[0]
                    // element_itr->GetGeometry().IntegrationPoints()[0].Weight() instead element_itr->GetValue(MP_VOLUME)
                    // pelem->pGetGeometry()

                    auto p_new_geometry = (is_pqmpm)
                        ? PartitionMasterMaterialPointsIntoSubPoints(
                            rBackgroundGridModelPart, xg[0], *element_itr, pelem->pGetGeometry(), MaxNumberOfResults, Tolerance)
                        : CreateQuadraturePointsUtility<Node<3>>::CreateFromCoordinates(
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
            for (int i = 0; i < static_cast<int>(rMPMModelPart.Conditions().size()); ++i) {

                auto condition_itr = rMPMModelPart.Conditions().begin() + i;
                std::vector<array_1d<double, 3>> xg;
                condition_itr->CalculateOnIntegrationPoints(MPC_COORD, xg, rMPMModelPart.GetProcessInfo());

                if (xg.size() == 1) {
                    // Only search for particle based BCs!
                    // Grid BCs are still applied on MP_model_part but we don't want to search for them.
                    typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();

                    Element::Pointer pelem;

                    // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                    bool is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);

                    if (is_found == true) {
                        pelem->Set(ACTIVE);
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

    typename Geometry<Node<3>>::Pointer PartitionMasterMaterialPointsIntoSubPoints(ModelPart& rBackgroundGridModelPart,
                                                    const array_1d<double, 3>& rCoordinates,
                                                    Element& rMasterMaterialPoint,
                                                    typename Geometry<Node<3>>::Pointer pGeometry,
                                                    const std::size_t MaxNumberOfResults,
                                                    const double Tolerance)
    {
        KRATOS_TRY;

        // TODO should we be attaching the each background grid element to each subpoint?

        const SizeType working_dim = pGeometry->WorkingSpaceDimension();
        KRATOS_ERROR_IF(working_dim > 2) << "PQMPM is currently limited to 2D!" << std::endl;
        const SizeType safety_buffer = 5;


        // Get volume and set up master domain bounding square
        std::vector<double> mp_volume;
        rMasterMaterialPoint.CalculateOnIntegrationPoints(MP_VOLUME, mp_volume, rBackgroundGridModelPart.GetProcessInfo());
        const double side_length = std::pow(mp_volume[0], double(1.0 / working_dim));
        const SizeType n_bounding_box_vertices = std::pow(2, working_dim);

        Point point_low(rCoordinates[0] - side_length, rCoordinates[1] - side_length, rCoordinates[2]);
        Point point_high(rCoordinates[0] + side_length, rCoordinates[1] + side_length, rCoordinates[2]);

        // Determine what elements we intersect
        // TODO try to reduce this search more with initial binning
        std::vector<Element&> intersected_elements;
        auto element_begin = rBackgroundGridModelPart.ElementsBegin();
        for (IndexType i = 0; i < rBackgroundGridModelPart.Elements().size(); ++i) {
            auto ele_it = element_begin + i;
            if (ele_it->GetGeometry().HasIntersection(point_low, point_high)) intersected_elements.push_back(*ele_it);
        }

        // Prepare containers
        const SizeType number_of_subpoints = intersected_elements.size();

        PointerVector<Node<3>> nodes_list(number_of_subpoints * element_begin->GetGeometry().PointsNumber());
        typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::IntegrationPointsArrayType ips(number_of_subpoints);
        Matrix N_matrix (number_of_subpoints, element_begin->GetGeometry().PointsNumber(),0.0);
        DenseVector<Matrix> DN_De_vector (number_of_subpoints, ZeroMatrix(number_of_subpoints,working_dim));

        // Temporary local containers
        Vector N(number_of_subpoints);
        Matrix DN_De (number_of_subpoints,working_dim);


        if (number_of_subpoints == 1) {
            // bounding box must be within entire element - same as normal case
            ips[0] = CreateSubPoint(rCoordinates, 1.0, intersected_elements[0].GetGeometry(),N,DN_De);

            for (size_t i = 0; i < N.size(); ++i) {
                N_matrix(0, i) = N[i];
                nodes_list[i] = intersected_elements[0].GetGeometry().Points()[i];
            }
            DN_De_vector[0] = DN_De;
        }
        else
        {
            // Setup master domain as a quad geom to use core geom tests
            Point::Pointer p1(new Point(rCoordinates[0] - side_length, rCoordinates[1] - side_length, rCoordinates[2]));
            Point::Pointer p2(new Point(rCoordinates[0] + side_length, rCoordinates[1] - side_length, rCoordinates[2]));
            Point::Pointer p3(new Point(rCoordinates[0] + side_length, rCoordinates[1] + side_length, rCoordinates[2]));
            Point::Pointer p4(new Point(rCoordinates[0] - side_length, rCoordinates[1] + side_length, rCoordinates[2]));
            Quadrilateral2D4<Point> master_domain (p1,p2,p3,p4);

            IndexType node_index = 0;

            for (size_t i = 0; i < number_of_subpoints; ++i)
            {
                if (CheckGeometryIsCompletelyWithinAnother(intersected_elements[i].GetGeometry(), master_domain)) { 
                    // whole element is completely inside bounding box
                    ips[i] = CreateSubPoint(intersected_elements[i].GetGeometry().Center(), 
                        intersected_elements[i].GetGeometry().DomainSize()/mp_volume[0], 
                        intersected_elements[i].GetGeometry(),N,DN_De);
                }
                else  {
                    // only some of the background element is within the bounding box

                }

                // Transfer local data to containers
                for (size_t j = 0; j < N.size(); ++j) {
                    N_matrix(i, node_index) = N[j];
                    nodes_list[node_index] = intersected_elements[0].GetGeometry().Points()[j];
                    node_index += 1;
                }
                DN_De_vector[i] = DN_De;
            }
        }
        // Set elements to active


        GeometryData::IntegrationMethod ThisDefaultMethod = pGeometry->GetDefaultIntegrationMethod();
        typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::IntegrationPointsContainerType ips_container;
        ips_container[ThisDefaultMethod] = ips;
        typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::ShapeFunctionsValuesContainerType shape_function_container;
        shape_function_container[ThisDefaultMethod] = N_matrix;
        typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::ShapeFunctionsLocalGradientsContainerType shape_function_derivatives_container;
        shape_function_derivatives_container[ThisDefaultMethod] = DN_De_vector;

        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
            ThisDefaultMethod,
            ips_container,
            shape_function_container,
            shape_function_derivatives_container);

        return CreateQuadraturePoint(
            pGeometry->WorkingSpaceDimension(),
            pGeometry->LocalSpaceDimension(),
            data_container,
            nodes_list);

        KRATOS_CATCH("");
    }

    //template<class TPointType>
    static typename Geometry<Node<3>>::Pointer CreateQuadraturePoint(
        SizeType WorkingSpaceDimension,
        SizeType LocalSpaceDimension,
        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainer,
        typename Geometry<Node<3>>::PointsArrayType rPoints)
    {
        if (WorkingSpaceDimension == 1 && LocalSpaceDimension == 1)
            return Kratos::make_shared<
            QuadraturePointPartitionedGeometry<Node<3>, 1>>(
                rPoints,
                rShapeFunctionContainer);
        else if (WorkingSpaceDimension == 2 && LocalSpaceDimension == 1)
            return Kratos::make_shared<
            QuadraturePointPartitionedGeometry<Node<3>, 2, 1>>(
                rPoints,
                rShapeFunctionContainer);
        else if (WorkingSpaceDimension == 2 && LocalSpaceDimension == 2)
            return Kratos::make_shared<
            QuadraturePointPartitionedGeometry<Node<3>, 2>>(
                rPoints,
                rShapeFunctionContainer);
        else if (WorkingSpaceDimension == 3 && LocalSpaceDimension == 2)
            return Kratos::make_shared<
            QuadraturePointPartitionedGeometry<Node<3>, 3, 2>>(
                rPoints,
                rShapeFunctionContainer);
        else if (WorkingSpaceDimension == 3 && LocalSpaceDimension == 3)
            return Kratos::make_shared<
            QuadraturePointPartitionedGeometry<Node<3>, 3>>(
                rPoints,
                rShapeFunctionContainer);
        else {
            KRATOS_ERROR << "Working/Local space dimension combinations are "
                << "not provided for QuadraturePointGeometry. WorkingSpaceDimension: "
                << WorkingSpaceDimension << ", LocalSpaceDimension: " << LocalSpaceDimension
                << std::endl;
        }
    }


    const bool CheckGeometryIsCompletelyWithinAnother(const GeometryType& rTestGeom, const GeometryType& rReferenceGeom)
    {
        array_1d<double, 3> local_coords;
        for (size_t i = 0; i < rTestGeom.PointsNumber(); i++) {
            if (!rReferenceGeom.IsInside(rTestGeom.GetPoint(i).Coordinates(), local_coords)) return false;
        }

        return true;
    }

    IntegrationPoint<3> CreateSubPoint(const array_1d<double, 3>& rGlobalCoords, const double rVolumeFraction, 
        const GeometryType& rBackgroundGridElementGeom, Vector& rN, Matrix& rDN_De)
    {
        array_1d<double, 3> local_coordinates;
        rBackgroundGridElementGeom.PointLocalCoordinates(local_coordinates, rGlobalCoords);
        rBackgroundGridElementGeom.ShapeFunctionsValues(rN, local_coordinates);
        rBackgroundGridElementGeom.ShapeFunctionsLocalGradients(rDN_De, local_coordinates);

        return IntegrationPoint<3>(local_coordinates, rVolumeFraction);
    }
} // end namespace MPMSearchElementUtility

} // end namespace Kratos

#endif // KRATOS_MPM_SEARCH_ELEMENT_UTILITY

