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

namespace Kratos
{
namespace MPMSearchElementUtility
{

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef Kratos::GeometricalObject::GeometryType::Pointer MPElement;

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
                                                    MPElement& rMasterMaterialPoint,
        typename Geometry<Node<3>>::Pointer pGeometry,
        const array_1d<double, 3>& rCoordinates,
                                                    const std::size_t MaxNumberOfResults,
                                                    const double Tolerance)
    {
        // Get volume

        // Find bounds

        // Do splitting algorithm

        // Add insert sub points as quadrature points in current MP

        // Set elements to active
        KRATOS_TRY;

        PointerVector<Node<3>> nodes_list; // all nodes

        GeometryData::IntegrationMethod ThisDefaultMethod = pGeometry->GetDefaultIntegrationMethod();
        typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::IntegrationPointsArrayType ips(size);

        IntegrationPoint<3> int_p(x, y, z, w);
        ips.push_back(int_p);
        //ips[whateverposition] = int_p;
        Matrix N_matrix;
        Matrix DN_De;
        DenseVector<Matrix> DN_De_vector;
        for (int i + ...)
        {
            array_1d<double, 3> local_coordinates;
            pGeometry->PointLocalCoordinates(local_coordinates, rCoordinates);

            IntegrationPoint<3> int_p(local_coordinates, integration_weight);

            Vector N;
            pGeometry->ShapeFunctionsValues(N, local_coordinates);
            for (IndexType i = 0; i < N.size(); ++i)
            {
                N_matrix(integration_point_index, i) = N[i];
            }

            Matrix DN_De;
            pGeometry->ShapeFunctionsLocalGradients(DN_De, local_coordinates);
        }
        ips[0] = ThisIntegrationPoint;

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
} // end namespace MPMSearchElementUtility

} // end namespace Kratos

#endif // KRATOS_MPM_SEARCH_ELEMENT_UTILITY

