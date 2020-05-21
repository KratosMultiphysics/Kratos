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

// System includes

// External includes

// Project includes
#include "custom_utilities/mpm_search_element_utility.h"

namespace Kratos
{
    void MPMSearchElementUtility::CreateBoundingBoxPoints(std::vector<array_1d<double, 3>>& rPointVector,
        const array_1d<double, 3>& rCenter, const double SideHalfLength, const SizeType WorkingDim)
    {
        KRATOS_TRY

            if (WorkingDim == 2)
            {
                if (rPointVector.size() != 4) rPointVector.resize(4);
                for (size_t i = 0; i < 4; ++i) {
                    rPointVector[i].clear();
                    rPointVector[i] += rCenter;
                }
                rPointVector[0][0] -= SideHalfLength;
                rPointVector[1][0] += SideHalfLength;
                rPointVector[2][0] += SideHalfLength;
                rPointVector[3][0] -= SideHalfLength;

                rPointVector[0][1] -= SideHalfLength;
                rPointVector[1][1] -= SideHalfLength;
                rPointVector[2][1] += SideHalfLength;
                rPointVector[3][1] += SideHalfLength;
            }
            else
            {
                if (rPointVector.size() != 8) rPointVector.resize(8);
                for (size_t i = 0; i < 8; ++i) {
                    rPointVector[i].clear();
                    rPointVector[i] += rCenter;
                }
                rPointVector[0][0] -= SideHalfLength;
                rPointVector[1][0] += SideHalfLength;
                rPointVector[2][0] += SideHalfLength;
                rPointVector[3][0] -= SideHalfLength;
                rPointVector[4][0] -= SideHalfLength;
                rPointVector[5][0] += SideHalfLength;
                rPointVector[6][0] += SideHalfLength;
                rPointVector[7][0] -= SideHalfLength;

                rPointVector[0][1] -= SideHalfLength;
                rPointVector[1][1] -= SideHalfLength;
                rPointVector[2][1] += SideHalfLength;
                rPointVector[3][1] += SideHalfLength;
                rPointVector[4][1] -= SideHalfLength;
                rPointVector[5][1] -= SideHalfLength;
                rPointVector[6][1] += SideHalfLength;
                rPointVector[7][1] += SideHalfLength;

                rPointVector[0][2] -= SideHalfLength;
                rPointVector[1][2] -= SideHalfLength;
                rPointVector[2][2] -= SideHalfLength;
                rPointVector[3][2] -= SideHalfLength;
                rPointVector[4][2] += SideHalfLength;
                rPointVector[5][2] += SideHalfLength;
                rPointVector[6][2] += SideHalfLength;
                rPointVector[7][2] += SideHalfLength;
            }

        KRATOS_CATCH("")
    }


    void MPMSearchElementUtility::Check(IntegrationPointsArrayType& rIntergrationSubPoints, const double Tolerance, const Matrix& rN, const DenseVector<Matrix>& rDN_De)
    {
        KRATOS_TRY

            double vol_frac_accum = 0.0;

        KRATOS_ERROR_IF(rIntergrationSubPoints.size() != rN.size1())
            << "Shape function rows must equal number of sub-points!";

        for (size_t i = 0; i < rIntergrationSubPoints.size(); i++)
        {
            KRATOS_ERROR_IF(rIntergrationSubPoints[i].Weight() < Tolerance)
                << "Volume fraction of sub-points is too small!";

            KRATOS_ERROR_IF(rIntergrationSubPoints[i].Weight() > 1.0)
                << "Volume fraction of sub-points is too large!";

            vol_frac_accum += rIntergrationSubPoints[i].Weight();
        }

        KRATOS_ERROR_IF(vol_frac_accum < (1.0 - rIntergrationSubPoints.size() * Tolerance))
            << "Volume fraction of sub-points does not approximately sum to 1.0."
            << " This probably means the background grid is not big enough";

        for (size_t j = 0; j < rN.size2(); j++)
        {
            SizeType nonzero_entries = 0;
            for (size_t i = 0; i < rIntergrationSubPoints.size(); i++) if (rN(i, j) != 0.0) nonzero_entries += 1;
            KRATOS_ERROR_IF(nonzero_entries != 1)
                << "There must be only one nonzero entry per shape function column!";
        }

        KRATOS_CATCH("")
    }


    const bool MPMSearchElementUtility::CheckAllPointsAreInGeom(
        const std::vector<array_1d<double, 3>>& rPoints,
        const GeometryType& rReferenceGeom,
        const double Tolerance)
    {
        KRATOS_TRY

            array_1d<double, 3> dummy_local_coords;
        bool is_coincident;
        for (size_t i = 0; i < rPoints.size(); ++i) {
            if (!rReferenceGeom.IsInside(rPoints[i], dummy_local_coords), Tolerance) {
                // the test point may directly lie on one of the ref geom nodes - test this
                is_coincident = false;
                for (size_t j = 0; j < rReferenceGeom.PointsNumber(); ++j) {
                    if (norm_2(rPoints[i] - rReferenceGeom.GetPoint(j).Coordinates()) < Tolerance)
                    {
                        is_coincident = true;
                        break;
                    }
                }
                if (!is_coincident) return false;
            }
        }
        return true;

        KRATOS_CATCH("")
    }


    void MPMSearchElementUtility::Check3DBackGroundMeshIsCubicAxisAligned(const std::vector<typename GeometryType::Pointer> rIntersectedGeometries)
    {
        KRATOS_TRY

            NodeType point_low, point_high;
        for (size_t i = 0; i < rIntersectedGeometries.size(); ++i) {
            if (rIntersectedGeometries[i]->GetGeometryType() != GeometryData::Kratos_Hexahedra3D8) {
                KRATOS_ERROR << "3D PQMPM CAN ONLY BE USED FOR AXIS-ALIGNED CUBIC BACKGROUND GRIDS";
            }
            rIntersectedGeometries[i]->BoundingBox(point_low, point_high);
            for (size_t j = 0; j < rIntersectedGeometries[i]->PointsNumber(); ++j) {
                for (size_t k = 0; k < 3; ++k) {
                    if (rIntersectedGeometries[i]->GetPoint(j).Coordinates()[k] != point_low[k]) {
                        if (rIntersectedGeometries[i]->GetPoint(j).Coordinates()[k] != point_high[k]) {
                            KRATOS_ERROR << "3D PQMPM CAN ONLY BE USED FOR AXIS-ALIGNED CUBIC BACKGROUND GRIDS";
                        }
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }


    const bool MPMSearchElementUtility::CheckNoPointsAreInGeom(
        const std::vector<array_1d<double, 3>>& rPoints,
        const GeometryType& rReferenceGeom,
        const double Tolerance)
    {
        KRATOS_TRY

            array_1d<double, 3> dummy_local_coords;
        for (size_t i = 0; i < rPoints.size(); ++i) {
            if (rReferenceGeom.IsInside(rPoints[i], dummy_local_coords), Tolerance) return false;
        }
        return true;

        KRATOS_CATCH("")
    }


    void MPMSearchElementUtility::Create2DPolygonBoundingSquareFromPoints(const std::vector<array_1d<double, 3>>& rPoints,
        std::vector<Boost2DPointType>& rPolygonPoints,
        Boost2DPolygonType& rPolygon,
        const bool XActive, const bool YActive, const bool ZActive)
    {
        KRATOS_TRY

        if (!XActive || !YActive || ZActive)  if (rPoints.size() != 8)
            KRATOS_ERROR << "ALL BOUNDING SQUARES SHOULD BE CONSTRUCTED IN XY SPACE EXCEPT FOR HEX BACKGROUND GRID";

        if (rPolygonPoints.size() != 4) rPolygonPoints.resize(4);

        rPolygon.clear();

        if (XActive && YActive && !ZActive)
        {
            for (size_t i = 0; i < rPolygonPoints.size(); ++i) {
                rPolygonPoints[i] = Boost2DPointType(rPoints[i][0], rPoints[i][1]);
            }
        }
        else if (!XActive && YActive && ZActive) // 3D case only!
        {
            rPolygonPoints[0] = Boost2DPointType(rPoints[0][1], rPoints[0][2]);
            rPolygonPoints[1] = Boost2DPointType(rPoints[4][1], rPoints[4][2]);
            rPolygonPoints[2] = Boost2DPointType(rPoints[7][1], rPoints[7][2]);
            rPolygonPoints[3] = Boost2DPointType(rPoints[3][1], rPoints[3][2]); // as per Hexahedra3D8 node ordering
        }
        else if (XActive && !YActive && ZActive)
        {
            rPolygonPoints[0] = Boost2DPointType(rPoints[0][0], rPoints[0][2]);
            rPolygonPoints[1] = Boost2DPointType(rPoints[1][0], rPoints[1][2]);
            rPolygonPoints[2] = Boost2DPointType(rPoints[5][0], rPoints[5][2]);
            rPolygonPoints[3] = Boost2DPointType(rPoints[4][0], rPoints[4][2]);
        }
        else
        {
            KRATOS_ERROR << "INVALID PLANE TO MAKE 2D POLYGON IN!";
        }
        rPolygon.outer().assign(rPolygonPoints.begin(), rPolygonPoints.end());
        boost::geometry::correct(rPolygon); // to close the polygon

        KRATOS_CATCH("")
    }

    void MPMSearchElementUtility::Create2DPolygonFromGeometry(const GeometryType& rGeom,
        std::vector<Boost2DPointType>& rPolygonPoints,
        Boost2DPolygonType& rPolygon,
        const bool XActive, const bool YActive, const bool ZActive)
    {
        KRATOS_TRY

        rPolygon.clear();

        if (rGeom.WorkingSpaceDimension() == 3)
        {
            if (rPolygonPoints.size() != 4) rPolygonPoints.resize(4);
            NodeType point_low, point_high;
            rGeom.BoundingBox(point_low, point_high);

            if (XActive && YActive && !ZActive)
            {
                rPolygonPoints[0] = Boost2DPointType(point_low[0], point_low[1]);
                rPolygonPoints[1] = Boost2DPointType(point_high[0], point_low[1]);
                rPolygonPoints[2] = Boost2DPointType(point_high[0], point_high[1]);
                rPolygonPoints[3] = Boost2DPointType(point_low[0], point_high[1]);
            }
            else if (XActive && !YActive && ZActive)
            {
                rPolygonPoints[0] = Boost2DPointType(point_low[0], point_low[2]);
                rPolygonPoints[1] = Boost2DPointType(point_high[0], point_low[2]);
                rPolygonPoints[2] = Boost2DPointType(point_high[0], point_high[2]);
                rPolygonPoints[3] = Boost2DPointType(point_low[0], point_high[2]);
            }
            else if (!XActive && YActive && ZActive)
            {
                rPolygonPoints[0] = Boost2DPointType(point_low[1], point_low[2]);
                rPolygonPoints[1] = Boost2DPointType(point_high[1], point_low[2]);
                rPolygonPoints[2] = Boost2DPointType(point_high[1], point_high[2]);
                rPolygonPoints[3] = Boost2DPointType(point_low[1], point_high[2]);
            }
            else
            {
                KRATOS_ERROR << "INVALID PLANE TO MAKE 2D POLYGON IN!";
            }
        }
        else
        {
            if (rPolygonPoints.size() != rGeom.PointsNumber()) rPolygonPoints.resize(rGeom.PointsNumber());
            for (size_t i = 0; i < rGeom.PointsNumber(); ++i) {
                rPolygonPoints[i] = Boost2DPointType(rGeom.GetPoint(i).X(), rGeom.GetPoint(i).Y());
            }
        }
        rPolygon.outer().assign(rPolygonPoints.begin(), rPolygonPoints.end());
        boost::geometry::correct(rPolygon); // to close the polygon

        KRATOS_CATCH("")
    }


    IntegrationPoint<3> MPMSearchElementUtility::CreateSubPoint(const array_1d<double, 3>& rGlobalCoords, const double rVolumeFraction,
        const GeometryType& rBackgroundGridElementGeom, Vector& rN, Matrix& rDN_De)
    {
        KRATOS_TRY

            array_1d<double, 3> local_coordinates;
        rBackgroundGridElementGeom.PointLocalCoordinates(local_coordinates, rGlobalCoords);
        rBackgroundGridElementGeom.ShapeFunctionsValues(rN, local_coordinates);
        rBackgroundGridElementGeom.ShapeFunctionsLocalGradients(rDN_De, local_coordinates);

        return IntegrationPoint<3>(local_coordinates, rVolumeFraction);

        KRATOS_CATCH("")
    }


    void MPMSearchElementUtility::Determine2DSubPoint(const GeometryType& rGridElement, const std::vector<array_1d<double, 3>>& rMasterDomainPoints,
        array_1d<double, 3>& rSubPointCoord, double& rSubPointVolume)
    {
        KRATOS_TRY

            // make boost polygon of current background element geometry
            std::vector<Boost2DPointType> polygon_grid_points(rGridElement.PointsNumber());
        Boost2DPolygonType polygon_grid;
        Create2DPolygonFromGeometry(rGridElement, polygon_grid_points, polygon_grid);

        // make boost polygon of bounding box
        std::vector<Boost2DPointType> polygon_box_points(4);
        Boost2DPolygonType polygon_box;
        Create2DPolygonBoundingSquareFromPoints(rMasterDomainPoints, polygon_box_points, polygon_box);

        // make boost polygon result container
        std::vector<Boost2DPolygonType> polygon_result_container;

        // reset accumulated quantities
        rSubPointVolume = 0.0;
        rSubPointCoord.clear();
        Boost2DPointType centroid_result;

        // accumulate result over intersected sub-polygons
        if (boost::geometry::intersection(polygon_grid, polygon_box, polygon_result_container)) {
            for (auto& polygon_result : polygon_result_container) {
                rSubPointVolume += boost::geometry::area(polygon_result);
                boost::geometry::centroid(polygon_result, centroid_result);
                rSubPointCoord[0] += centroid_result.get<0>();
                rSubPointCoord[1] += centroid_result.get<1>();
            }
        }
        else
            KRATOS_ERROR << "BOOST INTERSECTION FAILED ALTHOUGH KRATOS INTERSECTION WORKED!";

        rSubPointCoord /= double(polygon_result_container.size());

        KRATOS_CATCH("")
    }


    void MPMSearchElementUtility::Determine3DSubPoint(const GeometryType& rGridElement, const std::vector<array_1d<double, 3>>& rMasterDomainPoints,
        array_1d<double, 3>& rSubPointCoord, double& rSubPointVolume)
    {
        KRATOS_TRY

        // NOTE: THIS FUNCTION ASSUMES THE BACKGROUND GRID ELEMENT IS PERFECTLY CUBIC AND THE RESULTING INTERSECTION VOLUME IS A RECTANGULAR PRISM

        // make boost xy polygon of current background element geometry
        std::vector<Boost2DPointType> polygon_grid_xy_points(4);
        Boost2DPolygonType polygon_grid_xy;
        Create2DPolygonFromGeometry(rGridElement, polygon_grid_xy_points, polygon_grid_xy);

        // make boost yz polygon of current background element geometry
        std::vector<Boost2DPointType> polygon_grid_yz_points(4);
        Boost2DPolygonType polygon_grid_yz;
        Create2DPolygonFromGeometry(rGridElement, polygon_grid_yz_points, polygon_grid_yz, false, true, true);

        // make boost xy polygon of bounding box
        std::vector<Boost2DPointType> polygon_box_xy_points(4);
        Boost2DPolygonType polygon_box_xy;
        Create2DPolygonBoundingSquareFromPoints(rMasterDomainPoints, polygon_box_xy_points, polygon_box_xy);

        // make boost yz polygon of bounding box
        std::vector<Boost2DPointType> polygon_box_yz_points(4);
        Boost2DPolygonType polygon_box_yz;
        Create2DPolygonBoundingSquareFromPoints(rMasterDomainPoints, polygon_box_yz_points, polygon_box_yz, false, true, true);

        // make boost polygon result container
        std::vector<Boost2DPolygonType> polygon_xy_result_container;
        std::vector<Boost2DPolygonType> polygon_yz_result_container;

        // reset accumulated quantities
        double sub_volume_area = 0.0;

        rSubPointCoord.clear();
        Boost2DPointType centroid_result;

        // Determine area and x y coordinates from xy polygons
        if (boost::geometry::intersection(polygon_grid_xy, polygon_box_xy, polygon_xy_result_container)) {
            for (auto& polygon_result : polygon_xy_result_container) {
                sub_volume_area += boost::geometry::area(polygon_result);
                boost::geometry::centroid(polygon_result, centroid_result);
                rSubPointCoord[0] += centroid_result.get<0>();
                rSubPointCoord[1] += centroid_result.get<1>();
            }
        }
        else
            KRATOS_ERROR << "BOOST INTERSECTION FAILED ALTHOUGH KRATOS INTERSECTION WORKED!";

        rSubPointCoord /= double(polygon_xy_result_container.size()); // at the moment this is just the xy coords!

        // Perform yz polygon intersection to determine depth and z-position of sub-point
        // local x = global y
        // local y = global z
        array_1d<double, 2> sub_point_z_coord = ZeroVector(2);
        bool is_initialized = false;
        double min_z = 0.0;
        double max_z = 0.0;
        if (boost::geometry::intersection(polygon_grid_yz, polygon_box_yz, polygon_yz_result_container)) {
            for (auto& polygon_result : polygon_yz_result_container) {
                for (auto& result_point : polygon_result.outer()) {
                    if (!is_initialized) {
                        min_z = result_point.get<1>();
                        max_z = result_point.get<1>();
                        is_initialized = true;
                    }
                    else if (result_point.get<1>() < min_z) min_z = result_point.get<1>();
                    else if (result_point.get<1>() > max_z) max_z = result_point.get<1>();
                }
            }
        }
        else
            KRATOS_ERROR << "BOOST INTERSECTION FAILED ALTHOUGH KRATOS INTERSECTION WORKED!";

        rSubPointCoord[2] = 0.5 * (min_z + max_z);
        rSubPointVolume = sub_volume_area * (max_z - min_z);

        KRATOS_CATCH("")
    }


    typename Geometry<Node<3>>::Pointer MPMSearchElementUtility::CreateCustomQuadraturePoint(
        SizeType WorkingSpaceDimension,
        SizeType LocalSpaceDimension,
        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainer,
        typename Geometry<Node<3>>::PointsArrayType rPoints)
    {
        KRATOS_TRY

            if (WorkingSpaceDimension == 1 && LocalSpaceDimension == 1)
                return Kratos::make_shared<
                QuadraturePointPartitionedGeometry<Node<3>, 1>>(
                    rPoints, rShapeFunctionContainer);
            else if (WorkingSpaceDimension == 2 && LocalSpaceDimension == 1)
                return Kratos::make_shared<
                QuadraturePointPartitionedGeometry<Node<3>, 2, 1>>(
                    rPoints, rShapeFunctionContainer);
            else if (WorkingSpaceDimension == 2 && LocalSpaceDimension == 2)
                return Kratos::make_shared<
                QuadraturePointPartitionedGeometry<Node<3>, 2>>(
                    rPoints, rShapeFunctionContainer);
            else if (WorkingSpaceDimension == 3 && LocalSpaceDimension == 2)
                return Kratos::make_shared<
                QuadraturePointPartitionedGeometry<Node<3>, 3, 2>>(
                    rPoints, rShapeFunctionContainer);
            else if (WorkingSpaceDimension == 3 && LocalSpaceDimension == 3)
                return Kratos::make_shared<
                QuadraturePointPartitionedGeometry<Node<3>, 3>>(
                    rPoints, rShapeFunctionContainer);
            else {
                KRATOS_ERROR << "Working/Local space dimension combinations are "
                    << "not provided for QuadraturePointGeometry. WorkingSpaceDimension: "
                    << WorkingSpaceDimension << ", LocalSpaceDimension: " << LocalSpaceDimension
                    << std::endl;
            }
        KRATOS_CATCH("")
    }


    typename Geometry<Node<3>>::Pointer MPMSearchElementUtility::PartitionMasterMaterialPointsIntoSubPoints(const ModelPart& rBackgroundGridModelPart,
                                                    const array_1d<double, 3>& rCoordinates,
                                                    Element& rMasterMaterialPoint,
                                                    const typename Geometry<Node<3>>::Pointer pGeometry,
                                                    const double Tolerance)
    {
        KRATOS_TRY;

        const SizeType working_dim = pGeometry->WorkingSpaceDimension();

        const bool is_axisymmetric = (rBackgroundGridModelPart.GetProcessInfo().Has(IS_AXISYMMETRIC))
            ? rBackgroundGridModelPart.GetProcessInfo().GetValue(IS_AXISYMMETRIC)
            : false;

        // Get volume and set up master domain bounding points
        std::vector<double> mp_volume_vec;
        rMasterMaterialPoint.CalculateOnIntegrationPoints(MP_VOLUME, mp_volume_vec, rBackgroundGridModelPart.GetProcessInfo());
        const double side_half_length = std::pow(mp_volume_vec[0], 1.0 / double(working_dim))/2.0;
        const SizeType n_bounding_box_vertices = std::pow(2.0, working_dim);
        std::vector<array_1d<double, 3>> master_domain_points(n_bounding_box_vertices);
        CreateBoundingBoxPoints(master_domain_points, rCoordinates, side_half_length,working_dim);

        // If axisymmetric, we can't make a sub-point with x<0.
        if (is_axisymmetric) { if ((rCoordinates[0] - side_half_length) < std::numeric_limits<double>::epsilon()) {
            return CreateQuadraturePointsUtility<Node<3>>::CreateFromCoordinates(
                pGeometry, rCoordinates, rMasterMaterialPoint.GetGeometry().IntegrationPoints()[0].Weight());
            }
        }

        // Initially check if the bounding box volume scalar is less than the element volume scalar
        if (mp_volume_vec[0] <= pGeometry->DomainSize() && CheckAllPointsAreInGeom(master_domain_points, *pGeometry, Tolerance))
            return CreateQuadraturePointsUtility<Node<3>>::CreateFromCoordinates(
                pGeometry, rCoordinates, rMasterMaterialPoint.GetGeometry().IntegrationPoints()[0].Weight());
        else
        { // we need to do splitting. Initially determine all grid elements we intersect with
            std::vector<typename GeometryType::Pointer> intersected_geometries;
            std::vector<Element::Pointer> intersected_elements;
            const double z_mod = (working_dim == 3) ? 1.0 : 0.0; 
            const Point point_low(rCoordinates[0] - side_half_length, rCoordinates[1] - side_half_length, rCoordinates[2] - z_mod* side_half_length);
            const Point point_high(rCoordinates[0] + side_half_length, rCoordinates[1] + side_half_length, rCoordinates[2] + z_mod* side_half_length);
            SizeType number_of_nodes = 0;
            const double range_factor = (working_dim == 3) ? 2.0 : 1.414214; // 45 deg for each dim
            double center_to_center, maximum_contact_range;
            NodeType ele_point_low, ele_point_high;
            for (auto ele_it : rBackgroundGridModelPart.Elements()) 
            {
                center_to_center = norm_2(ele_it.GetGeometry().Center() - rCoordinates);
                ele_it.GetGeometry().BoundingBox(ele_point_low, ele_point_high);
                maximum_contact_range = range_factor * side_half_length + norm_2(ele_point_high - ele_point_low) / 2.0;
                if (center_to_center <= maximum_contact_range)
                {
                    if (ele_it.GetGeometry().HasIntersection(point_low, point_high)) {
                        number_of_nodes += ele_it.GetGeometry().PointsNumber();
                        intersected_geometries.push_back(ele_it.pGetGeometry());
                        intersected_elements.push_back(&ele_it);
                    }
                }
            }

            // If we are 3D, check background mesh are axis-aligned perfect cubes
            if (working_dim == 3)  Check3DBackGroundMeshIsCubicAxisAligned(intersected_geometries);

            // Prepare containers to hold all sub-points
            const SizeType number_of_sub_material_points = intersected_geometries.size();
            PointerVector<Node<3>> nodes_list(number_of_nodes);
            IntegrationPointsArrayType ips(number_of_sub_material_points);
            Matrix N_matrix(number_of_sub_material_points, number_of_nodes, 0.0);
            DenseVector<Matrix> DN_De_vector(number_of_sub_material_points);

            // Temporary local containers
            double sub_point_volume;
            array_1d<double, 3> sub_point_position;
            IndexType active_node_index = 0;
            IndexType active_subpoint_index = 0;

            // Loop over all intersected grid elements and make subpoints in each
            for (size_t i = 0; i < number_of_sub_material_points; ++i) {
                Matrix DN_De(intersected_geometries[i]->PointsNumber(), working_dim);
                Vector N(intersected_geometries[i]->PointsNumber());
                sub_point_position.clear();
                sub_point_volume = 0.0;
                IntegrationPoint<3> trial_subpoint;

                if (CheckNoPointsAreInGeom(master_domain_points, *intersected_geometries[i],Tolerance))  {
                    // whole element is completely inside bounding box

                    trial_subpoint = CreateSubPoint(intersected_geometries[i]->Center(),
                        intersected_geometries[i]->DomainSize() / mp_volume_vec[0],
                        *intersected_geometries[i], N, DN_De);
                }
                else  {
                    // only some of the background element is within the bounding box - most expensive check

                    if (working_dim == 2) {
                        Determine2DSubPoint(*intersected_geometries[i], master_domain_points, sub_point_position, sub_point_volume);
                        sub_point_position[2] = rCoordinates[2]; // set z coord of sub point to that of the master
                    }
                    else
                        Determine3DSubPoint(*intersected_geometries[i], master_domain_points, sub_point_position, sub_point_volume);
                    trial_subpoint = CreateSubPoint(sub_point_position, sub_point_volume / mp_volume_vec[0],
                        *intersected_geometries[i], N, DN_De);
                }

                // Transfer local data to containers
                if (trial_subpoint.Weight() > std::numeric_limits<double>::epsilon())
                {
                    intersected_elements[i]->Set(ACTIVE);
                    ips[active_subpoint_index] = trial_subpoint;
                    //ips_good.push_back(trial_subpoint);
                    DN_De_vector[active_subpoint_index] = DN_De;
                    for (size_t j = 0; j < N.size(); ++j) {
                        //nodes_list_good.push_back(intersected_geometries[i]->pGetPoint(j));
                        N_matrix(active_subpoint_index, active_node_index) = N[j];
                        nodes_list(active_node_index) = intersected_geometries[i]->pGetPoint(j);

                        active_node_index += 1;
                    }
                    active_subpoint_index += 1;
                }
            }

            if (active_subpoint_index == 1) return CreateQuadraturePointsUtility<Node<3>>::CreateFromCoordinates(
                    pGeometry, rCoordinates, rMasterMaterialPoint.GetGeometry().IntegrationPoints()[0].Weight());

            IntegrationPointsArrayType ips_good(active_subpoint_index);
            PointerVector<Node<3>> nodes_list_good(active_node_index);
            if (ips_good.size() == ips.size())
            {
                ips_good = ips;
                nodes_list_good = nodes_list;
            }
            else
            {
                N_matrix.resize(active_subpoint_index, active_node_index, true);
                DN_De_vector.resize(active_subpoint_index, true);
                for (size_t i = 0; i < active_subpoint_index; i++) ips_good[i] = ips[i];
                for (size_t i = 0; i < active_node_index; i++) nodes_list_good(i) = nodes_list(i);
            }

            Check(ips_good, std::numeric_limits<double>::epsilon(), N_matrix, DN_De_vector);

            GeometryData::IntegrationMethod ThisDefaultMethod = pGeometry->GetDefaultIntegrationMethod();
            typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::IntegrationPointsContainerType ips_container;
            ips_container[ThisDefaultMethod] = ips_good;
            typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::ShapeFunctionsValuesContainerType shape_function_container;
            shape_function_container[ThisDefaultMethod] = N_matrix;
            typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::ShapeFunctionsLocalGradientsContainerType shape_function_derivatives_container;
            shape_function_derivatives_container[ThisDefaultMethod] = DN_De_vector;

            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                ThisDefaultMethod,
                ips_container,
                shape_function_container,
                shape_function_derivatives_container);

            return CreateCustomQuadraturePoint( working_dim, pGeometry->LocalSpaceDimension(), data_container, nodes_list_good);
        }
        KRATOS_CATCH("");
    }
} // end namespace Kratos

