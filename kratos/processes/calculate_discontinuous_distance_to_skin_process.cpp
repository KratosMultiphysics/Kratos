//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Ruben Zorrilla
//
//  Collaborators:   Franziska Wahl
//

// System includes

// External includes

// Project includes
#include "geometries/plane_3d.h"
#include "processes/calculate_discontinuous_distance_to_skin_process.h"
#include "utilities/geometry_utilities.h"
#include "utilities/intersection_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/plane_approximation_utility.h"

namespace Kratos
{
    KRATOS_CREATE_LOCAL_FLAG(CalculateDiscontinuousDistanceToSkinProcessFlags, CALCULATE_ELEMENTAL_EDGE_DISTANCES, 0);
    KRATOS_CREATE_LOCAL_FLAG(CalculateDiscontinuousDistanceToSkinProcessFlags, CALCULATE_ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED, 1);

    template<std::size_t TDim>
    CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateDiscontinuousDistanceToSkinProcess(
        ModelPart& rVolumePart,
        ModelPart& rSkinPart)
        : mFindIntersectedObjectsProcess(rVolumePart, rSkinPart)
        , mrSkinPart(rSkinPart)
        , mrVolumePart(rVolumePart)
        , mOptions(CalculateDiscontinuousDistanceToSkinProcessFlags::CALCULATE_ELEMENTAL_EDGE_DISTANCES.AsFalse()
             | CalculateDiscontinuousDistanceToSkinProcessFlags::CALCULATE_ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED.AsFalse())
    {
    }

    template<std::size_t TDim>
    CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateDiscontinuousDistanceToSkinProcess(
        ModelPart& rVolumePart,
        ModelPart& rSkinPart,
        const Flags rOptions)
        : mFindIntersectedObjectsProcess(rVolumePart, rSkinPart)
        , mrSkinPart(rSkinPart)
        , mrVolumePart(rVolumePart)
        , mOptions(rOptions)
    {
    }

    template<std::size_t TDim>
    CalculateDiscontinuousDistanceToSkinProcess<TDim>::~CalculateDiscontinuousDistanceToSkinProcess()
    {
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::Initialize()
    {
        // Initialize the intersected objects process
        mFindIntersectedObjectsProcess.Initialize();

        // Initialize the elemental distances to the domain characteristic length
        const double initial_distance = this->CalculateCharacteristicLength();
        array_1d<double,mNumNodes> init_dist_vect;
        for (unsigned int i_node = 0; i_node < mNumNodes; ++i_node) {
            init_dist_vect[i_node] = initial_distance;
        }

        // Also initialize the embedded velocity of the fluid element and the TO_SPLIT flag.
        if (mOptions.Is(CalculateDiscontinuousDistanceToSkinProcessFlags::CALCULATE_ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED)) {
            // Initialize the edge distances vector
            array_1d<double, mNumEdges> init_edge_dist_vect;
            for (double& r_val : init_edge_dist_vect) {
                r_val = -1.0;
            }

            block_for_each(mrVolumePart.Elements(), [&](Element& rElement){
                rElement.Set(TO_SPLIT, false);
                rElement.SetValue(EMBEDDED_VELOCITY, ZeroVector(3));
                rElement.SetValue(ELEMENTAL_DISTANCES,init_dist_vect);
                rElement.SetValue(ELEMENTAL_EDGE_DISTANCES, init_edge_dist_vect);
                rElement.SetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED, init_edge_dist_vect);
            });
        } else if (mOptions.Is(CalculateDiscontinuousDistanceToSkinProcessFlags::CALCULATE_ELEMENTAL_EDGE_DISTANCES)) {
            // Initialize the edge distances vector
            array_1d<double, mNumEdges> init_edge_dist_vect;
            for (double& r_val : init_edge_dist_vect) {
                r_val = -1.0;
            }

            block_for_each(mrVolumePart.Elements(), [&](Element& rElement){
                rElement.Set(TO_SPLIT, false);
                rElement.SetValue(EMBEDDED_VELOCITY, ZeroVector(3));
                rElement.SetValue(ELEMENTAL_DISTANCES,init_dist_vect);
                rElement.SetValue(ELEMENTAL_EDGE_DISTANCES, init_edge_dist_vect);
            });
        } else {
            block_for_each(mrVolumePart.Elements(), [&](Element& rElement){
                rElement.Set(TO_SPLIT, false);
                rElement.SetValue(EMBEDDED_VELOCITY, ZeroVector(3));
                rElement.SetValue(ELEMENTAL_DISTANCES, init_dist_vect);
            });
        }
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::FindIntersections()
    {
        mFindIntersectedObjectsProcess.FindIntersections();
    }

    template<std::size_t TDim>
    std::vector<PointerVector<GeometricalObject>>& CalculateDiscontinuousDistanceToSkinProcess<TDim>::GetIntersections()
    {
        return mFindIntersectedObjectsProcess.GetIntersections();
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateDistances(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects)
    {
        const int number_of_elements = (mFindIntersectedObjectsProcess.GetModelPart1()).NumberOfElements();
        auto& r_elements = (mFindIntersectedObjectsProcess.GetModelPart1()).ElementsArray();

        if (mOptions.Is(CalculateDiscontinuousDistanceToSkinProcessFlags::CALCULATE_ELEMENTAL_EDGE_DISTANCES)) {
            #pragma omp parallel for schedule(dynamic)
            for (int i = 0; i < number_of_elements; ++i) {
                CalculateElementalAndEdgeDistances(*(r_elements[i]), rIntersectedObjects[i]);
            }
        } else {
            #pragma omp parallel for schedule(dynamic)
            for (int i = 0; i < number_of_elements; ++i) {
                CalculateElementalDistances(*(r_elements[i]), rIntersectedObjects[i]);
            }
        }
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::Clear()
    {
        mFindIntersectedObjectsProcess.Clear();
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::Execute()
    {
        this->Clear();
        this->Initialize();
        this->FindIntersections();
        this->CalculateDistances(this->GetIntersections());
    }

    /// Turn back information as a string.
    template<std::size_t TDim>
    std::string CalculateDiscontinuousDistanceToSkinProcess<TDim>::Info() const {
        return "CalculateDiscontinuousDistanceToSkinProcess";
    }

    /// Print information about this object.
    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::PrintData(std::ostream& rOStream) const
    {
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateElementalDistances(
        Element& rElement1,
        PointerVector<GeometricalObject>& rIntersectedObjects)
    {
        if (rIntersectedObjects.empty()) {
            rElement1.Set(TO_SPLIT, false);
            return;
        }

        // Get edges container
        const auto r_edges_container = rElement1.GetGeometry().GenerateEdges();

        // Compute the number of intersected edges
        array_1d<double, mNumEdges> cut_edges_ratio_vector;
        array_1d<double, mNumEdges> cut_extra_edges_ratio_vector;
        std::vector<array_1d <double,3> > int_pts_vector;
        const unsigned int n_cut_edges = ComputeEdgesIntersections(rElement1, rIntersectedObjects, r_edges_container,
            cut_edges_ratio_vector, cut_extra_edges_ratio_vector, int_pts_vector);

        // Check if there is intersection: 3 or more intersected edges for a tetrahedron
        // If there is only 1 or 2 intersected edges, intersection is not considered
        // If there is intersection, calculate the elemental distances
        const bool is_intersection = (n_cut_edges < rElement1.GetGeometry().WorkingSpaceDimension()) ? false : true;
        if (is_intersection){
            ComputeIntersectionPlaneElementalDistances(rElement1, rIntersectedObjects, int_pts_vector);
        }

        // Check if the element is split and set the TO_SPLIT flag accordingly
        const double epsilon = std::numeric_limits<double>::epsilon();
        SetToSplitFlag(rElement1, epsilon);
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateElementalAndEdgeDistances(
        Element& rElement1,
        PointerVector<GeometricalObject>& rIntersectedObjects)
    {
        if (rIntersectedObjects.empty()) {
            rElement1.Set(TO_SPLIT, false);
            return;
        }

        // Get edges container
        const auto r_edges_container = rElement1.GetGeometry().GenerateEdges();

        // Compute the number of intersected edges
        array_1d<double, mNumEdges> cut_edges_ratio_vector;
        array_1d<double, mNumEdges> cut_extra_edges_ratio_vector;
        std::vector<array_1d <double,3> > int_pts_vector;
        const unsigned int n_cut_edges = ComputeEdgesIntersections(rElement1, rIntersectedObjects, r_edges_container,
            cut_edges_ratio_vector, cut_extra_edges_ratio_vector, int_pts_vector);

        // Save the cut edges ratios in the ELEMENTAL_EDGE_DISTANCES variable
        rElement1.GetValue(ELEMENTAL_EDGE_DISTANCES) = cut_edges_ratio_vector;

        // Check if there is an intersection
        bool is_intersection = false;
        // Extrapolated edge distances were calculated (for Ausas incised elements)
        if (mOptions.Is(CalculateDiscontinuousDistanceToSkinProcessFlags::CALCULATE_ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED)) {
            // Save the cut edges ratios of the extrapolated geometry in the ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED variable
            rElement1.GetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED) = cut_extra_edges_ratio_vector;

            // Check whether element is incised (this includes case, in which three edges of tetrahedron are intersected)
            bool is_incised = false;
            for (std::size_t i = 0; i < cut_extra_edges_ratio_vector.size(); i++) {
                double tolerance = std::numeric_limits<double>::epsilon();
                if ( std::abs(cut_extra_edges_ratio_vector[i] - (-1.0)) > tolerance ) {
                    is_incised = true;
                }
            }

            // Calculate and save elemental (node) distances based on both: edge ratios of original and extrapolated geometry
            if (is_incised) {
                ComputeElementalDistancesFromEdgeRatios(rElement1, rIntersectedObjects, r_edges_container,
                    cut_edges_ratio_vector, cut_extra_edges_ratio_vector);
            }
            // If element is not incised, it is either not cut at all or completely intersected
            else {
                is_intersection = (n_cut_edges < rElement1.GetGeometry().WorkingSpaceDimension()) ? false : true;
            }
        } else {
            // 3D: 3 or more intersected edges for a tetrahedron
            // 2D: 2 or more intersected edges for a triangle
            is_intersection = (n_cut_edges < rElement1.GetGeometry().WorkingSpaceDimension()) ? false : true;
        }

        // If there is an intersection, calculate the elemental distances (node-based)
        if (is_intersection){
            ComputeIntersectionPlaneElementalDistances(rElement1, rIntersectedObjects, int_pts_vector);
        }

        // Check if the element is split and set the TO_SPLIT flag accordingly
        const double epsilon = std::numeric_limits<double>::epsilon();
        SetToSplitFlag(rElement1, epsilon);
    }

    template<std::size_t TDim>
    unsigned int CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputeEdgesIntersections(
        Element& rElement1,
        const PointerVector<GeometricalObject>& rIntersectedObjects,
        const Element::GeometryType::GeometriesArrayType& rEdgesContainer,
        array_1d<double,mNumEdges> &rCutEdgesRatioVector,
        array_1d<double,mNumEdges> &rCutExtraEdgesRatioVector,
        std::vector<array_1d <double,3> > &rIntersectionPointsArray)
    {
        // Initialize cut edges vectors and points arrays
        unsigned int n_cut_edges = 0;
        rIntersectionPointsArray.clear();
        array_1d<unsigned int, mNumEdges> cut_edges_vector(mNumEdges, 0);
        rCutEdgesRatioVector = array_1d<double, mNumEdges>(mNumEdges, -1.0);
        rCutExtraEdgesRatioVector = array_1d<double, mNumEdges>(mNumEdges, -1.0);

        // Initialize intersecting segments normal for extrapolated edge calculation
        array_1d<double,3> extra_geom_normal = ZeroVector(3);

        // Initialize average points and average normal
        array_1d<double,3> avg_pt;
        array_1d<double,3> avg_extra_geom_normal;
        std::vector<array_1d<double,3> > aux_pts;

        // Check which edges are intersected
        for (std::size_t i_edge = 0; i_edge < mNumEdges; ++i_edge){
            avg_pt = ZeroVector(3);
            avg_extra_geom_normal = ZeroVector(3);
            aux_pts.clear();

            // Check against all candidates to count the number of current edge intersections
            for (const auto &r_int_obj : rIntersectedObjects){
                // Call the compute intersection method
                Point int_pt;
                const auto &r_int_obj_geom = r_int_obj.GetGeometry();
                const int int_id = ComputeEdgeIntersection(r_int_obj_geom, rEdgesContainer[i_edge][0], rEdgesContainer[i_edge][1], int_pt);

                // There is intersection
                if (int_id == 1){
                    // Check if there is a close intersection (repeated intersection point)
                    bool is_repeated = false;
                    for (auto aux_pt : aux_pts){
                        const double aux_dist = norm_2(int_pt - aux_pt);
                        const double tol_edge = 1e-2*norm_2(rEdgesContainer[i_edge][0] - rEdgesContainer[i_edge][1]);
                        if (aux_dist < tol_edge){
                            is_repeated = true;
                            break;
                        }
                    }

                    // If the intersection pt. is not repeated, consider it
                    if (!is_repeated){
                        // Add the intersection pt. to the aux array pts.
                        aux_pts.push_back(int_pt);
                        // Increase the edge intersections counter
                        cut_edges_vector[i_edge] += 1;
                        // Save the intersection point for computing the average
                        avg_pt += int_pt;
                        // Get normal of intersecting segment for extrapolated cut edges calculation
                        if (mOptions.Is(CalculateDiscontinuousDistanceToSkinProcessFlags::CALCULATE_ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED)) {
                            array_1d<double,3> int_extra_geom_normal;
                            ComputeIntersectionNormalFromGeometry(r_int_obj_geom, int_extra_geom_normal);
                            avg_extra_geom_normal += int_extra_geom_normal;
                        }
                    }
                }
            }

            // Collect the current edge information
            if (cut_edges_vector[i_edge] != 0){
                // Average the edge intersection point and save it
                avg_pt /= cut_edges_vector[i_edge];
                rIntersectionPointsArray.push_back(avg_pt);
                // Save the ratio location of the average intersection point
                rCutEdgesRatioVector[i_edge] = ConvertIntersectionPointToEdgeRatio(rEdgesContainer[i_edge], avg_pt);
                // Increase the total intersected edges counter
                n_cut_edges++;
                // Get average normal of intersecting segments for the edge (for extrapolated cut edges calculation)
                if (mOptions.Is(CalculateDiscontinuousDistanceToSkinProcessFlags::CALCULATE_ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED)) {
                    avg_extra_geom_normal /= cut_edges_vector[i_edge];
                    extra_geom_normal += avg_extra_geom_normal;
                }
            }
        }

        // Calculate extrapolated edge distances (for extrapolated cut edges calculation)
        if (mOptions.Is(CalculateDiscontinuousDistanceToSkinProcessFlags::CALCULATE_ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED) && n_cut_edges > 0) {
            // Get average normal of intersecting segments for all cut edges
            extra_geom_normal /= n_cut_edges;
            // Compute the intersections of the element's edges with the extrapolated averaged geometry
            ComputeExtrapolatedEdgesIntersectionsIfIncised(rElement1, rEdgesContainer, n_cut_edges, rCutEdgesRatioVector, extra_geom_normal, rCutExtraEdgesRatioVector);
        }

        return n_cut_edges;
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputeIntersectionPlaneElementalDistances(
        Element& rElement,
        const PointerVector<GeometricalObject>& rIntersectedObjects,
        const std::vector<array_1d<double,3>>& rIntersectionPointsCoordinates)
    {
        const auto &r_geometry = rElement.GetGeometry();
        const unsigned int n_cut_edges = rIntersectionPointsCoordinates.size();

        // Get reference to ELEMENTAL_DISTANCES
        Vector& r_elemental_distances = rElement.GetValue(ELEMENTAL_DISTANCES);

        // If there are more than 3 (3D) or 2 (2D) intersected edges, compute the least squares plane approximation
        // using the ComputePlaneApproximation utility.
        // Otherwise, the distance is computed using the plane defined by the 3 (3D) or 2 (2D) intersection points.
        const bool do_plane_approx = (n_cut_edges == TDim) ? false : true;

        if (do_plane_approx){
            // Call the plane optimization utility
            array_1d<double,3> base_pt, normal;
            ComputePlaneApproximation(rElement, rIntersectionPointsCoordinates, base_pt, normal);

            // Compute the distance to the approximation plane
            Plane3D approximation_plane(normal, Point{base_pt});
            for (std::size_t i = 0; i < mNumNodes; i++) {
                r_elemental_distances[i] = approximation_plane.CalculateSignedDistance(r_geometry[i]);
            }
        } else {
            // Create a plane with the 3 intersection points (or 2 in 2D)
            Plane3D plane = SetIntersectionPlane(rIntersectionPointsCoordinates);

            // Compute the distance to the intersection plane
            for (std::size_t i = 0; i < mNumNodes; i++) {
                r_elemental_distances[i] = plane.CalculateSignedDistance(r_geometry[i]);
            }
        }

        // Correct the distance values orientation
        CorrectDistanceOrientation(r_geometry, rIntersectedObjects, r_elemental_distances);
    }

    template<std::size_t TDim>
    int CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputeEdgeIntersection(
        const Element::GeometryType& rIntObjGeometry,
        const Element::NodeType& rEdgePoint1,
        const Element::NodeType& rEdgePoint2,
        Point& rIntersectionPoint)
    {
        int intersection_flag = 0;
        const auto work_dim = rIntObjGeometry.WorkingSpaceDimension();
        if (work_dim == 2){
            intersection_flag = IntersectionUtilities::ComputeLineLineIntersection<Element::GeometryType>(
                rIntObjGeometry, rEdgePoint1.Coordinates(), rEdgePoint2.Coordinates(), rIntersectionPoint.Coordinates());
        } else if (work_dim == 3){
            intersection_flag = IntersectionUtilities::ComputeTriangleLineIntersection<Element::GeometryType>(
                rIntObjGeometry, rEdgePoint1.Coordinates(), rEdgePoint2.Coordinates(), rIntersectionPoint.Coordinates());
        } else {
            KRATOS_ERROR << "Working space dimension value equal to " << work_dim << ". Check your skin geometry implementation." << std::endl;
        }

        return intersection_flag;
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputeIntersectionNormal(
        const Element::GeometryType& rGeometry,
        const Vector& rElementalDistances,
        array_1d<double,3>& rNormal)
    {
        double volume;
        array_1d<double,mNumNodes> N;
        BoundedMatrix<double,mNumNodes,TDim> DN_DX;
        GeometryUtils::CalculateGeometryData(rGeometry, DN_DX, N, volume);

        rNormal = ZeroVector(3);
        for (std::size_t comp = 0; comp < TDim; ++comp){
            for (std::size_t i_node = 0; i_node < rGeometry.PointsNumber(); ++i_node){
                rNormal(comp) += DN_DX(i_node,comp)*rElementalDistances[i_node];
            }
        }
        rNormal /= norm_2(rNormal);
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputePlaneApproximation(
        const Element& rElement1,
        const std::vector< array_1d<double,3> >& rPointsCoord,
        array_1d<double,3>& rPlaneBasePointCoords,
        array_1d<double,3>& rPlaneNormal)
    {
        const auto work_dim = rElement1.GetGeometry().WorkingSpaceDimension();
        if (work_dim == 2){
            PlaneApproximationUtility<2>::ComputePlaneApproximation(rPointsCoord, rPlaneBasePointCoords, rPlaneNormal);
        } else if (work_dim == 3){
            PlaneApproximationUtility<3>::ComputePlaneApproximation(rPointsCoord, rPlaneBasePointCoords, rPlaneNormal);
        } else {
            KRATOS_ERROR << "Working space dimension value equal to " << work_dim << ". Check your skin geometry implementation." << std::endl;
        }
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::CorrectDistanceOrientation(
        const Element::GeometryType& rGeometry,
        const PointerVector<GeometricalObject>& rIntersectedObjects,
        Vector& rElementalDistances)
    {
        // Check the obtained intersection orientation (normal as distance gradient)
        array_1d<double,3> distance_normal;
        ComputeIntersectionNormal(rGeometry, rElementalDistances, distance_normal);

        // Vote the intersection orientation using the intersecting entities normals
        unsigned int n_pos = 0;
        unsigned int n_neg = 0;

        for (const auto &r_int_obj : rIntersectedObjects){
            const auto &r_int_obj_geom = r_int_obj.GetGeometry();

            array_1d<double, 3> r_int_obj_normal;
            ComputeIntersectionNormalFromGeometry(r_int_obj_geom, r_int_obj_normal);
            r_int_obj_normal /= norm_2(r_int_obj_normal);

            if (inner_prod(r_int_obj_normal, distance_normal) < 0.0){
                n_neg++;
            } else {
                n_pos++;
            }
        }

        // Negative votes win. Switch the distance values
        if (n_neg > n_pos){
            for (std::size_t i_node = 0; i_node < mNumNodes; ++i_node){
                rElementalDistances[i_node] *= -1.0;
            }
        }
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::SetToSplitFlag(
        Element& rElement,
        const double ZeroTolerance)
    {
        bool has_positive_distance = false;
        bool has_negative_distance = false;
        const auto& r_elemental_distances = rElement.GetValue(ELEMENTAL_DISTANCES);
        for (const double& r_dist : r_elemental_distances)
            if (r_dist > ZeroTolerance) {
                has_positive_distance = true;
            } else {
                has_negative_distance = true;
            }

        rElement.Set(TO_SPLIT, has_positive_distance && has_negative_distance);
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateEmbeddedVariableFromSkin(
        const Variable<double> &rVariable,
        const Variable<double> &rEmbeddedVariable)
    {
        this->CalculateEmbeddedVariableFromSkinSpecialization<double>(rVariable, rEmbeddedVariable);
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateEmbeddedVariableFromSkin(
        const Variable<array_1d<double,3>> &rVariable,
        const Variable<array_1d<double,3>> &rEmbeddedVariable)
    {
        this->CalculateEmbeddedVariableFromSkinSpecialization<array_1d<double,3>>(rVariable, rEmbeddedVariable);
    }

    template<>
    Plane3D CalculateDiscontinuousDistanceToSkinProcess<2>::SetIntersectionPlane(
        const std::vector<array_1d<double,3>> &rIntPtsVector)
    {
        // Since the Plane3D object only works in 3D, in 2D we set the intersection
        // plane by extruding the intersection point 0 in the z-direction.
        array_1d<double,3> z_coord_pt = rIntPtsVector[0];
        z_coord_pt[2] = 1.0;
        return Plane3D(Point{rIntPtsVector[0]}, Point{rIntPtsVector[1]}, Point{z_coord_pt});
    }

    template<>
    Plane3D CalculateDiscontinuousDistanceToSkinProcess<3>::SetIntersectionPlane(
        const std::vector<array_1d<double,3>> &rIntPtsVector)
    {
        return Plane3D(Point{rIntPtsVector[0]}, Point{rIntPtsVector[1]}, Point{rIntPtsVector[2]});
    }

    template<std::size_t TDim>
    double CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateCharacteristicLength()
    {
        // Get the background mesh model part
        auto &r_model_part = mFindIntersectedObjectsProcess.GetModelPart1();
        KRATOS_ERROR_IF(r_model_part.NumberOfNodes() == 0)
            << "Background mesh model part has no nodes." << std::endl;

        // Compute the domain characteristic length
        typedef CombinedReduction<MaxReduction<double>,MaxReduction<double>,MaxReduction<double>,MinReduction<double>,MinReduction<double>,MinReduction<double>> CustomReduction;
        double max_x, max_y, max_z, min_x, min_y, min_z;
        std::tie(max_x,max_y,max_z,min_x,min_y,min_z) = block_for_each<CustomReduction>(r_model_part.Nodes(),[](const Node<3>& rNode){
            return std::make_tuple(rNode[0],rNode[1],rNode[2],rNode[0],rNode[1],rNode[2]);}
        );

        const double char_length = std::sqrt(std::pow(max_x - min_x, 2) + std::pow(max_y - min_y, 2) + std::pow(max_z - min_z, 2));
        KRATOS_ERROR_IF(char_length < std::numeric_limits<double>::epsilon())
            << "Domain characteristic length is close to zero. Check if there is any node in the model part." << std::endl;

        return char_length;
    }

    template<>
    void inline CalculateDiscontinuousDistanceToSkinProcess<2>::ComputeIntersectionNormalFromGeometry(
        const Element::GeometryType &rGeometry,
        array_1d<double,3> &rIntObjNormal)
    {
        rIntObjNormal[0] = rGeometry[0].Y() - rGeometry[1].Y();
        rIntObjNormal[1] = rGeometry[1].X() - rGeometry[0].X();
        rIntObjNormal[2] = 0.0;
    }

    template<>
    void inline CalculateDiscontinuousDistanceToSkinProcess<3>::ComputeIntersectionNormalFromGeometry(
        const Element::GeometryType &rGeometry,
        array_1d<double,3> &rIntObjNormal)
    {
        MathUtils<double>::CrossProduct(rIntObjNormal, rGeometry[1]-rGeometry[0], rGeometry[2]-rGeometry[0]);
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputeExtrapolatedEdgesIntersectionsIfIncised(
        const Element& rElement,
        const Element::GeometryType::GeometriesArrayType& rEdgesContainer,
        unsigned int &rNumCutEdges,
        array_1d<double,mNumEdges>& rCutEdgesRatioVector,
        array_1d<double,3> &rExtraGeomNormal,
        array_1d<double,mNumEdges>& rCutExtraEdgesRatioVector)
    {
        if (TDim == 2 && rNumCutEdges == 1) {
            ComputeExtrapolatedGeometryIntersections(rElement, rEdgesContainer, rNumCutEdges, rCutEdgesRatioVector, rExtraGeomNormal, rCutExtraEdgesRatioVector);
        } else if (TDim == 3) {
            if (rNumCutEdges == 3) {
                // if all three cut edges share one node, then the element is intersected and not incised
                if (CheckIfCutEdgesShareNode(rElement, rEdgesContainer, rCutEdgesRatioVector)) {
                    return;
                }
                ComputeExtrapolatedGeometryIntersections(rElement, rEdgesContainer, rNumCutEdges, rCutEdgesRatioVector, rExtraGeomNormal, rCutExtraEdgesRatioVector);
            } else if (rNumCutEdges == 1 || rNumCutEdges == 2) {
                ComputeExtrapolatedGeometryIntersections(rElement, rEdgesContainer, rNumCutEdges, rCutEdgesRatioVector, rExtraGeomNormal, rCutExtraEdgesRatioVector);
            }
        }
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputeExtrapolatedGeometryIntersections(
        const Element& rElement,
        const Element::GeometryType::GeometriesArrayType& rEdgesContainer,
        unsigned int& rNumCutEdges,
        array_1d<double,mNumEdges>& rCutEdgesRatioVector,
        array_1d<double,3>& rExtraGeomNormal,
        array_1d<double,mNumEdges>& rCutExtraEdgesRatioVector)
    {
        // Calculate average point of intersection points from rCutEdgesRatioVector for intersection plane definition
        array_1d<double,3> avg_base_point = ZeroVector(3);
        for (std::size_t i_edge = 0; i_edge < mNumEdges; i_edge++) {
            // Calculate point coordinates and add to avg_point
            if (rCutEdgesRatioVector[i_edge] >= 0.0) {
                avg_base_point += ConvertEdgeRatioToIntersectionPoint(rEdgesContainer[i_edge], rCutEdgesRatioVector[i_edge]);
            }
        }
        avg_base_point /= rNumCutEdges;

        // Calculate intersections of each edge of element, which is not cut already, with the intersection plane
        for (std::size_t i_edge = 0; i_edge < mNumEdges; i_edge++) {
            if (rCutEdgesRatioVector[i_edge] == -1.0) {
                array_1d<double,3> extra_int_pt;
                const auto& edge_point_0 = rEdgesContainer[i_edge][0];
                const auto& edge_point_1 = rEdgesContainer[i_edge][1];
                int is_intersection = IntersectionUtilities::ComputePlaneLineIntersection(
                    avg_base_point, rExtraGeomNormal, edge_point_0.Coordinates(), edge_point_1.Coordinates(), extra_int_pt);

                // Calculate intersection ratio of edge and save it
                if (is_intersection == 1) {
                    rCutExtraEdgesRatioVector[i_edge] = ConvertIntersectionPointToEdgeRatio(rEdgesContainer[i_edge], extra_int_pt);
                }
            }
        }
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputeElementalDistancesFromEdgeRatios(
        Element& rElement,
        const PointerVector<GeometricalObject>& rIntersectedObjects,
        const Element::GeometryType::GeometriesArrayType& rEdgesContainer,
        const array_1d<double,mNumEdges> &rCutEdgesRatioVector,
        const array_1d<double,mNumEdges> &rCutExtraEdgesRatioVector)
    {
        // Get combined vector of intersection and extrapolated edge ratios and count total amount of cut edges
        array_1d<double, mNumEdges> combined_edge_ratios = array_1d<double, mNumEdges>(mNumEdges, -1.0);
        std::size_t n_cut_edges = 0;
        for (std::size_t i_edge = 0; i_edge < mNumEdges; i_edge++) {
            if (rCutEdgesRatioVector[i_edge] > 0.0){
                combined_edge_ratios[i_edge] = rCutEdgesRatioVector[i_edge];
                n_cut_edges++;
            } else if (rCutExtraEdgesRatioVector[i_edge] > 0.0) {
                combined_edge_ratios[i_edge] = rCutExtraEdgesRatioVector[i_edge];
                n_cut_edges++;
            } else {
                combined_edge_ratios[i_edge] = -1.0;
            }
        }

        // Check if element qualifies as being intersected after calculating extrapolated intersections and calculate elemental distances
        if (n_cut_edges >= TDim) {
            // Calculate points from nodes of edges and length ratio of intersections
            std::vector<array_1d <double,3> > intsect_pts_vector;
            ConvertRatiosToIntersectionPoints(rElement.GetGeometry(), rEdgesContainer, combined_edge_ratios, intsect_pts_vector);

            // Create plane from intersection points, calculate and correct node distances to plane
            ComputeIntersectionPlaneElementalDistances(rElement, rIntersectedObjects, intsect_pts_vector);
        }
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ConvertRatiosToIntersectionPoints(
        const Element::GeometryType& rGeometry,
        const Element::GeometryType::GeometriesArrayType& rEdgesContainer,
        const array_1d<double,mNumEdges> &rEdgeRatiosVector,
        std::vector<array_1d <double,3> > &rIntersectionPointsVector)
    {
        // Clear intersection points vector and allocate averaged point
        rIntersectionPointsVector.clear();
        array_1d<double,3> avg_pt;

        // Calculate intersection point of each edge that is intersected
        for (std::size_t i_edge = 0; i_edge < mNumEdges; ++i_edge){
            if (rEdgeRatiosVector[i_edge] >= 0.0){
                avg_pt = ConvertEdgeRatioToIntersectionPoint(rEdgesContainer[i_edge], rEdgeRatiosVector[i_edge]);
                rIntersectionPointsVector.push_back(avg_pt);
            }
        }
    }

    template<std::size_t TDim>
    double CalculateDiscontinuousDistanceToSkinProcess<TDim>::ConvertIntersectionPointToEdgeRatio(
        const Geometry<Node<3> >& rEdge,
        const array_1d<double,3>& rIntersectionPoint)
    {
        const double edge_length = rEdge.Length();
        KRATOS_ERROR_IF(edge_length < std::numeric_limits<double>::epsilon())
            << "Edge length of incised element is close to zero." << std::endl;
        const double dist_avg_pt = norm_2(rEdge[0] - rIntersectionPoint);
        return dist_avg_pt / edge_length;
    }

    template<std::size_t TDim>
    array_1d<double,3> CalculateDiscontinuousDistanceToSkinProcess<TDim>::ConvertEdgeRatioToIntersectionPoint(
        const Geometry<Node<3> >& rEdge,
        const double& rEdgeRatio)
    {
        return rEdge[0] + rEdgeRatio * (rEdge[1] - rEdge[0]);
    }

    template<std::size_t TDim>
    bool CalculateDiscontinuousDistanceToSkinProcess<TDim>::CheckIfCutEdgesShareNode(
        const Element& rElement,
        const Element::GeometryType::GeometriesArrayType& rEdgesContainer,
        const array_1d<double,mNumEdges>& rCutEdgesRatioVector) const
    {
        // Get nodes of cut edges (Point necessary to be able to use operator '=='!)
        std::vector<Point> nodes_0;
        std::vector<Point> nodes_1;
        for (std::size_t i_edge = 0; i_edge < mNumEdges; i_edge++) {
            if (rCutEdgesRatioVector[i_edge] > -1) {
                nodes_0.push_back(rEdgesContainer[i_edge][0]);
                nodes_1.push_back(rEdgesContainer[i_edge][1]);
            }
        }

        // Check if cut edges share a node - operator==
        bool is_shared = true;
        for (std::size_t i = 1; i < nodes_0.size(); i++) {
            if (!(nodes_0[0] == nodes_0[i] || nodes_0[0] == nodes_1[i])) {
                is_shared = false;
            }
        }
        if (!is_shared) {
            is_shared = true;
            for (std::size_t i = 1; i < nodes_0.size(); i++) {
                if (!(nodes_1[0] == nodes_0[i] || nodes_1[0] == nodes_1[i])) {
                    is_shared = false;
                }
            }
        }

        return is_shared;
    }

    template class Kratos::CalculateDiscontinuousDistanceToSkinProcess<2>;
    template class Kratos::CalculateDiscontinuousDistanceToSkinProcess<3>;

}  // namespace Kratos.
