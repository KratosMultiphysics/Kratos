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
#include <unordered_set>

// External includes

// Project includes
#include "geometries/plane_3d.h"
#include "processes/calculate_discontinuous_distance_to_skin_process.h"
#include "processes/find_global_nodal_entity_neighbours_process.h"
#include "utilities/geometry_utilities.h"
#include "utilities/intersection_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/plane_approximation_utility.h"
#include "utilities/global_pointer_utilities.h"

namespace Kratos
{
    KRATOS_CREATE_LOCAL_FLAG(CalculateDiscontinuousDistanceToSkinProcessFlags, CALCULATE_ELEMENTAL_EDGE_DISTANCES, 0);
    KRATOS_CREATE_LOCAL_FLAG(CalculateDiscontinuousDistanceToSkinProcessFlags, CALCULATE_ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED, 1);
    KRATOS_CREATE_LOCAL_FLAG(CalculateDiscontinuousDistanceToSkinProcessFlags, USE_POSITIVE_EPSILON_FOR_ZERO_VALUES, 2);

    template<std::size_t TDim>
    CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateDiscontinuousDistanceToSkinProcess(
        ModelPart& rVolumePart,
        ModelPart& rSkinPart)
        : mFindIntersectedObjectsProcess(rVolumePart, rSkinPart)
        , mrSkinPart(rSkinPart)
        , mrVolumePart(rVolumePart)
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
        KRATOS_WARNING("DEPRECATION") << "Please use the parameter constructor instead of the flag constructor." << std::endl;
        mCalculateElementalEdgeDistances = mOptions.Is(CalculateDiscontinuousDistanceToSkinProcessFlags::CALCULATE_ELEMENTAL_EDGE_DISTANCES);
        mCalculateElementalEdgeDistancesExtrapolated = mOptions.Is(CalculateDiscontinuousDistanceToSkinProcessFlags::CALCULATE_ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED);
        mUsePositiveEpsilonForZeroValues = mOptions.Is(CalculateDiscontinuousDistanceToSkinProcessFlags::USE_POSITIVE_EPSILON_FOR_ZERO_VALUES);
    }

    template<std::size_t TDim>
    CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateDiscontinuousDistanceToSkinProcess(
        ModelPart& rVolumePart,
        ModelPart& rSkinPart,
        Parameters rParameters )
        : mFindIntersectedObjectsProcess(rVolumePart, rSkinPart)
        , mrSkinPart(rSkinPart)
        , mrVolumePart(rVolumePart)
    {
        rParameters.RecursivelyValidateAndAssignDefaults(GetDefaultParameters());
        mCalculateElementalEdgeDistances = rParameters["calculate_elemental_edge_distances"].GetBool();
        mCalculateElementalEdgeDistancesExtrapolated = rParameters["calculate_elemental_edge_distances_extrapolated"].GetBool();
        mUsePositiveEpsilonForZeroValues = rParameters["use_positive_epsilon_for_zero_values"].GetBool();
        mpElementalDistancesVariable = &KratosComponents<Variable<Vector>>::Get(rParameters["elemental_distances_variable"].GetString());;
        mpElementalEdgeDistancesVariable = &KratosComponents<Variable<Vector>>::Get(rParameters["elemental_edge_distances_variable"].GetString());
        mpElementalEdgeDistancesExtrapolatedVariable = &KratosComponents<Variable<Vector>>::Get(rParameters["elemental_edge_distances_extrapolated_variable"].GetString());
        mpEmbeddedVelocityVariable = &KratosComponents<Variable<array_1d<double, 3>>>::Get(rParameters["embedded_velocity_variable"].GetString());
    }

    template<std::size_t TDim>
    const Parameters CalculateDiscontinuousDistanceToSkinProcess<TDim>::GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "elemental_distances_variable"                          :  "ELEMENTAL_DISTANCES",
            "elemental_edge_distances_variable"                     :  "ELEMENTAL_EDGE_DISTANCES",
            "elemental_edge_distances_extrapolated_variable"        :  "ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED",
            "embedded_velocity_variable"                            :  "EMBEDDED_VELOCITY",
            "calculate_elemental_edge_distances"                    : false,
            "calculate_elemental_edge_distances_extrapolated"       : false,
            "use_positive_epsilon_for_zero_values"                  : true
        })" );

        return default_parameters;
    }

    template<std::size_t TDim>
    CalculateDiscontinuousDistanceToSkinProcess<TDim>::~CalculateDiscontinuousDistanceToSkinProcess()
    {
    }

    template<std::size_t TDim>
    Process::Pointer CalculateDiscontinuousDistanceToSkinProcess<TDim>::Create(
        Model& rModel,
        Parameters ThisParameters
        )
    {
        const std::string& r_volume_part_name = ThisParameters["volume_model_part_name"].GetString();
        const std::string& r_skin_part_name = ThisParameters["skin_model_part_name"].GetString();
        return Kratos::make_shared<CalculateDiscontinuousDistanceToSkinProcess<TDim>>(
            rModel.GetModelPart(r_volume_part_name),
            rModel.GetModelPart(r_skin_part_name),
            ThisParameters);
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::Initialize()
    {
        // Initialize the intersected objects process
        mFindIntersectedObjectsProcess.ExecuteInitialize();

        // Initialize the elemental distances to the domain characteristic length
        const double initial_distance = this->CalculateCharacteristicLength();
        array_1d<double,mNumNodes> init_dist_vect;
        for (unsigned int i_node = 0; i_node < mNumNodes; ++i_node) {
            init_dist_vect[i_node] = initial_distance;
        }

        // Also initialize the embedded velocity of the fluid element and the TO_SPLIT flag.
        if (mCalculateElementalEdgeDistancesExtrapolated) {
            // Initialize the edge distances vector
            array_1d<double, mNumEdges> init_edge_dist_vect;
            for (double& r_val : init_edge_dist_vect) {
                r_val = -1.0;
            }

            block_for_each(mrVolumePart.Elements(), [&](Element& rElement){
                rElement.Set(TO_SPLIT, false);
                rElement.SetValue(*mpEmbeddedVelocityVariable, ZeroVector(3));
                rElement.SetValue(*mpElementalDistancesVariable,init_dist_vect);
                rElement.SetValue(*mpElementalEdgeDistancesVariable, init_edge_dist_vect);
                rElement.SetValue(*mpElementalEdgeDistancesExtrapolatedVariable, init_edge_dist_vect);
            });
        } else if (mCalculateElementalEdgeDistances) {
            // Initialize the edge distances vector
            array_1d<double, mNumEdges> init_edge_dist_vect;
            for (double& r_val : init_edge_dist_vect) {
                r_val = -1.0;
            }

            block_for_each(mrVolumePart.Elements(), [&](Element& rElement){
                rElement.Set(TO_SPLIT, false);
                rElement.SetValue(*mpEmbeddedVelocityVariable, ZeroVector(3));
                rElement.SetValue(*mpElementalDistancesVariable,init_dist_vect);
                rElement.SetValue(*mpElementalEdgeDistancesVariable, init_edge_dist_vect);
            });
        } else {
            block_for_each(mrVolumePart.Elements(), [&](Element& rElement){
                rElement.Set(TO_SPLIT, false);
                rElement.SetValue(*mpEmbeddedVelocityVariable, ZeroVector(3));
                rElement.SetValue(*mpElementalDistancesVariable, init_dist_vect);
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

        if (mCalculateElementalEdgeDistances) {
            #pragma omp parallel for schedule(dynamic)
            for (int i = 0; i < number_of_elements; ++i) {
                CalculateElementalAndEdgeDistances(*(r_elements[i]), rIntersectedObjects[i]);
            }
            if (mrVolumePart.GetCommunicator().GetDataCommunicator().MaxAll(mDetectedZeroDistanceValues)) {
                CheckAndCorrectEdgeDistances();
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
        rElement1.GetValue(*mpElementalEdgeDistancesVariable) = cut_edges_ratio_vector;

        // Check if there is an intersection
        bool is_intersection = false;
        // Extrapolated edge distances were calculated (for Ausas incised elements)
        if (mCalculateElementalEdgeDistancesExtrapolated) {
            // Save the cut edges ratios of the extrapolated geometry in the ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED variable
            rElement1.GetValue(*mpElementalEdgeDistancesExtrapolatedVariable) = cut_extra_edges_ratio_vector;

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
        std::vector<array_1d<double,3> > aux_avg_pts;

        // Check which edges are intersected
        for (std::size_t i_edge = 0; i_edge < mNumEdges; ++i_edge){
            avg_pt = ZeroVector(3);
            avg_extra_geom_normal = ZeroVector(3);
            aux_pts.clear();

            // Check against all candidates to count the number of current edge intersections
            const double edge_tolerance = 1e-6* rEdgesContainer[i_edge][0].Distance(rEdgesContainer[i_edge][1]);
            for (const auto &r_int_obj : rIntersectedObjects){
                // Call the compute intersection method
                Point int_pt;
                const auto &r_int_obj_geom = r_int_obj.GetGeometry();
                const int int_id = ComputeEdgeIntersection(r_int_obj_geom, rEdgesContainer[i_edge][0], rEdgesContainer[i_edge][1], int_pt);

                // There is intersection
                if (int_id == 1 || int_id == 3){
                    // If the intersection pt. is not repeated, consider it

                    if (!CheckIfPointIsRepeated(int_pt, aux_pts, edge_tolerance)) {
                        // Add the intersection pt. to the aux array pts.
                        aux_pts.push_back(int_pt);
                        // Increase the edge intersections counter
                        cut_edges_vector[i_edge] += 1;
                        // Save the intersection point for computing the average
                        avg_pt += int_pt;
                        // Get normal of intersecting segment for extrapolated cut edges calculation
                        if (mCalculateElementalEdgeDistancesExtrapolated) {
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
                // Save the ratio location of the average intersection point
                rCutEdgesRatioVector[i_edge] = ConvertIntersectionPointToEdgeRatio(rEdgesContainer[i_edge], avg_pt);
                // Increase the total intersected edges counter
                n_cut_edges++;
                // Check if the avg_pt is already present
                if (!CheckIfPointIsRepeated(avg_pt, aux_avg_pts, edge_tolerance)) {
                    rIntersectionPointsArray.push_back(avg_pt);
                    aux_avg_pts.push_back(avg_pt);
                }
                // Get average normal of intersecting segments for the edge (for extrapolated cut edges calculation)
                if (mCalculateElementalEdgeDistancesExtrapolated) {
                    avg_extra_geom_normal /= cut_edges_vector[i_edge];
                    extra_geom_normal += avg_extra_geom_normal;
                }
            }
        }

        // Calculate extrapolated edge distances (for extrapolated cut edges calculation)
        if (mCalculateElementalEdgeDistancesExtrapolated && n_cut_edges > 0) {
            // Get average normal of intersecting segments for all cut edges
            extra_geom_normal /= n_cut_edges;
            // Compute the intersections of the element's edges with the extrapolated averaged geometry
            ComputeExtrapolatedEdgesIntersectionsIfIncised(rElement1, rEdgesContainer, n_cut_edges, rCutEdgesRatioVector, extra_geom_normal, rCutExtraEdgesRatioVector);
        }

        return n_cut_edges;
    }

    template<std::size_t TDim>
    bool CalculateDiscontinuousDistanceToSkinProcess<TDim>::CheckIfPointIsRepeated(
        const array_1d<double,3>& rIntersectionPoint,
        const std::vector<array_1d<double,3>>&  rIntersectionPointsVector,
        const double& rEdgeTolerance)
    {
        // Check if there is a close intersection (repeated intersection point)
        for (const auto& aux_pt : rIntersectionPointsVector){
                const double aux_dist = norm_2(rIntersectionPoint - aux_pt);
                if (aux_dist < rEdgeTolerance){
                    return true;
            }
        }
        return false;
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputeElementalDistancesFromPlaneApproximation(
        Element& rElement,
        Vector& rElementalDistances,
        const std::vector<array_1d<double,3>>& rPointVector)

    {
        array_1d<double,3> base_pt, normal;
        ComputePlaneApproximation(rElement, rPointVector, base_pt, normal);

        // Compute the distance to the approximation plane
        Plane3D approximation_plane(normal, Point{base_pt});
        const auto &r_geometry = rElement.GetGeometry();
        for (std::size_t i = 0; i < mNumNodes; i++) {
            rElementalDistances[i] = approximation_plane.CalculateSignedDistance(r_geometry[i]);
        }
    }

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputeIntersectionPlaneElementalDistances(
        Element& rElement,
        const PointerVector<GeometricalObject>& rIntersectedObjects,
        const std::vector<array_1d<double,3>>& rIntersectionPointsCoordinates)
    {
        const auto &r_geometry = rElement.GetGeometry();
        const unsigned int n_cut_points = rIntersectionPointsCoordinates.size();

        // Get reference to ELEMENTAL_DISTANCES
        Vector& r_elemental_distances = rElement.GetValue(*mpElementalDistancesVariable);

        // If there are more than 3 (3D) or 2 (2D) intersected edges, compute the least squares plane approximation
        // using the ComputePlaneApproximation utility.
        // Otherwise, the distance is computed using the plane defined by the 3 (3D) or 2 (2D) intersection points.
        const bool do_plane_approx = (n_cut_points == TDim) ? false : true;

        if (do_plane_approx){
            if (n_cut_points > TDim) {
                // Call the plane optimization utility
                ComputeElementalDistancesFromPlaneApproximation(rElement, r_elemental_distances, rIntersectionPointsCoordinates);
            }
            else {
                // Not enough intersection points to build a plane
                // Use the intersected objects to create an approximation plane
                std::vector<array_1d<double,3>> int_pts_vector;
                for (const auto &r_int_obj : rIntersectedObjects) {
                    for (std::size_t i_int = 0; i_int < r_int_obj.GetGeometry().size(); i_int++) {
                        int_pts_vector.push_back(r_int_obj.GetGeometry()[i_int].Coordinates());
                    }
                }
                ComputeElementalDistancesFromPlaneApproximation(rElement, r_elemental_distances, int_pts_vector);
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

        ReplaceZeroDistances(r_elemental_distances);
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
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ReplaceZeroDistances(
        Vector& rElementalDistances)
    {
        const double multiplier = (mUsePositiveEpsilonForZeroValues)  ? mZeroToleranceMultiplier : -mZeroToleranceMultiplier;
        const double corrected_distance = multiplier*std::numeric_limits<double>::epsilon();
        for (double& r_distance : rElementalDistances) {
            if (std::abs(r_distance) < std::numeric_limits<double>::epsilon()) {
                r_distance = corrected_distance;
                if (!mDetectedZeroDistanceValues) {
                    mDetectedZeroDistanceValues = true;
                }
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
        const auto& r_elemental_distances = rElement.GetValue(*mpElementalDistancesVariable);
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
        KRATOS_ERROR_IF(r_model_part.GetCommunicator().GlobalNumberOfNodes() == 0)
            << "Background mesh model part has no nodes." << std::endl;

        // Compute the domain characteristic length
        typedef CombinedReduction<MaxReduction<double>,MaxReduction<double>,MaxReduction<double>,MinReduction<double>,MinReduction<double>,MinReduction<double>> CustomReduction;
        auto [max_x,max_y,max_z,min_x,min_y,min_z] = block_for_each<CustomReduction>(r_model_part.Nodes(),[](const Node& rNode){
            return std::make_tuple(rNode[0],rNode[1],rNode[2],rNode[0],rNode[1],rNode[2]);}
        );
        auto max_vector = r_model_part.GetCommunicator().GetDataCommunicator().MaxAll(std::vector<double>{max_x, max_y, max_z});
        auto min_vector = r_model_part.GetCommunicator().GetDataCommunicator().MinAll(std::vector<double>{min_x, min_y, min_z});

        const double char_length = std::sqrt(std::pow(max_vector[0] - min_vector[0], 2) + std::pow(max_vector[1] - min_vector[1], 2) + std::pow(max_vector[2] - min_vector[2], 2));
        KRATOS_ERROR_IF(char_length < std::numeric_limits<double>::epsilon())
            << "Domain characteristic length is close to zero. Check if there is any node in the model part." << std::endl;

        return char_length;
    }

    template<>
    void inline CalculateDiscontinuousDistanceToSkinProcess<2>::ComputeIntersectionNormalFromGeometry(
        const Element::GeometryType &rGeometry,
        array_1d<double,3> &rIntObjNormal)
    {
        rIntObjNormal[0] = rGeometry[1].Y() - rGeometry[0].Y();
        rIntObjNormal[1] = rGeometry[0].X() - rGeometry[1].X();
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
        } else if constexpr (TDim == 3) {
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
        const Geometry<Node >& rEdge,
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
        const Geometry<Node >& rEdge,
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

    template<std::size_t TDim>
    void CalculateDiscontinuousDistanceToSkinProcess<TDim>::CheckAndCorrectEdgeDistances()
    {

        KRATOS_TRY;

        if (!mAreNeighboursComputed) {
            FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType> find_nodal_elems_process(mrVolumePart);
            find_nodal_elems_process.Execute();
            mAreNeighboursComputed = true;
        }

        auto p_global_ptr_comm = CreatePointerCommunicator();

        // Proxy to retrieve cut edges from neighbours in other partitions. Works as a lambda in serial runs.
        auto get_neighbour_cut_edges = p_global_ptr_comm->Apply([this](GlobalPointer<Element> pElem) {

            std::unordered_map<IndexType, IndexType> global_id_to_local;
            for (IndexType i = 0; i < mNumNodes; ++i) {
                global_id_to_local[pElem->GetGeometry()[i].Id()] = i;
            }

            auto elem_dist = pElem->GetValue(*mpElementalDistancesVariable);
            std::vector<std::vector<IndexType>> cut_edges_vector;
            const auto edges_container = pElem->GetGeometry().GenerateEdges();
            for (std::size_t i_edge = 0; i_edge < mNumEdges; ++i_edge) {
                IndexType edge_i_id = edges_container[i_edge][0].Id();
                IndexType edge_j_id = edges_container[i_edge][1].Id();
                double i_distance = elem_dist[global_id_to_local[edge_i_id]];
                double j_distance = elem_dist[global_id_to_local[edge_j_id]];
                bool is_edge_cut = (i_distance * j_distance) < 0.0;
                std::vector<IndexType> edge_vector;
                edge_vector.push_back(std::min(edge_i_id, edge_j_id));
                edge_vector.push_back(std::max(edge_i_id, edge_j_id));
                if (is_edge_cut) {
                    cut_edges_vector.push_back(edge_vector);
                }
            }
            return cut_edges_vector;
        });

        auto check_edge_against_distances = [&, get_neighbour_cut_edges] (const Element::GeometryType::GeometriesArrayType& rEdgeContainer, Vector& rCuteEdgeVector) {
            for (std::size_t i_edge = 0; i_edge < mNumEdges; ++i_edge) {
                // Edge initially cut with distances exactly 0
                if (rCuteEdgeVector[i_edge] >= 0) {
                    IndexType edge_i_id = std::min(rEdgeContainer[i_edge][0].Id(), rEdgeContainer[i_edge][1].Id());
                    IndexType edge_j_id = std::max(rEdgeContainer[i_edge][0].Id(), rEdgeContainer[i_edge][1].Id());

                    // Get neighbour elements in a single vector without repetitions.
                    std::unordered_set<Element::WeakPointer, SharedPointerHasher<Element::WeakPointer>, SharedPointerComparator<Element::WeakPointer>> all_neighbour_set;
                    const auto& r_edge_containter_0  = rEdgeContainer[i_edge][0].GetValue(NEIGHBOUR_ELEMENTS).GetContainer();
                    for(const auto& p_neigh : r_edge_containter_0) {
                        all_neighbour_set.insert(p_neigh);
                    }
                    const auto& r_edge_containter_1  = rEdgeContainer[i_edge][1].GetValue(NEIGHBOUR_ELEMENTS).GetContainer();
                    for(const auto& p_neigh : r_edge_containter_1) {
                        all_neighbour_set.insert(p_neigh);
                    }

                    bool is_edge_cut = false;
                    for (const auto& g_ptr_neigh : all_neighbour_set) {
                        if (!is_edge_cut) {
                            auto cut_edges = get_neighbour_cut_edges.Get(g_ptr_neigh);
                            for (auto cut_edge: cut_edges) {
                                if (cut_edge[0]==edge_i_id && cut_edge[1]==edge_j_id) {
                                    is_edge_cut = true;
                                    break;
                                }
                            }
                        }
                    }

                    // Correcting edge if not intersected after applying an epsilon
                    if (!is_edge_cut) {
                        rCuteEdgeVector[i_edge] = -1;
                    }
                }
            }
        };

        for (auto &r_elem : mrVolumePart.Elements())  {
            const auto edges_container = r_elem.GetGeometry().GenerateEdges();
            auto &r_cut_edge_vector = r_elem.GetValue(*mpElementalEdgeDistancesVariable);
            check_edge_against_distances(edges_container, r_cut_edge_vector);

        }
        if (mCalculateElementalEdgeDistancesExtrapolated) {
            for (auto &r_elem : mrVolumePart.Elements())  {
                const auto edges_container = r_elem.GetGeometry().GenerateEdges();
                auto &r_cut_edge_extra_vector = r_elem.GetValue(*mpElementalEdgeDistancesExtrapolatedVariable);
                check_edge_against_distances(edges_container, r_cut_edge_extra_vector);
            }
        }

        KRATOS_CATCH(" ");
    }

    template<std::size_t TDim>
    GlobalPointerCommunicator<Element>::Pointer CalculateDiscontinuousDistanceToSkinProcess<TDim>::CreatePointerCommunicator()
    {
        KRATOS_TRY;


        const auto &r_comm = mrVolumePart.GetCommunicator().GetDataCommunicator();

        std::vector<int> elem_indices;
        elem_indices.reserve(mrVolumePart.NumberOfElements());
        for (const auto& r_elem : mrVolumePart.Elements()) {
            elem_indices.push_back(r_elem.Id());
        }

        const auto g_ptr_elem_map = GlobalPointerUtilities::RetrieveGlobalIndexedPointersMap(mrVolumePart.Elements(), elem_indices, r_comm);

        GlobalPointersVector<Element> g_ptr_elem_list;
        for(auto& r_item: g_ptr_elem_map) {
            g_ptr_elem_list.push_back(r_item.second);
        }

        for (auto& r_node : mrVolumePart.Nodes()) {
            for (const auto& g_ptr : r_node.GetValue(NEIGHBOUR_ELEMENTS).GetContainer()) {
                g_ptr_elem_list.push_back(g_ptr);
            }
        }

        g_ptr_elem_list.Unique();

        return Kratos::make_shared<GlobalPointerCommunicator<Element>>(r_comm, g_ptr_elem_list);

        KRATOS_CATCH(" ");
    }

    template class Kratos::CalculateDiscontinuousDistanceToSkinProcess<2>;
    template class Kratos::CalculateDiscontinuousDistanceToSkinProcess<3>;

}  // namespace Kratos.
