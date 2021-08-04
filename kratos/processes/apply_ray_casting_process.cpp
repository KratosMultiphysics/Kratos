//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "geometries/plane_3d.h"
#include "processes/apply_ray_casting_process.h"
#include "utilities/geometry_utilities.h"
#include "utilities/intersection_utilities.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{

	template<std::size_t TDim>
	ApplyRayCastingProcess<TDim>::ApplyRayCastingProcess(
		ModelPart& rVolumePart,
		ModelPart& rSkinPart)
		: mpFindIntersectedObjectsProcess(new FindIntersectedGeometricalObjectsProcess(rVolumePart, rSkinPart)),
		  mIsSearchStructureAllocated(true)
	{
	}

	template <std::size_t TDim>
	ApplyRayCastingProcess<TDim>::ApplyRayCastingProcess(
		ModelPart &rVolumePart,
		ModelPart &rSkinPart,
		const double RelativeTolerance)
		: mRelativeTolerance(RelativeTolerance),
		  mpFindIntersectedObjectsProcess(new FindIntersectedGeometricalObjectsProcess(rVolumePart, rSkinPart)),
		  mIsSearchStructureAllocated(true)
	{
	}

	template <std::size_t TDim>
	ApplyRayCastingProcess<TDim>::ApplyRayCastingProcess(
		FindIntersectedGeometricalObjectsProcess &TheFindIntersectedObjectsProcess,
		const double RelativeTolerance)
		: mRelativeTolerance(RelativeTolerance),
		  mpFindIntersectedObjectsProcess(&TheFindIntersectedObjectsProcess),
		  mIsSearchStructureAllocated(false)
	{
	}

	template<std::size_t TDim>
	ApplyRayCastingProcess<TDim>::~ApplyRayCastingProcess()
	{
        if(mIsSearchStructureAllocated)
            delete mpFindIntersectedObjectsProcess;
	}

	template<std::size_t TDim>
	void ApplyRayCastingProcess<TDim>::Execute()
	{
        if(mIsSearchStructureAllocated) // we have not initialized it yet
            mpFindIntersectedObjectsProcess->ExecuteInitialize();

		this->SetRayCastingTolerances();

		ModelPart& ModelPart1 = mpFindIntersectedObjectsProcess->GetModelPart1();

		block_for_each(ModelPart1.Nodes(), [&](Node<3>& rNode){
			double &r_node_distance = rNode.FastGetSolutionStepValue(DISTANCE);
			const double ray_distance = this->DistancePositionInSpace(rNode);
			if (ray_distance * r_node_distance < 0.0) {
				r_node_distance = -r_node_distance;
			}
		});
	}

	template<std::size_t TDim>
	double ApplyRayCastingProcess<TDim>::DistancePositionInSpace(const Node<3> &rNode)
	{
		array_1d<double,TDim> distances;
		unsigned int n_ray_pos(0), n_ray_neg(0);
        IntersectionsContainerType intersections;
		const auto &r_coords = rNode.Coordinates();

		// Loop the x,y and z (3D) ray directions
        for (unsigned int i_direction = 0; i_direction < TDim; i_direction++){
			// Initialize the current direction distance
			distances[i_direction] = 1.0;

            // Creating the ray
            double ray[3] = {r_coords[0], r_coords[1], r_coords[2]};
			auto &rp_octree = mpFindIntersectedObjectsProcess->GetOctreePointer();
			rp_octree->NormalizeCoordinates(ray);
			ray[i_direction] = 0; // Starting from the lower extreme

            this->GetRayIntersections(ray, i_direction, intersections);

            int ray_color = 1;
            auto i_intersection = intersections.begin();
            while (i_intersection != intersections.end()) {
                double int_d = r_coords[i_direction] - i_intersection->first; // Octree ray intersection distance
                if (int_d > mEpsilon) {
                    ray_color = -ray_color;
                    distances[i_direction] = int_d;
                } else if (int_d > -mEpsilon) {
                    distances[i_direction] = 0.0;
                    break;
                } else {
                    if (distances[i_direction] > -int_d) {
                        distances[i_direction] = -int_d;
					}
                    break;
                }

                i_intersection++;
            }

            distances[i_direction] *= ray_color;
			if (ray_color == -1) {
				n_ray_neg++;
			} else {
				n_ray_pos++;
			}
        }

		// Check the obtained cartesian ray colors
		// If this situation happens, do the "evolved Predator" raycasting to vote
		if (n_ray_neg != 0 && n_ray_pos != 0) {
			this->ComputeExtraRayColors(r_coords, distances);
		}

        double distance = (std::abs(distances[0]) > std::abs(distances[1])) ? distances[1] : distances[0];
		if (TDim == 3){
        	distance = (std::abs(distance) > std::abs(distances[2])) ? distances[2] : distance;
		}

        return distance;
	}

	template<std::size_t TDim>
	void ApplyRayCastingProcess<TDim>::ComputeExtraRayColors(
		const array_1d<double,3> &rCoords,
        array_1d<double,TDim> &rDistances)
	{
		// Set the extra ray origins
		array_1d<array_1d<double,3>, (TDim == 3) ? 9 : 5> extra_ray_origs;
		this->GetExtraRayOrigins(rCoords, extra_ray_origs);

		// Get the pointer to the base Octree binary
		auto &rp_octree = mpFindIntersectedObjectsProcess->GetOctreePointer();

		// Loop the extra rays to compute its color
		unsigned int n_ray_pos = 0; // Positive rays counter
		unsigned int n_ray_neg = 0; // Negative rays counter
		IntersectionsContainerType intersections; // Ray intersections container initialization
		for (unsigned int i_direction = 0; i_direction < TDim; ++i_direction) {
			for (unsigned int i_ray = 0; i_ray < extra_ray_origs.size(); ++i_ray) {
				// Creating the ray
				const auto aux_ray = extra_ray_origs[i_ray];
				double ray[3] = {aux_ray[0], aux_ray[1], aux_ray[2]};
				rp_octree->NormalizeCoordinates(ray);
				ray[i_direction] = 0; // Starting from the lower extreme
				this->CorrectExtraRayOrigin(ray); // Avoid extra ray normalized coordinates to be larger than 1 or 0
				this->GetRayIntersections(ray, i_direction, intersections);

				// Compute the extra rays intersections
				int ray_color = 1;
				auto i_intersection = intersections.begin();
				while (i_intersection != intersections.end()) {
					const double int_d = extra_ray_origs[i_ray][i_direction] - i_intersection->first; // Octree ray intersection distance
					if (int_d > mEpsilon) {
						ray_color = -ray_color;
					} else if (int_d > - mEpsilon) {
						break;
					} else {
						break;
					}
					i_intersection++;
				}

				// Update the extra rays color counters
				if (ray_color == -1) {
					n_ray_neg++;
				} else {
					n_ray_pos++;
				}
			}

		}

		// Do the extra rays voting
		int ray_color = n_ray_neg > n_ray_pos ? -1 : 1;
		for (unsigned int i_direction = 0; i_direction < TDim; ++i_direction) {
			rDistances[i_direction] = std::abs(rDistances[i_direction]) * ray_color;
		}
	}

	template<>
	void ApplyRayCastingProcess<2>::GetExtraRayOrigins(
        const array_1d<double,3> &rCoords,
		array_1d<array_1d<double,3>, 5> &rExtraRayOrigs)
	{
		if (rExtraRayOrigs.size() != 5) {
			rExtraRayOrigs.resize(5);
		}

		array_1d<double,3> aux_1; aux_1[0] = rCoords[0]; aux_1[1] = rCoords[1]; aux_1[2] = rCoords[2];
		array_1d<double,3> aux_2; aux_2[0] = rCoords[0] + 2*mExtraRayOffset; aux_2[1] = rCoords[1] + mExtraRayOffset; aux_2[2] = rCoords[2];
		array_1d<double,3> aux_3; aux_3[0] = rCoords[0] + mExtraRayOffset; aux_3[1] = rCoords[1] + 2*mExtraRayOffset; aux_3[2] = rCoords[2];
		array_1d<double,3> aux_4; aux_4[0] = rCoords[0] - 2*mExtraRayOffset; aux_4[1] = rCoords[1] - mExtraRayOffset; aux_4[2] = rCoords[2];
		array_1d<double,3> aux_5; aux_5[0] = rCoords[0] - mExtraRayOffset; aux_5[1] = rCoords[1] - 2*mExtraRayOffset; aux_5[2] = rCoords[2];

		rExtraRayOrigs[0] = aux_1;
		rExtraRayOrigs[1] = aux_2;
		rExtraRayOrigs[2] = aux_3;
		rExtraRayOrigs[3] = aux_4;
		rExtraRayOrigs[4] = aux_5;
	}

	template<>
	void ApplyRayCastingProcess<3>::GetExtraRayOrigins(
        const array_1d<double,3> &rCoords,
		array_1d<array_1d<double,3>,9> &rExtraRayOrigs)
	{
		if (rExtraRayOrigs.size() != 9) {
			rExtraRayOrigs.resize(9);
		}

		array_1d<double,3> aux_1; aux_1[0] = rCoords[0]; aux_1[1] = rCoords[1]; aux_1[2] = rCoords[2];
		array_1d<double,3> aux_2; aux_2[0] = rCoords[0] + 2*mExtraRayOffset; aux_2[1] = rCoords[1] + mExtraRayOffset; aux_2[2] = rCoords[2] - mExtraRayOffset;
		array_1d<double,3> aux_3; aux_3[0] = rCoords[0] + mExtraRayOffset; aux_3[1] = rCoords[1] + 2*mExtraRayOffset; aux_3[2] = rCoords[2] - mExtraRayOffset;
		array_1d<double,3> aux_4; aux_4[0] = rCoords[0] - 2*mExtraRayOffset; aux_4[1] = rCoords[1] - mExtraRayOffset; aux_4[2] = rCoords[2] - mExtraRayOffset;
		array_1d<double,3> aux_5; aux_5[0] = rCoords[0] - mExtraRayOffset; aux_5[1] = rCoords[1] - 2*mExtraRayOffset; aux_5[2] = rCoords[2] - mExtraRayOffset;
		array_1d<double,3> aux_6; aux_6[0] = rCoords[0] + 2*mExtraRayOffset; aux_6[1] = rCoords[1] + mExtraRayOffset; aux_6[2] = rCoords[2] + mExtraRayOffset;
		array_1d<double,3> aux_7; aux_7[0] = rCoords[0] + mExtraRayOffset; aux_7[1] = rCoords[1] + 2*mExtraRayOffset; aux_7[2] = rCoords[2] + mExtraRayOffset;
		array_1d<double,3> aux_8; aux_8[0] = rCoords[0] - 2*mExtraRayOffset; aux_8[1] = rCoords[1] - mExtraRayOffset; aux_8[2] = rCoords[2] + mExtraRayOffset;
		array_1d<double,3> aux_9; aux_9[0] = rCoords[0] - mExtraRayOffset; aux_9[1] = rCoords[1] - 2*mExtraRayOffset; aux_9[2] = rCoords[2] + mExtraRayOffset;

		rExtraRayOrigs[0] = aux_1;
		rExtraRayOrigs[1] = aux_2;
		rExtraRayOrigs[2] = aux_3;
		rExtraRayOrigs[3] = aux_4;
		rExtraRayOrigs[4] = aux_5;
		rExtraRayOrigs[5] = aux_6;
		rExtraRayOrigs[6] = aux_7;
		rExtraRayOrigs[7] = aux_8;
		rExtraRayOrigs[8] = aux_9;
	}

	template<std::size_t TDim>
	void ApplyRayCastingProcess<TDim>::CorrectExtraRayOrigin(double* ExtraRayCoords)
	{
		for (unsigned int d = 0; d < 3; ++d) {
			if (ExtraRayCoords[d] > 1.0) {
				ExtraRayCoords[d] = 1.0;
			} else if (ExtraRayCoords[d] < 0.0) {
				ExtraRayCoords[d] = 0.0;
			}
		}
	}

	template<std::size_t TDim>
	void ApplyRayCastingProcess<TDim>::GetRayIntersections(
		const double* ray,
		const unsigned int direction,
		std::vector<std::pair<double,Element::GeometryType*> >& rIntersections)
	{
		// This function passes the ray through the model and gives the hit point to all objects in its way
        // Ray is of dimension (3) normalized in (0,1)^3 space
        // Direction can be 0,1,2 which are x,y and z respectively

        // First clearing the intersections points vector
        rIntersections.clear();

		// Get the octree from the parent discontinuous distance process
        auto &rp_octree = mpFindIntersectedObjectsProcess->GetOctreePointer();

		// Compute the normalized ray key
		OctreeType::key_type ray_key[3] = {rp_octree->CalcKeyNormalized(ray[0]), rp_octree->CalcKeyNormalized(ray[1]), rp_octree->CalcKeyNormalized(ray[2])};

		// Getting the entrance cell from lower extreme
        OctreeType::key_type cell_key[3];
        auto cell = rp_octree->pGetCell(ray_key);
        while (cell) {
			// Get the current cell intersections
            const int cell_int = this->GetCellIntersections(cell, ray, ray_key, direction, rIntersections);
			KRATOS_ERROR_IF(cell_int != 0)
				<< "Error in GetCellIntersections for ray [" << ray[0] << "," << ray[1] << "," << ray[2] << "] with direction " << direction << std::endl;
			// And if it exists, go to the next cell
            if (cell->GetNeighbourKey(1 + direction * 2, cell_key)) {
                ray_key[direction] = cell_key[direction];
                cell = rp_octree->pGetCell(ray_key);
                ray_key[direction] -= 1 ; // The key returned by GetNeighbourKey is inside the cell (minkey +1), to ensure that the corresponding cell get in pGetCell is the right one.
            } else {
                cell = NULL;
			}
        }

        // Eliminate the repeated intersecting objects
        if (!rIntersections.empty()) {
            // Sort
            std::sort(rIntersections.begin(), rIntersections.end());
            // Unique
            auto i_begin = rIntersections.begin();
            auto i_intersection = rIntersections.begin();
            while (++i_begin != rIntersections.end()) {
                // considering the very near points as the same points
                if (std::abs(i_begin->first - i_intersection->first) > mEpsilon) // if the hit points are far enough they are not the same
                    *(++i_intersection) = *i_begin;
            }
            rIntersections.resize((++i_intersection) - rIntersections.begin());
        }
	}

	template<std::size_t TDim>
	int  ApplyRayCastingProcess<TDim>::GetCellIntersections(
		OctreeType::cell_type* cell,
		const double* ray,
		OctreeType::key_type* ray_key,
		const unsigned int direction,
		std::vector<std::pair<double, Element::GeometryType*> > &rIntersections)
	{
		//This function passes the ray through the cell and gives the hit point to all objects in its way
		//ray is of dimension (3) normalized in (0,1)^3 space
		// direction can be 0,1,2 which are x,y and z respectively
		typedef OctreeType::cell_type::object_container_type object_container_type;
		auto objects = (cell->pGetObjects());

		// There are no intersection in empty cells
		if (objects->empty()){
			return 0;
		}

		// Calculating the two extreme of the ray segment inside the cell
		double ray_point1[3] = {ray[0], ray[1], ray[2]};
		double ray_point2[3] = {ray[0], ray[1], ray[2]};
		double normalized_coordinate;

		auto &rp_octree = mpFindIntersectedObjectsProcess->GetOctreePointer();
		rp_octree->CalculateCoordinateNormalized(ray_key[direction], normalized_coordinate);
		ray_point1[direction] = normalized_coordinate;
		ray_point2[direction] = ray_point1[direction] + rp_octree->CalcSizeNormalized(cell);
		rp_octree->ScaleBackToOriginalCoordinate(ray_point1);
		rp_octree->ScaleBackToOriginalCoordinate(ray_point2);

		// This is a workaround to avoid the z-component in 2D
		if (TDim == 2){
			ray_point1[2] = 0.0;
			ray_point2[2] = 0.0;
		}

		for (object_container_type::iterator i_object = objects->begin(); i_object != objects->end(); ++i_object){
			double intersection[3]={0.0, 0.0, 0.0};
			const int is_intersected = ComputeRayIntersection((*i_object)->GetGeometry(), ray_point1, ray_point2, intersection);
			if (is_intersected == 1){ // There is an intersection but not coplanar
				rIntersections.push_back(std::pair<double, Element::GeometryType*>(intersection[direction], &((*i_object)->GetGeometry())));
			}
		}

		return 0;
	}

	template<>
	int ApplyRayCastingProcess<2>::ComputeRayIntersection(
		Element::GeometryType& rGeometry,
        const double* pRayPoint1,
        const double* pRayPoint2,
        double* pIntersectionPoint)
	{
		// Auxiliar arrays
		array_1d<double,3> ray_pt_1;
		array_1d<double,3> ray_pt_2;
		for (unsigned int i = 0; i < 3; ++i){
			ray_pt_1[i] = pRayPoint1[i];
			ray_pt_2[i] = pRayPoint2[i];
		}

		// Call the line - line intersection util
		array_1d<double,3> int_pt = ZeroVector(3);
		const double tolerance = mEpsilon;
		const int is_intersected =  IntersectionUtilities::ComputeLineLineIntersection(
			rGeometry,
			ray_pt_1,
			ray_pt_2,
			int_pt,
			tolerance);

		// Convert the auxiliar intersection point to the original type
		for (unsigned int i = 0; i < 3; ++i){
			pIntersectionPoint[i] = int_pt[i];
		}

		return is_intersected;
	}

	template<>
	int ApplyRayCastingProcess<3>::ComputeRayIntersection(
		Element::GeometryType& rGeometry,
        const double* pRayPoint1,
        const double* pRayPoint2,
        double* pIntersectionPoint)
	{
		// Auxiliar arrays
		array_1d<double,3> ray_pt_1;
		array_1d<double,3> ray_pt_2;
		for (unsigned int i = 0; i < 3; ++i){
			ray_pt_1[i] = pRayPoint1[i];
			ray_pt_2[i] = pRayPoint2[i];
		}

		// Call the line - triangle intersection util
		array_1d<double,3> int_pt = ZeroVector(3);
		const double tolerance = mEpsilon;
		const int is_intersected = IntersectionUtilities::ComputeTriangleLineIntersection(
			rGeometry,
			ray_pt_1,
			ray_pt_2,
			int_pt,
			tolerance);

		// Convert the auxiliar intersection point to the original type
		for (unsigned int i = 0; i < 3; ++i){
			pIntersectionPoint[i] = int_pt[i];
		}

		return is_intersected;
	}

	/// Turn back information as a string.
	template<std::size_t TDim>
	std::string ApplyRayCastingProcess<TDim>::Info() const
	{
		return "ApplyRayCastingProcess";
	}

	/// Print information about this object.
	template<std::size_t TDim>
	void ApplyRayCastingProcess<TDim>::PrintInfo(std::ostream& rOStream) const
	{
		rOStream << Info();
	}

	/// Print object's data.
	template<std::size_t TDim>
	void ApplyRayCastingProcess<TDim>::PrintData(std::ostream& rOStream) const
	{
	}

	template<std::size_t TDim>
	void ApplyRayCastingProcess<TDim>::CalculateCharacteristicLength()
	{
		ModelPart& model_part1 = mpFindIntersectedObjectsProcess->GetModelPart1();
        Point min_point(0.00, 0.00, 0.00);
        Point max_point(0.00, 0.00, 0.00);

        for(auto& node : model_part1.Nodes()) {
            for(std::size_t i = 0; i<3; i++) {
                min_point[i]  =  (min_point[i]  >  node[i] ) ?  node[i] : min_point[i];
                max_point[i] =  (max_point[i] <  node[i] ) ?  node[i] : max_point[i];
            }
        }

		mCharacteristicLength = norm_2(max_point - min_point);
		KRATOS_ERROR_IF(mCharacteristicLength < std::numeric_limits<double>::epsilon()) << "Domain characteristic length is close to zero. Check if there is any node in the model part." << std::endl;
	}

	template <std::size_t TDim>
	void ApplyRayCastingProcess<TDim>::SetRayCastingTolerances()
	{
		this->CalculateCharacteristicLength();
		mEpsilon = mRelativeTolerance * mCharacteristicLength;
		mExtraRayOffset = 2.0 * mRelativeTolerance * mCharacteristicLength;
	}



	template class Kratos::ApplyRayCastingProcess<2>;
	template class Kratos::ApplyRayCastingProcess<3>;

}  // namespace Kratos.
