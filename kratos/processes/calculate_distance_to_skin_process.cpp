//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//


// System includes


// External includes


// Project includes
#include "geometries/plane_3d.h"
#include "processes/calculate_distance_to_skin_process.h"
#include "utilities/geometry_utilities.h"
#include "utilities/intersection_utilities.h"


namespace Kratos
{

	template<std::size_t TDim>
	CalculateDistanceToSkinProcess<TDim>::CalculateDistanceToSkinProcess(ModelPart& rVolumePart, ModelPart& rSkinPart)
		: CalculateDiscontinuousDistanceToSkinProcess<TDim>(rVolumePart, rSkinPart)
	{
	}

	template<std::size_t TDim>
	CalculateDistanceToSkinProcess<TDim>::~CalculateDistanceToSkinProcess()
	{
	}

	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::Initialize()
	{
		CalculateDiscontinuousDistanceToSkinProcess<TDim>::Initialize();
		this->InitializeNodalDistances();
	}

	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::InitializeNodalDistances()
	{
		ModelPart& ModelPart1 = (CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess).GetModelPart1();

		for (auto& node : ModelPart1.Nodes())
		{
			node.GetSolutionStepValue(DISTANCE) = std::numeric_limits<double>::max();
		}
	}

	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::CalculateDistances(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects)
	{
		this->CalculateElementalDistances(rIntersectedObjects);
		this->CalculateNodalDistances();
		this->CalculateNodesDistances();
	}

	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::CalculateElementalDistances(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects)
	{
		const int number_of_elements = (CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess.GetModelPart1()).NumberOfElements();
		auto& r_elements = (CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess.GetModelPart1()).ElementsArray();

		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < number_of_elements; ++i)
		{
			Element &r_element = *(r_elements[i]);
			PointerVector<GeometricalObject>& r_element_intersections = rIntersectedObjects[i]; 

			// Check if the element has intersections
			if (r_element_intersections.empty())
			{
				r_element.Set(TO_SPLIT, false);
			} 
			else 
			{
				// This function assumes tetrahedra element and triangle intersected object as input at this moment
				constexpr int number_of_tetrahedra_points = TDim + 1;
				constexpr double epsilon = std::numeric_limits<double>::epsilon();
				Vector &elemental_distances = r_element.GetValue(ELEMENTAL_DISTANCES);

				if (elemental_distances.size() != number_of_tetrahedra_points)
					elemental_distances.resize(number_of_tetrahedra_points, false);

				for (int i = 0; i < number_of_tetrahedra_points; i++)
				{
					elemental_distances[i] = this->CalculateDistanceToNode(r_element, i, r_element_intersections, epsilon);
				}

				bool has_positive_distance = false;
				bool has_negative_distance = false;
				for (int i = 0; i < number_of_tetrahedra_points; i++)
					if (elemental_distances[i] > epsilon)
						has_positive_distance = true;
					else
						has_negative_distance = true;

				r_element.Set(TO_SPLIT, has_positive_distance && has_negative_distance);
			}
		}
	}

	template<std::size_t TDim>
	double CalculateDistanceToSkinProcess<TDim>::CalculateDistanceToNode(
		Element& rElement1,
		const int NodeIndex,
		PointerVector<GeometricalObject>& rIntersectedObjects,
		const double Epsilon)
	{
		double result_distance = std::numeric_limits<double>::max();
		for (auto it_int_obj : rIntersectedObjects.GetContainer()) {
			const auto &r_int_obj_geom = it_int_obj->GetGeometry();
			const double distance = this->CalculatePointDistance(r_int_obj_geom, rElement1.GetGeometry()[NodeIndex]);
			if (std::abs(result_distance) > distance)
			{
				if (distance < Epsilon) {
					result_distance = -Epsilon;
				} else {
					result_distance = distance;
					std::vector<array_1d<double,3>> plane_pts;
					for (unsigned int i_node = 0; i_node < r_int_obj_geom.PointsNumber(); ++i_node){
						plane_pts.push_back(r_int_obj_geom[i_node]);
					}
					Plane3D plane = this->SetIntersectionPlane(plane_pts);

					if (plane.CalculateSignedDistance(rElement1.GetGeometry()[NodeIndex]) < 0){
						result_distance = -result_distance;
					}
				}
			}
		}
		return result_distance;
	}

	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::CalculateNodalDistances()
	{
		ModelPart& ModelPart1 = (CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess).GetModelPart1();

		constexpr int number_of_tetrahedra_points = TDim + 1;
		for (auto& element : ModelPart1.Elements())
		{
			if (element.Is(TO_SPLIT))
			{
				const auto& r_elemental_distances = element.GetValue(ELEMENTAL_DISTANCES);
				for (int i = 0; i < number_of_tetrahedra_points; i++)
				{
					Node<3>& r_node = element.GetGeometry()[i];
					double& r_distance = r_node.GetSolutionStepValue(DISTANCE);
					if (std::abs(r_distance) > std::abs(r_elemental_distances[i])){
						r_distance = r_elemental_distances[i];
					}
				}
			}
		}
	}

	//TODO: This method has been adapted from the previous implementation. It is still pending to update it.
	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::CalculateNodesDistances()
	{
		ModelPart& ModelPart1 = (CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess).GetModelPart1();

		#pragma omp parallel for
		for(int k = 0 ; k < static_cast<int>(ModelPart1.NumberOfNodes()); ++k)
        {
			ModelPart::NodesContainerType::iterator itNode = ModelPart1.NodesBegin() + k;
            this->CalculateNodeDistance(*itNode);
        }

	}

	//TODO: This method has been adapted from the previous implementation. It is still pending to update it.
	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::CalculateNodeDistance(Node<3>& rNode)
	{
		double distance = DistancePositionInSpace(rNode);
		double& node_distance =  rNode.GetSolutionStepValue(DISTANCE);

		//const double epsilon = 1.00e-12;
		//if(fabs(node_distance) > fabs(distance))
		//    node_distance = distance;
		/*else*/ 
		
		if (distance * node_distance < 0.0) { // assigning the correct sign
			node_distance = -node_distance;
		}
	}

	//TODO: This method has been adapted from the previous implementation. It is still pending to update it.
	template<std::size_t TDim>
	double CalculateDistanceToSkinProcess<TDim>::DistancePositionInSpace(const Node<3> &rNode)
	{

		typedef Element::GeometryType triangle_type;
        typedef std::vector<std::pair<double, triangle_type*> > intersections_container_type;

        intersections_container_type intersections;

        const double epsilon = 1e-12;
		array_1d<double,TDim> distances;

        for (int i_direction = 0; i_direction < TDim; i_direction++){
			// Initialize the current direction distance
			// distances[i_direction] = std::numeric_limits<double>::max();
			distances[i_direction] = 1.0;

            // Creating the ray
			const array_1d<double,3> coords = rNode.Coordinates();
            double ray[3] = {coords[0], coords[1], coords[2]};
			if (rNode.Id() == 22){
				std::cout << "Ray: [" << ray[0] << " , " << ray[1] << " , " << ray[2] << "]" << std::endl;
			}

			OctreeType* pOctree = CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateDiscontinuousDistanceToSkinProcess::mFindIntersectedObjectsProcess.GetOctreePointer();
            pOctree->NormalizeCoordinates(ray);
			// This is a workaround to avoid the z-component in 2D
			// The octree sets it when normalizing the coordinates
			if (TDim == 2){
				ray[2] = 0.0;
			}
            ray[i_direction] = 0; // starting from the lower extreme

			if (rNode.Id() == 22){
				std::cout << "Ray: [" << ray[0] << " , " << ray[1] << " , " << ray[2] << "]" << std::endl;
			}
            this->GetRayIntersections(ray, i_direction, intersections);

            int ray_color = 1;
            std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = intersections.begin();
            while (i_intersection != intersections.end())
			{
				if (rNode.Id() == 22){
					KRATOS_WATCH(*(i_intersection->second))
				}
                double d = coords[i_direction] - i_intersection->first;
                if (d > epsilon) {
                    ray_color = -ray_color;
                    distances[i_direction] = d;
                } else if (d > -epsilon) {
                    distances[i_direction] = 0.0;
                    break;
                } else {
                    if (distances[i_direction] > -d) {
                        distances[i_direction] = -d;
					}
                    break;
                }

                i_intersection++;
            }

			if (intersections.size() == 1){
				std::cout << "Node: " << rNode.Id() << " has " << intersections.size() << " interserctions. d = " << distances[i_direction] << " direction: " << i_direction << " ray: [" << ray[0] << " , " << ray[1] << "]" << std::endl;
			}

            distances[i_direction] *= ray_color;
        }

        double distance = (std::abs(distances[0]) > std::abs(distances[1])) ? distances[1] : distances[0];
		if (TDim == 3){
        	distance = (std::abs(distance) > std::abs(distances[2])) ? distances[2] : distance;
		}

        return distance;
	}

	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::GetRayIntersections(
		double* ray,
		int direction,
		std::vector<std::pair<double,Element::GeometryType*> >& intersections)
	{
		//This function passes the ray through the model and gives the hit point to all objects in its way
        //ray is of dimension (3) normalized in (0,1)^3 space
        // direction can be 0,1,2 which are x,y and z respectively

        const double epsilon = 1.00e-12;

        // first clearing the intersections points vector
        intersections.clear();

        //OctreeType* octree = &mOctree;
        OctreeType* pOctree = CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess.GetOctreePointer();

        OctreeType::key_type ray_key[3] = {pOctree->CalcKeyNormalized(ray[0]), pOctree->CalcKeyNormalized(ray[1]), pOctree->CalcKeyNormalized(ray[2])};
        OctreeType::key_type cell_key[3];

        // getting the entrance cell from lower extreme
        OctreeType::cell_type* cell = pOctree->pGetCell(ray_key);

		std::cout  << "ray_key: [" << ray_key[0] << "," << ray_key[1] << "," << ray_key[2] << "]" << std::endl;

		unsigned int i_cell = 0;
        while (cell)
		{
			i_cell++;
            this->GetCellIntersections(cell, ray, ray_key, direction, intersections);
            std::cout << intersections.size() << std::endl;
			// go to the next cell
            if (cell->GetNeighbourKey(1 + direction * 2, cell_key))
			{
                ray_key[direction] = cell_key[direction];
                cell = pOctree->pGetCell(ray_key);
                ray_key[direction] -= 1 ;//the key returned by GetNeighbourKey is inside the cell (minkey +1), to ensure that the corresponding
                std::cout  << "\tray_key: [" << ray_key[0] << "," << ray_key[1] << "," << ray_key[2] << "]" << std::endl;
				//cell get in pGetCell is the right one.
            } else
                cell = NULL;
        }

		// KRATOS_WATCH(i_cell)

        // now eliminating the repeated objects
        if (!intersections.empty())
		{
            //sort
            std::sort(intersections.begin(), intersections.end());
            // unique
            std::vector<std::pair<double, Element::GeometryType*> >::iterator i_begin = intersections.begin();
            std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = intersections.begin();
            while (++i_begin != intersections.end())
			{
                // considering the very near points as the same points
                if (std::abs(i_begin->first - i_intersection->first) > epsilon) // if the hit points are far enough they are not the same
                    *(++i_intersection) = *i_begin;
            }
            intersections.resize((++i_intersection) - intersections.begin());
        }
	}

	template<std::size_t TDim>
	int  CalculateDistanceToSkinProcess<TDim>::GetCellIntersections(
		OctreeType::cell_type* cell,
		double* ray,
		OctreeType::key_type* ray_key,
		int direction,
		std::vector<std::pair<double, Element::GeometryType*> >& intersections)
	{
		//This function passes the ray through the cell and gives the hit point to all objects in its way
		//ray is of dimension (3) normalized in (0,1)^3 space
		// direction can be 0,1,2 which are x,y and z respectively

		typedef OctreeType::cell_type::object_container_type object_container_type;

		object_container_type* objects = (cell->pGetObjects());

		// There are no intersection in empty cells
		if (objects->empty())
			return 0;

		// calculating the two extreme of the ray segment inside the cell
		double ray_point1[3] = {ray[0], ray[1], ray[2]};
		double ray_point2[3] = {ray[0], ray[1], ray[2]};
		double normalized_coordinate;

		// OctreeType* pOctree = (CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess.GetOctree()).get();
		OctreeType* pOctree = CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess.GetOctreePointer();

		pOctree->CalculateCoordinateNormalized(ray_key[direction], normalized_coordinate);
		ray_point1[direction] = normalized_coordinate;
		ray_point2[direction] = ray_point1[direction] + pOctree->CalcSizeNormalized(cell);

		pOctree->ScaleBackToOriginalCoordinate(ray_point1);
		pOctree->ScaleBackToOriginalCoordinate(ray_point2);

		// This is a workaround to avoid the z-component in 2D
		if (TDim == 2){
			ray_point1[2] = 0.0;
			ray_point2[2] = 0.0;
		}

		for (object_container_type::iterator i_object = objects->begin(); i_object != objects->end(); i_object++)
		{
			double intersection[3]={0.00,0.00,0.00};

			int is_intersected = ComputeRayIntersection((*i_object)->GetGeometry(), ray_point1, ray_point2, intersection);
			if (is_intersected == 1){ // There is an intersection but not coplanar
				intersections.push_back(std::pair<double, Element::GeometryType*>(intersection[direction], &((*i_object)->GetGeometry())));
			}
		}

		return 0;
	}

	//TODO: This method has been adapted from the previous implementation. It is still pending to update it.
	template<std::size_t TDim>
	int CalculateDistanceToSkinProcess<TDim>::IntersectionTriangleSegment(
		Element::GeometryType& rGeometry,
		const double* RayPoint1,
		const double* RayPoint2,
		double* IntersectionPoint)
	{
		const double epsilon = 1.00e-12;

        array_1d<double,3>    u, v, n;             // triangle vectors
        array_1d<double,3>    dir, w0, w;          // ray vectors
        double     r, a, b;             // params to calc ray-plane intersect


        // get triangle edge vectors and plane normal
        u = rGeometry[1] - rGeometry[0];
        v = rGeometry[2] - rGeometry[0];

        MathUtils<double>::CrossProduct(n, u, v);             // cross product

        if (norm_2(n) == 0)            // triangle is degenerate
            return -1;                 // do not deal with this case

		double triangle_origin_distance = -inner_prod(n, rGeometry[0]);
		Point ray_point_1, ray_point_2;

		for(int i = 0 ; i < 3 ; i++)
        {
            dir[i] = RayPoint2[i] - RayPoint1[i];             // ray direction vector
            w0[i] = RayPoint1[i] - rGeometry[0][i];
			ray_point_1[i] = RayPoint1[i];
			ray_point_2[i] = RayPoint2[i];
		}

		double sign_distance_1 = inner_prod(n, ray_point_1) + triangle_origin_distance;
		double sign_distance_2 = inner_prod(n, ray_point_2) + triangle_origin_distance;

		if (sign_distance_1*sign_distance_2 > epsilon) // segment line point on the same side of plane
			return 0;
		a = -inner_prod(n,w0);
        b = inner_prod(n,dir);

        if (fabs(b) < epsilon) // ray is parallel to triangle plane
		{
            if (a == 0)                // ray lies in triangle plane
                return 2;
            else return 0;             // ray disjoint from plane
        }

        // get intersect point of ray with triangle plane
        r = a / b;
        if (r < 0.0)                   // ray goes away from triangle
            return 0;                  // => no intersect
        // for a segment, also test if (r > 1.0) => no intersect

        for(int i = 0 ; i < 3 ; i++)
            IntersectionPoint[i]  = RayPoint1[i] + r * dir[i];           // intersect point of ray and plane

        // is I inside T?
        double    uu, uv, vv, wu, wv, D;
        uu = inner_prod(u,u);
        uv = inner_prod(u,v);
        vv = inner_prod(v,v);

        for(int i = 0 ; i < 3 ; i++)
            w[i] = IntersectionPoint[i] - rGeometry[0][i];

        wu = inner_prod(w,u);
        wv = inner_prod(w,v);
        D = uv * uv - uu * vv;

        // get and test parametric coords
        double s, t;
        s = (uv * wv - vv * wu) / D;
        if (s < 0.0 - epsilon || s > 1.0 + epsilon)        // I is outside T
            return 0;
        t = (uv * wu - uu * wv) / D;
        if (t < 0.0 - epsilon || (s + t) > 1.0 + epsilon)  // I is outside T
            return 0;

        return 1;                      // I is in T

	}

	template<>
	int CalculateDistanceToSkinProcess<2>::ComputeRayIntersection(
		Element::GeometryType& rGeometry,
        const double* pRayPoint1,
        const double* pRayPoint2,
        double* pIntersectionPoint)
	{
		// Auxiliar arrays 
		array_1d<double,3> int_pt;
		array_1d<double,3> ray_pt_1;
		array_1d<double,3> ray_pt_2;
		for (unsigned int i = 0; i < 3; ++i){
			ray_pt_1[i] = pRayPoint1[i];
			ray_pt_2[i] = pRayPoint2[i];
		}

		// Call the line - line intersection util
		const double tolerance = 1.0e-6*rGeometry.Length();
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
	int CalculateDistanceToSkinProcess<3>::ComputeRayIntersection(
		Element::GeometryType& rGeometry,
        const double* pRayPoint1,
        const double* pRayPoint2,
        double* pIntersectionPoint)
	{
		return IntersectionTriangleSegment(rGeometry, pRayPoint1, pRayPoint2, pIntersectionPoint); 
	}

	template<>
	double inline CalculateDistanceToSkinProcess<2>::CalculatePointDistance(
		const Element::GeometryType &rIntObjGeom,
		const Point &rDistancePoint)
	{
		return GeometryUtils::PointDistanceToLineSegment3D(
			rIntObjGeom[0],
			rIntObjGeom[1],
			rDistancePoint);
	}

	template<>
	double inline CalculateDistanceToSkinProcess<3>::CalculatePointDistance(
		const Element::GeometryType &rIntObjGeom,
		const Point &rDistancePoint)
	{
		return GeometryUtils::PointDistanceToTriangle3D(
			rIntObjGeom[0],
			rIntObjGeom[1],
			rIntObjGeom[2],
			rDistancePoint);
	}

	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::Execute()
	{
		this->Initialize();
		this->FindIntersections();
		this->CalculateDistances(this->GetIntersections());
	}

	/// Turn back information as a string.
	template<std::size_t TDim>
	std::string CalculateDistanceToSkinProcess<TDim>::Info() const
	{
		return "CalculateDistanceToSkinProcess";
	}

	/// Print information about this object.
	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::PrintInfo(std::ostream& rOStream) const
	{
		rOStream << Info();
	}

	/// Print object's data.
	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::PrintData(std::ostream& rOStream) const
	{
	}

	template class Kratos::CalculateDistanceToSkinProcess<2>;
	template class Kratos::CalculateDistanceToSkinProcess<3>;

}  // namespace Kratos.
