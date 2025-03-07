#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/eigen/dense.h>
#include <fcpw/fcpw.h>
#ifdef FCPW_USE_GPU
    #include <nanobind/stl/string.h>
    #include <fcpw/fcpw_gpu.h>
#endif

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(py, m) {
    m.doc() = "FCPW Python bindings";

    nb::enum_<fcpw::PrimitiveType>(m, "primitive_type")
        .value("line_segment", fcpw::PrimitiveType::LineSegment)
        .value("triangle", fcpw::PrimitiveType::Triangle);

    nb::enum_<fcpw::AggregateType>(m, "aggregate_type")
        .value("baseline", fcpw::AggregateType::Baseline)
        .value("bvh_longest_axis_center", fcpw::AggregateType::Bvh_LongestAxisCenter)
        .value("bvh_overlap_surface_area", fcpw::AggregateType::Bvh_OverlapSurfaceArea)
        .value("bvh_surface_area", fcpw::AggregateType::Bvh_SurfaceArea)
        .value("bvh_overlap_volume", fcpw::AggregateType::Bvh_OverlapVolume)
        .value("bvh_volume", fcpw::AggregateType::Bvh_Volume);

    nb::enum_<fcpw::DistanceInfo>(m, "distance_info")
        .value("exact", fcpw::DistanceInfo::Exact)
        .value("bounded", fcpw::DistanceInfo::Bounded);

    nb::enum_<fcpw::BooleanOperation>(m, "boolean_operation")
        .value("union", fcpw::BooleanOperation::Union)
        .value("intersection", fcpw::BooleanOperation::Intersection)
        .value("difference", fcpw::BooleanOperation::Difference)
        .value("none", fcpw::BooleanOperation::None);

    nb::class_<fcpw::CsgTreeNode>(m, "csg_tree_node")
        .def(nb::init<int, int, bool, bool, fcpw::BooleanOperation>(),
            "child1"_a, "child2"_a, "is_leaf_child1"_a, "is_leaf_child2"_a, "operation"_a)
        .def_rw("child1", &fcpw::CsgTreeNode::child1)
        .def_rw("child2", &fcpw::CsgTreeNode::child2)
        .def_rw("is_leaf_child1", &fcpw::CsgTreeNode::isLeafChild1)
        .def_rw("is_leaf_child2", &fcpw::CsgTreeNode::isLeafChild2)
        .def_rw("operation", &fcpw::CsgTreeNode::operation);

    nb::class_<fcpw::Ray<2>>(m, "ray_2D")
        .def(nb::init<const fcpw::Vector<2>&, const fcpw::Vector<2>&, float>(),
            "o"_a, "d"_a, "t_max"_a=fcpw::maxFloat)
        .def_rw("o", &fcpw::Ray<2>::o)
        .def_rw("d", &fcpw::Ray<2>::d)
        .def_rw("inv_d", &fcpw::Ray<2>::invD)
        .def_rw("t_max", &fcpw::Ray<2>::tMax);

    nb::class_<fcpw::BoundingSphere<2>>(m, "bounding_sphere_2D")
        .def(nb::init<const fcpw::Vector<2>&, float>(),
            "c"_a, "r2"_a)
        .def_rw("c", &fcpw::BoundingSphere<2>::c)
        .def_rw("r2", &fcpw::BoundingSphere<2>::r2);

    nb::class_<fcpw::Interaction<2>>(m, "interaction_2D")
        .def(nb::init<>())
        .def("signed_distance", &fcpw::Interaction<2>::signedDistance, "x"_a)
        .def_rw("d", &fcpw::Interaction<2>::d)
        .def_rw("sign", &fcpw::Interaction<2>::sign)
        .def_rw("primitive_index", &fcpw::Interaction<2>::primitiveIndex)
        .def_rw("node_index", &fcpw::Interaction<2>::nodeIndex)
        .def_rw("reference_index", &fcpw::Interaction<2>::referenceIndex)
        .def_rw("object_index", &fcpw::Interaction<2>::objectIndex)
        .def_rw("p", &fcpw::Interaction<2>::p)
        .def_rw("n", &fcpw::Interaction<2>::n)
        .def_rw("uv", &fcpw::Interaction<2>::uv)
        .def_rw("distance_info", &fcpw::Interaction<2>::distanceInfo);

    nb::class_<fcpw::Ray<3>>(m, "ray_3D")
        .def(nb::init<const fcpw::Vector<3>&, const fcpw::Vector<3>&, float>(),
            "o"_a, "d"_a, "t_max"_a=fcpw::maxFloat)
        .def_rw("o", &fcpw::Ray<3>::o)
        .def_rw("d", &fcpw::Ray<3>::d)
        .def_rw("inv_d", &fcpw::Ray<3>::invD)
        .def_rw("t_max", &fcpw::Ray<3>::tMax);

    nb::class_<fcpw::BoundingSphere<3>>(m, "bounding_sphere_3D")
        .def(nb::init<const fcpw::Vector<3>&, float>(),
            "c"_a, "r2"_a)
        .def_rw("c", &fcpw::BoundingSphere<3>::c)
        .def_rw("r2", &fcpw::BoundingSphere<3>::r2);

    nb::class_<fcpw::Interaction<3>>(m, "interaction_3D")
        .def(nb::init<>())
        .def("signed_distance", &fcpw::Interaction<3>::signedDistance, "x"_a)
        .def_rw("d", &fcpw::Interaction<3>::d)
        .def_rw("sign", &fcpw::Interaction<3>::sign)
        .def_rw("primitive_index", &fcpw::Interaction<3>::primitiveIndex)
        .def_rw("node_index", &fcpw::Interaction<3>::nodeIndex)
        .def_rw("reference_index", &fcpw::Interaction<3>::referenceIndex)
        .def_rw("object_index", &fcpw::Interaction<3>::objectIndex)
        .def_rw("p", &fcpw::Interaction<3>::p)
        .def_rw("n", &fcpw::Interaction<3>::n)
        .def_rw("uv", &fcpw::Interaction<3>::uv)
        .def_rw("distance_info", &fcpw::Interaction<3>::distanceInfo);

    using UInt32List = std::vector<uint32_t>;
    nb::bind_vector<UInt32List>(m, "uint32_list");

    using Vector2fList = std::vector<fcpw::Vector<2>>;
    nb::bind_vector<Vector2fList>(m, "vector_2f_list");

    using Vector3fList = std::vector<fcpw::Vector<3>>;
    nb::bind_vector<Vector3fList>(m, "vector_3f_list");

    using Vector2iList = std::vector<fcpw::Vector2i>;
    nb::bind_vector<Vector2iList>(m, "vector_2i_list");

    using Vector3iList = std::vector<fcpw::Vector3i>;
    nb::bind_vector<Vector3iList>(m, "vector_3i_list");

    using Interaction2DList = std::vector<fcpw::Interaction<2>>;
    nb::bind_vector<Interaction2DList>(m, "interaction_2D_list");

    using Interaction3DList = std::vector<fcpw::Interaction<3>>;
    nb::bind_vector<Interaction3DList>(m, "interaction_3D_list");

    using PrimitiveTypeList = std::vector<fcpw::PrimitiveType>;
    nb::bind_vector<PrimitiveTypeList>(m, "primitive_type_list");

    using ObjectTypeList = std::vector<PrimitiveTypeList>;
    nb::bind_vector<ObjectTypeList>(m, "object_type_list");

    using Transform2DList = std::vector<Eigen::Matrix3f>;
    nb::bind_vector<Transform2DList>(m, "transform_2D_list");

    using Transform3DList = std::vector<Eigen::Matrix4f>;
    nb::bind_vector<Transform3DList>(m, "transform_3D_list");

    using Ray2DList = std::vector<fcpw::Ray<2>>;
    nb::bind_vector<Ray2DList>(m, "ray_2D_list");

    using Ray3DList = std::vector<fcpw::Ray<3>>;
    nb::bind_vector<Ray3DList>(m, "ray_3D_list");

    using BoundingSphere2DList = std::vector<fcpw::BoundingSphere<2>>;
    nb::bind_vector<BoundingSphere2DList>(m, "bounding_sphere_2D_list");

    using BoundingSphere3DList = std::vector<fcpw::BoundingSphere<3>>;
    nb::bind_vector<BoundingSphere3DList>(m, "bounding_sphere_3D_list");

    nb::class_<fcpw::Scene<2>>(m, "scene_2D")
        .def(nb::init<>())
        .def("set_object_count", &fcpw::Scene<2>::setObjectCount,
            "Sets the number of objects in the scene. Each call to this function resets the scene data.",
            "n_objects"_a)
        .def("set_object_vertices", &fcpw::Scene<2>::setObjectVertices,
            "Sets the vertex positions for an object.",
            "positions"_a, "object_index"_a)
        .def("set_object_line_segments", &fcpw::Scene<2>::setObjectLineSegments,
            "Sets the vertex indices of line segments for an object.",
            "indices"_a, "object_index"_a)
        .def("Set_object_types", &fcpw::Scene<2>::setObjectTypes,
            "Sets the primitive type for each object in the scene.\nE.g., [[line_segment], [line_segment], [line_segment]] specifies a scene with 3 objects made of line segments in 2D.\nObjects with multiple primitive types, or primitives embedded in different dimensions (e.g., [[triangle, line_segment], [triangle]]) are not currently supported.\nEach call to this function resets the scene data.",
            "object_types"_a)
        .def("set_object_vertex_count", &fcpw::Scene<2>::setObjectVertexCount,
            "Sets the number of vertices for an object.",
            "n_vertices"_a, "object_index"_a)
        .def("set_object_line_segment_count", &fcpw::Scene<2>::setObjectLineSegmentCount,
            "Sets the number of line segments for an object.",
            "n_line_segments"_a, "object_index"_a)
        .def("set_object_vertex", &fcpw::Scene<2>::setObjectVertex,
            "Sets the position of a vertex for an object.",
            "position"_a, "vertex_index"_a, "object_index"_a)
        .def("set_object_line_segment", &fcpw::Scene<2>::setObjectLineSegment,
            "Sets the vertex indices of a line segment for an object.",
            "indices"_a, "line_segment_index"_a, "object_index"_a)
        .def("set_object_instance_transforms",
            [](fcpw::Scene<2>& self, const Transform2DList& matrixTransforms, int objectIndex) {
                std::vector<fcpw::Transform<2>> transforms;
                for (int i = 0; i < (int)matrixTransforms.size(); i++) {
                    fcpw::Transform<2> t(matrixTransforms[i]);
                    transforms.emplace_back(t);
                }
                self.setObjectInstanceTransforms(transforms, objectIndex);
            },
            "Sets the instance transforms for an object.",
            "transforms"_a, "object_index"_a)
        .def("set_csg_tree_node", &fcpw::Scene<2>::setCsgTreeNode,
            "Sets the data for a node in the csg tree.\nNOTE: the root node of the csg tree must have index 0.",
            "csg_tree_node"_a, "node_index"_a)
        .def("compute_silhouettes", &fcpw::Scene<2>::computeSilhouettes,
            "Precomputes silhouette information for primitives in a scene to perform closest silhouette point queries.\nThe optional ignore_silhouette callback allows the user to specify which interior vertices in the line segment geometry\nto ignore for silhouette tests (arguments: vertex dihedral angle, index of an adjacent line segment).\nNOTE: does not currently support non-manifold geometry.",
            "ignore_silhouette"_a.none())
        .def("compute_object_normals", &fcpw::Scene<2>::computeObjectNormals,
            "Precomputes vertex normals for an object made of line segments.",
            "object_index"_a, "compute_weighted"_a=false)
        .def("build", &fcpw::Scene<2>::build,
            "Builds a (possibly vectorized) aggregate/accelerator for the scene.\nEach call to this function rebuilds the aggregate/accelerator for the scene from the specified geometry\n(except when reduce_memory_footprint is set to true which results in undefined behavior).\nIt is recommended to set vectorize to false for primitives that do not implement vectorized intersection and closest point queries.\nSet reduce_memory_footprint to true to reduce the memory footprint of fcpw when constructing an aggregate,\nhowever if you plan to access the scene data let it remain false.",
            "aggregate_type"_a, "vectorize"_a, "print_stats"_a=false, "reduce_memory_footprint"_a=false)
        .def("update_object_vertex", &fcpw::Scene<2>::updateObjectVertex,
            "Updates the position of a vertex for an object.",
            "position"_a, "vertex_index"_a, "object_index"_a)
        .def("update_object_vertices", &fcpw::Scene<2>::updateObjectVertices,
            "Updates the vertex positions for an object.",
            "positions"_a, "object_index"_a)
        .def("refit", &fcpw::Scene<2>::refit,
            "Refits the scene aggregate hierarchy after updating the geometry, via calls to update_object_vertex and compute_object_normals.\nNOTE: refitting of instanced aggregates is currently quite inefficient, since the shared aggregate is refit for each instance.",
            "print_stats"_a=false)
        .def("intersect", nb::overload_cast<fcpw::Ray<2>&, fcpw::Interaction<2>&, bool>(
            &fcpw::Scene<2>::intersect, nb::const_),
            "Intersects the scene with the given ray and returns whether there is a hit.\nIf check_for_occlusion is enabled, the interaction is not populated.",
            "r"_a, "i"_a, "check_for_occlusion"_a=false)
        .def("intersect", nb::overload_cast<fcpw::Ray<2>&, Interaction2DList&, bool, bool>(
            &fcpw::Scene<2>::intersect, nb::const_),
            "Intersects the scene with the given ray and returns the number of hits.\nBy default, returns the closest interaction if it exists.\nIf check_for_occlusion is enabled, the interactions vector is not populated.\nIf record_all_hits is enabled, sorts interactions by distance to the ray origin.",
            "r"_a, "is"_a, "check_for_occlusion"_a=false, "record_all_hits"_a=false)
        .def("intersect", nb::overload_cast<const fcpw::BoundingSphere<2>&, Interaction2DList&, bool>(
            &fcpw::Scene<2>::intersect, nb::const_),
            "Intersects the scene with the given sphere and returns the number of primitives inside the sphere: interactions contain the primitive indices.\nIf record_one_hit is set to true, randomly selects one geometric primitive inside the sphere (one for each aggregate in the hierarchy)\nand writes the selection pdf value to interaction_2D.d along with the primitive index.",
            "s"_a, "is"_a, "record_one_hit"_a=false)
        .def("intersect", nb::overload_cast<const fcpw::BoundingSphere<2>&, fcpw::Interaction<2>&, const fcpw::Vector<2>&, const std::function<float(float)>&>(
            &fcpw::Scene<2>::intersect, nb::const_),
            "Intersects the scene with the given sphere.\nThis method does not visit all primitives inside the sphere during traversal--the primitives visited are chosen stochastically.\nIt randomly selects one geometric primitive inside the sphere using the user specified weight function (function argument is the squared distance\nbetween the sphere and box/primitive centers) and samples a random point on that primitive (written to interaction_2D.p) using the random numbers rand_nums[2].\nThe selection pdf value is written to interaction_2D.d along with the primitive index.",
            "s"_a, "i"_a, "rand_nums"_a, "branch_traversal_weight"_a.none())
        .def("contains", nb::overload_cast<const fcpw::Vector<2>&>(
            &fcpw::Scene<2>::contains, nb::const_),
            "Checks whether a point is contained inside a scene.\nNOTE: the scene must be watertight.",
            "x"_a)
        .def("has_line_of_sight", nb::overload_cast<const fcpw::Vector<2>&, const fcpw::Vector<2>&>(
            &fcpw::Scene<2>::hasLineOfSight, nb::const_),
            "Checks whether there is a line of sight between between two points in the scene.",
            "xi"_a, "xj"_a)
        .def("find_closest_point", &fcpw::Scene<2>::findClosestPoint,
            "Finds the closest point in the scene to a query point.\nOptionally specify a conservative radius guess around the query point inside which the search is performed.",
            "x"_a, "i"_a, "squared_radius"_a=fcpw::maxFloat, "record_normal"_a=false)
        .def("find_closest_silhouette_point", &fcpw::Scene<2>::findClosestSilhouettePoint,
            "Finds the closest point on the visibility silhouette in the scene to a query point.\nOptionally specify a minimum radius to stop the closest silhouette search, a conservative maximum radius guess\naround the query point inside which the search is performed, as well as a precision parameter to help classify\nsilhouettes when the query point lies on the scene geometry.",
            "x"_a, "i"_a, "flip_normal_orientation"_a=false, "squared_min_radius"_a=0.0f,
            "squared_max_radius"_a=fcpw::maxFloat, "precision"_a=1e-3f, "record_normal"_a=false)
        .def("intersect", nb::overload_cast<Ray2DList&, Interaction2DList&, bool>(
            &fcpw::Scene<2>::intersect, nb::const_),
            "Intersects the scene with the given rays, returning the closest interaction if it exists.",
            "rays"_a, "interactions"_a, "check_for_occlusion"_a=false)
        .def("intersect", nb::overload_cast<const BoundingSphere2DList&, Interaction2DList&, const Vector2fList&, const std::function<float(float)>&>(
            &fcpw::Scene<2>::intersect, nb::const_),
            "Intersects the scene with the given spheres, randomly selecting one geometric primitive contained inside each sphere and sampling\na random point on that primitive (written to interaction_2D.p) using the random numbers randNums[2].\nThe selection pdf value is written to interaction_2D.d along with the primitive index.",
            "bounding_spheres"_a, "interactions"_a, "rand_nums"_a, "branch_traversal_weight"_a.none())
        .def("contains", nb::overload_cast<const Vector2fList&, UInt32List&>(
            &fcpw::Scene<2>::contains, nb::const_),
            "Checks whether points are contained inside a scene. NOTE: the scene must be watertight.",
            "points"_a, "result"_a)
        .def("has_line_of_sight", nb::overload_cast<const Vector2fList&, const Vector2fList&, UInt32List&>(
            &fcpw::Scene<2>::hasLineOfSight, nb::const_),
            "Checks whether there is a line of sight between between two sets of points in the scene.",
            "points_i"_a, "points_j"_a, "result"_a)
        .def("find_closest_points", &fcpw::Scene<2>::findClosestPoints,
            "Finds the closest points in the scene to the given query points, encoded as bounding spheres.\nThe radius of each bounding sphere specifies the conservative radius guess around the query point inside which the search is performed.",
            "bounding_spheres"_a, "interactions"_a, "record_normal"_a=false)
        .def("find_closest_silhouette_points", &fcpw::Scene<2>::findClosestSilhouettePoints,
            "Finds the closest points on the visibility silhouette in the scene to the given query points, encoded as bounding spheres.\nOptionally specify a minimum radius to stop the closest silhouette search, as well as a precision parameter to help classify silhouettes.",
            "bounding_spheres"_a, "interactions"_a, "flip_normal_orientation"_a, "squared_min_radius"_a=0.0f, "precision"_a=1e-3f, "record_normal"_a=false);

    nb::class_<fcpw::Scene<3>>(m, "scene_3D")
        .def(nb::init<>())
        .def("set_object_count", &fcpw::Scene<3>::setObjectCount,
            "Sets the number of objects in the scene. Each call to this function resets the scene data.",
            "n_objects"_a)
        .def("set_object_vertices", &fcpw::Scene<3>::setObjectVertices,
            "Sets the vertex positions for an object.",
            "positions"_a, "object_index"_a)
        .def("set_object_triangles", &fcpw::Scene<3>::setObjectTriangles,
            "Sets the vertex indices of triangles for an object.",
            "indices"_a, "object_index"_a)
        .def("set_object_types", &fcpw::Scene<3>::setObjectTypes,
            "Sets the primitive type for each object in the scene.\nE.g., [[triangle], [triangle]] specifies a scene with 2 objects made of triangles in 3D.\nObjects with multiple primitive types, or primitives embedded in different dimensions (e.g., [[triangle, line_segment], [triangle]]) are not currently supported.\nEach call to this function resets the scene data.",
            "object_types"_a)
        .def("set_object_vertex_count", &fcpw::Scene<3>::setObjectVertexCount,
            "Sets the number of vertices for an object.",
            "n_vertices"_a, "object_index"_a)
        .def("set_object_triangle_count", &fcpw::Scene<3>::setObjectTriangleCount,
            "Sets the number of triangles for an object.",
            "n_triangles"_a, "object_index"_a)
        .def("set_object_vertex", &fcpw::Scene<3>::setObjectVertex,
            "Sets the position of a vertex for an object.",
            "position"_a, "vertex_index"_a, "object_index"_a)
        .def("set_object_triangle", &fcpw::Scene<3>::setObjectTriangle,
            "Sets the vertex indices of a triangle for an object.",
            "indices"_a, "triangle_index"_a, "object_index"_a)
        .def("set_object_instance_transforms",
            [](fcpw::Scene<3>& self, const Transform3DList& matrixTransforms, int objectIndex) {
                std::vector<fcpw::Transform<3>> transforms;
                for (int i = 0; i < (int)matrixTransforms.size(); i++) {
                    fcpw::Transform<3> t(matrixTransforms[i]);
                    transforms.emplace_back(t);
                }
                self.setObjectInstanceTransforms(transforms, objectIndex);
            },
            "Sets the instance transforms for an object.",
            "transforms"_a, "object_index"_a)
        .def("set_csg_tree_node", &fcpw::Scene<3>::setCsgTreeNode,
            "Sets the data for a node in the csg tree.\nNOTE: the root node of the csg tree must have index 0.",
            "csg_tree_node"_a, "node_index"_a)
        .def("compute_silhouettes", &fcpw::Scene<3>::computeSilhouettes,
            "Precomputes silhouette information for primitives in a scene to perform closest silhouette point queries.\nThe optional ignore_silhouette callback allows the user to specify which interior edges in the triangle geometry\nto ignore for silhouette tests (arguments: edge dihedral angle, index of an adjacent triangle).\nNOTE: does not currently support non-manifold geometry.",
            "ignore_silhouette"_a.none())
        .def("compute_object_normals", &fcpw::Scene<3>::computeObjectNormals,
            "Precomputes edge normals for an object made of triangles.",
            "object_index"_a, "compute_weighted"_a=false)
        .def("build", &fcpw::Scene<3>::build,
            "Builds a (possibly vectorized) aggregate/accelerator for the scene.\nEach call to this function rebuilds the aggregate/accelerator for the scene from the specified geometry\n(except when reduce_memory_footprint is set to true which results in undefined behavior).\nIt is recommended to set vectorize to false for primitives that do not implement vectorized intersection and closest point queries.\nSet reduce_memory_footprint to true to reduce the memory footprint of fcpw when constructing an aggregate,\nhowever if you plan to access the scene data let it remain false.",
            "aggregate_type"_a, "vectorize"_a, "print_stats"_a=false, "reduce_memory_footprint"_a=false)
        .def("update_object_vertex", &fcpw::Scene<3>::updateObjectVertex,
            "Updates the position of a vertex for an object.",
            "position"_a, "vertex_index"_a, "object_index"_a)
        .def("update_object_vertices", &fcpw::Scene<3>::updateObjectVertices,
            "Updates the vertex positions for an object.",
            "positions"_a, "object_index"_a)
        .def("refit", &fcpw::Scene<3>::refit,
            "Refits the scene aggregate hierarchy after updating the geometry, via calls to update_object_vertex and compute_object_normals.\nNOTE: refitting of instanced aggregates is currently quite inefficient, since the shared aggregate is refit for each instance.",
            "print_stats"_a=false)
        .def("intersect", nb::overload_cast<fcpw::Ray<3>&, fcpw::Interaction<3>&, bool>(
            &fcpw::Scene<3>::intersect, nb::const_),
            "Intersects the scene with the given ray and returns whether there is a hit.\nIf check_for_occlusion is enabled, the interaction is not populated.",
            "r"_a, "i"_a, "check_for_occlusion"_a=false)
        .def("intersect", nb::overload_cast<fcpw::Ray<3>&, Interaction3DList&, bool, bool>(
            &fcpw::Scene<3>::intersect, nb::const_),
            "Intersects the scene with the given ray and returns the number of hits.\nBy default, returns the closest interaction if it exists.\nIf check_for_occlusion is enabled, the interactions vector is not populated.\nIf record_all_hits is enabled, sorts interactions by distance to the ray origin.",
            "r"_a, "is"_a, "check_for_occlusion"_a=false, "record_all_hits"_a=false)
        .def("intersect", nb::overload_cast<const fcpw::BoundingSphere<3>&, Interaction3DList&, bool>(
            &fcpw::Scene<3>::intersect, nb::const_),
            "Intersects the scene with the given sphere and returns the number of primitives inside the sphere: interactions contain the primitive indices.\nIf record_one_hit is set to true, randomly selects one geometric primitive inside the sphere (one for each aggregate in the hierarchy)\nand writes the selection pdf value to interaction_3D.d along with the primitive index.",
            "s"_a, "is"_a, "record_one_hit"_a=false)
        .def("intersect", nb::overload_cast<const fcpw::BoundingSphere<3>&, fcpw::Interaction<3>&, const fcpw::Vector<3>&, const std::function<float(float)>&>(
            &fcpw::Scene<3>::intersect, nb::const_),
            "Intersects the scene with the given sphere.\nThis method does not visit all primitives inside the sphere during traversal--the primitives visited are chosen stochastically.\nIt randomly selects one geometric primitive inside the sphere using the user specified weight function (function argument is the squared distance\nbetween the sphere and box/primitive centers) and samples a random point on that primitive (written to interaction_3D.p) using the random numbers rand_nums[3].\nThe selection pdf value is written to interaction_3D.d along with the primitive index.",
            "s"_a, "i"_a, "rand_nums"_a, "branch_traversal_weight"_a.none())
        .def("contains", nb::overload_cast<const fcpw::Vector<3>&>(
            &fcpw::Scene<3>::contains, nb::const_),
            "Checks whether a point is contained inside a scene.\nNOTE: the scene must be watertight.",
            "x"_a)
        .def("has_line_of_sight", nb::overload_cast<const fcpw::Vector<3>&, const fcpw::Vector<3>&>(
            &fcpw::Scene<3>::hasLineOfSight, nb::const_),
            "Checks whether there is a line of sight between between two points in the scene.",
            "xi"_a, "xj"_a)
        .def("find_closest_point", &fcpw::Scene<3>::findClosestPoint,
            "Finds the closest point in the scene to a query point.\nOptionally specify a conservative radius guess around the query point inside which the search is performed.",
            "x"_a, "i"_a, "squared_radius"_a=fcpw::maxFloat, "record_normal"_a=false)
        .def("find_closest_silhouette_point", &fcpw::Scene<3>::findClosestSilhouettePoint,
            "Finds the closest point on the visibility silhouette in the scene to a query point.\nOptionally specify a minimum radius to stop the closest silhouette search, a conservative maximum radius guess\naround the query point inside which the search is performed, as well as a precision parameter to help classify\nsilhouettes when the query point lies on the scene geometry.",
            "x"_a, "i"_a, "flipNormalOrientation"_a=false, "squaredMinRadius"_a=0.0f,
            "squared_max_radius"_a=fcpw::maxFloat, "precision"_a=1e-3f, "record_normal"_a=false)
        .def("intersect", nb::overload_cast<Ray3DList&, Interaction3DList&, bool>(
            &fcpw::Scene<3>::intersect, nb::const_),
            "Intersects the scene with the given rays, returning the closest interaction if it exists.",
            "rays"_a, "interactions"_a, "check_for_occlusion"_a=false)
        .def("intersect", nb::overload_cast<const BoundingSphere3DList&, Interaction3DList&, const Vector3fList&, const std::function<float(float)>&>(
            &fcpw::Scene<3>::intersect, nb::const_),
            "Intersects the scene with the given spheres, randomly selecting one geometric primitive contained inside each sphere and sampling\na random point on that primitive (written to interaction_3D.p) using the random numbers randNums[3].\nThe selection pdf value is written to interaction_3D.d along with the primitive index.",
            "bounding_spheres"_a, "interactions"_a, "rand_nums"_a, "branch_traversal_weight"_a.none())
        .def("contains", nb::overload_cast<const Vector3fList&, UInt32List&>(
            &fcpw::Scene<3>::contains, nb::const_),
            "Checks whether points are contained inside a scene. NOTE: the scene must be watertight.",
            "points"_a, "result"_a)
        .def("has_line_of_sight", nb::overload_cast<const Vector3fList&, const Vector3fList&, UInt32List&>(
            &fcpw::Scene<3>::hasLineOfSight, nb::const_),
            "Checks whether there is a line of sight between between two sets of points in the scene.",
            "points_i"_a, "points_j"_a, "result"_a)
        .def("find_closest_points", &fcpw::Scene<3>::findClosestPoints,
            "Finds the closest points in the scene to the given query points, encoded as bounding spheres.\nThe radius of each bounding sphere specifies the conservative radius guess around the query point inside which the search is performed.",
            "bounding_spheres"_a, "interactions"_a, "record_normal"_a=false)
        .def("find_closest_silhouette_points", &fcpw::Scene<3>::findClosestSilhouettePoints,
            "Finds the closest points on the visibility silhouette in the scene to the given query points, encoded as bounding spheres.\nOptionally specify a minimum radius to stop the closest silhouette search, as well as a precision parameter to help classify silhouettes.",
            "bounding_spheres"_a, "interactions"_a, "flip_normal_orientation"_a, "squared_min_radius"_a=0.0f, "precision"_a=1e-3f, "record_normal"_a=false);

#ifdef FCPW_USE_GPU
    nb::class_<fcpw::float2>(m, "float_2D")
        .def(nb::init<>())
        .def(nb::init<float, float>(),
            "x"_a, "y"_a)
        .def_rw("x", &fcpw::float2::x)
        .def_rw("y", &fcpw::float2::y);

    nb::class_<fcpw::float3>(m, "float_3D")
        .def(nb::init<>())
        .def(nb::init<float, float, float>(),
            "x"_a, "y"_a, "z"_a)
        .def_rw("x", &fcpw::float3::x)
        .def_rw("y", &fcpw::float3::y)
        .def_rw("z", &fcpw::float3::z);

    nb::class_<fcpw::GPURay>(m, "gpu_ray")
        .def(nb::init<>())
        .def(nb::init<const fcpw::float3&, const fcpw::float3&, float>(),
            "o"_a, "d"_a, "tMax"_a=fcpw::maxFloat)
        .def_rw("o", &fcpw::GPURay::o)
        .def_rw("d", &fcpw::GPURay::d)
        .def_rw("inv_d", &fcpw::GPURay::dInv)
        .def_rw("t_max", &fcpw::GPURay::tMax);

    nb::class_<fcpw::GPUBoundingSphere>(m, "gpu_bounding_sphere")
        .def(nb::init<>())
        .def(nb::init<const fcpw::float3&, float>(),
            "c"_a, "r2"_a)
        .def_rw("c", &fcpw::GPUBoundingSphere::c)
        .def_rw("r2", &fcpw::GPUBoundingSphere::r2);

    nb::class_<fcpw::GPUInteraction>(m, "gpu_interaction")
        .def(nb::init<>())
        .def_rw("p", &fcpw::GPUInteraction::p)
        .def_rw("n", &fcpw::GPUInteraction::n)
        .def_rw("uv", &fcpw::GPUInteraction::uv)
        .def_rw("d", &fcpw::GPUInteraction::d)
        .def_rw("index", &fcpw::GPUInteraction::index);

    using Float3DList = std::vector<fcpw::float3>;
    nb::bind_vector<Float3DList>(m, "float_3D_list");

    using GPURayList = std::vector<fcpw::GPURay>;
    nb::bind_vector<GPURayList>(m, "gpu_ray_list");

    using GPUBoundingSphereList = std::vector<fcpw::GPUBoundingSphere>;
    nb::bind_vector<GPUBoundingSphereList>(m, "gpu_bounding_sphere_list");

    using GPUInteractionList = std::vector<fcpw::GPUInteraction>;
    nb::bind_vector<GPUInteractionList>(m, "gpu_interaction_list");

    nb::class_<fcpw::GPUScene<2>>(m, "gpu_scene_2D")
        .def(nb::init<const std::string&, bool>(),
            "fcpw_directory_path"_a, "print_logs"_a=false)
        .def("transfer_to_gpu", &fcpw::GPUScene<2>::transferToGPU,
            "Transfers a binary (non-vectorized) BVH aggregate, constructed on the CPU using the 'build' function in the Scene class, to the GPU.\nNOTE: Currently only supports scenes with a single object, i.e., no CSG trees, instanced or transformed aggregates, or nested hierarchies of aggregates.\nWhen using 'build', set 'vectorize' to false.",
            "scene"_a)
        .def("refit", &fcpw::GPUScene<2>::refit,
            "Refits the BVH on the GPU after updating the geometry, either via calls to update_object_vertex in the Scene class, or directly in GPU code\nin the user's slang shaders (set updateGeometry to false if the geometry is updated directly on the GPU).\nNOTE: Before calling this function, the BVH must already have been transferred to the GPU.",
            "scene"_a, "update_geometry"_a=true)
        .def("intersect", nb::overload_cast<GPURayList&, GPUInteractionList&, bool>(
            &fcpw::GPUScene<2>::intersect),
            "Intersects the scene with the given rays, returning the closest interaction if it exists.",
            "rays"_a, "interactions"_a, "check_for_occlusion"_a=false)
        .def("intersect", nb::overload_cast<GPUBoundingSphereList&, Float3DList&, GPUInteractionList&>(
            &fcpw::GPUScene<2>::intersect),
            "Intersects the scene with the given spheres, randomly selecting one geometric primitive contained inside each sphere and sampling\na random point on that primitive (written to interaction.p) using the random numbers rand_nums[3] (float_3D.z is ignored).\nThe selection pdf value is written to interaction.d along with the primitive index.",
            "bounding_spheres"_a, "rand_nums"_a, "interactions"_a)
        .def("find_closest_points", &fcpw::GPUScene<2>::findClosestPoints,
            "Finds the closest points in the scene to the given query points, encoded as bounding spheres.\nThe radius of each bounding sphere specifies the conservative radius guess around the query point inside which the search is performed.",
            "bounding_spheres"_a, "interactions"_a)
        .def("find_closest_silhouette_points", &fcpw::GPUScene<2>::findClosestSilhouettePoints,
            "Finds the closest points on the visibility silhouette in the scene to the given query points, encoded as bounding spheres.\nOptionally specify a minimum radius to stop the closest silhouette search, as well as a precision parameter to help classify silhouettes.",
            "bounding_spheres"_a, "flip_normal_orientation"_a, "interactions"_a,
            "squared_min_radius"_a=0.0f, "precision"_a=1e-3f);

    nb::class_<fcpw::GPUScene<3>>(m, "gpu_scene_3D")
        .def(nb::init<const std::string&, bool>(),
            "fcpw_directory_path"_a, "print_logs"_a=false)
        .def("transfer_to_gpu", &fcpw::GPUScene<3>::transferToGPU,
            "Transfers a binary (non-vectorized) BVH aggregate, constructed on the CPU using the 'build' function in the Scene class, to the GPU.\nNOTE: Currently only supports scenes with a single object, i.e., no CSG trees, instanced or transformed aggregates, or nested hierarchies of aggregates.\nWhen using 'build', set 'vectorize' to false.",
            "scene"_a)
        .def("refit", &fcpw::GPUScene<3>::refit,
            "Refits the BVH on the GPU after updating the geometry, either via calls to update_object_vertex in the Scene class, or directly in GPU code\nin the user's slang shaders (set updateGeometry to false if the geometry is updated directly on the GPU).\nNOTE: Before calling this function, the BVH must already have been transferred to the GPU.",
            "scene"_a, "update_geometry"_a=true)
        .def("intersect", nb::overload_cast<GPURayList&, GPUInteractionList&, bool>(
            &fcpw::GPUScene<3>::intersect),
            "Intersects the scene with the given rays, returning the closest interaction if it exists.",
            "rays"_a, "interactions"_a, "check_for_occlusion"_a=false)
        .def("intersect", nb::overload_cast<GPUBoundingSphereList&, Float3DList&, GPUInteractionList&>(
            &fcpw::GPUScene<3>::intersect),
            "Intersects the scene with the given spheres, randomly selecting one geometric primitive contained inside each sphere and sampling\na random point on that primitive (written to interaction.p) using the random numbers rand_nums[3].\nThe selection pdf value is written to interaction.d along with the primitive index.",
            "bounding_spheres"_a, "rand_nums"_a, "interactions"_a)
        .def("find_closest_points", &fcpw::GPUScene<3>::findClosestPoints,
            "Finds the closest points in the scene to the given query points, encoded as bounding spheres.\nThe radius of each bounding sphere specifies the conservative radius guess around the query point inside which the search is performed.",
            "bounding_spheres"_a, "interactions"_a)
        .def("find_closest_silhouette_points", &fcpw::GPUScene<3>::findClosestSilhouettePoints,
            "Finds the closest points on the visibility silhouette in the scene to the given query points, encoded as bounding spheres.\nOptionally specify a minimum radius to stop the closest silhouette search, as well as a precision parameter to help classify silhouettes.",
            "bounding_spheres"_a, "flip_normal_orientation"_a, "interactions"_a,
            "squared_min_radius"_a=0.0f, "precision"_a=1e-3f);
#endif
}