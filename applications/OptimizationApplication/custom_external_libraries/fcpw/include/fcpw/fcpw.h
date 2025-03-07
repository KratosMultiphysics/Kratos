#pragma once

#include <fcpw/utilities/scene_data.h>

namespace fcpw {

enum class PrimitiveType {
    LineSegment,
    Triangle
};

template<size_t DIM>
class Scene {
public:
    // constructor
    Scene();

    /////////////////////////////////////////////////////////////////////////////////////////////
    // API to specify scene geometry

    // sets the number of objects in the scene; each call to this function resets the scene data
    void setObjectCount(int nObjects);

    // sets the vertex positions for an object
    void setObjectVertices(const std::vector<Vector<DIM>>& positions, int objectIndex);

    // sets the vertex indices of line segments for an object
    void setObjectLineSegments(const std::vector<Vector2i>& indices, int objectIndex);

    // sets the vertex indices of triangles for an object
    void setObjectTriangles(const std::vector<Vector3i>& indices, int objectIndex);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // alternate lower-level API to specify scene geometry

    // sets the primitive type for each object in the scene;
    // e.g., {{Triangle}, {Triangle}} specifies a scene with 2 objects made of triangles in 3D;
    // {{LineSegment}, {LineSegment}, {LineSegment}} specifies a scene with 3 objects made of
    // line segments in 2D---objects with multiple primitive types, or primitives embedded in
    // different dimensions (e.g., {{Triangle, LineSegment}, {Triangle}}) are not currently supported;
    // each call to this function resets the scene data
    void setObjectTypes(const std::vector<std::vector<PrimitiveType>>& objectTypes);

    // sets the number of vertices for an object
    void setObjectVertexCount(int nVertices, int objectIndex);

    // sets the number of line segments for an object
    void setObjectLineSegmentCount(int nLineSegments, int objectIndex);

    // sets the number of triangles for an object
    void setObjectTriangleCount(int nTriangles, int objectIndex);

    // sets the position of a vertex for an object
    void setObjectVertex(const Vector<DIM>& position, int vertexIndex, int objectIndex);

    // sets the vertex indices of a line segment for an object
    void setObjectLineSegment(const Vector2i& indices, int lineSegmentIndex, int objectIndex);

    // sets the vertex indices of a triangle for an object
    void setObjectTriangle(const Vector3i& indices, int triangleIndex, int objectIndex);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // API to specify instance transforms and csg tree

    // sets the instance transforms for an object
    void setObjectInstanceTransforms(const std::vector<Transform<DIM>>& transforms, int objectIndex);

    // sets the data for a node in the csg tree; NOTE: the root node of the csg tree must have index 0
    void setCsgTreeNode(const CsgTreeNode& csgTreeNode, int nodeIndex);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // API to build and refit the scene aggregate/accelerator

    // precomputes silhouette information for primitives in a scene to perform closest silhouette
    // point queries; the optional ignoreSilhouette callback allows the user to specify which
    // interior vertices/edges in the line segment/triangle geometry to ignore for silhouette tests
    // (arguments: vertex/edge dihedral angle, index of an adjacent line segment/triangle)
    // NOTE: does not currently support non-manifold geometry
    void computeSilhouettes(const std::function<bool(float, int)>& ignoreSilhouette={});

    // precomputes vertex or edge normals for an object made of line segments or triangles (resp.)
    void computeObjectNormals(int objectIndex, bool computeWeighted=false);

    // builds a (possibly vectorized) aggregate/accelerator for the scene; each call to this
    // function rebuilds the aggregate/accelerator for the scene from the specified geometry
    // (except when reduceMemoryFootprint is set to true which results in undefined behavior);
    // it is recommended to set vectorize to false for primitives that do not implement
    // vectorized intersection and closest point queries; set reduceMemoryFootprint to true
    // to reduce the memory footprint of fcpw when constructing an aggregate, however if you
    // plan to access the scene data let it remain false.
    void build(const AggregateType& aggregateType, bool vectorize,
               bool printStats=false, bool reduceMemoryFootprint=false);

    // updates the position of a vertex for an object
    void updateObjectVertex(const Vector<DIM>& position, int vertexIndex, int objectIndex);

    // updates the vertex positions for an object
    void updateObjectVertices(const std::vector<Vector<DIM>>& positions, int objectIndex);

    // refits the scene aggregate hierarchy after updating the geometry, via calls to
    // updateObjectVertex and computeObjectNormals; NOTE: refitting of instanced aggregates
    // is currently quite inefficient, since the shared aggregate is refit for each instance
    void refit(bool printStats=false);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // API to perform ray intersection and distance queries to the scene, among others

    // intersects the scene with the given ray and returns whether there is a hit;
    // if checkForOcclusion is enabled, the interaction is not populated
    bool intersect(Ray<DIM>& r, Interaction<DIM>& i, bool checkForOcclusion=false) const;

    // intersects the scene with the given ray and returns the number of hits;
    // by default, returns the closest interaction if it exists;
    // if checkForOcclusion is enabled, the interactions vector is not populated;
    // if recordAllHits is enabled, sorts interactions by distance to the ray origin
    int intersect(Ray<DIM>& r, std::vector<Interaction<DIM>>& is,
                  bool checkForOcclusion=false, bool recordAllHits=false) const;

    // intersects the scene with the given sphere and returns the number of primitives
    // inside the sphere; interactions contain the primitive indices; if recordOneHit
    // is set to true, randomly selects one geometric primitive inside the sphere
    // (one for each aggregate in the hierarchy) and writes the selection pdf value to
    // Interaction<DIM>::d along with the primitive index
    int intersect(const BoundingSphere<DIM>& s, std::vector<Interaction<DIM>>& is,
                  bool recordOneHit=false) const;

    // intersects the scene with the given sphere; this method is faster than the one
    // above, since it does not visit all primitives inside the sphere during traversal--
    // the primitives visited are chosen stochastically; this method randomly selects one
    // geometric primitive inside the sphere using the user specified weight function
    // (function argument is the squared distance between the sphere and box/primitive centers)
    // and samples a random point on that primitive (written to Interaction<DIM>::p) using the
    // random numbers randNums[DIM]; the selection pdf value is written to Interaction<DIM>::d
    // along with the primitive index
    int intersect(const BoundingSphere<DIM>& s, Interaction<DIM>& i, const Vector<DIM>& randNums,
                  const std::function<float(float)>& branchTraversalWeight={}) const;

    // checks whether a point is contained inside a scene; NOTE: the scene must be watertight
    bool contains(const Vector<DIM>& x) const;

    // checks whether there is a line of sight between between two points in the scene
    bool hasLineOfSight(const Vector<DIM>& xi, const Vector<DIM>& xj) const;

    // finds the closest point in the scene to a query point; optionally specify a conservative
    // radius guess around the query point inside which the search is performed
    bool findClosestPoint(const Vector<DIM>& x, Interaction<DIM>& i,
                          float squaredRadius=maxFloat, bool recordNormal=false) const;

    // finds the closest point on the visibility silhouette in the scene to a query point;
    // optionally specify a minimum radius to stop the closest silhouette search,
    // a conservative maximum radius guess around the query point inside which the
    // search is performed, as well as a precision parameter to help classify silhouettes
    // when the query point lies on the scene geometry
    bool findClosestSilhouettePoint(const Vector<DIM>& x, Interaction<DIM>& i,
                                    bool flipNormalOrientation=false, float squaredMinRadius=0.0f,
                                    float squaredMaxRadius=maxFloat, float precision=1e-3f,
                                    bool recordNormal=false) const;

    /////////////////////////////////////////////////////////////////////////////////////////////
    // API to perform bundled ray intersection and distance queries to the scene
    // NOTE: the functions below use std::thread internally to parallelize queries, however
    // using Intel TBB can give up to a 2x speedup if applied directly to the query API above

    // intersects the scene with the given rays, returning the closest interaction if it exists.
    void intersect(std::vector<Ray<DIM>>& rays,
                   std::vector<Interaction<DIM>>& interactions,
                   bool checkForOcclusion=false) const;

    // intersects the scene with the given spheres, randomly selecting one geometric primitive
    // contained inside each sphere and sampling a random point on that primitive (written to 
    // Interaction<DIM>::p) using the random numbers randNums[DIM]; the selection pdf value is
    // written to Interaction<DIM>::d along with the primitive index
    void intersect(const std::vector<BoundingSphere<DIM>>& boundingSpheres,
                   std::vector<Interaction<DIM>>& interactions,
                   const std::vector<Vector<DIM>>& randNums,
                   const std::function<float(float)>& branchTraversalWeight={}) const;

    // checks whether points are contained inside a scene; NOTE: the scene must be watertight
    void contains(const std::vector<Vector<DIM>>& points,
                  std::vector<uint32_t>& result) const;

    // checks whether there is a line of sight between between two sets of points in the scene
    void hasLineOfSight(const std::vector<Vector<DIM>>& pointsI,
                        const std::vector<Vector<DIM>>& pointsJ,
                        std::vector<uint32_t>& result) const;

    // finds the closest points in the scene to the given query points, encoded as bounding spheres.
    // The radius of each bounding sphere specifies the conservative radius guess around the query
    // point inside which the search is performed.
    void findClosestPoints(std::vector<BoundingSphere<DIM>>& boundingSpheres,
                           std::vector<Interaction<DIM>>& interactions,
                           bool recordNormal=false) const;

    // finds the closest points on the visibility silhouette in the scene to the given query points,
    // encoded as bounding spheres. Optionally specify a minimum radius to stop the closest silhouette
    // search, as well as a precision parameter to help classify silhouettes.
    void findClosestSilhouettePoints(std::vector<BoundingSphere<DIM>>& boundingSpheres,
                                     std::vector<Interaction<DIM>>& interactions,
                                     const std::vector<uint32_t>& flipNormalOrientation,
                                     float squaredMinRadius=0.0f, float precision=1e-3f,
                                     bool recordNormal=false) const;

    /////////////////////////////////////////////////////////////////////////////////////////////

    // returns a pointer to the underlying scene data; use at your own risk...
    SceneData<DIM>* getSceneData();

private:
    // member
    std::unique_ptr<SceneData<DIM>> sceneData;
};

} // namespace fcpw

#include "fcpw.inl"