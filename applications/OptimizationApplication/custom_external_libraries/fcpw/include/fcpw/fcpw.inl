#include <fcpw/aggregates/baseline.h>
#include <fcpw/aggregates/bvh.h>
#ifdef FCPW_USE_ENOKI
    #include <fcpw/aggregates/mbvh.h>
#else
// make dummy types to make it compile
template <size_t DIM>
struct MsnchNode {};

template <size_t DIM>
struct MbvhNode {};
#endif
#include <map>
#include <thread>

namespace fcpw {

template<size_t DIM>
inline Scene<DIM>::Scene():
sceneData(new SceneData<DIM>())
{

}

template<size_t DIM>
inline void Scene<DIM>::setObjectTypes(const std::vector<std::vector<PrimitiveType>>& objectTypes)
{
    std::cerr << "setObjectTypes(): DIM: " << DIM << std::endl;
    exit(EXIT_FAILURE);
}

template<>
inline void Scene<2>::setObjectTypes(const std::vector<std::vector<PrimitiveType>>& objectTypes)
{
    // clear old data
    sceneData->clearAggregateData();
    sceneData->clearObjectData();

    // initialize soup and object vectors
    int nObjects = (int)objectTypes.size();
    int nLineSegmentObjects = 0;
    sceneData->soups.resize(nObjects);
    sceneData->instanceTransforms.resize(nObjects);

    for (int i = 0; i < nObjects; i++) {
        if (objectTypes[i].size() != 1) {
            std::cerr << "Mixed or empty primitive types are not supported!" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (objectTypes[i][0] == PrimitiveType::LineSegment) {
            sceneData->soupToObjectsMap[i].emplace_back(std::make_pair(ObjectType::LineSegments,
                                                                       nLineSegmentObjects));
            nLineSegmentObjects++;
        }
    }

    sceneData->lineSegmentObjects.resize(nLineSegmentObjects);
}

template<>
inline void Scene<3>::setObjectTypes(const std::vector<std::vector<PrimitiveType>>& objectTypes)
{
    // clear old data
    sceneData->clearAggregateData();
    sceneData->clearObjectData();

    // initialize soup and object vectors
    int nObjects = (int)objectTypes.size();
    int nTriangleObjects = 0;
    sceneData->soups.resize(nObjects);
    sceneData->instanceTransforms.resize(nObjects);

    for (int i = 0; i < nObjects; i++) {
        if (objectTypes[i].size() != 1) {
            std::cerr << "Mixed or empty primitive types are not supported!" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (objectTypes[i][0] == PrimitiveType::Triangle) {
            sceneData->soupToObjectsMap[i].emplace_back(std::make_pair(ObjectType::Triangles,
                                                                       nTriangleObjects));
            nTriangleObjects++;
        }
    }

    sceneData->triangleObjects.resize(nTriangleObjects);
}

template<size_t DIM>
inline void Scene<DIM>::setObjectVertexCount(int nVertices, int objectIndex)
{
    sceneData->soups[objectIndex].positions.resize(nVertices);
}

template<size_t DIM>
inline void Scene<DIM>::setObjectLineSegmentCount(int nLineSegments, int objectIndex)
{
    std::cerr << "setObjectLineSegmentCount(): DIM: " << DIM << std::endl;
    exit(EXIT_FAILURE);
}

template<>
inline void Scene<2>::setObjectLineSegmentCount(int nLineSegments, int objectIndex)
{
    // resize soup indices
    PolygonSoup<2>& soup = sceneData->soups[objectIndex];
    int nIndices = (int)soup.indices.size();
    soup.indices.resize(nIndices + 2*nLineSegments);

    // allocate line segments
    const std::vector<std::pair<ObjectType, int>>& objectsMap = sceneData->soupToObjectsMap[objectIndex];
    for (int i = 0; i < (int)objectsMap.size(); i++) {
        if (objectsMap[i].first == ObjectType::LineSegments) {
            int lineSegmentObjectIndex = objectsMap[i].second;
            sceneData->lineSegmentObjects[lineSegmentObjectIndex] =
                    std::unique_ptr<std::vector<LineSegment>>(new std::vector<LineSegment>(nLineSegments));
            break;
        }
    }
}

template<size_t DIM>
inline void Scene<DIM>::setObjectTriangleCount(int nTriangles, int objectIndex)
{
    std::cerr << "setObjectTriangleCount(): DIM: " << DIM << std::endl;
    exit(EXIT_FAILURE);
}

template<>
inline void Scene<3>::setObjectTriangleCount(int nTriangles, int objectIndex)
{
    // resize soup indices
    PolygonSoup<3>& soup = sceneData->soups[objectIndex];
    int nIndices = (int)soup.indices.size();
    soup.indices.resize(nIndices + 3*nTriangles);

    // allocate triangles
    const std::vector<std::pair<ObjectType, int>>& objectsMap = sceneData->soupToObjectsMap[objectIndex];
    for (int i = 0; i < (int)objectsMap.size(); i++) {
        if (objectsMap[i].first == ObjectType::Triangles) {
            int triangleObjectIndex = objectsMap[i].second;
            sceneData->triangleObjects[triangleObjectIndex] =
                    std::unique_ptr<std::vector<Triangle>>(new std::vector<Triangle>(nTriangles));
            break;
        }
    }
}

template<size_t DIM>
inline void Scene<DIM>::setObjectVertex(const Vector<DIM>& position, int vertexIndex, int objectIndex)
{
    sceneData->soups[objectIndex].positions[vertexIndex] = position;
}

template<size_t DIM>
inline void Scene<DIM>::setObjectLineSegment(const Vector2i& indices, int lineSegmentIndex, int objectIndex)
{
    std::cerr << "setObjectLineSegment(): DIM: " << DIM << std::endl;
    exit(EXIT_FAILURE);
}

template<>
inline void Scene<2>::setObjectLineSegment(const Vector2i& indices, int lineSegmentIndex, int objectIndex)
{
    // update soup indices
    PolygonSoup<2>& soup = sceneData->soups[objectIndex];
    soup.indices[2*lineSegmentIndex + 0] = indices[0];
    soup.indices[2*lineSegmentIndex + 1] = indices[1];

    // update line segment indices
    int lineSegmentObjectIndex = sceneData->soupToObjectsMap[objectIndex][0].second;
    LineSegment& lineSegment = (*sceneData->lineSegmentObjects[lineSegmentObjectIndex])[lineSegmentIndex];
    lineSegment.soup = &soup;
    lineSegment.indices[0] = indices[0];
    lineSegment.indices[1] = indices[1];
    lineSegment.setIndex(lineSegmentIndex);
}

template<size_t DIM>
inline void Scene<DIM>::setObjectTriangle(const Vector3i& indices, int triangleIndex, int objectIndex)
{
    std::cerr << "setObjectTriangle(): DIM: " << DIM << std::endl;
    exit(EXIT_FAILURE);
}

template<>
inline void Scene<3>::setObjectTriangle(const Vector3i& indices, int triangleIndex, int objectIndex)
{
    // update soup indices
    PolygonSoup<3>& soup = sceneData->soups[objectIndex];
    soup.indices[3*triangleIndex + 0] = indices[0];
    soup.indices[3*triangleIndex + 1] = indices[1];
    soup.indices[3*triangleIndex + 2] = indices[2];

    // update triangle indices
    int triangleObjectIndex = sceneData->soupToObjectsMap[objectIndex][0].second;
    Triangle& triangle = (*sceneData->triangleObjects[triangleObjectIndex])[triangleIndex];
    triangle.soup = &soup;
    triangle.indices[0] = indices[0];
    triangle.indices[1] = indices[1];
    triangle.indices[2] = indices[2];
    triangle.setIndex(triangleIndex);
}

template<size_t DIM>
inline void Scene<DIM>::setObjectCount(int nObjects)
{
    // clear old data
    sceneData->clearAggregateData();
    sceneData->clearObjectData();

    // initialize soup
    sceneData->soups.resize(nObjects);
    sceneData->instanceTransforms.resize(nObjects);
}

template<size_t DIM>
inline void Scene<DIM>::setObjectVertices(const std::vector<Vector<DIM>>& positions, int objectIndex)
{
    int nVertices = (int)positions.size();
    setObjectVertexCount(nVertices, objectIndex);

    for (int i = 0; i < nVertices; i++) {
        setObjectVertex(positions[i], i, objectIndex);
    }
}

template<size_t DIM>
inline void Scene<DIM>::setObjectLineSegments(const std::vector<Vector2i>& indices, int objectIndex)
{
    std::cerr << "setObjectLineSegments(): DIM: " << DIM << std::endl;
    exit(EXIT_FAILURE);
}

template<>
inline void Scene<2>::setObjectLineSegments(const std::vector<Vector2i>& indices, int objectIndex)
{
    const std::vector<std::pair<ObjectType, int>>& objectsMap = sceneData->soupToObjectsMap[objectIndex];
    for (int i = 0; i < (int)objectsMap.size(); i++) {
        if (objectsMap[i].first == ObjectType::LineSegments) {
            std::cerr << "Already set line segments!" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    // create line segments
    int nLineSegments = (int)indices.size();
    int nLineSegmentObjects = (int)sceneData->lineSegmentObjects.size();
    sceneData->soupToObjectsMap[objectIndex].emplace_back(std::make_pair(ObjectType::LineSegments,
                                                                         nLineSegmentObjects));
    sceneData->lineSegmentObjects.resize(nLineSegmentObjects + 1);
    sceneData->lineSegmentObjects[nLineSegmentObjects] =
        std::unique_ptr<std::vector<LineSegment>>(new std::vector<LineSegment>(nLineSegments));

    // resize soup indices
    PolygonSoup<2>& soup = sceneData->soups[objectIndex];
    int nIndices = (int)soup.indices.size();
    soup.indices.resize(nIndices + 2*nLineSegments);

    // set line segments
    for (int i = 0; i < nLineSegments; i++) {
        setObjectLineSegment(indices[i], i, objectIndex);
    }
}

template<size_t DIM>
inline void Scene<DIM>::setObjectTriangles(const std::vector<Vector3i>& indices, int objectIndex)
{
    std::cerr << "setObjectTriangles(): DIM: " << DIM << std::endl;
    exit(EXIT_FAILURE);
}

template<>
inline void Scene<3>::setObjectTriangles(const std::vector<Vector3i>& indices, int objectIndex)
{
    const std::vector<std::pair<ObjectType, int>>& objectsMap = sceneData->soupToObjectsMap[objectIndex];
    for (int i = 0; i < (int)objectsMap.size(); i++) {
        if (objectsMap[i].first == ObjectType::Triangles) {
            std::cerr << "Already set triangles!" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    // create line segment object
    int nTriangles = (int)indices.size();
    int nTriangleObjects = (int)sceneData->triangleObjects.size();
    sceneData->soupToObjectsMap[objectIndex].emplace_back(std::make_pair(ObjectType::Triangles,
                                                                         nTriangleObjects));
    sceneData->triangleObjects.resize(nTriangleObjects + 1);
    sceneData->triangleObjects[nTriangleObjects] =
        std::unique_ptr<std::vector<Triangle>>(new std::vector<Triangle>(nTriangles));

    // resize soup indices
    PolygonSoup<3>& soup = sceneData->soups[objectIndex];
    int nIndices = (int)soup.indices.size();
    soup.indices.resize(nIndices + 3*nTriangles);

    // set triangles
    for (int i = 0; i < nTriangles; i++) {
        setObjectTriangle(indices[i], i, objectIndex);
    }
}

template<size_t DIM>
inline void Scene<DIM>::setObjectInstanceTransforms(const std::vector<Transform<DIM>>& transforms, int objectIndex)
{
    std::vector<Transform<DIM>>& objectTransforms = sceneData->instanceTransforms[objectIndex];
    objectTransforms.insert(objectTransforms.end(), transforms.begin(), transforms.end());
}

template<size_t DIM>
inline void Scene<DIM>::setCsgTreeNode(const CsgTreeNode& csgTreeNode, int nodeIndex)
{
    sceneData->csgTree[nodeIndex] = csgTreeNode;
}

inline int assignEdgeIndices(const std::vector<Triangle>& triangles, PolygonSoup<3>& soup)
{
    int E = 0;
    int N = (int)triangles.size();
    std::map<std::pair<int, int>, int> indexMap;
    soup.eIndices.clear();

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 3; j++) {
            int k = (j + 1)%3;
            int I = triangles[i].indices[j];
            int J = triangles[i].indices[k];
            if (I > J) std::swap(I, J);
            std::pair<int, int> e(I, J);

            if (indexMap.find(e) == indexMap.end()) indexMap[e] = E++;
            soup.eIndices.emplace_back(indexMap[e]);
        }
    }

    return E;
}

template<size_t DIM>
inline void Scene<DIM>::computeSilhouettes(const std::function<bool(float, int)>& ignoreSilhouette)
{
    std::cerr << "computeSilhouettes(): DIM: " << DIM << std::endl;
    exit(EXIT_FAILURE);
}

template<>
inline void Scene<2>::computeSilhouettes(const std::function<bool(float, int)>& ignoreSilhouette)
{
    int nLineSegmentObjects = (int)sceneData->lineSegmentObjects.size();
    sceneData->silhouetteVertexObjects.resize(nLineSegmentObjects);
    sceneData->ignoreSilhouette = ignoreSilhouette;

    for (auto& kv: sceneData->soupToObjectsMap) {
        int objectIndex = kv.first;
        const std::vector<std::pair<ObjectType, int>>& objectsMap = kv.second;
        PolygonSoup<2>& soup = sceneData->soups[objectIndex];

        if (objectsMap[0].first == ObjectType::LineSegments) {
            // allocate silhouette vertices
            int lineSegmentObjectIndex = objectsMap[0].second;
            int nLineSegments = (int)soup.indices.size()/2;
            int nSilhouetteVertices = (int)soup.positions.size();
            sceneData->silhouetteVertexObjects[lineSegmentObjectIndex] =
                std::unique_ptr<std::vector<SilhouetteVertex>>(new std::vector<SilhouetteVertex>(nSilhouetteVertices));

            // assign soup and indices to silhouette vertices
            for (int i = 0; i < nLineSegments; i++) {
                const LineSegment& lineSegment = (*sceneData->lineSegmentObjects[lineSegmentObjectIndex])[i];

                SilhouetteVertex& silhouetteVertex1 = (*sceneData->silhouetteVertexObjects[lineSegmentObjectIndex])[lineSegment.indices[0]];
                silhouetteVertex1.soup = &soup;
                silhouetteVertex1.indices[1] = lineSegment.indices[0];
                silhouetteVertex1.indices[2] = lineSegment.indices[1];
                silhouetteVertex1.setIndex(lineSegment.indices[0]);

                SilhouetteVertex& silhouetteVertex2 = (*sceneData->silhouetteVertexObjects[lineSegmentObjectIndex])[lineSegment.indices[1]];
                silhouetteVertex2.soup = &soup;
                silhouetteVertex2.indices[0] = lineSegment.indices[0];
                silhouetteVertex2.indices[1] = lineSegment.indices[1];
                silhouetteVertex2.setIndex(lineSegment.indices[1]);
            }
        }
    }
}

template<>
inline void Scene<3>::computeSilhouettes(const std::function<bool(float, int)>& ignoreSilhouette)
{
    int nTriangleObjects = (int)sceneData->triangleObjects.size();
    sceneData->silhouetteEdgeObjects.resize(nTriangleObjects);
    sceneData->ignoreSilhouette = ignoreSilhouette;

    for (auto& kv: sceneData->soupToObjectsMap) {
        int objectIndex = kv.first;
        const std::vector<std::pair<ObjectType, int>>& objectsMap = kv.second;
        PolygonSoup<3>& soup = sceneData->soups[objectIndex];

        if (objectsMap[0].first == ObjectType::Triangles) {
            // allocate silhouette edges
            int triangleObjectIndex = objectsMap[0].second;
            int nTriangles = (int)soup.indices.size()/3;
            int nSilhouetteEdges = assignEdgeIndices(*sceneData->triangleObjects[triangleObjectIndex], soup);
            sceneData->silhouetteEdgeObjects[triangleObjectIndex] =
                std::unique_ptr<std::vector<SilhouetteEdge>>(new std::vector<SilhouetteEdge>(nSilhouetteEdges));

            // assign soup and indices to silhouette edges
            for (int i = 0; i < nTriangles; i++) {
                const Triangle& triangle = (*sceneData->triangleObjects[triangleObjectIndex])[i];

                for (int j = 0; j < 3; j++) {
                    int I = j - 1 < 0 ? 2 : j - 1;
                    int J = j + 0;
                    int K = j + 1 > 2 ? 0 : j + 1;
                    int eIndex = soup.eIndices[3*triangle.getIndex() + j];

                    float orientation = 1;
                    if (triangle.indices[J] > triangle.indices[K]) {
                        std::swap(J, K);
                        orientation *= -1;
                    }

                    SilhouetteEdge& silhouetteEdge = (*sceneData->silhouetteEdgeObjects[triangleObjectIndex])[eIndex];
                    silhouetteEdge.soup = &soup;
                    silhouetteEdge.indices[orientation == 1 ? 0 : 3] = triangle.indices[I];
                    silhouetteEdge.indices[1] = triangle.indices[J];
                    silhouetteEdge.indices[2] = triangle.indices[K];
                    silhouetteEdge.setIndex(eIndex);
                }
            }
        }
    }
}

template<size_t DIM, typename PrimitiveType>
inline void computeNormals(const std::vector<PrimitiveType>& primitives,
                           PolygonSoup<DIM>& soup, bool computeWeighted)
{
    // do nothing
}

template<>
inline void computeNormals<2, LineSegment>(const std::vector<LineSegment>& lineSegments,
                                           PolygonSoup<2>& soup, bool computeWeighted)
{
    int N = (int)lineSegments.size();
    int V = (int)soup.positions.size();
    soup.vNormals.clear();
    soup.vNormals.resize(V, Vector2::Zero());

    for (int i = 0; i < N; i++) {
        Vector2 n = lineSegments[i].normal(true);
        float a = computeWeighted ? lineSegments[i].surfaceArea() : 1.0f;

        soup.vNormals[lineSegments[i].indices[0]] += a*n;
        soup.vNormals[lineSegments[i].indices[1]] += a*n;
    }

    for (int i = 0; i < V; i++) {
        soup.vNormals[i].normalize();
    }
}

template<>
inline void computeNormals<3, Triangle>(const std::vector<Triangle>& triangles,
                                        PolygonSoup<3>& soup, bool computeWeighted)
{
    int N = (int)triangles.size();
    int V = (int)soup.positions.size();
    int E = assignEdgeIndices(triangles, soup);
    soup.vNormals.clear();
    soup.eNormals.clear();
    soup.vNormals.resize(V, Vector3::Zero());
    soup.eNormals.resize(E, Vector3::Zero());

    for (int i = 0; i < N; i++) {
        Vector3 n = triangles[i].normal(true);
        float area = triangles[i].surfaceArea();

        for (int j = 0; j < 3; j++) {
            float angle = computeWeighted ? triangles[i].angle(j) : 1.0f;

            soup.vNormals[triangles[i].indices[j]] += angle*n;
            soup.eNormals[soup.eIndices[3*triangles[i].getIndex() + j]] += area*n;
        }
    }

    for (int i = 0; i < V; i++) soup.vNormals[i].normalize();
    for (int i = 0; i < E; i++) soup.eNormals[i].normalize();
}

template<size_t DIM>
inline void Scene<DIM>::computeObjectNormals(int objectIndex, bool computeWeighted)
{
    std::cerr << "computeObjectNormals(): DIM: " << DIM << std::endl;
    exit(EXIT_FAILURE);
}

template<>
inline void Scene<2>::computeObjectNormals(int objectIndex, bool computeWeighted)
{
    const std::vector<std::pair<ObjectType, int>>& objectsMap = sceneData->soupToObjectsMap[objectIndex];
    PolygonSoup<2>& soup = sceneData->soups[objectIndex];

    if (objectsMap[0].first == ObjectType::LineSegments) {
        int lineSegmentObjectIndex = objectsMap[0].second;
        computeNormals<2, LineSegment>(*sceneData->lineSegmentObjects[lineSegmentObjectIndex],
                                       soup, computeWeighted);
    }
}

template<>
inline void Scene<3>::computeObjectNormals(int objectIndex, bool computeWeighted)
{
    const std::vector<std::pair<ObjectType, int>>& objectsMap = sceneData->soupToObjectsMap[objectIndex];
    PolygonSoup<3>& soup = sceneData->soups[objectIndex];

    if (objectsMap[0].first == ObjectType::Triangles) {
        int triangleObjectIndex = objectsMap[0].second;
        computeNormals<3, Triangle>(*sceneData->triangleObjects[triangleObjectIndex],
                                    soup, computeWeighted);
    }
}

template<typename NodeType, typename PrimitiveType>
inline void sortLineSegmentSoupPositions(const std::vector<NodeType>& flatTree,
                                         std::vector<PrimitiveType *>& lineSegments,
                                         PolygonSoup<2>& soup)
{
    int V = (int)soup.positions.size();
    std::vector<Vector2> sortedPositions(V), sortedVertexNormals(V);
    soup.indexMap.resize(V);
    std::fill(soup.indexMap.begin(), soup.indexMap.end(), -1);
    int v = 0;

    // collect sorted positions, updating line segment and soup indices
    for (int i = 0; i < (int)flatTree.size(); i++) {
        const NodeType& node(flatTree[i]);

        for (int j = 0; j < node.nReferences; j++) { // leaf node if nReferences > 0
            int referenceIndex = node.referenceOffset + j;
            PrimitiveType *lineSegment = lineSegments[referenceIndex];

            for (int k = 0; k < 2; k++) {
                int vIndex = lineSegment->indices[k];

                if (soup.indexMap[vIndex] == -1) {
                    sortedPositions[v] = soup.positions[vIndex];
                    if (soup.vNormals.size() > 0) sortedVertexNormals[v] = soup.vNormals[vIndex];
                    soup.indexMap[vIndex] = v++;
                }

                soup.indices[2*lineSegment->getIndex() + k] = soup.indexMap[vIndex];
                lineSegment->indices[k] = soup.indexMap[vIndex];
            }
        }
    }

    // update to sorted positions
    soup.positions = std::move(sortedPositions);
    if (soup.vNormals.size() > 0) soup.vNormals = std::move(sortedVertexNormals);
}

template<typename NodeType, typename PrimitiveType>
inline void sortTriangleSoupPositions(const std::vector<NodeType>& flatTree,
                                      std::vector<PrimitiveType *>& triangles,
                                      PolygonSoup<3>& soup)
{
    int V = (int)soup.positions.size();
    std::vector<Vector3> sortedPositions(V), sortedVertexNormals(V);
    soup.indexMap.resize(V);
    std::fill(soup.indexMap.begin(), soup.indexMap.end(), -1);
    int v = 0;

    // collect sorted positions, updating triangle and soup indices
    for (int i = 0; i < (int)flatTree.size(); i++) {
        const NodeType& node(flatTree[i]);

        for (int j = 0; j < node.nReferences; j++) { // leaf node if nReferences > 0
            int referenceIndex = node.referenceOffset + j;
            PrimitiveType *triangle = triangles[referenceIndex];

            for (int k = 0; k < 3; k++) {
                int vIndex = triangle->indices[k];

                if (soup.indexMap[vIndex] == -1) {
                    sortedPositions[v] = soup.positions[vIndex];
                    if (soup.vNormals.size() > 0) sortedVertexNormals[v] = soup.vNormals[vIndex];
                    soup.indexMap[vIndex] = v++;
                }

                soup.indices[3*triangle->getIndex() + k] = soup.indexMap[vIndex];
                triangle->indices[k] = soup.indexMap[vIndex];
            }
        }
    }

    // update to sorted positions
    soup.positions = std::move(sortedPositions);
    if (soup.vNormals.size() > 0) soup.vNormals = std::move(sortedVertexNormals);
}

template<size_t DIM, typename NodeType, typename PrimitiveType, typename SilhouetteType>
inline void sortSoupPositions(const std::vector<NodeType>& flatTree,
                              std::vector<PrimitiveType *>& primitives,
                              std::vector<SilhouetteType *>& silhouettes,
                              PolygonSoup<DIM>& soup)
{
    // do nothing
}

template<>
inline void sortSoupPositions<2, BvhNode<2>, LineSegment, SilhouettePrimitive<2>>(const std::vector<BvhNode<2>>& flatTree,
                                                                                  std::vector<LineSegment *>& lineSegments,
                                                                                  std::vector<SilhouettePrimitive<2> *>& silhouettes,
                                                                                  PolygonSoup<2>& soup)
{
    sortLineSegmentSoupPositions<BvhNode<2>, LineSegment>(flatTree, lineSegments, soup);
}

template<>
inline void sortSoupPositions<2, SnchNode<2>, LineSegment, SilhouetteVertex>(const std::vector<SnchNode<2>>& flatTree,
                                                                             std::vector<LineSegment *>& lineSegments,
                                                                             std::vector<SilhouetteVertex *>& silhouetteVertices,
                                                                             PolygonSoup<2>& soup)
{
    sortLineSegmentSoupPositions<SnchNode<2>, LineSegment>(flatTree, lineSegments, soup);

    for (int i = 0; i < (int)lineSegments.size(); i++) {
        int index1 = lineSegments[i]->indices[0];
        int index2 = lineSegments[i]->indices[1];

        SilhouetteVertex *silhouetteVertex1 = silhouetteVertices[index1];
        silhouetteVertex1->indices[1] = index1;
        silhouetteVertex1->indices[2] = index2;
        silhouetteVertex1->setIndex(index1);

        SilhouetteVertex *silhouetteVertex2 = silhouetteVertices[index2];
        silhouetteVertex2->indices[0] = index1;
        silhouetteVertex2->indices[1] = index2;
        silhouetteVertex2->setIndex(index2);
    }
}

template<>
inline void sortSoupPositions<3, BvhNode<3>, Triangle, SilhouettePrimitive<3>>(const std::vector<BvhNode<3>>& flatTree,
                                                                               std::vector<Triangle *>& triangles,
                                                                               std::vector<SilhouettePrimitive<3> *>& silhouettes,
                                                                               PolygonSoup<3>& soup)
{
    sortTriangleSoupPositions<BvhNode<3>, Triangle>(flatTree, triangles, soup);
}

template<>
inline void sortSoupPositions<3, SnchNode<3>, Triangle, SilhouetteEdge>(const std::vector<SnchNode<3>>& flatTree,
                                                                        std::vector<Triangle *>& triangles,
                                                                        std::vector<SilhouetteEdge *>& silhouetteEdges,
                                                                        PolygonSoup<3>& soup)
{
    sortTriangleSoupPositions<SnchNode<3>, Triangle>(flatTree, triangles, soup);

    for (int i = 0; i < silhouetteEdges.size(); i++) {
        SilhouetteEdge *silhouetteEdge = silhouetteEdges[i];

        for (int j = 0; j < 4; j++) {
            int vIndex = silhouetteEdge->indices[j];
            if (vIndex != -1) silhouetteEdge->indices[j] = soup.indexMap[vIndex];
        }
    }
}

#ifdef FCPW_USE_ENOKI
template<size_t DIM,
         typename BvhNodeType,
         typename PrimitiveType,
         typename SilhouetteType,
         typename NodeType,
         typename LeafNodeType,
         typename SilhouetteLeafNodeType>
inline std::unique_ptr<Aggregate<DIM>> makeVectorizedAggregate(std::vector<PrimitiveType *>& primitives,
                                                               std::vector<SilhouetteType *>& silhouettes,
                                                               const Bvh<DIM, BvhNodeType, PrimitiveType, SilhouetteType> *bvh,
                                                               bool printStats)
{
    using MbvhType = Mbvh<FCPW_SIMD_WIDTH, DIM,
                          PrimitiveType,
                          SilhouetteType,
                          NodeType,
                          LeafNodeType,
                          SilhouetteLeafNodeType>;
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    std::unique_ptr<MbvhType> mbvh = std::unique_ptr<MbvhType>(new MbvhType(primitives, silhouettes));
    mbvh->template initialize<BvhNodeType>(bvh);

    if (printStats) {
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        duration<double> timeSpan = duration_cast<duration<double>>(t2 - t1);
        std::cout << FCPW_MBVH_BRANCHING_FACTOR << "-BVH construction time: " << timeSpan.count() << " seconds" << std::endl;
        mbvh->printStats();
    }

    return mbvh;
}
#endif

template<size_t DIM,
         typename NodeType,
         typename PrimitiveType,
         typename SilhouetteType,
         typename VectorizedNodeType>
inline std::unique_ptr<Aggregate<DIM>> makeAggregate(const AggregateType& aggregateType,
                                                     std::vector<PrimitiveType *>& primitives,
                                                     std::vector<SilhouetteType *>& silhouettes,
                                                     bool vectorize, bool printStats,
                                                     SortPositionsFunc<DIM, NodeType, PrimitiveType, SilhouetteType> sortPositions={},
                                                     const std::function<bool(float, int)>& ignoreSilhouette={})
{
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    std::unique_ptr<Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>> bvh = nullptr;
    bool packLeaves = false;
    int leafSize = 4;

#ifdef FCPW_USE_ENOKI
    if (vectorize) {
        packLeaves = true;
        leafSize = FCPW_SIMD_WIDTH;
    }
#endif

    if (aggregateType == AggregateType::Bvh_LongestAxisCenter) {
        bvh = std::unique_ptr<Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>>(new Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>(
                CostHeuristic::LongestAxisCenter, primitives, silhouettes, sortPositions, ignoreSilhouette, false, leafSize));

    } else if (aggregateType == AggregateType::Bvh_SurfaceArea) {
        bvh = std::unique_ptr<Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>>(new Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>(
                CostHeuristic::SurfaceArea, primitives, silhouettes, sortPositions, ignoreSilhouette, packLeaves, leafSize));

    } else if (aggregateType == AggregateType::Bvh_OverlapSurfaceArea) {
        bvh = std::unique_ptr<Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>>(new Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>(
                CostHeuristic::OverlapSurfaceArea, primitives, silhouettes, sortPositions, ignoreSilhouette, packLeaves, leafSize));

    } else if (aggregateType == AggregateType::Bvh_Volume) {
        bvh = std::unique_ptr<Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>>(new Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>(
                CostHeuristic::Volume, primitives, silhouettes, sortPositions, ignoreSilhouette, packLeaves, leafSize));

    } else if (aggregateType == AggregateType::Bvh_OverlapVolume) {
        bvh = std::unique_ptr<Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>>(new Bvh<DIM, NodeType, PrimitiveType, SilhouetteType>(
                CostHeuristic::OverlapVolume, primitives, silhouettes, sortPositions, ignoreSilhouette, packLeaves, leafSize));

    } else {
        return std::unique_ptr<Baseline<DIM, PrimitiveType, SilhouetteType>>(
                new Baseline<DIM, PrimitiveType, SilhouetteType>(primitives, silhouettes));
    }

    if (printStats) {
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        duration<double> timeSpan = duration_cast<duration<double>>(t2 - t1);
        std::cout << "BVH construction time: " << timeSpan.count() << " seconds" << std::endl;
        bvh->printStats();
    }

#ifdef FCPW_USE_ENOKI
    if (vectorize) {
        return makeVectorizedAggregate<DIM, NodeType,
                                       PrimitiveType,
                                       SilhouetteType,
                                       VectorizedNodeType,
                                       MbvhLeafNode<FCPW_SIMD_WIDTH, DIM>,
                                       MbvhSilhouetteLeafNode<FCPW_SIMD_WIDTH, DIM>>(
                                        primitives, silhouettes, bvh.get(), printStats);
    }
#endif

    return bvh;
}

template<size_t DIM>
inline void buildGeometricAggregates(const AggregateType& aggregateType, bool vectorize, bool printStats,
                                     const std::function<bool(float, int)>& ignoreSilhouette,
                                     std::unique_ptr<SceneData<DIM>>& sceneData,
                                     std::vector<std::unique_ptr<Aggregate<DIM>>>& objectAggregates)
{
    std::cerr << "buildGeometricAggregates(): DIM: " << DIM << std::endl;
    exit(EXIT_FAILURE);
}

template<>
inline void buildGeometricAggregates<2>(const AggregateType& aggregateType, bool vectorize, bool printStats,
                                        const std::function<bool(float, int)>& ignoreSilhouette,
                                        std::unique_ptr<SceneData<2>>& sceneData,
                                        std::vector<std::unique_ptr<Aggregate<2>>>& objectAggregates)
{
    // allocate space for line segment object ptrs
    int nObjects = (int)sceneData->soups.size();
    int nLineSegmentObjectPtrs = 0;

    for (int i = 0; i < nObjects; i++) {
        const std::vector<std::pair<ObjectType, int>>& objectsMap = sceneData->soupToObjectsMap[i];
        if (objectsMap[0].first == ObjectType::LineSegments) nLineSegmentObjectPtrs++;
    }

    objectAggregates.resize(nObjects);
    sceneData->lineSegmentObjectPtrs.resize(nLineSegmentObjectPtrs);
    sceneData->silhouetteVertexObjectPtrs.resize(nLineSegmentObjectPtrs);

    // populate the object ptrs and make their aggregates
    int nAggregates = 0;
    nLineSegmentObjectPtrs = 0;

    for (int i = 0; i < nObjects; i++) {
        const std::vector<std::pair<ObjectType, int>>& objectsMap = sceneData->soupToObjectsMap[i];

        if (objectsMap[0].first == ObjectType::LineSegments) {
            // soup contains line segments, set line segment object ptrs
            int lineSegmentObjectIndex = objectsMap[0].second;
            std::vector<LineSegment>& lineSegmentObject = *sceneData->lineSegmentObjects[lineSegmentObjectIndex];
            std::vector<LineSegment *>& lineSegmentObjectPtr = sceneData->lineSegmentObjectPtrs[nLineSegmentObjectPtrs];

            for (int j = 0; j < (int)lineSegmentObject.size(); j++) {
                lineSegmentObjectPtr.emplace_back(&lineSegmentObject[j]);
            }

            if (sceneData->silhouetteVertexObjects.size() > 0) {
                // soup contains silhouette vertices, set silhouette vertex object ptrs
                std::vector<SilhouetteVertex>& silhouetteVertexObject = *sceneData->silhouetteVertexObjects[lineSegmentObjectIndex];
                std::vector<SilhouetteVertex *>& silhouetteVertexObjectPtr = sceneData->silhouetteVertexObjectPtrs[nLineSegmentObjectPtrs];

                for (int j = 0; j < (int)silhouetteVertexObject.size(); j++) {
                    silhouetteVertexObjectPtr.emplace_back(&silhouetteVertexObject[j]);
                }

                using SortLineSegmentPositionsFunc = std::function<void(const std::vector<SnchNode<2>>&,
                                                                        std::vector<LineSegment *>&,
                                                                        std::vector<SilhouetteVertex *>&)>;
                SortLineSegmentPositionsFunc sortLineSegmentPositions = std::bind(&sortSoupPositions<2, SnchNode<2>, LineSegment, SilhouetteVertex>,
                                                                                  std::placeholders::_1, std::placeholders::_2,
                                                                                  std::placeholders::_3, std::ref(sceneData->soups[i]));
                objectAggregates[i] = makeAggregate<2, SnchNode<2>, LineSegment, SilhouetteVertex, MsnchNode<2>>(
                    aggregateType, lineSegmentObjectPtr, silhouetteVertexObjectPtr, vectorize, printStats, sortLineSegmentPositions, ignoreSilhouette);

            } else {
                using SortLineSegmentPositionsFunc = std::function<void(const std::vector<BvhNode<2>>&,
                                                                        std::vector<LineSegment *>&,
                                                                        std::vector<SilhouettePrimitive<2> *>&)>;
                SortLineSegmentPositionsFunc sortLineSegmentPositions = std::bind(&sortSoupPositions<2, BvhNode<2>, LineSegment, SilhouettePrimitive<2>>,
                                                                                  std::placeholders::_1, std::placeholders::_2,
                                                                                  std::placeholders::_3, std::ref(sceneData->soups[i]));
                objectAggregates[i] = makeAggregate<2, BvhNode<2>, LineSegment, SilhouettePrimitive<2>, MbvhNode<2>>(
                    aggregateType, lineSegmentObjectPtr, sceneData->silhouetteObjectPtrStub, vectorize, printStats, sortLineSegmentPositions);
            }

            nLineSegmentObjectPtrs++;
        }

        objectAggregates[i]->setIndex(nAggregates++);
    }
}

template<>
inline void buildGeometricAggregates<3>(const AggregateType& aggregateType, bool vectorize, bool printStats,
                                        const std::function<bool(float, int)>& ignoreSilhouette,
                                        std::unique_ptr<SceneData<3>>& sceneData,
                                        std::vector<std::unique_ptr<Aggregate<3>>>& objectAggregates)
{
    // allocate space for line segment and triangle object ptrs
    int nObjects = (int)sceneData->soups.size();
    int nTriangleObjectPtrs = 0;

    for (int i = 0; i < nObjects; i++) {
        const std::vector<std::pair<ObjectType, int>>& objectsMap = sceneData->soupToObjectsMap[i];
        if (objectsMap[0].first == ObjectType::Triangles) nTriangleObjectPtrs++;
    }

    objectAggregates.resize(nObjects);
    sceneData->triangleObjectPtrs.resize(nTriangleObjectPtrs);
    sceneData->silhouetteEdgeObjectPtrs.resize(nTriangleObjectPtrs);

    // populate the object ptrs and make their aggregates
    int nAggregates = 0;
    nTriangleObjectPtrs = 0;

    for (int i = 0; i < nObjects; i++) {
        const std::vector<std::pair<ObjectType, int>>& objectsMap = sceneData->soupToObjectsMap[i];

        if (objectsMap[0].first == ObjectType::Triangles) {
            // soup contains triangles, set triangle object ptrs
            int triangleObjectIndex = objectsMap[0].second;
            std::vector<Triangle>& triangleObject = *sceneData->triangleObjects[triangleObjectIndex];
            std::vector<Triangle *>& triangleObjectPtr = sceneData->triangleObjectPtrs[nTriangleObjectPtrs];

            for (int j = 0; j < (int)triangleObject.size(); j++) {
                triangleObjectPtr.emplace_back(&triangleObject[j]);
            }

            if (sceneData->silhouetteEdgeObjects.size() > 0) {
                // soup contains silhouette edges, set silhouette edge object ptrs
                std::vector<SilhouetteEdge>& silhouetteEdgeObject = *sceneData->silhouetteEdgeObjects[triangleObjectIndex];
                std::vector<SilhouetteEdge *>& silhouetteEdgeObjectPtr = sceneData->silhouetteEdgeObjectPtrs[nTriangleObjectPtrs];

                for (int j = 0; j < (int)silhouetteEdgeObject.size(); j++) {
                    silhouetteEdgeObjectPtr.emplace_back(&silhouetteEdgeObject[j]);
                }

                using SortTrianglePositionsFunc = std::function<void(const std::vector<SnchNode<3>>&,
                                                                     std::vector<Triangle *>&,
                                                                     std::vector<SilhouetteEdge *>&)>;
                SortTrianglePositionsFunc sortTrianglePositions = std::bind(&sortSoupPositions<3, SnchNode<3>, Triangle, SilhouetteEdge>,
                                                                            std::placeholders::_1, std::placeholders::_2,
                                                                            std::placeholders::_3, std::ref(sceneData->soups[i]));
                objectAggregates[i] = makeAggregate<3, SnchNode<3>, Triangle, SilhouetteEdge, MsnchNode<3>>(
                    aggregateType, triangleObjectPtr, silhouetteEdgeObjectPtr, vectorize, printStats, sortTrianglePositions, ignoreSilhouette);

            } else {
                using SortTrianglePositionsFunc = std::function<void(const std::vector<BvhNode<3>>&,
                                                                     std::vector<Triangle *>&,
                                                                     std::vector<SilhouettePrimitive<3> *>&)>;
                SortTrianglePositionsFunc sortTrianglePositions = std::bind(&sortSoupPositions<3, BvhNode<3>, Triangle, SilhouettePrimitive<3>>,
                                                                            std::placeholders::_1, std::placeholders::_2,
                                                                            std::placeholders::_3, std::ref(sceneData->soups[i]));
                objectAggregates[i] = makeAggregate<3, BvhNode<3>, Triangle, SilhouettePrimitive<3>, MbvhNode<3>>(
                    aggregateType, triangleObjectPtr, sceneData->silhouetteObjectPtrStub, vectorize, printStats, sortTrianglePositions);
            }

            nTriangleObjectPtrs++;
        }

        objectAggregates[i]->setIndex(nAggregates++);
    }
}

template<size_t DIM>
inline std::unique_ptr<Aggregate<DIM>> buildCsgAggregateRecursive(
                                        int nodeIndex, std::unordered_map<int, CsgTreeNode>& csgTree,
                                        std::vector<std::unique_ptr<Aggregate<DIM>>>& aggregateInstances,
                                        int& nAggregates)
{
    const CsgTreeNode& node = csgTree[nodeIndex];
    std::unique_ptr<Aggregate<DIM>> instance1 = nullptr;
    std::unique_ptr<Aggregate<DIM>> instance2 = nullptr;

    if (node.isLeafChild1) {
        instance1 = std::move(aggregateInstances[node.child1]);

    } else {
        instance1 = buildCsgAggregateRecursive<DIM>(node.child1, csgTree, aggregateInstances, nAggregates);
        instance1->setIndex(nAggregates++);
    }

    if (node.isLeafChild2) {
        instance2 = std::move(aggregateInstances[node.child2]);

    } else {
        instance2 = buildCsgAggregateRecursive<DIM>(node.child2, csgTree, aggregateInstances, nAggregates);
        instance2->setIndex(nAggregates++);
    }

    return std::unique_ptr<CsgNode<DIM, Aggregate<DIM>, Aggregate<DIM>>>(
        new CsgNode<DIM, Aggregate<DIM>, Aggregate<DIM>>(std::move(instance1), std::move(instance2), node.operation));
}

template<size_t DIM>
inline void Scene<DIM>::build(const AggregateType& aggregateType, bool vectorize,
                              bool printStats, bool reduceMemoryFootprint)
{
    // clear old aggregate data
    sceneData->clearAggregateData();

    // build geometric aggregates
    std::vector<std::unique_ptr<Aggregate<DIM>>> objectAggregates;
    buildGeometricAggregates<DIM>(aggregateType, vectorize, printStats,
                                  sceneData->ignoreSilhouette,
                                  sceneData, objectAggregates);
    int nAggregates = (int)objectAggregates.size();

    // build aggregate instances and instance ptrs
    for (int i = 0; i < (int)sceneData->soups.size(); i++) {
        int nObjectInstances = (int)sceneData->instanceTransforms[i].size();

        if (nObjectInstances == 0) {
            sceneData->aggregateInstancePtrs.emplace_back(objectAggregates[i].get());
            sceneData->aggregateInstances.emplace_back(std::move(objectAggregates[i]));

        } else {
            std::shared_ptr<Aggregate<DIM>> aggregate = std::move(objectAggregates[i]);
            for (int j = 0; j < nObjectInstances; j++) {
                std::unique_ptr<TransformedAggregate<DIM>> transformedAggregate(
                    new TransformedAggregate<DIM>(aggregate, sceneData->instanceTransforms[i][j]));
                transformedAggregate->setIndex(nAggregates++);

                sceneData->aggregateInstancePtrs.emplace_back(transformedAggregate.get());
                sceneData->aggregateInstances.emplace_back(std::move(transformedAggregate));
            }
        }
    }

    // build root aggregate
    if (sceneData->aggregateInstances.size() == 1) {
        // clear the vectors of aggregate instances if there is only a single aggregate
        sceneData->aggregate = std::move(sceneData->aggregateInstances[0]);
        sceneData->aggregateInstancePtrs.clear();
        sceneData->aggregateInstances.clear();

    } else if (sceneData->csgTree.size() > 0) {
        // build csg tree
        sceneData->aggregate = buildCsgAggregateRecursive<DIM>(0, sceneData->csgTree,
                                                               sceneData->aggregateInstances,
                                                               nAggregates);
        sceneData->aggregate->setIndex(nAggregates++);
        sceneData->aggregateInstancePtrs.clear();
        sceneData->aggregateInstances.clear();

    } else {
        // make aggregate of aggregates
        if (sceneData->silhouetteVertexObjects.size() > 0 || sceneData->silhouetteEdgeObjects.size() > 0) {
            sceneData->aggregate = makeAggregate<DIM, SnchNode<DIM>, Aggregate<DIM>, SilhouettePrimitive<DIM>, MsnchNode<DIM>>(
                aggregateType, sceneData->aggregateInstancePtrs, sceneData->silhouetteObjectPtrStub, false, printStats);

        } else {
            sceneData->aggregate = makeAggregate<DIM, BvhNode<DIM>, Aggregate<DIM>, SilhouettePrimitive<DIM>, MbvhNode<DIM>>(
                aggregateType, sceneData->aggregateInstancePtrs, sceneData->silhouetteObjectPtrStub, false, printStats);
        }

        sceneData->aggregate->setIndex(nAggregates++);
    }

    // reduce memory footprint of aggregate
    if (reduceMemoryFootprint) {
        sceneData->instanceTransforms.clear();

        for (int i = 0; i < (int)sceneData->soups.size(); i++) {
            PolygonSoup<DIM>& soup = sceneData->soups[i];
            soup.indices.clear();
        }
    }
}

template<size_t DIM>
inline void Scene<DIM>::updateObjectVertex(const Vector<DIM>& position, int vertexIndex, int objectIndex)
{
    PolygonSoup<DIM>& soup = sceneData->soups[objectIndex];
    soup.positions[soup.indexMap[vertexIndex]] = position;
}

template<size_t DIM>
inline void Scene<DIM>::updateObjectVertices(const std::vector<Vector<DIM>>& positions, int objectIndex)
{
    int nVertices = (int)positions.size();
    for (int i = 0; i < nVertices; i++) {
        updateObjectVertex(positions[i], i, objectIndex);
    }
}

template<size_t DIM>
inline void Scene<DIM>::refit(bool printStats)
{
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    sceneData->aggregate->refit(); // refit is applied recursively to all children

    if (printStats) {
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        duration<double> timeSpan = duration_cast<duration<double>>(t2 - t1);
        std::cout << "refit time: " << timeSpan.count() << " seconds" << std::endl;
    }
}

template<size_t DIM>
inline bool Scene<DIM>::intersect(Ray<DIM>& r, Interaction<DIM>& i, bool checkForOcclusion) const
{
    return sceneData->aggregate->intersect(r, i, checkForOcclusion);
}

template<size_t DIM>
inline int Scene<DIM>::intersect(Ray<DIM>& r, std::vector<Interaction<DIM>>& is,
                                 bool checkForOcclusion, bool recordAllHits) const
{
    return sceneData->aggregate->intersect(r, is, checkForOcclusion, recordAllHits);
}

template<size_t DIM>
inline int Scene<DIM>::intersect(const BoundingSphere<DIM>& s, std::vector<Interaction<DIM>>& is,
                                 bool recordOneHit) const
{
    return sceneData->aggregate->intersect(s, is, recordOneHit);
}

template<size_t DIM>
inline int Scene<DIM>::intersect(const BoundingSphere<DIM>& s,
                                 Interaction<DIM>& i, const Vector<DIM>& randNums,
                                 const std::function<float(float)>& branchTraversalWeight) const
{
    return sceneData->aggregate->intersect(s, i, randNums, branchTraversalWeight);
}

template<size_t DIM>
inline bool Scene<DIM>::contains(const Vector<DIM>& x) const
{
    return sceneData->aggregate->contains(x);
}

template<size_t DIM>
inline bool Scene<DIM>::hasLineOfSight(const Vector<DIM>& xi, const Vector<DIM>& xj) const
{
    return sceneData->aggregate->hasLineOfSight(xi, xj);
}

template<size_t DIM>
inline bool Scene<DIM>::findClosestPoint(const Vector<DIM>& x, Interaction<DIM>& i,
                                         float squaredRadius, bool recordNormal) const
{
    BoundingSphere<DIM> s(x, squaredRadius);
    return sceneData->aggregate->findClosestPoint(s, i, recordNormal);
}

template<size_t DIM>
inline bool Scene<DIM>::findClosestSilhouettePoint(const Vector<DIM>& x, Interaction<DIM>& i,
                                                   bool flipNormalOrientation, float squaredMinRadius,
                                                   float squaredMaxRadius, float precision, bool recordNormal) const
{
    BoundingSphere<DIM> s(x, squaredMaxRadius);
    return sceneData->aggregate->findClosestSilhouettePoint(s, i, flipNormalOrientation, squaredMinRadius,
                                                            precision, recordNormal);
}

template<size_t DIM>
inline void Scene<DIM>::intersect(std::vector<Ray<DIM>>& rays,
                                  std::vector<Interaction<DIM>>& interactions,
                                  bool checkForOcclusion) const
{
    int nQueries = (int)rays.size();
    interactions.clear();
    interactions.resize(nQueries);

    auto callback = [&](int start, int end) {
        for (int i = start; i < end; i++) {
            sceneData->aggregate->intersect(rays[i], interactions[i], checkForOcclusion);
        }
    };

    int nThreads = std::thread::hardware_concurrency();
    int nQueriesPerThread = nQueries/nThreads;
    std::vector<std::thread> threads;

    for (int i = 0; i < nThreads; i++) {
        int start = i*nQueriesPerThread;
        int end = (i == nThreads - 1) ? nQueries : (i + 1)*nQueriesPerThread;
        threads.emplace_back(callback, start, end);
    }

    for (auto& t: threads) {
        t.join();
    }
}

template<size_t DIM>
inline void Scene<DIM>::intersect(const std::vector<BoundingSphere<DIM>>& boundingSpheres,
                                  std::vector<Interaction<DIM>>& interactions,
                                  const std::vector<Vector<DIM>>& randNums,
                                  const std::function<float(float)>& branchTraversalWeight) const
{
    int nQueries = (int)boundingSpheres.size();
    interactions.clear();
    interactions.resize(nQueries);

    auto callback = [&](int start, int end) {
        for (int i = start; i < end; i++) {
            sceneData->aggregate->intersect(boundingSpheres[i], interactions[i],
                                            randNums[i], branchTraversalWeight);
        }
    };

    int nThreads = std::thread::hardware_concurrency();
    int nQueriesPerThread = nQueries/nThreads;
    std::vector<std::thread> threads;

    for (int i = 0; i < nThreads; i++) {
        int start = i*nQueriesPerThread;
        int end = (i == nThreads - 1) ? nQueries : (i + 1)*nQueriesPerThread;
        threads.emplace_back(callback, start, end);
    }

    for (auto& t: threads) {
        t.join();
    }
}

template<size_t DIM>
inline void Scene<DIM>::contains(const std::vector<Vector<DIM>>& points,
                                 std::vector<uint32_t>& result) const
{
    int nQueries = (int)points.size();
    result.clear();
    result.resize(nQueries);

    auto callback = [&](int start, int end) {
        for (int i = start; i < end; i++) {
            result[i] = sceneData->aggregate->contains(points[i]) ? 1 : 0;
        }
    };

    int nThreads = std::thread::hardware_concurrency();
    int nQueriesPerThread = nQueries/nThreads;
    std::vector<std::thread> threads;

    for (int i = 0; i < nThreads; i++) {
        int start = i*nQueriesPerThread;
        int end = (i == nThreads - 1) ? nQueries : (i + 1)*nQueriesPerThread;
        threads.emplace_back(callback, start, end);
    }

    for (auto& t: threads) {
        t.join();
    }
}

template<size_t DIM>
inline void Scene<DIM>::hasLineOfSight(const std::vector<Vector<DIM>>& pointsI,
                                       const std::vector<Vector<DIM>>& pointsJ,
                                       std::vector<uint32_t>& result) const
{
    int nQueries = (int)pointsI.size();
    result.clear();
    result.resize(nQueries);

    auto callback = [&](int start, int end) {
        for (int i = start; i < end; i++) {
            result[i] = sceneData->aggregate->hasLineOfSight(pointsI[i], pointsJ[i]) ? 1 : 0;
        }
    };

    int nThreads = std::thread::hardware_concurrency();
    int nQueriesPerThread = nQueries/nThreads;
    std::vector<std::thread> threads;

    for (int i = 0; i < nThreads; i++) {
        int start = i*nQueriesPerThread;
        int end = (i == nThreads - 1) ? nQueries : (i + 1)*nQueriesPerThread;
        threads.emplace_back(callback, start, end);
    }

    for (auto& t: threads) {
        t.join();
    }
}

template<size_t DIM>
inline void Scene<DIM>::findClosestPoints(std::vector<BoundingSphere<DIM>>& boundingSpheres,
                                          std::vector<Interaction<DIM>>& interactions,
                                          bool recordNormal) const
{
    int nQueries = (int)boundingSpheres.size();
    interactions.clear();
    interactions.resize(nQueries);

    auto callback = [&](int start, int end) {
        for (int i = start; i < end; i++) {
            sceneData->aggregate->findClosestPoint(boundingSpheres[i], interactions[i], recordNormal);
        }
    };

    int nThreads = std::thread::hardware_concurrency();
    int nQueriesPerThread = nQueries/nThreads;
    std::vector<std::thread> threads;

    for (int i = 0; i < nThreads; i++) {
        int start = i*nQueriesPerThread;
        int end = (i == nThreads - 1) ? nQueries : (i + 1)*nQueriesPerThread;
        threads.emplace_back(callback, start, end);
    }

    for (auto& t: threads) {
        t.join();
    }
}

template<size_t DIM>
inline void Scene<DIM>::findClosestSilhouettePoints(std::vector<BoundingSphere<DIM>>& boundingSpheres,
                                                    std::vector<Interaction<DIM>>& interactions,
                                                    const std::vector<uint32_t>& flipNormalOrientation,
                                                    float squaredMinRadius, float precision,
                                                    bool recordNormal) const
{
    int nQueries = (int)boundingSpheres.size();
    interactions.clear();
    interactions.resize(nQueries);

    auto callback = [&](int start, int end) {
        for (int i = start; i < end; i++) {
            sceneData->aggregate->findClosestSilhouettePoint(boundingSpheres[i], interactions[i],
                                                             flipNormalOrientation[i] == 1,
                                                             squaredMinRadius, precision, recordNormal);
        }
    };

    int nThreads = std::thread::hardware_concurrency();
    int nQueriesPerThread = nQueries/nThreads;
    std::vector<std::thread> threads;

    for (int i = 0; i < nThreads; i++) {
        int start = i*nQueriesPerThread;
        int end = (i == nThreads - 1) ? nQueries : (i + 1)*nQueriesPerThread;
        threads.emplace_back(callback, start, end);
    }

    for (auto& t: threads) {
        t.join();
    }
}

template<size_t DIM>
inline SceneData<DIM>* Scene<DIM>::getSceneData()
{
    return sceneData.get();
}

} // namespace fcpw