#include <fcpw/utilities/scene_loader.h>
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include <atomic>

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

#include "args/args.hxx"

using namespace fcpw;
using namespace std::chrono;

static bool vizScene = false;
static bool plotInteriorPoints = false;
static bool checkCorrectness = false;
static bool checkPerformance = false;
static bool computeSilhouettes = false;
static int nQueries = 10000;

template<size_t DIM>
void splitBoxRecursive(BoundingBox<DIM> boundingBox,
                       std::vector<BoundingBox<DIM>>& boxes, int depth)
{
    if (depth == 0) {
        boxes.emplace_back(boundingBox);

    } else {
        int splitDim = boundingBox.maxDimension();
        float splitCoord = (boundingBox.pMin[splitDim] + boundingBox.pMax[splitDim])*0.5f;

        BoundingBox<DIM> boxLeft = boundingBox;
        boxLeft.pMax[splitDim] = splitCoord;
        splitBoxRecursive<DIM>(boxLeft, boxes, depth - 1);

        BoundingBox<DIM> boxRight = boundingBox;
        boxRight.pMin[splitDim] = splitCoord;
        splitBoxRecursive<DIM>(boxRight, boxes, depth - 1);
    }
}

template<size_t DIM>
void generateScatteredPointsAndRays(std::vector<Vector<DIM>>& scatteredPoints,
                                    std::vector<Vector<DIM>>& randomDirections,
                                    std::vector<float>& randomSquaredRadii,
                                    const BoundingBox<DIM>& boundingBox)
{
    // parition the scene bounding box into boxes
    std::vector<BoundingBox<DIM>> boxes;
    splitBoxRecursive<DIM>(boundingBox, boxes, 6);

    // generate queries in each box
    int nBoxes = (int)boxes.size();
    int nQueriesPerBox = std::ceil((float)nQueries/nBoxes);

    for (int i = 0; i < nBoxes; i++) {
        Vector<DIM> e = boxes[i].extent();

        for (int j = 0; j < nQueriesPerBox; j++) {
            Vector<DIM> o = boxes[i].pMin + e.cwiseProduct(uniformRealRandomVector<DIM>());
            Vector<DIM> d = uniformRealRandomVector<DIM>(-1.0f, 1.0f);
            float r2 = 0.1f*uniformRealRandomNumber()*e.norm();
            if (std::fabs(e[DIM - 1]) < 5*epsilon) {
                o[DIM - 1] = 0.0f;
                d[DIM - 1] = 0.0f;
            }
            d.normalize();

            scatteredPoints.emplace_back(o);
            randomDirections.emplace_back(d);
            randomSquaredRadii.emplace_back(r2);
        }
    }

    // resize if required
    scatteredPoints.resize(nQueries);
    randomDirections.resize(nQueries);
    randomSquaredRadii.resize(nQueries);
}

template<size_t DIM>
void tagInteriorPoints(const std::unique_ptr<Aggregate<DIM>>& aggregate,
                       const std::vector<Vector<DIM>>& queryPoints,
                       std::vector<bool>& isInterior)
{
    isInterior.resize(nQueries, false);
    auto run = [&](const tbb::blocked_range<int>& range) {
        for (int i = range.begin(); i < range.end(); ++i) {
            if (aggregate->contains(queryPoints[i])) {
                isInterior[i] = true;
            }
        }
    };

    tbb::blocked_range<int> range(0, nQueries);
    tbb::parallel_for(range, run);
}

template<size_t DIM>
void isolateInteriorPoints(const std::unique_ptr<Aggregate<DIM>>& aggregate,
                           const std::vector<Vector<DIM>>& queryPoints,
                           std::vector<Vector<DIM>>& interiorPoints)
{
    std::vector<bool> isInterior;
    tagInteriorPoints<DIM>(aggregate, queryPoints, isInterior);

    for (int i = 0; i < nQueries; i++) {
        if (isInterior[i]) interiorPoints.emplace_back(queryPoints[i]);
    }
}

template<size_t DIM>
void timeRayIntersectionQueries(const std::unique_ptr<Aggregate<DIM>>& aggregate,
                                const std::vector<Vector<DIM>>& rayOrigins,
                                const std::vector<Vector<DIM>>& rayDirections,
                                const std::vector<int>& indices,
                                const std::string& aggregateType,
                                bool queriesCoherent=false)
{
    std::atomic<int> totalNodesVisited(0);
    std::atomic<int> maxNodesVisited(0);
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    auto time = [&](const tbb::blocked_range<int>& range) {
        int nodesVisitedByThread = 0;
        int maxNodesVisitedByThread = 0;
        Interaction<DIM> cPrev;

        for (int i = range.begin(); i < range.end(); ++i) {
            int I = indices[i];

            int nodesVisited = 0;
            Interaction<DIM> c;
            Ray<DIM> r(rayOrigins[I], rayDirections[I]);
            bool hit = aggregate->intersectFromNode(r, c, 0, aggregate->getIndex(), nodesVisited);
            nodesVisitedByThread += nodesVisited;
            maxNodesVisitedByThread = std::max(maxNodesVisitedByThread, nodesVisited);

            if (hit) cPrev = c;
        }

        totalNodesVisited += nodesVisitedByThread;
        if (maxNodesVisited < maxNodesVisitedByThread) {
            maxNodesVisited = maxNodesVisitedByThread; // not thread-safe, but ok for test
        }
    };

    tbb::blocked_range<int> range(0, nQueries);
    tbb::parallel_for(range, time);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> timeSpan = duration_cast<duration<double>>(t2 - t1);
    std::cout << rayOrigins.size() << " ray intersection queries"
              << (queriesCoherent ? " with backtracking search" : "")
              << " took " << timeSpan.count() << " seconds with "
              << aggregateType << "; "
              << (totalNodesVisited/nQueries) << " nodes visited on avg and max "
              << maxNodesVisited << " nodes visited"
              << std::endl;
}

template<size_t DIM>
void timeSphereIntersectionQueries(const std::unique_ptr<Aggregate<DIM>>& aggregate,
                                   const std::vector<Vector<DIM>>& sphereCenters,
                                   const std::vector<float>& sphereSquaredRadii,
                                   const std::vector<int>& indices,
                                   const std::string& aggregateType)
{
    std::atomic<int> totalNodesVisited(0);
    std::atomic<int> maxNodesVisited(0);
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    auto time = [&](const tbb::blocked_range<int>& range) {
        int nodesVisitedByThread = 0;
        int maxNodesVisitedByThread = 0;

        for (int i = range.begin(); i < range.end(); ++i) {
            int I = indices[i];

            int nodesVisited = 0;
            std::vector<Interaction<DIM>> cs;
            BoundingSphere<DIM> s(sphereCenters[I], sphereSquaredRadii[I]);
            int hit = aggregate->intersectFromNode(s, cs, 0, aggregate->getIndex(), nodesVisited);
            nodesVisitedByThread += nodesVisited;
            maxNodesVisitedByThread = std::max(maxNodesVisitedByThread, nodesVisited);
        }

        totalNodesVisited += nodesVisitedByThread;
        if (maxNodesVisited < maxNodesVisitedByThread) {
            maxNodesVisited = maxNodesVisitedByThread; // not thread-safe, but ok for test
        }
    };

    tbb::blocked_range<int> range(0, nQueries);
    tbb::parallel_for(range, time);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> timeSpan = duration_cast<duration<double>>(t2 - t1);
    std::cout << sphereCenters.size() << " sphere intersection queries"
              << " took " << timeSpan.count() << " seconds with "
              << aggregateType << "; "
              << (totalNodesVisited/nQueries) << " nodes visited on avg and max "
              << maxNodesVisited << " nodes visited"
              << std::endl;
}

template<size_t DIM>
void timeClosestPointQueries(const std::unique_ptr<Aggregate<DIM>>& aggregate,
                             const std::vector<Vector<DIM>>& queryPoints,
                             const std::vector<int>& indices,
                             const std::string& aggregateType,
                             bool queriesCoherent=false)
{
    std::atomic<int> totalNodesVisited(0);
    std::atomic<int> maxNodesVisited(0);
    std::atomic<bool> stopQueries(false);
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    auto time = [&](const tbb::blocked_range<int>& range) {
        int nodesVisitedByThread = 0;
        int maxNodesVisitedByThread = 0;
        Interaction<DIM> cPrev;
        Vector<DIM> queryPrev = Vector<DIM>::Zero();

        for (int i = range.begin(); i < range.end(); ++i) {
            if (stopQueries) break;
            int I = indices[i];
            float distPrev = (queryPoints[I] - queryPrev).norm();
            float r2 = cPrev.nodeIndex == -1 ? maxFloat : std::pow((cPrev.d + distPrev)*1.25, 2);

            int nodesVisited = 0;
            Interaction<DIM> c;
            BoundingSphere<DIM> s(queryPoints[I], r2);
            bool found = aggregate->findClosestPointFromNode(s, c, 0, aggregate->getIndex(), nodesVisited);
            nodesVisitedByThread += nodesVisited;
            maxNodesVisitedByThread = std::max(maxNodesVisitedByThread, nodesVisited);

            if (found) cPrev = c;
            else {
                std::cerr << "Closest points not found!" << std::endl;
                stopQueries = true;
            }

            queryPrev = queryPoints[I];
        }

        totalNodesVisited += nodesVisitedByThread;
        if (maxNodesVisited < maxNodesVisitedByThread) {
            maxNodesVisited = maxNodesVisitedByThread; // not thread-safe, but ok for test
        }
    };

    tbb::blocked_range<int> range(0, nQueries);
    tbb::parallel_for(range, time);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> timeSpan = duration_cast<duration<double>>(t2 - t1);
    std::cout << queryPoints.size() << " closest point queries"
              << (queriesCoherent ? " with backtracking search" : "")
              << " took " << timeSpan.count() << " seconds with "
              << aggregateType << "; "
              << (totalNodesVisited/nQueries) << " nodes visited on avg and max "
              << maxNodesVisited << " nodes visited"
              << std::endl;
}

template<size_t DIM>
void timeClosestSilhouettePointQueries(const std::unique_ptr<Aggregate<DIM>>& aggregate,
                                       const std::vector<Vector<DIM>>& sphereCenters,
                                       const std::vector<int>& indices,
                                       const std::string& aggregateType)
{
    if (!computeSilhouettes) return;
    std::atomic<int> totalNodesVisited(0);
    std::atomic<int> maxNodesVisited(0);
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    std::vector<bool> isInterior;
    tagInteriorPoints<DIM>(aggregate, sphereCenters, isInterior);

    auto time = [&](const tbb::blocked_range<int>& range) {
        int nodesVisitedByThread = 0;
        int maxNodesVisitedByThread = 0;

        for (int i = range.begin(); i < range.end(); ++i) {
            int I = indices[i];

            int nodesVisited = 0;
            Interaction<DIM> c;
            bool flipNormalOrientation = !isInterior[I];
            BoundingSphere<DIM> s(sphereCenters[I], maxFloat);
            bool found = aggregate->findClosestSilhouettePointFromNode(s, c, 0, aggregate->getIndex(), nodesVisited,
                                                                       flipNormalOrientation, 1e-6f, 1e-3f);
            nodesVisitedByThread += nodesVisited;
            maxNodesVisitedByThread = std::max(maxNodesVisitedByThread, nodesVisited);
        }

        totalNodesVisited += nodesVisitedByThread;
        if (maxNodesVisited < maxNodesVisitedByThread) {
            maxNodesVisited = maxNodesVisitedByThread; // not thread-safe, but ok for test
        }
    };

    tbb::blocked_range<int> range(0, nQueries);
    tbb::parallel_for(range, time);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> timeSpan = duration_cast<duration<double>>(t2 - t1);
    std::cout << sphereCenters.size() << " closest silhouette point queries"
              << " took " << timeSpan.count() << " seconds with "
              << aggregateType << "; "
              << (totalNodesVisited/nQueries) << " nodes visited on avg and max "
              << maxNodesVisited << " nodes visited"
              << std::endl;
}

template<size_t DIM>
void testRayIntersectionQueries(const std::unique_ptr<Aggregate<DIM>>& aggregate1,
                                const std::unique_ptr<Aggregate<DIM>>& aggregate2,
                                const std::vector<Vector<DIM>>& rayOrigins,
                                const std::vector<Vector<DIM>>& rayDirections,
                                const std::vector<int>& indices)
{
    std::atomic<bool> stopQueries(false);
    auto test = [&](const tbb::blocked_range<int>& range) {
        Interaction<DIM> cPrev;

        for (int i = range.begin(); i < range.end(); ++i) {
            if (stopQueries) break;
            int I = indices[i];

            Interaction<DIM> c1;
            Ray<DIM> r1(rayOrigins[I], rayDirections[I]);
            bool hit1 = aggregate1->intersect(r1, c1);

            int nodesVisited = 0;
            Interaction<DIM> c2;
            Ray<DIM> r2(rayOrigins[I], rayDirections[I]);
            bool hit2 = aggregate2->intersectFromNode(r2, c2, 0, aggregate2->getIndex(), nodesVisited);

            if ((hit1 != hit2) || (hit1 && hit2 && c1 != c2)) {
                std::cerr << "d1: " << c1.d << " d2: " << c2.d
                          << "\np1: " << c1.p << " p2: " << c2.p
                          << "\nIntersections do not match!" << std::endl;
                stopQueries = true;
            }

            if (hit2) cPrev = c2;

            std::vector<Interaction<DIM>> c3;
            Ray<DIM> r3(rayOrigins[I], rayDirections[I]);
            int hit3 = aggregate1->intersect(r3, c3, false, true);

            nodesVisited = 0;
            std::vector<Interaction<DIM>> c4;
            Ray<DIM> r4(rayOrigins[I], rayDirections[I]);
            int hit4 = aggregate2->intersectFromNode(r4, c4, 0, aggregate2->getIndex(), nodesVisited, false, true);

            if (hit3 != hit4) {
                std::cerr << "Number of intersections do not match!"
                          << " hits1: " << hit3
                          << " hits2: " << hit4
                          << std::endl;
                stopQueries = true;
            }
        }
    };

    tbb::blocked_range<int> range(0, nQueries);
    tbb::parallel_for(range, test);
}

template<size_t DIM>
void testSphereIntersectionQueries(const std::unique_ptr<Aggregate<DIM>>& aggregate1,
                                   const std::unique_ptr<Aggregate<DIM>>& aggregate2,
                                   const std::vector<Vector<DIM>>& sphereCenters,
                                   const std::vector<float>& sphereSquaredRadii,
                                   const std::vector<int>& indices)
{
    std::atomic<bool> stopQueries(false);

    auto test = [&](const tbb::blocked_range<int>& range) {
        for (int i = range.begin(); i < range.end(); ++i) {
            if (stopQueries) break;
            int I = indices[i];

            BoundingSphere<DIM> s(sphereCenters[I], sphereSquaredRadii[I]);
            std::vector<Interaction<DIM>> c1;
            int hit1 = aggregate1->intersect(s, c1);

            int nodesVisited = 0;
            std::vector<Interaction<DIM>> c2;
            int hit2 = aggregate2->intersectFromNode(s, c2, 0, aggregate2->getIndex(), nodesVisited);

            if ((hit1 != hit2) || (hit1 && hit2 && c1.size() != c2.size())) {
                std::cerr << "Number of intersections do not match!"
                          << " hit1: " << hit1
                          << " hit2: " << hit2
                          << std::endl;
                for (int j = 0; j < (int)c1.size(); j++) {
                    std::cerr << "c1 pdf: " << c1[j].d
                              << " pIndex: " << c1[j].primitiveIndex
                              << std::endl;
                }
                for (int j = 0; j < (int)c2.size(); j++) {
                    std::cerr << "c2 pdf: " << c2[j].d
                              << " pIndex: " << c2[j].primitiveIndex
                              << std::endl;
                }
                stopQueries = true;
            }

            std::vector<Interaction<DIM>> c3;
            int hit3 = aggregate1->intersect(s, c3, true);

            nodesVisited = 0;
            std::vector<Interaction<DIM>> c4;
            int hit4 = aggregate2->intersectFromNode(s, c4, 0, aggregate2->getIndex(), nodesVisited, true);

            if ((hit3 != hit4) || (hit3 && hit4 && c3.size() != c4.size())) {
                std::cerr << "Number of intersections do not match!"
                          << " hit1: " << hit3
                          << " hit2: " << hit4
                          << std::endl;
                for (int j = 0; j < (int)c3.size(); j++) {
                    std::cerr << "c3 pdf: " << c3[j].d
                              << " pIndex: " << c3[j].primitiveIndex
                              << std::endl;
                }
                for (int j = 0; j < (int)c4.size(); j++) {
                    std::cerr << "c4 pdf: " << c4[j].d
                              << " pIndex: " << c4[j].primitiveIndex
                              << std::endl;
                }
                stopQueries = true;
            }
        }
    };

    tbb::blocked_range<int> range(0, nQueries);
    tbb::parallel_for(range, test);
}

template<size_t DIM>
void testClosestPointQueries(const std::unique_ptr<Aggregate<DIM>>& aggregate1,
                             const std::unique_ptr<Aggregate<DIM>>& aggregate2,
                             const std::vector<Vector<DIM>>& queryPoints,
                             const std::vector<int>& indices)
{
    std::atomic<bool> stopQueries(false);
    auto test = [&](const tbb::blocked_range<int>& range) {
        Interaction<DIM> cPrev;
        Vector<DIM> queryPrev = Vector<DIM>::Zero();

        for (int i = range.begin(); i < range.end(); ++i) {
            if (stopQueries) break;

            int I = indices[i];
            float distPrev = (queryPoints[I] - queryPrev).norm();
            float r2 = cPrev.nodeIndex == -1 ? maxFloat : std::pow((cPrev.d + distPrev)*1.25, 2);

            Interaction<DIM> c1;
            BoundingSphere<DIM> s1(queryPoints[I], maxFloat);
            bool found1 = aggregate1->findClosestPoint(s1, c1);

            int nodesVisited = 0;
            Interaction<DIM> c2;
            BoundingSphere<DIM> s2(queryPoints[I], r2);
            bool found2 = aggregate2->findClosestPointFromNode(s2, c2, 0, aggregate2->getIndex(), nodesVisited);

            if (found1 != found2 || std::fabs(c1.d - c2.d) > 1e-6) {
                std::cerr << "d1: " << c1.d << " d2: " << c2.d
                          << "\np1: " << c1.p << " p2: " << c2.p
                          << "\nClosest points do not match!" << std::endl;
                stopQueries = true;
            }

            if (found2) cPrev = c2;
            else {
                std::cerr << "Closest points not found!" << std::endl;
                stopQueries = true;
            }

            queryPrev = queryPoints[I];
        }
    };

    tbb::blocked_range<int> range(0, nQueries);
    tbb::parallel_for(range, test);
}

template<size_t DIM>
void testClosestSilhouettePointQueries(const std::unique_ptr<Aggregate<DIM>>& aggregate1,
                                       const std::unique_ptr<Aggregate<DIM>>& aggregate2,
                                       const std::vector<Vector<DIM>>& sphereCenters,
                                       const std::vector<int>& indices)
{
    if (!computeSilhouettes) return;
    std::atomic<bool> stopQueries(false);
    std::vector<bool> isInterior;
    tagInteriorPoints<DIM>(aggregate2, sphereCenters, isInterior);

    auto test = [&](const tbb::blocked_range<int>& range) {
        for (int i = range.begin(); i < range.end(); ++i) {
            if (stopQueries) break;
            int I = indices[i];
            bool flipNormalOrientation = !isInterior[I];

            Interaction<DIM> c1;
            BoundingSphere<DIM> s1(sphereCenters[I], maxFloat);
            bool found1 = aggregate1->findClosestSilhouettePoint(s1, c1, flipNormalOrientation, 1e-6f, 1e-3f);

            int nodesVisited = 0;
            Interaction<DIM> c2;
            BoundingSphere<DIM> s2(sphereCenters[I], maxFloat);
            bool found2 = aggregate2->findClosestSilhouettePointFromNode(s2, c2, 0, aggregate2->getIndex(), nodesVisited,
                                                                         flipNormalOrientation, 1e-6f, 1e-3f);

            if (found1 != found2 || std::fabs(c1.d - c2.d) > 1e-3f) {
                std::cerr << "d1: " << c1.d << " d2: " << c2.d
                          << "\np1: " << c1.p << " p2: " << c2.p
                          << "\npIndex1: " << c1.primitiveIndex << " pIndex2: " << c2.primitiveIndex
                          << "\nClosest silhouette points do not match!" << std::endl;
                stopQueries = true;
            }
        }
    };

    tbb::blocked_range<int> range(0, nQueries);
    tbb::parallel_for(range, test);
}

template<size_t DIM>
void visualizeScene(SceneData<DIM> *sceneData,
                    const std::vector<Vector<DIM>>& queryPoints,
                    const std::vector<Vector<DIM>>& randomDirections,
                    const std::vector<float>& randomSquaredRadii,
                    const std::vector<Vector<DIM>>& interiorPoints)
{
    // set a few options
    polyscope::options::programName = "Aggregate Tests";
    polyscope::options::verbosity = 0;
    polyscope::options::usePrefsFile = false;
    polyscope::options::autocenterStructures = false;

    // initialize polyscope
    polyscope::init();

    if (DIM == 2) {
        // register query points and interior points
        polyscope::registerPointCloud2D("Query_Points", queryPoints);
        if (plotInteriorPoints) polyscope::registerPointCloud2D("Interior_Points", interiorPoints);

        for (int i = 0; i < (int)sceneData->soups.size(); i++) {
            std::string meshName = "Polygon_Soup_" + std::to_string(i);

            if (sceneData->soupToObjectsMap[i][0].first == ObjectType::LineSegments) {
                // register curve network
                int N = (int)sceneData->soups[i].indices.size()/2;
                std::vector<std::vector<int>> indices(N, std::vector<int>(2));
                for (int j = 0; j < N; j++) {
                    for (int k = 0; k < 2; k++) {
                        indices[j][k] = sceneData->soups[i].indices[2*j + k];
                    }
                }

                polyscope::registerCurveNetwork2D(meshName, sceneData->soups[i].positions, indices);
            }
        }

        // add direction vectors
        polyscope::getPointCloud("Query_Points")->addVectorQuantity2D("Random_Directions", randomDirections);

    } else if (DIM == 3) {
        // register query points and interior points
        polyscope::registerPointCloud("Query_Points", queryPoints);
        if (plotInteriorPoints) polyscope::registerPointCloud("Interior_Points", interiorPoints);

        for (int i = 0; i < (int)sceneData->soups.size(); i++) {
            std::string meshName = "Polygon_Soup_" + std::to_string(i);

            if (sceneData->soupToObjectsMap[i][0].first == ObjectType::Triangles) {
                // register surface mesh
                int N = (int)sceneData->soups[i].indices.size()/3;
                std::vector<std::vector<int>> indices(N, std::vector<int>(3));
                for (int j = 0; j < N; j++) {
                    for (int k = 0; k < 3; k++) {
                        indices[j][k] = sceneData->soups[i].indices[3*j + k];
                    }
                }

                polyscope::registerSurfaceMesh(meshName, sceneData->soups[i].positions, indices);
            }
        }

        // add direction vectors
        polyscope::getPointCloud("Query_Points")->addVectorQuantity("Random_Directions", randomDirections);
    }

    polyscope::getPointCloud("Query_Points")->addScalarQuantity("Random_Squared_Radii", randomSquaredRadii);
    polyscope::getPointCloud("Query_Points")->setPointRadiusQuantity("Random_Squared_Radii");

    // give control to polyscope gui
    polyscope::show();
}

template<size_t DIM>
void run()
{
    // build baseline scene
    Scene<DIM> scene;
    SceneLoader<DIM> sceneLoader;
    sceneLoader.loadFiles(scene, false);
    if (computeSilhouettes) scene.computeSilhouettes();
    scene.build(AggregateType::Baseline, false, true);
    SceneData<DIM> *sceneData = scene.getSceneData();

    // generate random points and rays used to visualize csg
    BoundingBox<DIM> boundingBox = sceneData->aggregate->boundingBox();
    std::vector<Vector<DIM>> queryPoints, randomDirections;
    std::vector<float> randomSquaredRadii;
    generateScatteredPointsAndRays<DIM>(queryPoints, randomDirections, randomSquaredRadii, boundingBox);

    // generate indices
    std::vector<int> indices, shuffledIndices;
    for (int i = 0; i < nQueries; i++) {
        indices.emplace_back(i);
        shuffledIndices.emplace_back(i);
    }

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(shuffledIndices.begin(), shuffledIndices.end(), std::default_random_engine(seed));

    std::vector<std::string> bvhTypes({"Bvh_LongestAxisCenter", "Bvh_OverlapSurfaceArea",
                                       "Bvh_SurfaceArea", "Bvh_OverlapVolume", "Bvh_Volume"});

    if (checkPerformance) {
        std::cout << "Running performance tests..." << std::endl;

        // benchmark baseline queries
        //timeRayIntersectionQueries<DIM>(sceneData->aggregate, queryPoints, randomDirections,
        //                              shuffledIndices, "Baseline");
        //timeSphereIntersectionQueries<DIM>(sceneData->aggregate, queryPoints, randomSquaredRadii,
        //                                 shuffledIndices, "Baseline");
        //timeClosestPointQueries<DIM>(sceneData->aggregate, queryPoints,
        //                           shuffledIndices, "Baseline");
        //timeClosestSilhouettePointQueries<DIM>(sceneData->aggregate, queryPoints,
        //                                     shuffledIndices, "Baseline");

        // build bvh aggregates and benchmark queries
        for (int bvh = 1; bvh < 6; bvh++) {
            for (int vec = 0; vec < 2; vec++) {
                scene.build(static_cast<AggregateType>(bvh), vec == 1, true);
                sceneData = scene.getSceneData();

                timeRayIntersectionQueries<DIM>(sceneData->aggregate, queryPoints, randomDirections,
                                                shuffledIndices, bvhTypes[bvh - 1]);
                timeRayIntersectionQueries<DIM>(sceneData->aggregate, queryPoints, randomDirections,
                                                indices, bvhTypes[bvh - 1], true);
                //timeSphereIntersectionQueries<DIM>(sceneData->aggregate, queryPoints, randomSquaredRadii,
                //                                 shuffledIndices, bvhTypes[bvh - 1]);
                timeClosestPointQueries<DIM>(sceneData->aggregate, queryPoints,
                                             shuffledIndices, bvhTypes[bvh - 1]);
                timeClosestPointQueries<DIM>(sceneData->aggregate, queryPoints,
                                             indices, bvhTypes[bvh - 1], true);
                timeClosestSilhouettePointQueries<DIM>(sceneData->aggregate, queryPoints,
                                                       shuffledIndices, bvhTypes[bvh - 1]);

#ifndef FCPW_USE_ENOKI
                break;
#endif
            }

            std::cout << std::endl;
        }
    }

    if (checkCorrectness) {
        std::cout << "Running correctness tests..." << std::endl;

        // build baseline aggregate
        scene.build(AggregateType::Baseline, false, true);
        sceneData = scene.getSceneData();

        // build bvh aggregates and compare results with baseline
        Scene<DIM> bvhScene;
        sceneLoader.loadFiles(bvhScene, false);
        if (computeSilhouettes) bvhScene.computeSilhouettes();

        for (int bvh = 1; bvh < 6; bvh++) {
            std::cout << "Testing " << bvhTypes[bvh - 1] << " results against Baseline" << std::endl;

            for (int vec = 0; vec < 2; vec++) {
                bvhScene.build(static_cast<AggregateType>(bvh), vec == 1, true);
                SceneData<DIM> *bvhSceneData = bvhScene.getSceneData();

                testRayIntersectionQueries<DIM>(sceneData->aggregate, bvhSceneData->aggregate,
                                                queryPoints, randomDirections, shuffledIndices);
                testRayIntersectionQueries<DIM>(sceneData->aggregate, bvhSceneData->aggregate,
                                                queryPoints, randomDirections, indices);
                //testSphereIntersectionQueries<DIM>(sceneData->aggregate, bvhSceneData->aggregate,
                //                                 queryPoints, randomSquaredRadii, shuffledIndices);
                testClosestPointQueries<DIM>(sceneData->aggregate, bvhSceneData->aggregate,
                                             queryPoints, shuffledIndices);
                testClosestPointQueries<DIM>(sceneData->aggregate, bvhSceneData->aggregate,
                                             queryPoints, indices);
                testClosestSilhouettePointQueries<DIM>(sceneData->aggregate, bvhSceneData->aggregate,
                                                       queryPoints, shuffledIndices);

#ifndef FCPW_USE_ENOKI
                break;
#endif
            }

            std::cout << std::endl;
        }
    }

    if (vizScene) {
        // build bvh aggregate
        scene.build(AggregateType::Bvh_LongestAxisCenter, false, true);
        sceneData = scene.getSceneData();

        // isolate interior points among query points
        std::vector<Vector<DIM>> interiorPoints;
        if (plotInteriorPoints) isolateInteriorPoints<DIM>(sceneData->aggregate, queryPoints, interiorPoints);

        // visualize scene
        visualizeScene<DIM>(sceneData, queryPoints, randomDirections, randomSquaredRadii, interiorPoints);
    }
}

int main(int argc, const char *argv[]) {
    // configure the argument parser
    args::ArgumentParser parser("aggregate tests");
    args::Group group(parser, "", args::Group::Validators::DontCare);
    args::Flag vizScene(group, "bool", "visualize scene", {"vizScene"});
    args::Flag plotInteriorPoints(group, "bool", "plot interior points", {"plotInteriorPoints"});
    args::Flag checkCorrectness(group, "bool", "check aggregate correctness", {"checkCorrectness"});
    args::Flag checkPerformance(group, "bool", "check aggregate performance", {"checkPerformance"});
    args::Flag computeSilhouettes(group, "bool", "compute silhouettes", {"computeSilhouettes"});
    args::ValueFlag<int> dim(parser, "integer", "scene dimension", {"dim"});
    args::ValueFlag<int> nQueries(parser, "integer", "number of queries", {"nQueries"});
    args::ValueFlagList<std::string> lineSegmentFilenames(parser, "string", "line segment soup filenames", {"lFile"});
    args::ValueFlagList<std::string> triangleFilenames(parser, "string", "triangle soup filenames", {"tFile"});

    // parse args
    try {
        parser.ParseCLI(argc, argv);

    } catch (const args::Help&) {
        std::cout << parser;
        return 0;

    } catch (const args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    int DIM = args::get(dim);
    if (!dim) {
        std::cerr << "Please specify dimension" << std::endl;
        return EXIT_FAILURE;
    }

    if (!lineSegmentFilenames && !triangleFilenames) {
        std::cerr << "Please specify either line segment or triangle soup filenames" << std::endl;
        return EXIT_FAILURE;
    }

    // set global flags
    if (vizScene) ::vizScene = args::get(vizScene);
    if (plotInteriorPoints) ::plotInteriorPoints = args::get(plotInteriorPoints);
    if (checkCorrectness) ::checkCorrectness = args::get(checkCorrectness);
    if (checkPerformance) ::checkPerformance = args::get(checkPerformance);
    if (computeSilhouettes) ::computeSilhouettes = args::get(computeSilhouettes);
    if (nQueries) ::nQueries = args::get(nQueries);
    if (lineSegmentFilenames) {
        for (const auto lsf: args::get(lineSegmentFilenames)) {
            files.emplace_back(std::make_pair(lsf, LoadingOption::ObjLineSegments));
        }
    }
    if (triangleFilenames) {
        for (const auto tsf: args::get(triangleFilenames)) {
            files.emplace_back(std::make_pair(tsf, LoadingOption::ObjTriangles));
        }
    }

    // run app
    if (DIM == 2) run<2>();
    else if (DIM == 3) run<3>();
    return 0;
}