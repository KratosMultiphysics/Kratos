#include <fcpw/utilities/scene_loader.h>
#include <fcpw/fcpw_gpu.h>
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

static bool refitBvh = false;
static bool vizScene = false;
static bool plotInteriorPoints = false;
static bool computeSilhouettes = false;
static int nQueries = 1048576;

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
    // partition the scene bounding box into boxes
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
void runCPURayIntersectionQueries(const std::unique_ptr<Aggregate<DIM>>& aggregate,
                                  const std::vector<Vector<DIM>>& rayOrigins,
                                  const std::vector<Vector<DIM>>& rayDirections,
                                  const std::vector<int>& indices,
                                  std::vector<Interaction<DIM>>& interactions)
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    auto run = [&](const tbb::blocked_range<int>& range) {
        for (int i = range.begin(); i < range.end(); ++i) {
            int I = indices[i];

            Interaction<DIM> c;
            Ray<DIM> r(rayOrigins[I], rayDirections[I]);
            bool hit = aggregate->intersect(r, c, false);

            if (hit) {
                interactions[I] = c;
            }
        }
    };

    tbb::blocked_range<int> range(0, nQueries);
    tbb::parallel_for(range, run);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> timeSpan = duration_cast<duration<double>>(t2 - t1);
    std::cout << nQueries << " ray intersection queries"
              << " took " << timeSpan.count() << " seconds"
              << std::endl;
}

template<size_t DIM>
void runCPUSphereIntersectionQueries(const std::unique_ptr<Aggregate<DIM>>& aggregate,
                                     const std::vector<Vector<DIM>>& sphereCenters,
                                     const std::vector<float>& sphereSquaredRadii,
                                     const std::vector<Vector<DIM>>& randNums,
                                     const std::vector<int>& indices,
                                     std::vector<Interaction<DIM>>& interactions)
{
    std::function<float(float)> branchTraversalWeight = [](float r2) -> float {
        return 1.0f;
    };

    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    auto run = [&](const tbb::blocked_range<int>& range) {
        for (int i = range.begin(); i < range.end(); ++i) {
            int I = indices[i];

            Interaction<DIM> c;
            BoundingSphere<DIM> s(sphereCenters[I], sphereSquaredRadii[I]);
            int hits = aggregate->intersect(s, c, randNums[I], branchTraversalWeight);

            if (hits > 0) {
                interactions[I] = c;
            }
        }
    };

    tbb::blocked_range<int> range(0, nQueries);
    tbb::parallel_for(range, run);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> timeSpan = duration_cast<duration<double>>(t2 - t1);
    std::cout << nQueries << " sphere intersection queries"
              << " took " << timeSpan.count() << " seconds"
              << std::endl;
}

template<size_t DIM>
void runCPUClosestPointQueries(const std::unique_ptr<Aggregate<DIM>>& aggregate,
                               const std::vector<Vector<DIM>>& queryPoints,
                               const std::vector<int>& indices,
                               std::vector<Interaction<DIM>>& interactions)
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    auto run = [&](const tbb::blocked_range<int>& range) {
        for (int i = range.begin(); i < range.end(); ++i) {
            int I = indices[i];

            Interaction<DIM> c;
            BoundingSphere<DIM> s(queryPoints[I], maxFloat);
            bool found = aggregate->findClosestPoint(s, c, false);

            if (found) {
                interactions[I] = c;
            }
        }
    };

    tbb::blocked_range<int> range(0, nQueries);
    tbb::parallel_for(range, run);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> timeSpan = duration_cast<duration<double>>(t2 - t1);
    std::cout << nQueries << " closest point queries"
              << " took " << timeSpan.count() << " seconds"
              << std::endl;
}

template<size_t DIM>
void runCPUClosestSilhouettePointQueries(const std::unique_ptr<Aggregate<DIM>>& aggregate,
                                         const std::vector<Vector<DIM>>& queryPoints,
                                         const std::vector<int>& indices,
                                         const std::vector<bool>& flipNormalOrientation,
                                         std::vector<Interaction<DIM>>& interactions)
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    auto run = [&](const tbb::blocked_range<int>& range) {
        for (int i = range.begin(); i < range.end(); ++i) {
            int I = indices[i];

            Interaction<DIM> c;
            BoundingSphere<DIM> s(queryPoints[I], maxFloat);
            bool found = aggregate->findClosestSilhouettePoint(s, c, flipNormalOrientation[i], 1e-6f, 1e-3f, false);

            if (found) {
                interactions[I] = c;
            }
        }
    };

    tbb::blocked_range<int> range(0, nQueries);
    tbb::parallel_for(range, run);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> timeSpan = duration_cast<duration<double>>(t2 - t1);
    std::cout << nQueries << " closest silhouette point queries"
              << " took " << timeSpan.count() << " seconds"
              << std::endl;
}

template<size_t DIM>
void compareInteractions(const std::vector<Interaction<DIM>>& cpuInteractions,
                         const std::vector<GPUInteraction>& gpuInteractions)
{
    for (int i = 0; i < nQueries; i++) {
        const Interaction<DIM>& cpuInteraction = cpuInteractions[i];
        const GPUInteraction& gpuInteraction = gpuInteractions[i];
        bool differentIndices = (cpuInteraction.primitiveIndex == -1 && gpuInteraction.index != FCPW_GPU_UINT_MAX) ||
                                (cpuInteraction.primitiveIndex != -1 && gpuInteraction.index == FCPW_GPU_UINT_MAX);

        if (differentIndices) {
            std::cout << "#" << i << "/" << nQueries << std::endl;
            std::cout << "CPU Interaction"
                      << "\n\tp: " << cpuInteraction.p[0] << " " << cpuInteraction.p[1] << " " << cpuInteraction.p[2]
                      << "\n\tn: " << cpuInteraction.n[0] << " " << cpuInteraction.n[1] << " " << cpuInteraction.n[2]
                      << "\n\tuv: " << cpuInteraction.uv[0] << " " << cpuInteraction.uv[1]
                      << "\n\td: " << cpuInteraction.d
                      << "\n\tindex: " << cpuInteraction.primitiveIndex
                      << std::endl;
            std::cout << "GPU Interaction"
                      << "\n\tp: " << gpuInteraction.p.x << " " << gpuInteraction.p.y << " " << gpuInteraction.p.z
                      << "\n\tn: " << gpuInteraction.n.x << " " << gpuInteraction.n.y << " " << gpuInteraction.n.z
                      << "\n\tuv: " << gpuInteraction.uv.x << " " << gpuInteraction.uv.y
                      << "\n\td: " << gpuInteraction.d
                      << "\n\tindex: " << gpuInteraction.index
                      << std::endl;
            std::cout << std::endl;
        }
    }
}

template<size_t DIM>
void visualizeScene(SceneData<DIM> *sceneData,
                    const std::vector<Vector<DIM>>& queryPoints,
                    const std::vector<Vector<DIM>>& randomDirections,
                    const std::vector<float>& randomSquaredRadii,
                    const std::vector<Vector<DIM>>& interiorPoints)
{
    // set a few options
    polyscope::options::programName = "GPU Bvh Traversal Tests";
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
    // build bvh on CPU
    Scene<DIM> scene;
    SceneLoader<DIM> sceneLoader;
    sceneLoader.loadFiles(scene, false);
    if (computeSilhouettes) scene.computeSilhouettes();
    scene.build(AggregateType::Bvh_OverlapSurfaceArea, false, true);
    SceneData<DIM> *sceneData = scene.getSceneData();
    bool isLineSegmentGeometry = sceneData->lineSegmentObjects.size() > 0;

    // generate random points and rays
    BoundingBox<DIM> boundingBox = sceneData->aggregate->boundingBox();
    std::vector<Vector<DIM>> queryPoints, randomDirections;
    std::vector<float> randomSquaredRadii;
    generateScatteredPointsAndRays<DIM>(queryPoints, randomDirections, randomSquaredRadii, boundingBox);

    // generate indices, random numbers and normal orientations for queries
    std::vector<int> indices;
    std::vector<Vector<DIM>> randNums;
    std::vector<bool> flipNormalOrientation;
    std::vector<GPURay> gpuRays(nQueries);
    std::vector<GPUBoundingSphere> gpuBoundingSpheres(nQueries);
    std::vector<GPUBoundingSphere> gpuInfiniteSpheres(nQueries);
    std::vector<float3> gpuRandNums(nQueries);
    std::vector<uint32_t> gpuFlipNormalOrientation(nQueries);

    std::vector<bool> isInterior;
    tagInteriorPoints<DIM>(sceneData->aggregate, queryPoints, isInterior);

    for (int i = 0; i < nQueries; i++) {
        float3 queryPoint = float3{queryPoints[i][0],
                                   queryPoints[i][1],
                                   DIM == 2 ? 0.0f : queryPoints[i][2]};
        float3 randomDirection = float3{randomDirections[i][0],
                                        randomDirections[i][1],
                                        DIM == 2 ? 0.0f : randomDirections[i][2]};
        float randomSquaredRadius = randomSquaredRadii[i];

        indices.emplace_back(i);
        randNums.emplace_back(uniformRealRandomVector<DIM>());
        flipNormalOrientation.emplace_back(!isInterior[i]);
        gpuRays[i] = GPURay(queryPoint, randomDirection);
        gpuBoundingSpheres[i] = GPUBoundingSphere(queryPoint, randomSquaredRadius);
        gpuInfiniteSpheres[i] = GPUBoundingSphere(queryPoint, maxFloat);
        gpuRandNums[i] = float3{randNums[i][0],
                                randNums[i][1],
                                DIM == 2 ? 0.0f : randNums[i][2]};
        gpuFlipNormalOrientation[i] = flipNormalOrientation[i] ? 1 : 0;
    }

    // transfer scene to GPU
    std::filesystem::path fpcwDirectoryPath = std::filesystem::current_path().parent_path();
    GPUScene<DIM> gpuScene(fpcwDirectoryPath.string(), true);
    gpuScene.transferToGPU(scene);

    // refit GPU BVH
    if (refitBvh) {
        std::cout << "Refitting GPU BVH..." << std::endl;
        gpuScene.refit(scene);
    }

    // run GPU queries
    std::cout << "Running GPU queries..." << std::endl;
    std::vector<GPUInteraction> gpuRayInteractions;
    gpuScene.intersect(gpuRays, gpuRayInteractions);

    std::vector<GPUInteraction> gpuSphereInteractions;
    gpuScene.intersect(gpuBoundingSpheres, gpuRandNums, gpuSphereInteractions);

    std::vector<GPUInteraction> gpuCPQInteractions;
    gpuScene.findClosestPoints(gpuInfiniteSpheres, gpuCPQInteractions);

    std::vector<GPUInteraction> gpuCSPQInteractions;
    if (computeSilhouettes) {
        gpuScene.findClosestSilhouettePoints(gpuInfiniteSpheres, gpuFlipNormalOrientation, gpuCSPQInteractions);
    }

    // run CPU queries
    std::cout << "\nRunning CPU queries..." << std::endl;
    std::vector<Interaction<DIM>> cpuSphereInteractions(nQueries);
    runCPUSphereIntersectionQueries<DIM>(sceneData->aggregate, queryPoints, randomSquaredRadii,
                                         randNums, indices, cpuSphereInteractions);

    // build vectorized bvh for remaining CPU queries
    scene.build(AggregateType::Bvh_OverlapSurfaceArea, true, true); 
    sceneData = scene.getSceneData();

    std::vector<Interaction<DIM>> cpuRayInteractions(nQueries);
    runCPURayIntersectionQueries<DIM>(sceneData->aggregate, queryPoints, randomDirections,
                                      indices, cpuRayInteractions);

    std::vector<Interaction<DIM>> cpuCPQInteractions(nQueries);
    runCPUClosestPointQueries<DIM>(sceneData->aggregate, queryPoints, indices, cpuCPQInteractions);

    std::vector<Interaction<DIM>> cpuCSPQInteractions(nQueries);
    if (computeSilhouettes) {
        runCPUClosestSilhouettePointQueries<DIM>(sceneData->aggregate, queryPoints, indices,
                                                 flipNormalOrientation, cpuCSPQInteractions);
    }

    // compare CPU & GPU results
    std::cout << "\nComparing CPU & GPU ray intersection query results..." << std::endl;
    compareInteractions<DIM>(cpuRayInteractions, gpuRayInteractions);
    std::cout << "\nComparing CPU & GPU sphere intersection query results..." << std::endl;
    compareInteractions<DIM>(cpuSphereInteractions, gpuSphereInteractions);
    std::cout << "\nComparing CPU & GPU closest point query results..." << std::endl;
    compareInteractions<DIM>(cpuCPQInteractions, gpuCPQInteractions);
    if (computeSilhouettes) {
        std::cout << "\nComparing CPU & GPU closest silhouette point query results..." << std::endl;
        compareInteractions<DIM>(cpuCSPQInteractions, gpuCSPQInteractions);
    }
    std::cout << "\ndone" << std::endl;

    if (vizScene) {
        // isolate interior points among query points
        std::vector<Vector<DIM>> interiorPoints;
        if (plotInteriorPoints) isolateInteriorPoints<DIM>(sceneData->aggregate, queryPoints, interiorPoints);

        // visualize scene
        visualizeScene<DIM>(sceneData, queryPoints, randomDirections, randomSquaredRadii, interiorPoints);
    }
}

int main(int argc, const char *argv[]) {
    // configure the argument parser
    args::ArgumentParser parser("gpu bvh traversal tests");
    args::Group group(parser, "", args::Group::Validators::DontCare);
    args::Flag refitBvh(group, "bool", "refit bvh", {"refitBvh"});
    args::Flag vizScene(group, "bool", "visualize scene", {"vizScene"});
    args::Flag plotInteriorPoints(group, "bool", "plot interior points", {"plotInteriorPoints"});
    args::Flag computeSilhouettes(group, "bool", "compute silhouettes", {"computeSilhouettes"});
    args::ValueFlag<int> dim(parser, "integer", "scene dimension", {"dim"});
    args::ValueFlag<int> nQueries(parser, "integer", "number of queries", {"nQueries"});
    args::ValueFlag<std::string> lineSegmentFilename(parser, "string", "line segment soup filename", {"lFile"});
    args::ValueFlag<std::string> triangleFilename(parser, "string", "triangle soup filename", {"tFile"});

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

    if (!lineSegmentFilename && !triangleFilename) {
        std::cerr << "Please specify either line segment or triangle soup filename" << std::endl;
        return EXIT_FAILURE;
    }

    // set global flags
    if (refitBvh) ::refitBvh = args::get(refitBvh);
    if (vizScene) ::vizScene = args::get(vizScene);
    if (plotInteriorPoints) ::plotInteriorPoints = args::get(plotInteriorPoints);
    if (computeSilhouettes) ::computeSilhouettes = args::get(computeSilhouettes);
    if (nQueries) ::nQueries = args::get(nQueries);
    if (lineSegmentFilename) {
        files.emplace_back(std::make_pair(args::get(lineSegmentFilename), LoadingOption::ObjLineSegments));
    }
    if (triangleFilename) {
        files.emplace_back(std::make_pair(args::get(triangleFilename), LoadingOption::ObjTriangles));
    }

    // run app
    if (DIM == 2) run<2>();
    else if (DIM == 3) run<3>();
    return 0;
}