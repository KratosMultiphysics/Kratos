#include <fcpw/utilities/scene_loader.h>
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

#include "args/args.hxx"

using namespace fcpw;

template<size_t DIM>
void generateScatteredPointsAndRays(int nPoints, std::vector<Vector<DIM>>& scatteredPoints,
                                    std::vector<Vector<DIM>>& randomDirections,
                                    const BoundingBox<DIM>& boundingBox)
{
    Vector<DIM> e = boundingBox.extent();

    for (int i = 0; i < nPoints; i++) {
        Vector<DIM> o = boundingBox.pMin + e.cwiseProduct(uniformRealRandomVector<DIM>());
        Vector<DIM> d = uniformRealRandomVector<DIM>(-1.0f, 1.0f);
        if (std::fabs(e[DIM - 1]) < 5*epsilon) {
            o[DIM - 1] = 0.0f;
            d[DIM - 1] = 0.0f;
        }
        d.normalize();

        scatteredPoints.emplace_back(o);
        randomDirections.emplace_back(d);
    }
}

template<size_t DIM>
bool raymarch(const std::unique_ptr<Aggregate<DIM>>& aggregate,
              const BoundingBox<DIM>& boundingBox,
              Ray<DIM> r, Interaction<DIM>& i)
{
    r.tMax = 0.0f;
    Vector<DIM> x = r(r.tMax);

    while (boundingBox.contains(x)) {
        Interaction<DIM> c;
        BoundingSphere<DIM> s(x, maxFloat);
        bool found = aggregate->findClosestPoint(s, c);
        if (!found) {
            std::cout << "No closest point found while raymarching!" << std::endl;
            break;
        }

        r.tMax += c.d;
        x = r(r.tMax);

        if (c.d < 1e-5) {
            i = c;
            return true;
        }
    }

    return false;
}

template<size_t DIM>
void clampToCsg(const std::string& method,
                const std::vector<Vector<DIM>>& scatteredPoints,
                const std::vector<Vector<DIM>>& randomDirections,
                const std::unique_ptr<Aggregate<DIM>>& aggregate,
                const BoundingBox<DIM>& boundingBox,
                std::vector<Vector<DIM>>& clampedPoints)
{
    int N = (int)scatteredPoints.size();
    std::vector<Vector<DIM>> hitPoints(N);
    std::vector<bool> didHit(N, false);

    auto clamp = [&](const tbb::blocked_range<int>& range) {
        for (int i = range.begin(); i < range.end(); ++i) {
            Ray<DIM> r(scatteredPoints[i], randomDirections[i]);

            if (method == "intersect") {
                std::vector<Interaction<DIM>> cs;
                int hits = aggregate->intersect(r, cs, false, true);
                if (hits > 0) {
                    hitPoints[i] = cs[0].p;
                    didHit[i] = true;
                }

            } else {
                Interaction<DIM> c;
                if (raymarch<DIM>(aggregate, boundingBox, r, c)) {
                    hitPoints[i] = c.p;
                    didHit[i] = true;
                }
            }
        }
    };

    tbb::blocked_range<int> range(0, N);
    tbb::parallel_for(range, clamp);

    // add points for which hit was found
    for (int i = 0; i < N; i++) {
        if (didHit[i]) clampedPoints.emplace_back(hitPoints[i]);
    }
}

template<size_t DIM>
void guiCallback(const std::unique_ptr<Aggregate<DIM>>& aggregate,
                 const BoundingBox<DIM>& boundingBox,
                 std::vector<Vector<DIM>>& intersectedPoints,
                 std::vector<Vector<DIM>>& raymarchedPoints)
{
    // make ui elements 100 pixels wide, instead of full width
    ImGui::PushItemWidth(100);

    if (ImGui::Button("Add Samples to Visualize CSG")) {
        // generate random points and rays used to visualize csg
        std::vector<Vector<DIM>> scatteredPoints, randomDirections;
        generateScatteredPointsAndRays<DIM>(1000, scatteredPoints, randomDirections, boundingBox);

        // intersect and raymarch points
        clampToCsg<DIM>("intersect", scatteredPoints, randomDirections,
                        aggregate, boundingBox, intersectedPoints);
        clampToCsg<DIM>("raymarch", scatteredPoints, randomDirections,
                        aggregate, boundingBox, raymarchedPoints);

        // register point clouds
        if (DIM == 2) {
            polyscope::registerPointCloud2D("Intersected_Points", intersectedPoints);
            polyscope::registerPointCloud2D("Raymarched_Points", raymarchedPoints);

        } else if (DIM == 3) {
            polyscope::registerPointCloud("Intersected_Points", intersectedPoints);
            polyscope::registerPointCloud("Raymarched_Points", raymarchedPoints);
        }
    }

    ImGui::PopItemWidth();
}

template<size_t DIM>
void visualizeScene(SceneData<DIM> *sceneData, const BoundingBox<DIM>& boundingBox,
                    std::vector<Vector<DIM>>& intersectedPoints,
                    std::vector<Vector<DIM>>& raymarchedPoints)
{
    // set a few options
    polyscope::options::programName = "CSG Tests";
    polyscope::options::verbosity = 0;
    polyscope::options::usePrefsFile = false;
    polyscope::options::autocenterStructures = false;

    // initialize polyscope
    polyscope::init();

    if (DIM == 2) {
        for (int i = 0; i < (int)sceneData->soups.size(); i++) {
            std::string meshName = "Polygon_Soup_" + std::to_string(i);

            if (sceneData->soupToObjectsMap[i][0].first == ObjectType::LineSegments) {
                // register curve network
                int N = (int)sceneData->soups[i].indices.size()/2;
                std::vector<std::vector<int>> indices(N, std::vector<int>(2));
                const std::vector<Vector<DIM>>& positions = sceneData->soups[i].positions;
                for (int j = 0; j < N; j++) {
                    for (int k = 0; k < 2; k++) {
                        indices[j][k] = sceneData->soups[i].indices[2*j + k];
                    }
                }

                if (sceneData->instanceTransforms[i].size() > 0) {
                    for (int j = 0; j < (int)sceneData->instanceTransforms[i].size(); j++) {
                        std::string transformedMeshName = meshName + "_" + std::to_string(j);
                        std::vector<Vector<DIM>> transformedPositions;

                        for (int k = 0; k < (int)positions.size(); k++) {
                            transformedPositions.emplace_back(sceneData->instanceTransforms[i][j]*positions[k]);
                        }

                        polyscope::registerCurveNetwork2D(transformedMeshName, transformedPositions, indices);
                    }

                } else {
                    polyscope::registerCurveNetwork2D(meshName, positions, indices);
                }
            }
        }

        // register point clouds
        polyscope::registerPointCloud2D("Intersected_Points", intersectedPoints);
        polyscope::registerPointCloud2D("Raymarched_Points", raymarchedPoints);

    } else if (DIM == 3) {
        for (int i = 0; i < (int)sceneData->soups.size(); i++) {
            std::string meshName = "Polygon_Soup_" + std::to_string(i);

            if (sceneData->soupToObjectsMap[i][0].first == ObjectType::Triangles) {
                // register surface mesh
                int N = (int)sceneData->soups[i].indices.size()/3;
                std::vector<std::vector<int>> indices(N, std::vector<int>(3));
                const std::vector<Vector<DIM>>& positions = sceneData->soups[i].positions;
                for (int j = 0; j < N; j++) {
                    for (int k = 0; k < 3; k++) {
                        indices[j][k] = sceneData->soups[i].indices[3*j + k];
                    }
                }

                if (sceneData->instanceTransforms[i].size() > 0) {
                    for (int j = 0; j < (int)sceneData->instanceTransforms[i].size(); j++) {
                        std::string transformedMeshName = meshName + "_" + std::to_string(j);
                        std::vector<Vector<DIM>> transformedPositions;

                        for (int k = 0; k < (int)positions.size(); k++) {
                            transformedPositions.emplace_back(sceneData->instanceTransforms[i][j]*positions[k]);
                        }

                        polyscope::registerSurfaceMesh(transformedMeshName, transformedPositions, indices);
                    }

                } else {
                    polyscope::registerSurfaceMesh(meshName, positions, indices);
                }
            }
        }

        // register point clouds
        polyscope::registerPointCloud("Intersected_Points", intersectedPoints);
        polyscope::registerPointCloud("Raymarched_Points", raymarchedPoints);
    }

    // register callback
    polyscope::state::userCallback = std::bind(&guiCallback<DIM>, std::cref(sceneData->aggregate),
                                               std::cref(boundingBox),
                                               std::ref(intersectedPoints),
                                               std::ref(raymarchedPoints));

    // give control to polyscope gui
    polyscope::show();
}

template<size_t DIM>
void run()
{
    // build scene and csg
    Scene<DIM> scene;
    SceneLoader<DIM> sceneLoader;
    sceneLoader.loadFiles(scene, true);
    scene.build(AggregateType::Bvh_LongestAxisCenter, false, true);
    SceneData<DIM> *sceneData = scene.getSceneData();

    // generate random points and rays used to visualize csg
    BoundingBox<DIM> boundingBox = sceneData->aggregate->boundingBox();
    std::vector<Vector<DIM>> scatteredPoints, randomDirections;
    generateScatteredPointsAndRays<DIM>(1000, scatteredPoints, randomDirections, boundingBox);

    // intersect and raymarch points
    std::vector<Vector<DIM>> intersectedPoints, raymarchedPoints;
    clampToCsg<DIM>("intersect", scatteredPoints, randomDirections,
                    sceneData->aggregate, boundingBox, intersectedPoints);
    clampToCsg<DIM>("raymarch", scatteredPoints, randomDirections,
                    sceneData->aggregate, boundingBox, raymarchedPoints);

    // visualize scene
    visualizeScene<DIM>(sceneData, boundingBox, intersectedPoints, raymarchedPoints);
}

int main(int argc, const char *argv[]) {
    // configure the argument parser
    args::ArgumentParser parser("csg tests");
    args::Group group(parser, "", args::Group::Validators::DontCare);
    args::ValueFlag<int> dim(parser, "integer", "scene dimension", {"dim"});
    args::ValueFlag<std::string> csgFilename(parser, "string", "csg filename", {"csgFile"});
    args::ValueFlag<std::string> instanceFilename(parser, "string", "instance filename", {"instanceFile"});
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

    if (!csgFilename || (!lineSegmentFilenames && !triangleFilenames)) {
        std::cerr << "Please specify csg and either line segment or triangle soup filenames" << std::endl;
        return EXIT_FAILURE;
    }

    // set global flags
    if (csgFilename) ::csgFilename = args::get(csgFilename);
    if (instanceFilename) ::instanceFilename = args::get(instanceFilename);
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

    if (files.size() <= 1) {
        std::cerr << "Specify atleast 2 soups" << std::endl;
        return EXIT_FAILURE;
    }

    // run app
    if (DIM == 2) run<2>();
    else if (DIM == 3) run<3>();
    return 0;
}