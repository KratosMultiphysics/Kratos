/*
    This demo demonstrates how to perform closest point queries using FCPW.
    Refer to the fcpw.h header file for the full API, including how to perform
    distance and ray intersection queries on a collection of query points as
    well as a single query point.

    The demo can be run from the command line using the following commands:
    > mkdir build
    > cd build
    > cmake -DFCPW_BUILD_DEMO=ON [-DFCPW_ENABLE_GPU_SUPPORT=ON] ..
    > make -j4
    > ./demos/demo [--useGpu]
*/

#include <fcpw/fcpw.h>
#include <fstream>
#include <sstream>
#include <filesystem>
#ifdef FCPW_USE_GPU
    #include <fcpw/fcpw_gpu.h>
#endif
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

#include "args/args.hxx"

using namespace fcpw;

// helper class for loading polygons from obj files
struct Index {
    Index() {}
    Index(int v) : position(v) {}

    bool operator<(const Index& i) const {
        if (position < i.position) return true;
        if (position > i.position) return false;

        return false;
    }

    int position;
};

// helper function for loading polygons from obj files
inline Index parseFaceIndex(const std::string& token) {
    std::stringstream in(token);
    std::string indexString;
    int indices[3] = {1, 1, 1};

    int i = 0;
    while (std::getline(in, indexString, '/')) {
        if (indexString != "\\") {
            std::stringstream ss(indexString);
            ss >> indices[i++];
        }
    }

    // decrement since indices in OBJ files are 1-based
    return Index(indices[0] - 1);
}

void loadObj(const std::string& objFilePath,
             std::vector<Vector<3>>& positions,
             std::vector<Vector3i>& indices)
{
    // initialize
    std::ifstream in(objFilePath);
    if (in.is_open() == false) {
        std::cerr << "Unable to open file: " << objFilePath << std::endl;
        exit(EXIT_FAILURE);
    }

    // parse
    std::string line;
    positions.clear();
    indices.clear();

    while (getline(in, line)) {
        std::stringstream ss(line);
        std::string token;
        ss >> token;

        if (token == "v") {
            float x, y, z;
            ss >> x >> y >> z;

            positions.emplace_back(Vector3(x, y, z));
        
        } else if (token == "f") {
            int j = 0;
            Vector3i face;

            while (ss >> token) {
                Index index = parseFaceIndex(token);

                if (index.position < 0) {
                    getline(in, line);
                    size_t i = line.find_first_not_of("\t\n\v\f\r ");
                    index = parseFaceIndex(line.substr(i));
                }

                face[j++] = index.position;
            }

            indices.emplace_back(face);
        }
    }

    // close
    in.close();
}

void loadFcpwScene(const std::vector<Vector<3>>& positions,
                   const std::vector<Vector3i>& indices,
                   bool buildVectorizedBvh, Scene<3>& scene)
{
    // load positions and indices
    scene.setObjectCount(1);
    scene.setObjectVertices(positions, 0);
    scene.setObjectTriangles(indices, 0);

    // build scene on CPU
    AggregateType aggregateType = AggregateType::Bvh_SurfaceArea;
    bool printStats = false;
    bool reduceMemoryFootprint = false;
    scene.build(aggregateType, buildVectorizedBvh, printStats, reduceMemoryFootprint);
}

template <typename T>
void performClosestPointQueries(const std::vector<Vector<3>>& queryPoints,
                                std::vector<Vector<3>>& closestPoints,
                                T& scene)
{
    // do nothing
}

template <>
void performClosestPointQueries(const std::vector<Vector<3>>& queryPoints,
                                std::vector<Vector<3>>& closestPoints,
                                Scene<3>& scene)
{
    // initialize bounding spheres
    std::vector<BoundingSphere<3>> boundingSpheres;
    for (const Vector<3>& q: queryPoints) {
        boundingSpheres.emplace_back(BoundingSphere<3>(q, maxFloat));
    }

    // perform cpqs
    std::vector<Interaction<3>> interactions;
    scene.findClosestPoints(boundingSpheres, interactions);

    // extract closest points
    closestPoints.clear();
    for (const Interaction<3>& i: interactions) {
        closestPoints.emplace_back(i.p);
    }
}

#ifdef FCPW_USE_GPU

template <>
void performClosestPointQueries(const std::vector<Vector<3>>& queryPoints,
                                std::vector<Vector<3>>& closestPoints,
                                GPUScene<3>& gpuScene)
{
    // initialize bounding spheres
    std::vector<GPUBoundingSphere> boundingSpheres;
    for (const Vector<3>& q: queryPoints) {
        float3 queryPoint = float3{q[0], q[1], q[2]};
        boundingSpheres.emplace_back(GPUBoundingSphere(queryPoint, maxFloat));
    }

    // perform cpqs on GPU
    std::vector<GPUInteraction> interactions;
    gpuScene.findClosestPoints(boundingSpheres, interactions);

    // extract closest points
    closestPoints.clear();
    for (const GPUInteraction& i: interactions) {
        closestPoints.emplace_back(Vector<3>(i.p.x, i.p.y, i.p.z));
    }
}

#endif

template <typename T>
void guiCallback(std::vector<Vector<3>>& queryPoints,
                 std::vector<Vector<3>>& closestPoints,
                 T& scene)
{
    // animate query points
    for (Vector<3>& q: queryPoints) {
        q[0] += 0.001 * std::sin(10.0 * q[1]);
        q[1] += 0.001 * std::cos(10.0 * q[0]);
    }

    // perform closest point queries
    performClosestPointQueries(queryPoints, closestPoints, scene);

    // plot results
    polyscope::registerPointCloud("query points", queryPoints);
    polyscope::registerPointCloud("closest points", closestPoints);

    std::vector<Vector2i> edgeIndices;
    std::vector<Vector<3>> edgePositions = queryPoints;
    edgePositions.insert(edgePositions.end(), closestPoints.begin(), closestPoints.end());
    for (int i = 0; i < (int)queryPoints.size(); i++) {
        edgeIndices.emplace_back(Vector2i(i, i + queryPoints.size()));
    }

    auto network = polyscope::registerCurveNetwork("edges", edgePositions, edgeIndices);
    network->setRadius(0.005, false);
}

template <typename T>
void visualize(const std::vector<Vector<3>>& positions,
               const std::vector<Vector3i>& indices,
               std::vector<Vector<3>>& queryPoints,
               std::vector<Vector<3>>& closestPoints,
               T& scene)
{
    // set a few options
    polyscope::options::programName = "FCPW Demo";
    polyscope::options::verbosity = 0;
    polyscope::options::usePrefsFile = false;
    polyscope::options::autocenterStructures = false;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;

    // initialize polyscope
    polyscope::init();

    // register mesh and callback
    polyscope::registerSurfaceMesh("mesh", positions, indices);
    polyscope::state::userCallback = std::bind(&guiCallback<T>, std::ref(queryPoints),
                                               std::ref(closestPoints), std::ref(scene));

    // give control to polyscope gui
    polyscope::show();
}

void run(bool useGpu)
{
    // load obj file
    std::filesystem::path currentDirectory = std::filesystem::current_path().parent_path();
    std::string objFilePath = (currentDirectory / "demos" / "dragon.obj").string();
    std::vector<Vector<3>> positions;
    std::vector<Vector3i> indices;
    loadObj(objFilePath, positions, indices);

    // generate random query points for closest point queries
    Vector<3> boxMin = Vector<3>::Constant(std::numeric_limits<float>::infinity());
    Vector<3> boxMax = -Vector<3>::Constant(std::numeric_limits<float>::infinity());
    for (const Vector<3>& p: positions) {
        boxMin = boxMin.cwiseMin(p);
        boxMax = boxMax.cwiseMax(p);
    }

    int numQueryPoints = 100;
    std::vector<Vector<3>> queryPoints;
    Vector<3> boxExtent = boxMax - boxMin;
    for (int i = 0; i < numQueryPoints; i++) {
        queryPoints.emplace_back(boxMin + boxExtent.cwiseProduct(uniformRealRandomVector<3>()));
    }

    if (useGpu) {
#ifdef FCPW_USE_GPU
        // load fcpw scene on CPU
        Scene<3> scene;
        loadFcpwScene(positions, indices, false, scene); // NOTE: must build non-vectorized CPU BVH

        // transfer scene to GPU
        bool printStats = false;
        GPUScene<3> gpuScene(currentDirectory.string(), printStats);
        gpuScene.transferToGPU(scene);

        // visualize results
        std::vector<Vector<3>> closestPoints;
        visualize(positions, indices, queryPoints, closestPoints, gpuScene);
#else
        std::cerr << "GPU support not enabled" << std::endl;
        exit(EXIT_FAILURE);
#endif

    } else {
        // load fcpw scene
        Scene<3> scene;
        loadFcpwScene(positions, indices, true, scene);

        // visualize results
        std::vector<Vector<3>> closestPoints;
        visualize(positions, indices, queryPoints, closestPoints, scene);
    }
}

int main(int argc, const char *argv[]) {
    // configure the argument parser
    args::ArgumentParser parser("fcpw demo");
    args::Group group(parser, "", args::Group::Validators::DontCare);
    args::Flag useGpu(group, "bool", "use GPU", {"useGpu"});

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

    run(args::get(useGpu));

    return 0;
}