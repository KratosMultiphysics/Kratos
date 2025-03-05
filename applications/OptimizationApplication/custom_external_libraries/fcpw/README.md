<p align="center">
  <img src="https://raw.githubusercontent.com/rohan-sawhney/fcpw/master/logo.png" height="140" width="100">
</p>
<h1 align="center"><em>FCPW: Fastest Closest Points in the West</em></h1>

*FCPW* is an acceleration library for performing fast geometric queries, such as closest points and ray intersections, on 3D triangle meshes and 2D line segment meshes. It is available in C++ and Python, and offers GPU acceleration via the [Slang shading language](https://shader-slang.com/slang/user-guide/) and CPU vectorization via the [Enoki library](https://github.com/mitsuba-renderer/enoki).

# C++ API

The easiest and most direct way to use *FCPW* is through the `Scene` class which provides methods to load geometry, build and refit the acceleration structure and perform geometric queries. Here is an example of performing a closest point query on a 3D triangle mesh:

```c++
#include <fcpw/fcpw.h>
using namespace fcpw;

// initialize a 3d scene
Scene<3> scene;

// load positions and indices of a single triangle mesh
scene.setObjectCount(1);
scene.setObjectVertices(positions, 0);
scene.setObjectTriangles(indices, 0);

// build acceleration structure
AggregateType aggregateType = AggregateType::Bvh_SurfaceArea;
bool buildVectorizedBvh = true;
scene.build(aggregateType, buildVectorizedBvh);

// perform a closest point query
Interaction<3> interaction;
bool found = scene.findClosestPoint(queryPoint, interaction);

// access distance and closest point via interaction.d and interaction.p (resp.)
```

The `Scene` class is templated on dimension, which enables it to work with geometric data in any dimension as long as primitives are defined for the dimension of interest. In addition to triangles in 3D, *FCPW* also provides support for line segments in 2D. The [Interaction](https://github.com/rohan-sawhney/fcpw/blob/master/include/fcpw/core/interaction.h) class stores information relevant to the query, such as the distance to and closest point on the primitive.

If your scene consists of multiple objects with the same type of primitives (e.g. triangles), it is better to "flatten" those objects into a single object (with a list of `positions` and `indices`). *FCPW* then builds a single acceleration structure over all primitives in the scene. If multiple objects are loaded, *FCPW* instead builds a hierarchy of acceleration structures, with a structure for each object in the scene.

Refer to [demo.cpp](https://github.com/rohan-sawhney/fcpw/blob/master/demos/demo.cpp) for a complete demo, and [fcpw.h](https://github.com/rohan-sawhney/fcpw/blob/master/include/fcpw/fcpw.h) for the full API, which includes the list of supported geometric queries, heuristrics for constructing the accerlation structure, as well as support for refitting, instancing and CSG operations. Details regarding GPU support and installation are provided below.

# Python API

*FCPW*'s C++ and Python APIs follow a similar structure. However, since Python does not support templates, Python classes and functions use explicit `_2D` or `_3D` dimension tags:

```python
import fcpw

# initialize a 3d scene
scene = fcpw.scene_3D()

# load positions and indices of a single triangle mesh
scene.set_object_count(1)
scene.set_object_vertices(positions, 0)
scene.set_object_triangles(indices, 0)

# build acceleration structure
aggregate_type = fcpw.aggregate_type.bvh_surface_area
build_vectorized_bvh = True
scene.build(aggregate_type, build_vectorized_bvh)

# initialize bounding spheres
bounding_spheres = fcpw.bounding_sphere_3D_list()
for q in query_points:
	bounding_spheres.append(fcpw.bounding_sphere_3D(q, np.inf))

# perform several closest point queries
interactions = fcpw.interaction_3D_list()
scene.find_closest_points(bounding_spheres, interactions)

# extract closest points
closest_points = [i.p for i in interactions]
```

Refer to [demo.py](https://github.com/rohan-sawhney/fcpw/blob/master/demos/demo.py) for a complete demo. The full API can be viewed in the Python console using `help(fcpw)`.

# GPU Support

GPU support for a large collection of query points is provided through the [GPUScene](https://github.com/rohan-sawhney/fcpw/blob/master/include/fcpw/fcpw_gpu.h) class in C++, and `gpu_scene_*D` classes in Python. *FCPW* currently requires the acceleration structure to be built on the CPU and then transferred to the GPU. For instance, in C++ we have:

```c++
#include <fcpw/fcpw_gpu.h>
using namespace fcpw;

// initialize a 3d scene and load geometry (same as above)
Scene<3> scene;
...

// build acceleration structure on CPU
bool buildVectorizedCPUBvh = false; // NOTE: must build non-vectorized structure
scene.build(AggregateType::Bvh_SurfaceArea, buildVectorizedCPUBvh);

// transfer scene to GPU
GPUScene<3> gpuScene("PATH_TO_FCPW_DIRECTORY");
gpuScene.transferToGPU(scene);

// initialize bounding spheres 
std::vector<GPUBoundingSphere> boundingSpheres;
for (auto q: queryPoints) {
	float3 queryPoint = float3{q[0], q[1], q[2]};
	boundingSpheres.emplace_back(GPUBoundingSphere(queryPoint, INFINITY));
}

// perform several closest point queries on GPU
std::vector<GPUInteraction> interactions;
gpuScene.findClosestPoints(boundingSpheres, interactions);
```

and in Python:

```python
import fcpw

# initialize a 3d scene and load geometry (same as above)
scene = fcpw.scene_3D()
...

# build acceleration structure on CPU
build_vectorized_cpu_bvh = False # NOTE: must build non-vectorized structure
scene.build(fcpw.aggregate_type.bvh_surface_area, build_vectorized_cpu_bvh)

# transfer scene to GPU
gpu_scene = fcpw.gpu_scene_3D("PATH_TO_FCPW_DIRECTORY")
gpu_scene.transfer_to_gpu(scene)

# initialize bounding spheres
bounding_spheres = fcpw.gpu_bounding_sphere_list()
for q in query_points:
	gpu_query_point = fcpw.float_3D(q[0], q[1], q[2])
	bounding_spheres.append(fcpw.gpu_bounding_sphere(gpu_query_point, np.inf))

# perform several closest point queries on GPU
interactions = fcpw.gpu_interaction_list()
gpu_scene.find_closest_points(bounding_spheres, interactions)
```

Refer to [demo.cpp](https://github.com/rohan-sawhney/fcpw/blob/master/demos/demo.cpp) and [demo.py](https://github.com/rohan-sawhney/fcpw/blob/master/demos/demo.py) for complete demos. GPU support is available on Linux and Windows ([Slang](https://shader-slang.com/slang/user-guide/) currently only has unofficial support for macOS). Installation details for Slang are provided below.

For those developing directly on the GPU, the acceleration structure can also be accessed through [bvh.slang](https://github.com/rohan-sawhney/fcpw/blob/master/include/fcpw/gpu/bvh.slang). Slang offers CUDA, HLSL, Vulkan, OpenGL and Metal (experimental) as compilation targets. Refer to its [user guide](https://shader-slang.com/slang/user-guide/get-started.html) for further details on compiling to these targets.

# C++ Installation

*FCPW* is developed as a header-only library. It can be downloaded and compiled using CMake as follows:

```
git clone https://github.com/rohan-sawhney/fcpw.git
cd fcpw && git submodule update --init --recursive
[...clone additional dependencies...] // instructions below
mkdir build && cd build
cmake [-DFCPW_BUILD_DEMO=ON] [-DFCPW_USE_ENOKI=OFF] [-DFCPW_ENABLE_GPU_SUPPORT=ON] ..
make -j8 // for Linux and Mac; cmake generates a visual studio project for Windows.
```

The `FCPW_BUILD_DEMO` option requires [polyscope](https://polyscope.run) as a dependency in the `deps` repo. Clone polyscope using:

```
git clone --recurse-submodules https://github.com/nmwsharp/polyscope.git deps/polyscope
```

The C++ demo can be run from the `build` directory with the command:
```
./demos/demo [--useGpu]
```

For CPU vectorization, Enoki is included by default as a submodule. It can be disabled with the command `-DFCPW_USE_ENOKI=OFF`, in which case *FCPW* falls back to [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) for non-vectorized CPU queries.

To include *FCPW* in your project without GPU support, add the following lines to your CMakeLists.txt file:

```
add_subdirectory(fcpw)
target_link_libraries(YOUR_TARGET fcpw)
target_include_directories(YOUR_TARGET PRIVATE ${FCPW_EIGEN_INCLUDES})
target_include_directories(YOUR_TARGET PRIVATE ${FCPW_ENOKI_INCLUDES})
```

If your prefer to directly include *FCPW* header files without CMake, then you'll have to define a few extra variables before including the library:

```
#define FCPW_USE_ENOKI
#define FCPW_SIMD_WIDTH 4 // change based on the SIMD width supported by your machine
#include <fcpw/fcpw.h>
```

Finally, for GPU support the `FCPW_ENABLE_GPU_SUPPORT` option assumes Slang binaries are available under `deps/slang`. Copy-paste the contents of the [latest release](https://github.com/shader-slang/slang/releases) under this directory for the relevant platform. Additionally include the following lines in your CMakeLists.txt file:

```
target_link_libraries(YOUR_TARGET ${FCPW_SLANG_LIBRARY})
target_link_libraries(YOUR_TARGET ${FCPW_SLANG_GLSLANG_LIBRARY})
target_link_libraries(YOUR_TARGET ${FCPW_GFX_LIBRARY})
target_include_directories(YOUR_TARGET PRIVATE ${FCPW_SLANG_INCLUDES})
```

On Windows, you may need to download necessary DLL files from the official DirectX Shader Compiler repository on GitHub [here](https://github.com/microsoft/DirectXShaderCompiler/releases). Copy `dxil.dll` and `dxcompiler.dll` to `C:\Windows\System32\` for 64-bit systems or to `C:\Windows\SysWOW64\` for 32-bit systems. This makes the DLLs available to all applications system-wide.

# Python Installation

Pre-built Python wheels are available for Linux, Windows and macOS on [PyPI](https://pypi.org/project/fcpw/):

```
pip install fcpw
```

These wheels can also be downloaded and installed from [Releases](https://github.com/rohan-sawhney/fcpw/releases):

```
pip install fcpw-*.whl
```

Alternatively, to build Python bindings on your local machine, first clone [nanobind](https://nanobind.readthedocs.io/en/latest/) using:

```
git clone --recurse-submodules https://github.com/wjakob/nanobind.git deps/nanobind
```

To build and install the bindings, run:

```
pip install . [--config-settings=cmake.define.FCPW_ENABLE_GPU_SUPPORT=ON]
```

For GPU support, you will need to follow the instructions in the previous section on fetching and placing Slang binaries under `deps/slang`. Finally, launch the Python demo from the `demos` folder using:

```
python -m pip install polyscope
python demo.py [--use_gpu]
```

# Citation
```
@software{FCPW,
author = {Sawhney, Rohan},
title = {FCPW: Fastest Closest Points in the West},
version = {1.0},
year = {2021}
}
```

# Author
[Rohan Sawhney](http://www.rohansawhney.io)

# License

Released under the [MIT License](https://opensource.org/licenses/MIT).
