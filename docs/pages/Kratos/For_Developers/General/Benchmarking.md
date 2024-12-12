---
title: Benchmarking in Kratos Multiphysics  
keywords: Benchmarking, Google Benchmark, Performance Measurement  
tags: [Benchmarking]  
sidebar: kratos_for_developers
summary:  Step-by-step guide on adding and executing benchmarks in Kratos Multiphysics using Google Benchmark.
---

## Introduction  

Benchmarking is now supported in **Kratos Multiphysics**, enabling developers to measure and analyze the performance of their code effectively. The system integrates seamlessly with [Google Benchmark](https://github.com/google/benchmark), which is included as an external dependency, similar to the **GTest** suite.

### Key Features:
- Automatic generation of a binary per benchmark source file located in the `benchmark` directory.
- Benchmark executables are compiled to the binary `benchmark` directory.
- Activated using the `KRATOS_BUILD_BENCHMARK` CMake option (disabled by default).

### What is Google Benchmark?

[Google Benchmark](https://github.com/google/benchmark) is a library developed by *Google* for microbenchmarking C++ code. It provides:

- Precision and reproducibility in measurements.
- Built-in support for statistical analysis of timing data.
- Customizable inputs to benchmarks for more flexible scenarios.
- Support for benchmarking in multi-threaded and single-threaded environments.

## Writing Benchmarks
​
### Overview of Google Benchmark API
​
A *Google Benchmark* consists of a **function** that performs repeated measurements of a specific code snippet. The framework automatically manages:
​
- Warm-up runs to stabilize performance.
- Statistical aggregation of results (mean, median, standard deviation, etc.).
- Handling of scaling issues with large inputs.
​
#### Essential Functions:
​
- `benchmark::State`: Represents the state of the benchmark, used to control iterations.
- `BENCHMARK(...)`: Registers a function as a benchmark.
- `BENCHMARK_MAIN()`: Main entry point for running benchmarks.
​
#### Benchmark Configuration Options:
​
1. **Repeats**: Specify the number of iterations with `--benchmark_repetitions=N`.
2. **Time Limit**: Set the time for each benchmark with `--benchmark_min_time=N`.
3. **Custom Arguments**: Pass custom inputs to benchmarks for parameterized testing.

### Example Benchmark (C++ Code)

Below is an example showcasing how to benchmark the nearest point calculation on a triangle in *Kratos*:

~~~cpp
// System includes

// External includes
#include <benchmark/benchmark.h>

// Project includes
#include "geometries/point.h"
#include "geometries/triangle_3d_3.h"
#include "utilities/geometry_utilities/nearest_point_utilities.h"

namespace Kratos
{

// Sample data for benchmarking
Point::Pointer p_point_1(make_shared<Point>( 0.0, 0.0, 0.0));
Point::Pointer p_point_2(make_shared<Point>( 1.0, 0.0, 0.0));
Point::Pointer p_point_3(make_shared<Point>( 0.0, 1.0, 0.0));

Triangle3D3<Point> triangle(p_point_1, p_point_3, p_point_2);
Point nearest_point(0.0, 0.0, 0.0);
Point inside_point(0.2, 0.1, 0.00);

static void BM_TriangleNearestPoint(benchmark::State& state) {
    for (auto _ : state) {
        NearestPointUtilities::TriangleNearestPoint(inside_point, triangle, nearest_point);
    }
}

// Register the function as a benchmark
BENCHMARK(BM_TriangleNearestPoint);

}  // namespace Kratos

BENCHMARK_MAIN();
~~~

### Advanced Features of Google Benchmark

#### Custom Arguments

Benchmarks can take custom arguments to test specific cases. Example:

~~~cpp
static void BM_CustomArguments(benchmark::State& state) {
    for (auto _ : state) {
        int input_size = state.range(0);
        // Your code here, e.g., filling arrays of size input_size
    }
}

BENCHMARK(BM_CustomArguments)->Arg(100)->Arg(1000)->Arg(10000);
~~~

#### Multithreading Support

Google Benchmark can run benchmarks in multi-threaded mode:

~~~cpp
static void BM_MultiThreaded(benchmark::State& state) {
    for (auto _ : state) {
        // Thread-safe code
    }
}

BENCHMARK(BM_MultiThreaded)->Threads(1)->Threads(2)->Threads(4);
~~~

## Running Benchmarks

### Execution

Once benchmarks are compiled, executables are created in the `benchmark` directory. Run them directly via the command line:

```
./benchmark_name
```

### Configuration Options

- **JSON Output**: Export results in JSON format:
  
```
./benchmark_name --benchmark_format=json --benchmark_out=output.json
```

- **Statistical Runs**: Run benchmarks multiple times for statistical accuracy:
  
```
./benchmark_name --benchmark_repetitions=10
```

- **Custom Time Limits**: Specify time per benchmark:
  
```
./benchmark_name --benchmark_min_time=2.0
```

## Comparing Benchmark Results  

The benchmark framework includes a Python utility script to compare results between multiple JSON outputs. This script, located in `kratos/python_scripts/benchmark/generate_plot_google_benchmark.py`, generates visualizations for easier comparison.

### Example Usage:  
To compare two JSON results:  
```
python generate_plot_google_benchmark.py --filenames result1.json result2.json
```

### Example Output:  
A plot similar to the one below will be generated, showing the performance (in terms of CPU time and real time) of different benchmarks across the compared files.  

![Benchmark Performance Graph](images/benchmark.png)

## Python Script Explanation  

### Purpose:
The provided Python script consolidates and visualizes benchmarking data from multiple JSON files generated by Google Benchmark.  

#### Script Functionality:  
1. **Load JSON Outputs:** Extract benchmark results from each JSON file.  
2. **Merge Data:** Combine results into a single pandas DataFrame for analysis.  
3. **Generate Plot:** Create a grouped bar chart showing `cpu_time` and `real_time` for each benchmark.  

## About the Google Benchmark Library

### Benefits

- Accurate Timing: Adjusts for background noise and other factors.
- Reproducibility: Consistent and repeatable measurements.
- Customizability: Configure inputs, threads, and scaling scenarios.

### Useful Links

- [Google Benchmark Documentation](https://github.com/google/benchmark)
- [Google Benchmark Usage Guide](https://github.com/google/benchmark/blob/main/docs/user_guide.md)