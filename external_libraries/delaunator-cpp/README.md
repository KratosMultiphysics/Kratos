# delaunator-cpp

A really fast C++ library for
[Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation) of 2D points.

delaunator-cpp is a C++ port from https://github.com/mapbox/delaunator a JavaScript implementation of very fast 2D Delaunay algorithm.

[![Build Status](https://travis-ci.org/delfrrr/delaunator-cpp.svg?branch=master)](https://travis-ci.org/delfrrr/delaunator-cpp)
[![badge](https://mapbox.s3.amazonaws.com/cpp-assets/hpp-skel-badge_blue.svg)](https://github.com/mapbox/hpp-skel)

## Features

* Probably the fastest C++ open source 2D Delaunay implementation
* Example showing triangulation of GeoJson points

## Usage

`examples/basic.cpp`

```CPP
#include <delaunator.hpp>
#include <cstdio>

int main() {
    /* x0, y0, x1, y1, ... */
    std::vector<double> coords = {-1, 1, 1, 1, 1, -1, -1, -1};

    //triangulation happens here
    delaunator::Delaunator d(coords);

    for(std::size_t i = 0; i < d.triangles.size(); i+=3) {
        printf(
            "Triangle points: [[%f, %f], [%f, %f], [%f, %f]]\n",
            d.coords[2 * d.triangles[i]],        //tx0
            d.coords[2 * d.triangles[i] + 1],    //ty0
            d.coords[2 * d.triangles[i + 1]],    //tx1
            d.coords[2 * d.triangles[i + 1] + 1],//ty1
            d.coords[2 * d.triangles[i + 2]],    //tx2
            d.coords[2 * d.triangles[i + 2] + 1] //ty2
        );
    }
}
```

[See more examples here](./examples)

## Benchmarks

```
Run on (4 X 2300 MHz CPU s)
2018-09-29 09:27:28
------------------------------------------------------------
Benchmark                     Time           CPU Iterations
------------------------------------------------------------
BM_45K_geojson_nodes         22 ms         22 ms         32
BM_uniform/2000               1 ms          1 ms        982
BM_uniform/100000            63 ms         62 ms          9
BM_uniform/200000           140 ms        140 ms          4
BM_uniform/500000           400 ms        399 ms          2
BM_uniform/1000000          994 ms        993 ms          1
```

Library is ~10% faster then JS version for 1M uniform points ([details](https://github.com/delfrrr/delaunator-cpp/pull/8#issuecomment-422690056))
