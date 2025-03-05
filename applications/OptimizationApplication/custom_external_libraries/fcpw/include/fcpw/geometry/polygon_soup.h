#pragma once

#include <fcpw/core/primitive.h>

namespace fcpw {

template<size_t DIM>
struct PolygonSoup {
    // constructor
    PolygonSoup() {}

    // constructor
    PolygonSoup(const std::vector<int>& indices_,
                const std::vector<Vector<DIM>>& positions_):
                indices(indices_), positions(positions_) {}

    // members
    std::vector<int> indices /* a.k.a. vIndices */, eIndices, tIndices;
    std::vector<int> indexMap;
    std::vector<Vector<DIM>> positions;
    std::vector<Vector<DIM - 1>> textureCoordinates;
    std::vector<Vector<DIM>> vNormals, eNormals; // normalized values
};

} // namespace fcpw