#pragma once

#include <vector>

namespace ANurbs {

template <typename TVector>
struct Polygon
{
    using Path = std::vector<TVector>;

    Path outer_path;
    std::vector<Path> inner_paths;

    int
    NbLoops() const
    {
        return inner_paths.size() + 1;
    }

    int
    NbVerticesOfLoop(
        const int index) const
    {
        if (index == 0) {
            return static_cast<int>(outer_path.size());
        } else {
            return static_cast<int>(inner_paths[index - 1].size());
        }
    }

    TVector
    VertexOfLoop(
        int loopIndex,
        int vertexIndex) const
    {
        if (loopIndex == 0) {
            return outer_path[vertexIndex];
        } else {
            return inner_paths[loopIndex - 1][vertexIndex];
        }
    }

    int
    NbVertices() const
    {
        int nbVertices = outer_path.size();

        for (const auto& path : inner_paths) {
            nbVertices += static_cast<int>(path.size());
        }

        return nbVertices;
    }

    TVector
    Vertex(
        int index) const
    {
        if (index < outer_path.size()) {
            return outer_path[index];
        }

        std::size_t offset = outer_path.size();

        for (const auto& path : inner_paths) {
            if (index < offset + path.size()) {
                return path[index - offset];
            }

            offset += path.size();
        }

        throw std::exception();
    }
};

} // namespace ANurbs
