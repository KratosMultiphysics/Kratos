//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/formats/mff.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

Mesh read_mff(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);

    std::vector<std::string> toks;
    std::string t;
    while (in >> t) {
        for (char& c : t)
            if (c == 'D' || c == 'd')
                c = 'E';
        toks.push_back(t);
    }

    Mesh mesh;
    if (toks.empty()) {
        mesh.AssignPoints(NDArray(DType::Float64, {0, 0}));
        return mesh;
    }
    std::size_t count = static_cast<std::size_t>(std::strtoll(toks[0].c_str(), nullptr, 10));
    if (count + 1 > toks.size())
        count = toks.size() - 1;
    NDArray values(DType::Float64, {count});
    for (std::size_t i = 0; i < count; ++i)
        values.As<double>()[i] = std::strtod(toks[i + 1].c_str(), nullptr);
    mesh.AssignPoints(NDArray(DType::Float64, {count, 0}));
    mesh.AddPointData("mff:field", std::move(values));
    return mesh;
}

void write_mff(const std::string& rPath, const Mesh& rMesh) {
    std::ofstream f(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);

    // pick the first field: first point_data, else first non-unv:pid cell_data
    std::vector<double> values;
    auto point_names = rMesh.PointDataNames();
    if (!point_names.empty()) {
        const NDArray& arr = rMesh.PointData(point_names.front());
        values.resize(arr.Size());
        for (std::size_t i = 0; i < arr.Size(); ++i)
            values[i] = detail::read_double(arr, i);
    } else {
        for (const auto& name : rMesh.CellDataNames()) {
            if (name == "unv:pid")
                continue;
            for (std::size_t bi = 0; bi < rMesh.CellDataNumBlocks(name); ++bi) {
                const NDArray& blk = rMesh.CellData(name, bi);
                for (std::size_t i = 0; i < blk.Size(); ++i)
                    values.push_back(detail::read_double(blk, i));
            }
            break;
        }
    }

    f << values.size() << "\n";
    char buf[64];
    for (double v : values) {
        std::snprintf(buf, sizeof(buf), "%.16E\n", v);
        f << buf;
    }
}

}  // namespace meshioplusplus
