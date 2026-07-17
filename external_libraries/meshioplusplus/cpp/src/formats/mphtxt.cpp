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
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/formats/mphtxt.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

std::string comsol_to_meshio(const std::string& rT) {
    static const std::unordered_map<std::string, std::string> m = {
        {"vtx", "vertex"},     {"edg", "line"},         {"tri", "triangle"}, {"quad", "quad"},
        {"tet", "tetra"},      {"prism", "wedge"},      {"pyr", "pyramid"},  {"hex", "hexahedron"},
        {"edg2", "line3"},     {"tri2", "triangle6"},   {"quad2", "quad9"},  {"tet2", "tetra10"},
        {"prism2", "wedge18"}, {"hex2", "hexahedron27"}};
    auto it = m.find(rT);
    return it == m.end() ? std::string() : it->second;
}

std::string meshio_to_comsol(const std::string& rT) {
    static const std::unordered_map<std::string, std::string> m = {
        {"vertex", "vtx"},     {"line", "edg"},         {"triangle", "tri"}, {"quad", "quad"},
        {"tetra", "tet"},      {"wedge", "prism"},      {"pyramid", "pyr"},  {"hexahedron", "hex"},
        {"line3", "edg2"},     {"triangle6", "tri2"},   {"quad9", "quad2"},  {"tetra10", "tet2"},
        {"wedge18", "prism2"}, {"hexahedron27", "hex2"}};
    auto it = m.find(rT);
    return it == m.end() ? std::string() : it->second;
}

const std::vector<int>* perm_of(const std::string& rT) {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"quad", {0, 1, 3, 2}}, {"hexahedron", {0, 1, 3, 2, 4, 5, 7, 6}}};
    auto it = m.find(rT);
    return it == m.end() ? nullptr : &it->second;
}

struct Cursor {
    std::vector<std::string> mT;
    std::size_t mI = 0;
    const std::string& Tok() {
        if (mI >= mT.size())
            throw ReadError("mphtxt: unexpected end of file");
        return mT[mI++];
    }
    long long Integer() { return std::strtoll(Tok().c_str(), nullptr, 10); }
    double Real() { return std::strtod(Tok().c_str(), nullptr); }
    std::string Str() {
        Integer();  // length prefix
        return Tok();
    }
};

}  // namespace

Mesh read_mphtxt(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    Cursor c;
    std::string line;
    while (std::getline(in, line)) {
        std::size_t h = line.find('#');
        if (h != std::string::npos)
            line = line.substr(0, h);
        std::istringstream iss(line);
        std::string w;
        while (iss >> w)
            c.mT.push_back(w);
    }

    c.Integer();  // version major
    c.Integer();  // version minor
    for (long long k = c.Integer(); k > 0; --k)
        c.Str();  // tags
    const long long n_types = c.Integer();
    for (long long k = 0; k < n_types; ++k)
        c.Str();  // type names

    Mesh mesh;
    std::vector<NDArray> geom;

    for (long long obj = 0; obj < n_types; ++obj) {
        c.Integer();
        c.Integer();
        c.Integer();  // object type indices
        c.Str();      // class name
        c.Integer();  // object version
        const long long sdim = c.Integer();
        const long long n_points = c.Integer();
        const long long lowest = c.Integer();
        NDArray pts(DType::Float64,
                    {static_cast<std::size_t>(n_points), static_cast<std::size_t>(sdim)});
        for (long long p = 0; p < n_points * sdim; ++p)
            pts.As<double>()[p] = c.Real();
        mesh.AssignPoints(std::move(pts));

        const long long n_eltypes = c.Integer();
        for (long long e = 0; e < n_eltypes; ++e) {
            std::string ctype = c.Str();
            std::string mtype = comsol_to_meshio(ctype);
            if (mtype.empty())
                throw ReadError("mphtxt: unknown element type " + ctype);
            const long long nn = c.Integer();
            const long long ne = c.Integer();
            NDArray conn(DType::Int64,
                         {static_cast<std::size_t>(ne), static_cast<std::size_t>(nn)});
            std::vector<std::int64_t> raw(ne * nn);
            for (long long v = 0; v < ne * nn; ++v)
                raw[v] = c.Integer() - lowest;
            const std::vector<int>* p = perm_of(mtype);
            for (long long r = 0; r < ne; ++r)
                for (long long j = 0; j < nn; ++j)
                    conn.As<std::int64_t>()[r * nn + j] =
                        p ? raw[r * nn + (*p)[j]] : raw[r * nn + j];

            const long long npar_per = c.Integer();
            const long long npar = c.Integer();
            for (long long v = 0; v < npar * npar_per; ++v)
                c.Tok();
            const long long ngeom = c.Integer();
            NDArray g(DType::Int64, {static_cast<std::size_t>(ngeom)});
            for (long long v = 0; v < ngeom; ++v)
                g.As<std::int64_t>()[v] = c.Integer();
            const long long nud = c.Integer();
            for (long long v = 0; v < nud * 2; ++v)
                c.Integer();

            mesh.AddCellBlock(mtype, std::move(conn));
            geom.push_back(std::move(g));
        }
        break;  // first mesh object only
    }

    if (!geom.empty())
        mesh.AddCellData("mphtxt:geom", std::move(geom));
    return mesh;
}

void write_mphtxt(const std::string& rPath, const Mesh& rMesh) {
    std::ofstream f(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t sdim = rMesh.PointDim();

    struct Blk {
        std::size_t mIdx;
        Mesh::CellView mCb;
    };
    std::vector<Blk> blocks;
    for (std::size_t k = 0; k < rMesh.NumCellBlocks(); ++k) {
        const auto cb = rMesh.Cells(k);
        if (!meshio_to_comsol(cb.Type()).empty())
            blocks.push_back({k, cb});
        else
            throw WriteError("mphtxt: unsupported cell type " + cb.Type());
    }

    const bool has_geom = rMesh.HasCellData("mphtxt:geom");

    f << "# Created by meshio++ (C++ core)\n\n";
    f << "0 1\n";
    f << "1 # number of tags\n5 mesh1\n";
    f << "1 # number of types\n3 obj\n\n";
    f << "0 0 1\n4 Mesh # class\n2 # version\n";
    f << sdim << " # sdim\n";
    f << rMesh.NumPoints() << " # number of mesh points\n";
    f << "1 # lowest mesh point index\n\n# Mesh point coordinates\n";
    const NDArray& points = rMesh.Points();
    char buf[32];
    for (std::size_t i = 0; i < rMesh.NumPoints(); ++i) {
        for (std::size_t cc = 0; cc < sdim; ++cc) {
            std::snprintf(buf, sizeof(buf), "%.16g", detail::read_double(points, i * sdim + cc));
            f << buf << (cc + 1 == sdim ? '\n' : ' ');
        }
    }
    f << "\n" << blocks.size() << " # number of element types\n\n";

    int ti = 0;
    for (const auto& b : blocks) {
        const auto cb = b.mCb;
        std::string ctype = meshio_to_comsol(cb.Type());
        const std::vector<int>* p = perm_of(cb.Type());
        const NDArray& conn = cb.Conn();
        std::size_t nn = detail::cols(conn);
        std::size_t ne = cb.NumCells();
        f << "# Type #" << (++ti) << "\n\n";
        f << ctype.size() << " " << ctype << " # type name\n\n";
        f << nn << " # number of nodes per element\n";
        f << ne << " # number of elements\n# Elements\n";
        for (std::size_t r = 0; r < ne; ++r) {
            for (std::size_t j = 0; j < nn; ++j) {
                std::size_t src = p ? (*p)[j] : j;
                f << (detail::read_int(conn, r * nn + src) + 1) << (j + 1 == nn ? '\n' : ' ');
            }
        }
        f << "\n" << nn << " # number of parameter values per element\n";
        f << "0 # number of parameters\n# Parameters\n\n";
        f << ne << " # number of geometric entity indices\n# Geometric entity indices\n";
        const NDArray* g = (has_geom && b.mIdx < rMesh.CellDataNumBlocks("mphtxt:geom"))
                               ? &rMesh.CellData("mphtxt:geom", b.mIdx)
                               : nullptr;
        for (std::size_t r = 0; r < ne; ++r)
            f << (g ? detail::read_int(*g, r) : 0) << "\n";
        f << "\n0 # number of up/down pairs\n# Up/down\n\n";
    }
}

}  // namespace meshioplusplus
