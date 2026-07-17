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
#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/types.hpp"

namespace meshioplusplus {

namespace {

// ---- type maps (subset; ported from gmsh/common.py) --------------------------
const std::unordered_map<int, std::string>& gmsh_to_meshio_type() {
    static const std::unordered_map<int, std::string> m = {
        {1, "line"},          {2, "triangle"},      {3, "quad"},           {4, "tetra"},
        {5, "hexahedron"},    {6, "wedge"},         {7, "pyramid"},        {8, "line3"},
        {9, "triangle6"},     {10, "quad9"},        {11, "tetra10"},       {12, "hexahedron27"},
        {13, "wedge18"},      {14, "pyramid14"},    {15, "vertex"},        {16, "quad8"},
        {17, "hexahedron20"}, {18, "wedge15"},      {19, "pyramid13"},     {21, "triangle10"},
        {23, "triangle15"},   {25, "triangle21"},   {26, "line4"},         {27, "line5"},
        {28, "line6"},        {29, "tetra20"},      {30, "tetra35"},       {31, "tetra56"},
        {36, "quad16"},       {37, "quad25"},       {38, "quad36"},        {62, "line7"},
        {63, "line8"},        {64, "line9"},        {65, "line10"},        {66, "line11"},
        {71, "tetra84"},      {72, "tetra120"},     {73, "tetra165"},      {74, "tetra220"},
        {75, "tetra286"},     {92, "hexahedron64"}, {93, "hexahedron125"},
    };
    return m;
}

const std::unordered_map<std::string, int>& meshio_to_gmsh_type() {
    static const std::unordered_map<std::string, int> m = [] {
        std::unordered_map<std::string, int> r;
        for (const auto& kv : gmsh_to_meshio_type())
            r[kv.second] = kv.first;
        return r;
    }();
    return m;
}

// Permutation P such that meshio_row[j] = gmsh_row[P[j]]; empty = identity.
const std::vector<int>& gmsh_to_meshio_perm(const std::string& rT) {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"tetra10", {0, 1, 2, 3, 4, 5, 6, 7, 9, 8}},
        {"hexahedron20", {0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 9, 16, 18, 19, 17, 10, 12, 14, 15}},
        {"hexahedron27", {0,  1,  2,  3,  4,  5,  6,  7,  8,  11, 13, 9,  16, 18,
                          19, 17, 10, 12, 14, 15, 22, 23, 21, 24, 20, 25, 26}},
        {"wedge15", {0, 1, 2, 3, 4, 5, 6, 9, 7, 12, 14, 13, 8, 10, 11}},
        {"pyramid13", {0, 1, 2, 3, 4, 5, 8, 10, 6, 7, 9, 11, 12}},
    };
    static const std::vector<int> empty;
    auto it = m.find(rT);
    return it == m.end() ? empty : it->second;
}

const std::vector<int>& meshio_to_gmsh_perm(const std::string& rT) {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"tetra10", {0, 1, 2, 3, 4, 5, 6, 7, 9, 8}},
        {"hexahedron20", {0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 16, 9, 17, 10, 18, 19, 12, 15, 13, 14}},
        {"hexahedron27", {0,  1,  2,  3,  4,  5,  6,  7,  8,  11, 16, 9,  17, 10,
                          18, 19, 12, 15, 13, 14, 24, 22, 20, 21, 23, 25, 26}},
        {"wedge15", {0, 1, 2, 3, 4, 5, 6, 8, 12, 7, 13, 14, 9, 11, 10}},
        {"pyramid13", {0, 1, 2, 3, 4, 5, 8, 9, 6, 10, 7, 11, 12}},
    };
    static const std::vector<int> empty;
    auto it = m.find(rT);
    return it == m.end() ? empty : it->second;
}

std::string trim(const std::string& rS) {
    std::size_t b = 0, e = rS.size();
    while (b < e && std::isspace(static_cast<unsigned char>(rS[b])))
        ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(rS[e - 1])))
        --e;
    return rS.substr(b, e - b);
}

struct Cursor {
    const std::string& mBuf;
    std::size_t mPos = 0;
    explicit Cursor(const std::string& rB) : mBuf(rB) {}
    bool eof() const { return mPos >= mBuf.size(); }

    std::string read_line() {
        std::size_t start = mPos;
        while (mPos < mBuf.size() && mBuf[mPos] != '\n')
            ++mPos;
        std::string line = mBuf.substr(start, mPos - start);
        if (mPos < mBuf.size())
            ++mPos;
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        return line;
    }
    std::string next_nonblank() {
        while (!eof()) {
            std::string l = read_line();
            if (!trim(l).empty())
                return l;
        }
        return "";
    }
    void skip_to_end(const std::string& rEnv) {
        std::string target = "$End" + rEnv;
        while (!eof()) {
            if (trim(read_line()) == target)
                return;
        }
    }
    double next_double() {
        const char* base = mBuf.c_str();
        char* endp = nullptr;
        double v = std::strtod(base + mPos, &endp);
        if (endp == base + mPos)
            throw ReadError("Gmsh: expected a number");
        mPos = static_cast<std::size_t>(endp - base);
        return v;
    }
    std::int64_t next_int() { return static_cast<std::int64_t>(next_double()); }

    std::int32_t read_i32() {
        std::int32_t v;
        std::memcpy(&v, mBuf.data() + mPos, 4);
        mPos += 4;
        return v;
    }
    double read_f64() {
        double v;
        std::memcpy(&v, mBuf.data() + mPos, 8);
        mPos += 8;
        return v;
    }
    // Read an unsigned integer of `sz` bytes (little-endian host).
    std::uint64_t read_uint(int sz) {
        std::uint64_t v = 0;
        std::memcpy(&v, mBuf.data() + mPos, static_cast<std::size_t>(sz));
        mPos += static_cast<std::size_t>(sz);
        return v;
    }
};

struct EBlock {
    std::string mType;
    std::size_t mN = 0;
    std::size_t mCount = 0;
    std::size_t mNumTags = 0;
    std::vector<std::int64_t> mConn;  // count*n, 0-based gmsh ids
    std::vector<std::int64_t> mTags;  // count*num_tags
};

void store_value(NDArray& rA, std::size_t i, double d) {
    switch (rA.Dtype()) {
        case DType::Float64:
            rA.As<double>()[i] = d;
            break;
        case DType::Float32:
            rA.As<float>()[i] = static_cast<float>(d);
            break;
        case DType::Int64:
            rA.As<std::int64_t>()[i] = static_cast<std::int64_t>(d);
            break;
        case DType::Int32:
            rA.As<std::int32_t>()[i] = static_cast<std::int32_t>(d);
            break;
        default:
            rA.As<double>()[i] = d;
            break;
    }
}

void read_physical_names(Cursor& rCur, std::unordered_map<std::string, NDArray>& rFieldData) {
    std::int64_t num = std::stoll(trim(rCur.read_line()));
    for (std::int64_t i = 0; i < num; ++i) {
        std::string line = rCur.read_line();
        std::istringstream iss(line);
        long long dim, tag;
        iss >> dim >> tag;
        std::size_t q1 = line.find('"');
        std::size_t q2 = line.rfind('"');
        std::string name =
            (q1 != std::string::npos && q2 > q1) ? line.substr(q1 + 1, q2 - q1 - 1) : "";
        NDArray v(DType::Int64, {2});
        v.As<std::int64_t>()[0] = tag;  // physical number
        v.As<std::int64_t>()[1] = dim;
        rFieldData.emplace(name, std::move(v));
    }
    rCur.skip_to_end("PhysicalNames");
}

void read_nodes(Cursor& rCur, bool is_ascii, NDArray& rPoints,
                std::vector<std::int64_t>& rPointTags) {
    std::int64_t num = std::stoll(trim(rCur.read_line()));
    rPoints = NDArray(DType::Float64, {static_cast<std::size_t>(num), 3});
    rPointTags.resize(num);
    double* pp = rPoints.As<double>();
    if (is_ascii) {
        for (std::int64_t i = 0; i < num; ++i) {
            rPointTags[i] = rCur.next_int();
            pp[i * 3 + 0] = rCur.next_double();
            pp[i * 3 + 1] = rCur.next_double();
            pp[i * 3 + 2] = rCur.next_double();
        }
    } else {
        for (std::int64_t i = 0; i < num; ++i) {
            rPointTags[i] = rCur.read_i32();
            pp[i * 3 + 0] = rCur.read_f64();
            pp[i * 3 + 1] = rCur.read_f64();
            pp[i * 3 + 2] = rCur.read_f64();
        }
    }
    rCur.skip_to_end("Nodes");
}

void append_element(std::vector<EBlock>& rBlocks, const std::string& rType, std::size_t n,
                    std::size_t num_tags, const std::int64_t* pTags, const std::int64_t* pNodes) {
    if (rBlocks.empty() || rBlocks.back().mType != rType || rBlocks.back().mNumTags != num_tags) {
        EBlock b;
        b.mType = rType;
        b.mN = n;
        b.mNumTags = num_tags;
        rBlocks.push_back(std::move(b));
    }
    EBlock& cur = rBlocks.back();
    for (std::size_t j = 0; j < num_tags; ++j)
        cur.mTags.push_back(pTags[j]);
    for (std::size_t j = 0; j < n; ++j)
        cur.mConn.push_back(pNodes[j] - 1);
    ++cur.mCount;
}

void read_elements(Cursor& rCur, bool is_ascii, std::vector<EBlock>& rBlocks) {
    std::int64_t total = std::stoll(trim(rCur.read_line()));
    const auto& g2m = gmsh_to_meshio_type();
    const auto& nnpc = num_nodes_per_cell();

    if (is_ascii) {
        for (std::int64_t e = 0; e < total; ++e) {
            std::string line = rCur.read_line();
            std::istringstream iss(line);
            std::vector<std::int64_t> v;
            long long x;
            while (iss >> x)
                v.push_back(x);
            int gtype = static_cast<int>(v[1]);
            std::size_t num_tags = static_cast<std::size_t>(v[2]);
            auto it = g2m.find(gtype);
            if (it == g2m.end())
                throw ReadError("Gmsh element type " + std::to_string(gtype) +
                                " not supported by the C++ reader");
            std::size_t n = static_cast<std::size_t>(nnpc.at(it->second));
            append_element(rBlocks, it->second, n, num_tags, v.data() + 3, v.data() + 3 + num_tags);
        }
    } else {
        std::int64_t done = 0;
        while (done < total) {
            int gtype = rCur.read_i32();
            std::int32_t nelem = rCur.read_i32();
            std::int32_t num_tags = rCur.read_i32();
            auto it = g2m.find(gtype);
            if (it == g2m.end())
                throw ReadError("Gmsh element type " + std::to_string(gtype) +
                                " not supported by the C++ reader");
            std::size_t n = static_cast<std::size_t>(nnpc.at(it->second));
            std::vector<std::int64_t> tags(num_tags), nodes(n);
            for (std::int32_t k = 0; k < nelem; ++k) {
                rCur.read_i32();  // element id
                for (std::int32_t j = 0; j < num_tags; ++j)
                    tags[j] = rCur.read_i32();
                for (std::size_t j = 0; j < n; ++j)
                    nodes[j] = rCur.read_i32();
                append_element(rBlocks, it->second, n, num_tags, tags.data(), nodes.data());
            }
            done += nelem;
        }
    }
    rCur.skip_to_end("Elements");
}

// NodeData / ElementData
void read_data(Cursor& rCur, const std::string& rTag, bool is_ascii,
               std::unordered_map<std::string, NDArray>& rOut) {
    std::int64_t num_str = std::stoll(trim(rCur.read_line()));
    std::string name;
    for (std::int64_t i = 0; i < num_str; ++i) {
        std::string s = trim(rCur.read_line());
        if (i == 0) {
            // strip quotes
            std::size_t q1 = s.find('"'), q2 = s.rfind('"');
            name = (q1 != std::string::npos && q2 > q1) ? s.substr(q1 + 1, q2 - q1 - 1) : s;
        }
    }
    std::int64_t num_real = std::stoll(trim(rCur.read_line()));
    for (std::int64_t i = 0; i < num_real; ++i)
        rCur.read_line();
    std::int64_t num_int = std::stoll(trim(rCur.read_line()));
    std::vector<std::int64_t> itags(num_int);
    for (std::int64_t i = 0; i < num_int; ++i)
        itags[i] = std::stoll(trim(rCur.read_line()));
    std::size_t ncomp = static_cast<std::size_t>(itags[1]);
    std::size_t nitems = static_cast<std::size_t>(itags[2]);

    NDArray data(DType::Float64, {nitems, ncomp});
    double* dp = data.As<double>();
    if (is_ascii) {
        for (std::size_t i = 0; i < nitems; ++i) {
            rCur.next_int();  // index
            for (std::size_t c = 0; c < ncomp; ++c)
                dp[i * ncomp + c] = rCur.next_double();
        }
    } else {
        for (std::size_t i = 0; i < nitems; ++i) {
            rCur.read_i32();  // index
            for (std::size_t c = 0; c < ncomp; ++c)
                dp[i * ncomp + c] = rCur.read_f64();
        }
    }
    rCur.skip_to_end(rTag);
    if (ncomp == 1)
        data.Reshape({nitems});
    rOut.emplace(name, std::move(data));
}

NDArray slice_rows(const NDArray& rA, std::size_t r0, std::size_t r1) {
    std::size_t nc = rA.Shape().size() >= 2 ? rA.Shape()[1] : 1;
    std::size_t isz = dtype_size(rA.Dtype());
    std::vector<std::size_t> shape = rA.Shape();
    shape[0] = r1 - r0;
    NDArray out(rA.Dtype(), shape);
    if (r1 > r0)
        std::memcpy(out.Data(), rA.Data() + r0 * nc * isz, (r1 - r0) * nc * isz);
    return out;
}

// ---- version 4.1 -------------------------------------------------------------

struct E41 {
    std::string mType;
    std::size_t mN = 0;
    std::size_t mCount = 0;
    int mEntityTag = 0;
    NDArray mConn;  // (count, n) Int64, 0-based gmsh node ids; moved into the
                    // cell block directly when the tag remap is the identity.
};

void read_nodes_41(Cursor& rCur, bool is_ascii, int data_size, NDArray& rPoints,
                   std::vector<std::int64_t>& rTags,
                   std::vector<std::array<std::int64_t, 2>>& rDimTags) {
    auto rd_size = [&]() -> std::int64_t {
        return is_ascii ? rCur.next_int() : static_cast<std::int64_t>(rCur.read_uint(data_size));
    };
    auto rd_int = [&]() -> int {
        return is_ascii ? static_cast<int>(rCur.next_int()) : rCur.read_i32();
    };
    auto rd_dbl = [&]() -> double { return is_ascii ? rCur.next_double() : rCur.read_f64(); };

    std::int64_t num_blocks = rd_size();
    std::int64_t num_nodes = rd_size();
    rd_size();  // min tag
    rd_size();  // max tag
    rPoints = NDArray(DType::Float64, {static_cast<std::size_t>(num_nodes), 3});
    rTags.resize(num_nodes);
    rDimTags.resize(num_nodes);
    double* pp = rPoints.As<double>();

    std::size_t idx = 0;
    for (std::int64_t b = 0; b < num_blocks; ++b) {
        int dim = rd_int();
        int entity_tag = rd_int();
        int parametric = rd_int();
        if (parametric != 0)
            throw ReadError("parametric Gmsh nodes not supported");
        std::int64_t nb = rd_size();
        const std::size_t nbz = static_cast<std::size_t>(nb);
        if (!is_ascii && data_size == 8) {
            // Native-endian, contiguous: bulk-copy tags (u64) and coords (3*f64).
            std::memcpy(&rTags[idx], rCur.mBuf.data() + rCur.mPos, nbz * 8);
            rCur.mPos += nbz * 8;
            for (std::size_t i = 0; i < nbz; ++i)
                rTags[idx + i] -= 1;
            std::memcpy(pp + idx * 3, rCur.mBuf.data() + rCur.mPos, nbz * 3 * 8);
            rCur.mPos += nbz * 3 * 8;
        } else {
            for (std::int64_t i = 0; i < nb; ++i)
                rTags[idx + i] = rd_size() - 1;
            for (std::int64_t i = 0; i < nb; ++i) {
                pp[(idx + i) * 3 + 0] = rd_dbl();
                pp[(idx + i) * 3 + 1] = rd_dbl();
                pp[(idx + i) * 3 + 2] = rd_dbl();
            }
        }
        for (std::int64_t i = 0; i < nb; ++i)
            rDimTags[idx + i] = {dim, entity_tag};
        idx += static_cast<std::size_t>(nb);
    }
    rCur.skip_to_end("Nodes");
}

void read_elements_41(Cursor& rCur, bool is_ascii, int data_size, std::vector<E41>& rBlocks) {
    auto rd_size = [&]() -> std::int64_t {
        return is_ascii ? rCur.next_int() : static_cast<std::int64_t>(rCur.read_uint(data_size));
    };
    auto rd_int = [&]() -> int {
        return is_ascii ? static_cast<int>(rCur.next_int()) : rCur.read_i32();
    };

    std::int64_t num_blocks = rd_size();
    rd_size();  // num elements
    rd_size();  // min tag
    rd_size();  // max tag
    const auto& g2m = gmsh_to_meshio_type();
    const auto& nnpc = num_nodes_per_cell();

    for (std::int64_t b = 0; b < num_blocks; ++b) {
        rd_int();  // entity dim
        int entity_tag = rd_int();
        int etype = rd_int();
        std::int64_t num_ele = rd_size();
        auto it = g2m.find(etype);
        if (it == g2m.end())
            throw ReadError("Gmsh element type " + std::to_string(etype) +
                            " not supported by the C++ reader");
        std::size_t n = static_cast<std::size_t>(nnpc.at(it->second));
        E41 blk;
        blk.mType = it->second;
        blk.mN = n;
        blk.mCount = static_cast<std::size_t>(num_ele);
        blk.mEntityTag = entity_tag;
        const std::size_t nez = static_cast<std::size_t>(num_ele);
        blk.mConn = NDArray(DType::Int64, {nez, n});
        std::int64_t* dst = blk.mConn.As<std::int64_t>();
        if (!is_ascii && data_size == 8) {
            // Each element is [tag, node0..node(n-1)] u64, native-endian and
            // contiguous. Decode the nodes straight from the slurped buffer into
            // the owning connectivity array (drop the tag), one parallel pass.
            const std::size_t stride = n + 1;
            const char* base = rCur.mBuf.data() + rCur.mPos;
            parallel_for_bw(nez, [&](std::size_t e) {
                const char* row = base + (e * stride + 1) * 8;  // skip element tag
                for (std::size_t j = 0; j < n; ++j) {
                    std::uint64_t v;
                    std::memcpy(&v, row + j * 8, 8);
                    dst[e * n + j] = static_cast<std::int64_t>(v) - 1;
                }
            });
            rCur.mPos += nez * stride * 8;
        } else {
            std::size_t p = 0;
            for (std::int64_t e = 0; e < num_ele; ++e) {
                rd_size();  // element tag
                for (std::size_t j = 0; j < n; ++j)
                    dst[p++] = rd_size() - 1;
            }
        }
        rBlocks.push_back(std::move(blk));
    }
    rCur.skip_to_end("Elements");
}

Mesh read_gmsh41_body(Cursor& rCur, bool is_ascii, int data_size) {
    NDArray points(DType::Float64, {0, 3});
    std::vector<std::int64_t> point_tags;
    std::vector<std::array<std::int64_t, 2>> dim_tags;
    std::vector<E41> eblocks;
    std::unordered_map<std::string, NDArray> field_data, point_data, cell_data_raw;

    while (!rCur.eof()) {
        std::string line = rCur.next_nonblank();
        if (line.empty())
            break;
        if (line[0] != '$')
            throw ReadError("Gmsh: unexpected line " + line);
        std::string env = trim(line.substr(1));
        if (env == "PhysicalNames")
            read_physical_names(rCur, field_data);
        else if (env == "Entities")
            throw ReadError("Gmsh $Entities not supported by the C++ reader");
        else if (env == "Nodes")
            read_nodes_41(rCur, is_ascii, data_size, points, point_tags, dim_tags);
        else if (env == "Elements")
            read_elements_41(rCur, is_ascii, data_size, eblocks);
        else if (env == "Periodic")
            throw ReadError("Gmsh $Periodic not supported by the C++ reader");
        else if (env == "NodeData")
            read_data(rCur, "NodeData", is_ascii, point_data);
        else if (env == "ElementData")
            read_data(rCur, "ElementData", is_ascii, cell_data_raw);
        else
            rCur.skip_to_end(env);
    }

    // When node tags are contiguous 0..N-1 (the common case) the tag->row remap
    // is the identity, so we can skip building it *and* skip the random-access
    // gather below (the connectivity is already the final mesh indexing).
    bool remap_identity = true;
    for (std::size_t i = 0; i < point_tags.size(); ++i)
        if (point_tags[i] != static_cast<std::int64_t>(i)) {
            remap_identity = false;
            break;
        }
    std::vector<std::int64_t> remap;
    if (!remap_identity) {
        std::int64_t max_tag = 0;
        for (auto t : point_tags)
            max_tag = std::max(max_tag, t);
        remap.assign(static_cast<std::size_t>(max_tag) + 1, -1);
        // Scatter: node tags are unique, so writes never alias -> parallel.
        parallel_for_bw(point_tags.size(), [&](std::size_t i) {
            remap[static_cast<std::size_t>(point_tags[i])] = static_cast<std::int64_t>(i);
        });
    }

    Mesh mesh;
    mesh.AssignPoints(std::move(points));
    for (auto& kv : point_data)
        mesh.AddPointData(kv.first, std::move(kv.second));
    for (auto& kv : field_data)
        mesh.AddFieldData(kv.first, std::move(kv.second));

    // Node entity (dim, tag) -> gmsh:dim_tags point data.
    NDArray dt(DType::Int64, {dim_tags.size(), 2});
    parallel_for_bw(dim_tags.size(), [&](std::size_t i) {
        dt.As<std::int64_t>()[i * 2 + 0] = dim_tags[i][0];
        dt.As<std::int64_t>()[i * 2 + 1] = dim_tags[i][1];
    });
    mesh.AddPointData("gmsh:dim_tags", std::move(dt));

    std::vector<NDArray> geom_blocks;
    for (auto& b : eblocks) {
        const std::vector<int>& perm = gmsh_to_meshio_perm(b.mType);
        const int* prm = perm.empty() ? nullptr : perm.data();
        if (remap_identity && !prm) {
            // Identity remap, no reorder -> the connectivity is already final:
            // move the owning (count, n) array straight into the cell block.
            mesh.AddCellBlock(b.mType, std::move(b.mConn));
        } else {
            NDArray data(DType::Int64, {b.mCount, b.mN});
            std::int64_t* dp = data.As<std::int64_t>();
            const std::int64_t* cn = b.mConn.As<std::int64_t>();
            if (remap_identity) {
                parallel_for_bw(b.mCount, [&](std::size_t r) {
                    for (std::size_t j = 0; j < b.mN; ++j)
                        dp[r * b.mN + j] = cn[r * b.mN + static_cast<std::size_t>(prm[j])];
                });
            } else {
                // Gather through the prebuilt read-only remap -> parallel by row.
                parallel_for_bw(b.mCount, [&](std::size_t r) {
                    for (std::size_t j = 0; j < b.mN; ++j) {
                        std::size_t src = prm ? static_cast<std::size_t>(prm[j]) : j;
                        dp[r * b.mN + j] = remap[static_cast<std::size_t>(cn[r * b.mN + src])];
                    }
                });
            }
            mesh.AddCellBlock(b.mType, std::move(data));
        }

        NDArray ge(DType::Int32, {b.mCount});
        std::int32_t* gep = ge.As<std::int32_t>();
        const std::int32_t etag = b.mEntityTag;
        parallel_for_bw(b.mCount, [&](std::size_t r) { gep[r] = etag; });
        geom_blocks.push_back(std::move(ge));
    }

    for (auto& kv : cell_data_raw) {
        std::vector<NDArray> per_block;
        std::size_t offset = 0;
        for (const auto& b : eblocks) {
            per_block.push_back(slice_rows(kv.second, offset, offset + b.mCount));
            offset += b.mCount;
        }
        mesh.AddCellData(kv.first, std::move(per_block));
    }
    if (!geom_blocks.empty())
        mesh.AddCellData("gmsh:geometrical", std::move(geom_blocks));

    return mesh;
}

}  // namespace

Mesh read_gmsh(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    // Bulk slurp (seek+read) rather than char-by-char istreambuf_iterator.
    in.seekg(0, std::ios::end);
    std::streamoff flen = in.tellg();
    in.seekg(0, std::ios::beg);
    std::string buf;
    if (flen > 0) {
        buf.resize(static_cast<std::size_t>(flen));
        in.read(buf.data(), flen);
    }
    Cursor cur(buf);

    if (trim(cur.read_line()) != "$MeshFormat")
        throw ReadError("Expected $MeshFormat");
    std::string fmt = cur.read_line();
    std::istringstream fss(fmt);
    std::string version;
    int file_type = 0, data_size = 8;
    fss >> version >> file_type >> data_size;
    bool is_ascii = (file_type == 0);
    if (!is_ascii) {
        cur.read_i32();  // endianness marker
        // consume trailing newline before $EndMeshFormat
        if (cur.mPos < buf.size() && buf[cur.mPos] == '\n')
            ++cur.mPos;
    }
    cur.skip_to_end("MeshFormat");

    if (version == "4.1" || version == "4")
        return read_gmsh41_body(cur, is_ascii, data_size);
    if (version.rfind("2", 0) != 0)
        throw ReadError("C++ Gmsh reader handles versions 2.2 and 4.1 only");

    NDArray points(DType::Float64, {0, 3});
    std::vector<std::int64_t> point_tags;
    std::vector<EBlock> eblocks;
    std::unordered_map<std::string, NDArray> field_data, point_data, cell_data_raw;

    while (!cur.eof()) {
        std::string line = cur.next_nonblank();
        if (line.empty())
            break;
        if (line[0] != '$')
            throw ReadError("Gmsh: unexpected line " + line);
        std::string env = trim(line.substr(1));
        if (env == "PhysicalNames")
            read_physical_names(cur, field_data);
        else if (env == "Nodes")
            read_nodes(cur, is_ascii, points, point_tags);
        else if (env == "Elements")
            read_elements(cur, is_ascii, eblocks);
        else if (env == "Periodic")
            throw ReadError("Gmsh $Periodic not supported by the C++ reader");
        else if (env == "NodeData")
            read_data(cur, "NodeData", is_ascii, point_data);
        else if (env == "ElementData")
            read_data(cur, "ElementData", is_ascii, cell_data_raw);
        else
            cur.skip_to_end(env);
    }

    // Build node-tag remap (gmsh ids are 1-based, possibly non-contiguous).
    std::int64_t max_tag = 0;
    for (auto t : point_tags)
        max_tag = std::max(max_tag, t - 1);
    std::vector<std::int64_t> remap(static_cast<std::size_t>(max_tag) + 1, -1);
    // Scatter: node tags are unique, so writes never alias -> parallel.
    parallel_for_bw(point_tags.size(), [&](std::size_t i) {
        remap[static_cast<std::size_t>(point_tags[i] - 1)] = static_cast<std::int64_t>(i);
    });

    Mesh mesh;
    mesh.AssignPoints(std::move(points));
    for (auto& kv : point_data)
        mesh.AddPointData(kv.first, std::move(kv.second));
    for (auto& kv : field_data)
        mesh.AddFieldData(kv.first, std::move(kv.second));

    // Determine which tag columns are present across all blocks.
    std::size_t min_tags = eblocks.empty() ? 0 : SIZE_MAX;
    for (const auto& b : eblocks)
        min_tags = std::min(min_tags, b.mNumTags);

    std::vector<NDArray> physical_blocks, geometrical_blocks;
    for (const auto& b : eblocks) {
        const std::vector<int>& perm = gmsh_to_meshio_perm(b.mType);
        NDArray data(DType::Int64, {b.mCount, b.mN});
        std::int64_t* dp = data.As<std::int64_t>();
        // Gather through the prebuilt read-only remap -> parallel over rows.
        parallel_for_bw(b.mCount, [&](std::size_t r) {
            for (std::size_t j = 0; j < b.mN; ++j) {
                std::size_t src = perm.empty() ? j : static_cast<std::size_t>(perm[j]);
                std::int64_t gid = b.mConn[r * b.mN + src];
                dp[r * b.mN + j] = remap[static_cast<std::size_t>(gid)];
            }
        });
        mesh.AddCellBlock(b.mType, std::move(data));

        if (min_tags >= 1) {
            NDArray ph(DType::Int32, {b.mCount});
            std::int32_t* php = ph.As<std::int32_t>();
            parallel_for_bw(b.mCount, [&](std::size_t r) {
                php[r] = static_cast<std::int32_t>(b.mTags[r * b.mNumTags + 0]);
            });
            physical_blocks.push_back(std::move(ph));
        }
        if (min_tags >= 2) {
            NDArray ge(DType::Int32, {b.mCount});
            std::int32_t* gep = ge.As<std::int32_t>();
            parallel_for_bw(b.mCount, [&](std::size_t r) {
                gep[r] = static_cast<std::int32_t>(b.mTags[r * b.mNumTags + 1]);
            });
            geometrical_blocks.push_back(std::move(ge));
        }
    }

    // Split ElementData (concatenated over blocks) back per block.
    for (auto& kv : cell_data_raw) {
        std::vector<NDArray> per_block;
        std::size_t offset = 0;
        for (const auto& b : eblocks) {
            per_block.push_back(slice_rows(kv.second, offset, offset + b.mCount));
            offset += b.mCount;
        }
        mesh.AddCellData(kv.first, std::move(per_block));
    }
    if (!physical_blocks.empty())
        mesh.AddCellData("gmsh:physical", std::move(physical_blocks));
    if (!geometrical_blocks.empty())
        mesh.AddCellData("gmsh:geometrical", std::move(geometrical_blocks));

    return mesh;
}

// ---- writer ------------------------------------------------------------------

namespace {

void write_physical_names(std::ostream& rOs, const Mesh& rMesh) {
    std::vector<std::tuple<long long, long long, std::string>> sortable;  // dim, num, name
    for (const auto& name : rMesh.FieldDataNames()) {
        const NDArray& d = rMesh.FieldData(name);
        if (d.Size() < 2)
            continue;
        long long num = detail::read_int(d, 0);
        long long dim = detail::read_int(d, 1);
        sortable.emplace_back(dim, num, name);
    }
    if (sortable.empty())
        return;
    std::sort(sortable.begin(), sortable.end());
    rOs << "$PhysicalNames\n" << sortable.size() << "\n";
    for (auto& e : sortable)
        rOs << std::get<0>(e) << ' ' << std::get<1>(e) << " \"" << std::get<2>(e) << "\"\n";
    rOs << "$EndPhysicalNames\n";
}

// Writes the cell-data array named `rName` as one $ElementData-style section,
// concatenated across cell blocks.
void write_data(std::ostream& rOs, const char* pTag, const std::string& rName, const Mesh& rMesh,
                bool binary) {
    // Concatenate blocks.
    const std::size_t nblocks = rMesh.CellDataNumBlocks(rName);
    std::size_t total = 0, ncomp = 1;
    for (std::size_t k = 0; k < nblocks; ++k) {
        const NDArray& b = rMesh.CellData(rName, k);
        total += b.Shape().empty() ? 0 : b.Shape()[0];
        ncomp = b.Shape().size() >= 2 ? b.Shape()[1] : 1;
    }
    rOs << "$" << pTag << "\n1\n\"" << rName << "\"\n1\n0\n3\n0\n"
        << ncomp << "\n"
        << total << "\n";
    std::int64_t idx = 1;
    for (std::size_t k = 0; k < nblocks; ++k) {
        const NDArray& b = rMesh.CellData(rName, k);
        std::size_t rows = b.Shape().empty() ? 0 : b.Shape()[0];
        for (std::size_t r = 0; r < rows; ++r) {
            if (binary) {
                std::int32_t id = static_cast<std::int32_t>(idx);
                rOs.write(reinterpret_cast<const char*>(&id), 4);
                for (std::size_t c = 0; c < ncomp; ++c) {
                    double v = detail::read_double(b, r * ncomp + c);
                    rOs.write(reinterpret_cast<const char*>(&v), 8);
                }
            } else {
                rOs << idx;
                char buf[32];
                for (std::size_t c = 0; c < ncomp; ++c) {
                    std::snprintf(buf, sizeof(buf), " %.17g",
                                  detail::read_double(b, r * ncomp + c));
                    rOs << buf;
                }
                rOs << '\n';
            }
            ++idx;
        }
    }
    if (binary)
        rOs << '\n';
    rOs << "$End" << pTag << "\n";
}

}  // namespace

void write_gmsh22(const std::string& rPath, const Mesh& rMesh, bool binary) {
    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t num_points = rMesh.NumPoints();
    const NDArray& points = rMesh.Points();
    const std::size_t dim = points.Shape().size() >= 2 ? points.Shape()[1] : 0;
    const std::size_t nblocks = rMesh.NumCellBlocks();

    // Tag cell data ("gmsh:physical"/"gmsh:geometrical") is written inline with
    // the elements; per-block zeros stand in when a tag column is absent.
    const bool has_physical = rMesh.HasCellData("gmsh:physical");
    const bool has_geometrical = rMesh.HasCellData("gmsh:geometrical");
    std::vector<NDArray> zeros_phys, zeros_geom;
    if (!has_physical)
        for (const auto cb : rMesh.CellRange())
            zeros_phys.emplace_back(DType::Int32, std::vector<std::size_t>{cb.NumCells()});
    if (!has_geometrical)
        for (const auto cb : rMesh.CellRange())
            zeros_geom.emplace_back(DType::Int32, std::vector<std::size_t>{cb.NumCells()});

    os << "$MeshFormat\n2.2 " << (binary ? 1 : 0) << " 8\n";
    if (binary) {
        std::int32_t one = 1;
        os.write(reinterpret_cast<const char*>(&one), 4);
        os << '\n';
    }
    os << "$EndMeshFormat\n";

    write_physical_names(os, rMesh);

    // Nodes.
    os << "$Nodes\n" << num_points << "\n";
    if (binary) {
        for (std::size_t i = 0; i < num_points; ++i) {
            std::int32_t id = static_cast<std::int32_t>(i + 1);
            os.write(reinterpret_cast<const char*>(&id), 4);
            for (std::size_t c = 0; c < 3; ++c) {
                double v = (c < dim) ? detail::read_double(points, i * dim + c) : 0.0;
                os.write(reinterpret_cast<const char*>(&v), 8);
            }
        }
        os << '\n';
    } else {
        char buf[80];
        for (std::size_t i = 0; i < num_points; ++i) {
            double x = (0 < dim) ? detail::read_double(points, i * dim + 0) : 0.0;
            double y = (1 < dim) ? detail::read_double(points, i * dim + 1) : 0.0;
            double z = (2 < dim) ? detail::read_double(points, i * dim + 2) : 0.0;
            std::snprintf(buf, sizeof(buf), "%zu %.16e %.16e %.16e\n", i + 1, x, y, z);
            os << buf;
        }
    }
    os << "$EndNodes\n";

    // Elements.
    std::size_t total_cells = 0;
    for (const auto cb : rMesh.CellRange())
        total_cells += cb.NumCells();
    os << "$Elements\n" << total_cells << "\n";
    const auto& m2g = meshio_to_gmsh_type();
    std::size_t consecutive = 0;
    for (std::size_t k = 0; k < nblocks; ++k) {
        const auto cb = rMesh.Cells(k);
        auto it = m2g.find(cb.Type());
        if (it == m2g.end())
            throw WriteError("Gmsh writer: unsupported cell type " + cb.Type());
        int gtype = it->second;
        const NDArray& conn = cb.Conn();
        std::size_t n = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
        const std::vector<int>& perm = meshio_to_gmsh_perm(cb.Type());
        std::size_t count = cb.NumCells();
        const NDArray& ph = has_physical ? rMesh.CellData("gmsh:physical", k) : zeros_phys[k];
        const NDArray& ge = has_geometrical ? rMesh.CellData("gmsh:geometrical", k) : zeros_geom[k];

        if (binary) {
            std::int32_t hdr[3] = {gtype, static_cast<std::int32_t>(count), 2};
            os.write(reinterpret_cast<const char*>(hdr), 12);
            for (std::size_t r = 0; r < count; ++r) {
                std::int32_t id = static_cast<std::int32_t>(consecutive + r + 1);
                std::int32_t t0 = static_cast<std::int32_t>(detail::read_int(ph, r));
                std::int32_t t1 = static_cast<std::int32_t>(detail::read_int(ge, r));
                os.write(reinterpret_cast<const char*>(&id), 4);
                os.write(reinterpret_cast<const char*>(&t0), 4);
                os.write(reinterpret_cast<const char*>(&t1), 4);
                for (std::size_t j = 0; j < n; ++j) {
                    std::size_t src = perm.empty() ? j : static_cast<std::size_t>(perm[j]);
                    std::int32_t node =
                        static_cast<std::int32_t>(detail::read_int(conn, r * n + src) + 1);
                    os.write(reinterpret_cast<const char*>(&node), 4);
                }
            }
        } else {
            for (std::size_t r = 0; r < count; ++r) {
                os << (consecutive + r + 1) << ' ' << gtype << " 2 " << detail::read_int(ph, r)
                   << ' ' << detail::read_int(ge, r);
                for (std::size_t j = 0; j < n; ++j) {
                    std::size_t src = perm.empty() ? j : static_cast<std::size_t>(perm[j]);
                    os << ' ' << (detail::read_int(conn, r * n + src) + 1);
                }
                os << '\n';
            }
        }
        consecutive += count;
    }
    if (binary)
        os << '\n';
    os << "$EndElements\n";

    for (const auto& name : rMesh.PointDataNames()) {
        if (name == "gmsh:dim_tags")
            continue;
        // Reusing write_data (cell-data-shaped) for point data is awkward; inline:
        const NDArray& d = rMesh.PointData(name);
        std::size_t ncomp = d.Shape().size() >= 2 ? d.Shape()[1] : 1;
        std::size_t rows = d.Shape().empty() ? 0 : d.Shape()[0];
        os << "$NodeData\n1\n\"" << name << "\"\n1\n0\n3\n0\n" << ncomp << "\n" << rows << "\n";
        char buf[32];
        for (std::size_t r = 0; r < rows; ++r) {
            if (binary) {
                std::int32_t id = static_cast<std::int32_t>(r + 1);
                os.write(reinterpret_cast<const char*>(&id), 4);
                for (std::size_t c = 0; c < ncomp; ++c) {
                    double v = detail::read_double(d, r * ncomp + c);
                    os.write(reinterpret_cast<const char*>(&v), 8);
                }
            } else {
                os << (r + 1);
                for (std::size_t c = 0; c < ncomp; ++c) {
                    std::snprintf(buf, sizeof(buf), " %.17g",
                                  detail::read_double(d, r * ncomp + c));
                    os << buf;
                }
                os << '\n';
            }
        }
        if (binary)
            os << '\n';
        os << "$EndNodeData\n";
    }

    for (const auto& name : rMesh.CellDataNames()) {
        if (name == "gmsh:physical" || name == "gmsh:geometrical" || name == "cell_tags")
            continue;
        write_data(os, "ElementData", name, rMesh, binary);
    }
}

void write_gmsh41(const std::string& rPath, const Mesh& rMesh, bool binary) {
    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t num_points = rMesh.NumPoints();
    const NDArray& points = rMesh.Points();
    const std::size_t dim = points.Shape().size() >= 2 ? points.Shape()[1] : 0;
    const int data_size = 8;

    auto put_u64 = [&](std::uint64_t v) { os.write(reinterpret_cast<const char*>(&v), 8); };
    auto put_i32 = [&](std::int32_t v) { os.write(reinterpret_cast<const char*>(&v), 4); };
    auto put_f64 = [&](double v) { os.write(reinterpret_cast<const char*>(&v), 8); };

    // "gmsh:geometrical" supplies the per-block entity tag below; the other
    // tag names are excluded from the $NodeData/$ElementData sections.
    const bool has_geometrical = rMesh.HasCellData("gmsh:geometrical");

    const auto& topo = topological_dimension();
    auto cell_dim = [&](const std::string& t) -> int {
        auto it = topo.find(t);
        return it == topo.end() ? 0 : it->second;
    };

    os << "$MeshFormat\n4.1 " << (binary ? 1 : 0) << " " << data_size << "\n";
    if (binary) {
        put_i32(1);
        os << '\n';
    }
    os << "$EndMeshFormat\n";

    write_physical_names(os, rMesh);

    // Nodes: a single entity block (no $Entities is emitted).
    int node_dim = rMesh.NumCellBlocks() == 0 ? 0 : cell_dim(rMesh.Cells(0).Type());
    os << "$Nodes\n";
    if (binary) {
        put_u64(1);
        put_u64(num_points);
        put_u64(1);
        put_u64(num_points);
        put_i32(node_dim);
        put_i32(0);
        put_i32(0);
        put_u64(num_points);
        // Node tags 1..num_points and the (3-padded) coords, each as one write
        // instead of a stream call per scalar (native endianness).
        std::vector<std::uint64_t> ntags(num_points);
        for (std::size_t i = 0; i < num_points; ++i)
            ntags[i] = i + 1;
        os.write(reinterpret_cast<const char*>(ntags.data()),
                 static_cast<std::streamsize>(num_points * 8));
        std::vector<double> cbuf(num_points * 3, 0.0);
        detail::dispatch_dtype(points.Dtype(), [&]<class T>() {
            const T* src = points.As<T>();
            parallel_for_bw(num_points, [&](std::size_t i) {
                for (std::size_t c = 0; c < dim && c < 3; ++c)
                    cbuf[i * 3 + c] = static_cast<double>(src[i * dim + c]);
            });
        });
        os.write(reinterpret_cast<const char*>(cbuf.data()),
                 static_cast<std::streamsize>(num_points * 3 * 8));
        os << '\n';
    } else {
        os << "1 " << num_points << " 1 " << num_points << "\n";
        os << node_dim << " 0 0 " << num_points << "\n";
        for (std::size_t i = 0; i < num_points; ++i)
            os << (i + 1) << "\n";
        char buf[80];
        for (std::size_t i = 0; i < num_points; ++i) {
            double x = (0 < dim) ? detail::read_double(points, i * dim + 0) : 0.0;
            double y = (1 < dim) ? detail::read_double(points, i * dim + 1) : 0.0;
            double z = (2 < dim) ? detail::read_double(points, i * dim + 2) : 0.0;
            std::snprintf(buf, sizeof(buf), "%.16e %.16e %.16e\n", x, y, z);
            os << buf;
        }
    }
    os << "$EndNodes\n";

    // Elements: one block per cell block.
    std::size_t total_cells = 0;
    for (const auto cb : rMesh.CellRange())
        total_cells += cb.NumCells();
    const auto& m2g = meshio_to_gmsh_type();
    os << "$Elements\n";
    if (binary) {
        put_u64(rMesh.NumCellBlocks());
        put_u64(total_cells);
        put_u64(1);
        put_u64(total_cells);
    } else {
        os << rMesh.NumCellBlocks() << " " << total_cells << " 1 " << total_cells << "\n";
    }
    std::size_t tag0 = 1;
    for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci) {
        const auto cb = rMesh.Cells(ci);
        auto it = m2g.find(cb.Type());
        if (it == m2g.end())
            throw WriteError("Gmsh writer: unsupported cell type " + cb.Type());
        int gtype = it->second;
        int bdim = cell_dim(cb.Type());
        int entity_tag =
            (has_geometrical && ci < rMesh.CellDataNumBlocks("gmsh:geometrical") &&
             rMesh.CellData("gmsh:geometrical", ci).Size() > 0)
                ? static_cast<int>(detail::read_int(rMesh.CellData("gmsh:geometrical", ci), 0))
                : 0;
        const NDArray& conn = cb.Conn();
        std::size_t n = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
        const std::vector<int>& perm = meshio_to_gmsh_perm(cb.Type());
        std::size_t count = cb.NumCells();
        if (binary) {
            put_i32(bdim);
            put_i32(entity_tag);
            put_i32(gtype);
            put_u64(count);
            // One buffer per block: [tag, node0..node(n-1)] u64, native, one write.
            const int* prm = perm.empty() ? nullptr : perm.data();
            const std::size_t stride = n + 1;
            std::vector<std::uint64_t> ebuf(count * stride);
            const std::uint64_t base = tag0;
            detail::dispatch_dtype(conn.Dtype(), [&]<class T>() {
                const T* src = conn.As<T>();
                parallel_for_bw(count, [&](std::size_t r) {
                    std::uint64_t* o = ebuf.data() + r * stride;
                    o[0] = base + r;
                    for (std::size_t j = 0; j < n; ++j) {
                        std::size_t sc = prm ? static_cast<std::size_t>(prm[j]) : j;
                        o[j + 1] = static_cast<std::uint64_t>(src[r * n + sc]) + 1;
                    }
                });
            });
            os.write(reinterpret_cast<const char*>(ebuf.data()),
                     static_cast<std::streamsize>(ebuf.size() * 8));
        } else {
            os << bdim << " " << entity_tag << " " << gtype << " " << count << "\n";
            for (std::size_t r = 0; r < count; ++r) {
                os << (tag0 + r);
                for (std::size_t j = 0; j < n; ++j) {
                    std::size_t src = perm.empty() ? j : static_cast<std::size_t>(perm[j]);
                    os << " " << (detail::read_int(conn, r * n + src) + 1);
                }
                os << "\n";
            }
        }
        tag0 += count;
    }
    if (binary)
        os << '\n';
    os << "$EndElements\n";

    for (const auto& name : rMesh.PointDataNames()) {
        if (name == "gmsh:dim_tags")
            continue;
        const NDArray& d = rMesh.PointData(name);
        std::size_t ncomp = d.Shape().size() >= 2 ? d.Shape()[1] : 1;
        std::size_t rows = d.Shape().empty() ? 0 : d.Shape()[0];
        os << "$NodeData\n1\n\"" << name << "\"\n1\n0\n3\n0\n" << ncomp << "\n" << rows << "\n";
        char buf[32];
        for (std::size_t r = 0; r < rows; ++r) {
            if (binary) {
                std::int32_t id = static_cast<std::int32_t>(r + 1);
                os.write(reinterpret_cast<const char*>(&id), 4);
                for (std::size_t c = 0; c < ncomp; ++c)
                    put_f64(detail::read_double(d, r * ncomp + c));
            } else {
                os << (r + 1);
                for (std::size_t c = 0; c < ncomp; ++c) {
                    std::snprintf(buf, sizeof(buf), " %.17g",
                                  detail::read_double(d, r * ncomp + c));
                    os << buf;
                }
                os << '\n';
            }
        }
        if (binary)
            os << '\n';
        os << "$EndNodeData\n";
    }

    for (const auto& name : rMesh.CellDataNames()) {
        if (name == "gmsh:physical" || name == "gmsh:geometrical" || name == "cell_tags")
            continue;
        write_data(os, "ElementData", name, rMesh, binary);
    }
}

}  // namespace meshioplusplus
