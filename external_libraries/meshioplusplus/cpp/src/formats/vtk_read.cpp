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
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vtk.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/types.hpp"

namespace meshioplusplus {

namespace {

DType dtype_from_vtk_token(std::string t) {
    for (auto& ch : t)
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    if (t == "float")
        return DType::Float32;
    if (t == "double")
        return DType::Float64;
    if (t == "int" || t == "vtktypeint64" || t == "long")
        return DType::Int64;
    if (t == "vtktypeint8" || t == "char")
        return DType::Int8;
    if (t == "vtktypeint16" || t == "short")
        return DType::Int16;
    if (t == "vtktypeint32")
        return DType::Int32;
    if (t == "vtktypeuint8" || t == "unsigned_char")
        return DType::UInt8;
    if (t == "vtktypeuint16")
        return DType::UInt16;
    if (t == "vtktypeuint32")
        return DType::UInt32;
    if (t == "vtktypeuint64")
        return DType::UInt64;
    throw ReadError("VTK data type '" + t + "' not supported by the C++ reader");
}

void store(NDArray& rA, std::size_t i, double d, std::int64_t v) {
    switch (rA.Dtype()) {
        case DType::Float32:
            rA.As<float>()[i] = static_cast<float>(d);
            break;
        case DType::Float64:
            rA.As<double>()[i] = d;
            break;
        case DType::Int8:
            rA.As<std::int8_t>()[i] = static_cast<std::int8_t>(v);
            break;
        case DType::Int16:
            rA.As<std::int16_t>()[i] = static_cast<std::int16_t>(v);
            break;
        case DType::Int32:
            rA.As<std::int32_t>()[i] = static_cast<std::int32_t>(v);
            break;
        case DType::Int64:
            rA.As<std::int64_t>()[i] = v;
            break;
        case DType::UInt8:
            rA.As<std::uint8_t>()[i] = static_cast<std::uint8_t>(v);
            break;
        case DType::UInt16:
            rA.As<std::uint16_t>()[i] = static_cast<std::uint16_t>(v);
            break;
        case DType::UInt32:
            rA.As<std::uint32_t>()[i] = static_cast<std::uint32_t>(v);
            break;
        case DType::UInt64:
            rA.As<std::uint64_t>()[i] = static_cast<std::uint64_t>(v);
            break;
    }
}

struct Cursor {
    const std::string& mBuf;
    std::size_t mPos = 0;

    explicit Cursor(const std::string& rB) : mBuf(rB) {}

    bool Eof() const { return mPos >= mBuf.size(); }

    std::string ReadLine() {
        std::size_t start = mPos;
        while (mPos < mBuf.size() && mBuf[mPos] != '\n')
            ++mPos;
        std::string line = mBuf.substr(start, mPos - start);
        if (mPos < mBuf.size())
            ++mPos;  // skip '\n'
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        return line;
    }

    void ConsumeEol() {
        while (mPos < mBuf.size() && mBuf[mPos] != '\n' &&
               std::isspace(static_cast<unsigned char>(mBuf[mPos])))
            ++mPos;
        if (mPos < mBuf.size() && mBuf[mPos] == '\n')
            ++mPos;
    }

    // Read `count` values of dtype `dt`, ascii or big-endian binary.
    NDArray ReadValues(DType dt, std::size_t count, bool is_ascii) {
        NDArray a = NDArray::Uninit(dt, {count});  // every element written below
        const std::size_t isz = dtype_size(dt);
        if (is_ascii) {
            const bool flt = detail::is_float_dtype(dt);
            const char* base = mBuf.c_str();
            for (std::size_t i = 0; i < count; ++i) {
                char* endp = nullptr;
                if (flt) {
                    double x = std::strtod(base + mPos, &endp);
                    if (endp == base + mPos)
                        throw ReadError("VTK ascii parse error");
                    store(a, i, x, 0);
                } else {
                    long long x = std::strtoll(base + mPos, &endp, 10);
                    if (endp == base + mPos)
                        throw ReadError("VTK ascii parse error");
                    store(a, i, 0.0, static_cast<std::int64_t>(x));
                }
                mPos = static_cast<std::size_t>(endp - base);
            }
        } else {
            if (mPos + count * isz > mBuf.size())
                throw ReadError("VTK binary truncated");
            char* out = reinterpret_cast<char*>(a.Data());
            // Element offsets are i*isz -> byte-swap in parallel (bswap intrinsic).
            const char* src = mBuf.data() + mPos;
            const int w = static_cast<int>(isz);
            parallel_for_bw(
                count, [&](std::size_t i) { detail::bswap_copy(out + i * isz, src + i * isz, w); });
            mPos += count * isz;
        }
        ConsumeEol();
        return a;
    }
};

std::vector<std::string> split(const std::string& rS) {
    std::vector<std::string> out;
    std::istringstream iss(rS);
    std::string tok;
    while (iss >> tok)
        out.push_back(tok);
    return out;
}

std::string upper(std::string s) {
    for (auto& c : s)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

std::vector<std::int64_t> to_int64(const NDArray& rA) {
    std::vector<std::int64_t> v(rA.Size());
    std::int64_t* dst = v.data();
    // Hoist the per-element dtype switch out of the loop, then bulk-convert.
    detail::dispatch_dtype(rA.Dtype(), [&]<class T>() {
        const T* src = rA.As<T>();
        parallel_for_bw(rA.Size(),
                        [&](std::size_t i) { dst[i] = static_cast<std::int64_t>(src[i]); });
    });
    return v;
}

}  // namespace

Mesh read_vtk(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    // Bulk slurp (seek+read) rather than char-by-char istreambuf_iterator.
    in.seekg(0, std::ios::end);
    std::streamoff len = in.tellg();
    in.seekg(0, std::ios::beg);
    std::string buf;
    if (len > 0) {
        buf.resize(static_cast<std::size_t>(len));
        in.read(buf.data(), len);
    }
    Cursor cur(buf);

    std::string header = cur.ReadLine();
    const bool is_v5 = header.find("Version 5") != std::string::npos;
    cur.ReadLine();  // title
    std::string dtype_line = upper(cur.ReadLine());
    bool is_ascii;
    if (dtype_line.find("ASCII") != std::string::npos)
        is_ascii = true;
    else if (dtype_line.find("BINARY") != std::string::npos)
        is_ascii = false;
    else
        throw ReadError("Unknown VTK data type line: " + dtype_line);

    Mesh mesh;
    std::vector<std::int64_t> conn, offsets, types;
    // Held alive so reconstruct_cells can read the int64 connectivity buffer
    // directly (VTK 5.1), skipping a to_int64 copy of the whole connectivity.
    NDArray conn_nd;
    const std::int64_t* conn_ptr = nullptr;
    bool conn_owned = false;  // conn_nd owns the int64 connectivity (VTK 5.1)
    std::unordered_map<std::string, NDArray> cell_data_raw;
    std::string active;  // POINT_DATA or CELL_DATA

    while (!cur.Eof()) {
        std::string line = cur.ReadLine();
        if (line.empty())
            continue;
        std::vector<std::string> tok = split(line);
        if (tok.empty())
            continue;
        std::string section = upper(tok[0]);

        if (section == "DATASET") {
            if (tok.size() < 2 || upper(tok[1]) != "UNSTRUCTURED_GRID")
                throw ReadError("C++ VTK reader only handles UNSTRUCTURED_GRID");
        } else if (section == "POINTS") {
            std::size_t n = std::stoull(tok[1]);
            DType dt = dtype_from_vtk_token(tok[2]);
            NDArray pts = cur.ReadValues(dt, n * 3, is_ascii);
            pts.Reshape({n, 3});
            mesh.AssignPoints(std::move(pts));
        } else if (section == "CELLS") {
            if (is_v5) {
                std::size_t num_off = std::stoull(tok[1]);
                std::size_t num_idx = std::stoull(tok[2]);
                std::string l = cur.ReadLine();
                if (upper(l).rfind("OFFSETS", 0) != 0)
                    throw ReadError("Expected OFFSETS (VTK 5.1 layout)");
                DType odt = dtype_from_vtk_token(split(l)[1]);
                std::vector<std::int64_t> off_all =
                    to_int64(cur.ReadValues(odt, num_off, is_ascii));
                l = cur.ReadLine();
                if (upper(l).rfind("CONNECTIVITY", 0) != 0)
                    throw ReadError("Expected CONNECTIVITY");
                DType cdt = dtype_from_vtk_token(split(l)[1]);
                conn_nd = cur.ReadValues(cdt, num_idx, is_ascii);
                if (conn_nd.Dtype() == DType::Int64) {
                    // Already int64 (vtktypeint64) -> read the buffer directly.
                    conn_ptr = conn_nd.As<std::int64_t>();
                    conn_owned = true;
                } else {
                    conn = to_int64(conn_nd);
                    conn_ptr = conn.data();
                }
                // off_all has a leading 0; end-offsets are the remainder.
                offsets.assign(off_all.begin() + 1, off_all.end());
            } else {
                // Version 4.2: interleaved [count, nodes...]; int32 values.
                std::size_t num_cells = std::stoull(tok[1]);
                std::size_t total = std::stoull(tok[2]);
                DType dt = is_ascii ? DType::Int64 : DType::Int32;
                std::vector<std::int64_t> raw = to_int64(cur.ReadValues(dt, total, is_ascii));
                conn.reserve(total - num_cells);
                offsets.reserve(num_cells);
                std::size_t p = 0;
                std::int64_t running = 0;
                for (std::size_t i = 0; i < num_cells; ++i) {
                    std::int64_t n = raw[p++];
                    for (std::int64_t j = 0; j < n; ++j)
                        conn.push_back(raw[p++]);
                    running += n;
                    offsets.push_back(running);
                }
                conn_ptr = conn.data();
            }
        } else if (section == "CELL_TYPES") {
            std::size_t n = std::stoull(tok[1]);
            DType dt = is_ascii ? DType::Int64 : DType::Int32;
            types = to_int64(cur.ReadValues(dt, n, is_ascii));
        } else if (section == "POINT_DATA") {
            active = "POINT_DATA";
        } else if (section == "CELL_DATA") {
            active = "CELL_DATA";
        } else if (section == "FIELD") {
            std::size_t k = std::stoull(tok[2]);
            for (std::size_t fi = 0; fi < k; ++fi) {
                std::vector<std::string> ft = split(cur.ReadLine());
                if (!ft.empty() && upper(ft[0]) == "METADATA") {
                    while (true) {
                        std::string ml = cur.ReadLine();
                        bool blank = true;
                        for (char c : ml)
                            if (!std::isspace(static_cast<unsigned char>(c)))
                                blank = false;
                        if (blank)
                            break;
                    }
                    ft = split(cur.ReadLine());
                }
                std::string name = ft[0];
                std::size_t ncomp = std::stoull(ft[1]);
                std::size_t ntuples = std::stoull(ft[2]);
                DType dt = dtype_from_vtk_token(ft[3]);
                NDArray arr = cur.ReadValues(dt, ncomp * ntuples, is_ascii);
                if (ncomp != 1)
                    arr.Reshape({ntuples, ncomp});
                if (active == "POINT_DATA")
                    mesh.AddPointData(name, std::move(arr));
                else
                    cell_data_raw.emplace(name, std::move(arr));
            }
        } else if (section == "METADATA") {
            while (true) {
                std::string ml = cur.ReadLine();
                bool blank = true;
                for (char c : ml)
                    if (!std::isspace(static_cast<unsigned char>(c)))
                        blank = false;
                if (blank || cur.Eof())
                    break;
            }
        } else {
            throw ReadError("VTK section '" + section + "' not supported by the C++ reader");
        }
    }

    // Fast path (zero copy): a single cell type spanning all cells, non-special,
    // with an identity VTK->meshio node order and regular end-offsets
    // (offsets[i] == (i+1)*n) means the owning int64 connectivity NDArray is
    // already the block data -> reshape and move it straight into the cell block
    // instead of gathering a fresh copy.
    bool moved = false;
    if (conn_owned && !types.empty()) {
        const int vt = static_cast<int>(types[0]);
        bool single = true;
        for (std::size_t i = 1; i < types.size(); ++i)
            if (types[i] != types[0]) {
                single = false;
                break;
            }
        const auto& tmap = vtk_to_meshio_type();
        auto it = tmap.find(vt);
        if (single && it != tmap.end() && !is_special_cell(it->second) &&
            vtk_to_meshio_order(vt).empty()) {
            auto nit = num_nodes_per_cell().find(it->second);
            if (nit != num_nodes_per_cell().end()) {
                const std::size_t n = static_cast<std::size_t>(nit->second);
                const std::size_t ncells = types.size();
                bool regular = conn_nd.Size() == ncells * n && offsets.size() == ncells;
                for (std::size_t i = 0; regular && i < ncells; ++i)
                    if (offsets[i] != static_cast<std::int64_t>((i + 1) * n))
                        regular = false;
                if (regular) {
                    conn_nd.Reshape({ncells, n});
                    mesh.AddCellBlock(it->second, std::move(conn_nd));
                    for (auto& kv : cell_data_raw)
                        mesh.AppendCellData(kv.first, std::move(kv.second));
                    moved = true;
                }
            }
        }
    }
    if (!moved)
        detail::reconstruct_cells(conn_ptr, offsets, types, cell_data_raw, mesh);
    return mesh;
}

}  // namespace meshioplusplus
