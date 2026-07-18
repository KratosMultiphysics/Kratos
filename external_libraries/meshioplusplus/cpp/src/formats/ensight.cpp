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
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/ensight.hpp"
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

// ---------------------------------------------------------------------------
// Shared tables
// ---------------------------------------------------------------------------

// EnSight Gold element keyword <-> meshio cell type (fixed-node types only;
// nsided/nfaced are handled separately as ragged blocks).
struct EnsightTypeEntry {
    const char* mKeyword;
    const char* mMeshioType;
    int mNumNodes;
};

const std::vector<EnsightTypeEntry>& ensight_type_table() {
    static const std::vector<EnsightTypeEntry> table = {
        {"point", "vertex", 1},         {"bar2", "line", 2},
        {"bar3", "line3", 3},           {"tria3", "triangle", 3},
        {"tria6", "triangle6", 6},      {"quad4", "quad", 4},
        {"quad8", "quad8", 8},          {"tetra4", "tetra", 4},
        {"tetra10", "tetra10", 10},     {"pyramid5", "pyramid", 5},
        {"pyramid13", "pyramid13", 13}, {"penta6", "wedge", 6},
        {"penta15", "wedge15", 15},     {"hexa8", "hexahedron", 8},
        {"hexa20", "hexahedron20", 20}};
    return table;
}

const EnsightTypeEntry* ensight_entry_from_keyword(const std::string& rKeyword) {
    static const std::unordered_map<std::string, const EnsightTypeEntry*> map = [] {
        std::unordered_map<std::string, const EnsightTypeEntry*> m;
        for (const auto& e : ensight_type_table())
            m.emplace(e.mKeyword, &e);
        return m;
    }();
    auto it = map.find(rKeyword);
    return it == map.end() ? nullptr : it->second;
}

const EnsightTypeEntry* ensight_entry_from_meshio(const std::string& rType) {
    static const std::unordered_map<std::string, const EnsightTypeEntry*> map = [] {
        std::unordered_map<std::string, const EnsightTypeEntry*> m;
        for (const auto& e : ensight_type_table())
            m.emplace(e.mMeshioType, &e);
        return m;
    }();
    auto it = map.find(rType);
    return it == map.end() ? nullptr : it->second;
}

// meshio <-> EnSight node-order permutation. meshio ordering is VTK ordering
// for every supported type; the only difference is the prism triangle
// winding of penta15 vs wedge15 — the same involution VTK's EnSight readers
// apply, so one table serves both directions (result[j] = source index for
// position j). Kratos's penta15/hexa20 group swaps fix Kratos-specific
// ordering and must NOT be ported here.
const std::vector<int>* ensight_permutation(const std::string& rMeshioType) {
    static const std::vector<int> wedge15 = {0, 2, 1, 3, 5, 4, 8, 7, 6, 11, 10, 9, 12, 14, 13};
    if (rMeshioType == "wedge15")
        return &wedge15;
    return nullptr;
}

// ---------------------------------------------------------------------------
// Small string / path helpers
// ---------------------------------------------------------------------------

bool ensight_has_suffix(const std::string& rPath, const std::string& rSuffix) {
    if (rPath.size() < rSuffix.size())
        return false;
    return rPath.compare(rPath.size() - rSuffix.size(), rSuffix.size(), rSuffix) == 0;
}

std::string ensight_trim(const std::string& rLine) {
    std::size_t b = 0, e = rLine.size();
    while (b < e && (std::isspace(static_cast<unsigned char>(rLine[b])) || rLine[b] == '\0'))
        ++b;
    while (e > b &&
           (std::isspace(static_cast<unsigned char>(rLine[e - 1])) || rLine[e - 1] == '\0'))
        --e;
    return rLine.substr(b, e - b);
}

bool ensight_starts_with(const std::string& rStr, const char* pPrefix) {
    return rStr.rfind(pPrefix, 0) == 0;
}

std::string ensight_dirname(const std::string& rPath) {
    std::size_t slash = rPath.find_last_of("/\\");
    return slash == std::string::npos ? std::string() : rPath.substr(0, slash + 1);
}

std::string ensight_basename(const std::string& rPath) {
    std::size_t slash = rPath.find_last_of("/\\");
    return slash == std::string::npos ? rPath : rPath.substr(slash + 1);
}

// "<stem>.case" or "<stem>.geo" -> {case path, geo path}.
std::pair<std::string, std::string> ensight_case_geo_paths(const std::string& rPath, bool& rOk) {
    rOk = true;
    if (ensight_has_suffix(rPath, ".case"))
        return {rPath, rPath.substr(0, rPath.size() - 5) + ".geo"};
    if (ensight_has_suffix(rPath, ".geo"))
        return {rPath.substr(0, rPath.size() - 4) + ".case", rPath};
    rOk = false;
    return {"", ""};
}

std::string ensight_read_whole_file(const std::string& rPath, const char* pWhat) {
    std::ifstream in(rPath, std::ios::binary);
    if (!in)
        throw ReadError(std::string("EnSight: could not open ") + pWhat + ": " + rPath);
    in.seekg(0, std::ios::end);
    const std::streamoff size = in.tellg();
    in.seekg(0, std::ios::beg);
    std::string data(static_cast<std::size_t>(size < 0 ? 0 : size), '\0');
    if (!data.empty())
        in.read(data.data(), static_cast<std::streamsize>(data.size()));
    return data;
}

// ---------------------------------------------------------------------------
// Geometry cursors: one record/number stream over ASCII or C-binary data
// ---------------------------------------------------------------------------

class EnsightCursor {
public:
    virtual ~EnsightCursor() = default;
    /// True once only trailing whitespace/padding remains.
    virtual bool AtEnd() = 0;
    /// Consume and return the next string record (trimmed line / 80-char record).
    virtual std::string NextRecord() = 0;
    /// Return the next string record without consuming it ("" at end).
    virtual std::string PeekRecord() = 0;
    virtual std::int64_t NextInt() = 0;
    virtual void ReadInts(std::size_t n, std::int64_t* pDst) = 0;
    virtual void ReadFloats(std::size_t n, double* pDst) = 0;
    virtual void SkipInts(std::size_t n) = 0;
    /// Binary-only hook: peek the next int32 and enable byte-swapping when it
    /// is implausible as-is but plausible swapped. With PreferSmaller (used at
    /// the tiny part-number record, where e.g. bswap(1) = 16777216 is still
    /// "plausible"), the smaller of two plausible interpretations wins.
    /// No-op for ASCII.
    virtual void CheckSwap(std::int64_t Max, bool PreferSmaller) {
        (void)Max;
        (void)PreferSmaller;
    }
};

class EnsightAsciiCursor final : public EnsightCursor {
public:
    explicit EnsightAsciiCursor(std::string Text) : mText(std::move(Text)) {}

    bool AtEnd() override {
        std::size_t p = mPos;
        while (p < mText.size() &&
               (std::isspace(static_cast<unsigned char>(mText[p])) || mText[p] == '\0'))
            ++p;
        return p >= mText.size();
    }

    std::string NextRecord() override {
        while (mPos < mText.size()) {
            std::size_t eol = mText.find('\n', mPos);
            if (eol == std::string::npos)
                eol = mText.size();
            std::string line = ensight_trim(mText.substr(mPos, eol - mPos));
            mPos = eol < mText.size() ? eol + 1 : eol;
            if (!line.empty())
                return line;
        }
        throw ReadError("EnSight: unexpected end of geometry file");
    }

    std::string PeekRecord() override {
        if (AtEnd())
            return "";
        const std::size_t saved = mPos;
        std::string rec = NextRecord();
        mPos = saved;
        return rec;
    }

    std::int64_t NextInt() override {
        const char* start = mText.c_str() + mPos;
        char* end = nullptr;
        const std::int64_t v = std::strtoll(start, &end, 10);
        if (end == start)
            throw ReadError("EnSight: expected an integer in geometry file");
        mPos = static_cast<std::size_t>(end - mText.c_str());
        return v;
    }

    void ReadInts(std::size_t n, std::int64_t* pDst) override {
        for (std::size_t i = 0; i < n; ++i)
            pDst[i] = NextInt();
    }

    void ReadFloats(std::size_t n, double* pDst) override {
        for (std::size_t i = 0; i < n; ++i) {
            const char* start = mText.c_str() + mPos;
            char* end = nullptr;
            pDst[i] = std::strtod(start, &end);
            if (end == start)
                throw ReadError("EnSight: expected a number in geometry file");
            mPos = static_cast<std::size_t>(end - mText.c_str());
        }
    }

    void SkipInts(std::size_t n) override {
        for (std::size_t i = 0; i < n; ++i)
            NextInt();
    }

private:
    std::string mText;
    std::size_t mPos = 0;
};

class EnsightBinaryCursor final : public EnsightCursor {
public:
    // Text is the whole file; the cursor starts after the leading
    // "C Binary" 80-char record.
    explicit EnsightBinaryCursor(std::string Data) : mData(std::move(Data)), mPos(80) {}

    bool AtEnd() override { return mPos >= mData.size(); }

    std::string NextRecord() override {
        if (mPos + 80 > mData.size())
            throw ReadError("EnSight: truncated binary geometry file");
        std::string rec(mData.data() + mPos, 80);
        mPos += 80;
        const std::size_t nul = rec.find('\0');
        if (nul != std::string::npos)
            rec.resize(nul);
        return ensight_trim(rec);
    }

    std::string PeekRecord() override {
        if (mPos + 80 > mData.size())
            return "";
        const std::size_t saved = mPos;
        std::string rec = NextRecord();
        mPos = saved;
        return rec;
    }

    std::int64_t NextInt() override {
        std::int32_t v;
        Require(4);
        std::memcpy(&v, mData.data() + mPos, 4);
        mPos += 4;
        if (mSwap)
            detail::bswap_inplace(reinterpret_cast<char*>(&v), 4);
        return v;
    }

    void ReadInts(std::size_t n, std::int64_t* pDst) override {
        Require(4 * n);
        const char* src = mData.data() + mPos;
        for (std::size_t i = 0; i < n; ++i) {
            std::int32_t v;
            std::memcpy(&v, src + 4 * i, 4);
            if (mSwap)
                detail::bswap_inplace(reinterpret_cast<char*>(&v), 4);
            pDst[i] = v;
        }
        mPos += 4 * n;
    }

    void ReadFloats(std::size_t n, double* pDst) override {
        Require(4 * n);
        const char* src = mData.data() + mPos;
        for (std::size_t i = 0; i < n; ++i) {
            float v;
            std::memcpy(&v, src + 4 * i, 4);
            if (mSwap)
                detail::bswap_inplace(reinterpret_cast<char*>(&v), 4);
            pDst[i] = v;
        }
        mPos += 4 * n;
    }

    void SkipInts(std::size_t n) override {
        Require(4 * n);
        mPos += 4 * n;
    }

    void CheckSwap(std::int64_t Max, bool PreferSmaller) override {
        if (mPos + 4 > mData.size())
            return;
        std::int32_t v;
        std::memcpy(&v, mData.data() + mPos, 4);
        if (mSwap)
            detail::bswap_inplace(reinterpret_cast<char*>(&v), 4);
        std::int32_t s = v;
        detail::bswap_inplace(reinterpret_cast<char*>(&s), 4);
        const bool v_ok = v >= 0 && v <= Max;
        const bool s_ok = s >= 0 && s <= Max;
        if ((!v_ok && s_ok) || (PreferSmaller && v_ok && s_ok && s < v))
            mSwap = !mSwap;
    }

private:
    void Require(std::size_t n) {
        if (mPos + n > mData.size())
            throw ReadError("EnSight: truncated binary geometry file");
    }

    std::string mData;
    std::size_t mPos;
    bool mSwap = false;
};

// ---------------------------------------------------------------------------
// .case parsing
// ---------------------------------------------------------------------------

// Parse the .case file and return the resolved geometry file path.
std::string ensight_parse_case(const std::string& rCasePath) {
    const std::string data = ensight_read_whole_file(rCasePath, "case file");

    std::string section;
    std::string format_type;
    std::string model_value;
    std::istringstream stream(data);
    std::string raw;
    while (std::getline(stream, raw)) {
        std::string line = ensight_trim(raw);
        if (line.empty() || line[0] == '#')
            continue;
        if (line == "FORMAT" || line == "GEOMETRY" || line == "VARIABLE" || line == "TIME" ||
            line == "FILE" || line == "MATERIAL" || line == "SCRIPTS") {
            section = line;
            continue;
        }
        if (section == "FORMAT" && ensight_starts_with(line, "type:"))
            format_type = ensight_trim(line.substr(5));
        else if (section == "GEOMETRY" && ensight_starts_with(line, "model:"))
            model_value = ensight_trim(line.substr(6));
    }

    if (format_type.find("ensight gold") == std::string::npos)
        throw ReadError("EnSight: case file is not 'type: ensight gold' (got '" + format_type +
                        "')");
    if (model_value.empty())
        throw ReadError("EnSight: case file has no GEOMETRY 'model:' entry");

    // model: [ts] [fs] filename [change_coords_only] — drop leading integer
    // timeset/fileset tokens, take the first remaining token as the filename.
    std::istringstream toks(model_value);
    std::vector<std::string> tokens;
    std::string tok;
    while (toks >> tok)
        tokens.push_back(tok);
    std::size_t first = 0;
    while (first < tokens.size()) {
        char* end = nullptr;
        (void)std::strtoll(tokens[first].c_str(), &end, 10);
        if (end == tokens[first].c_str() || *end != '\0')
            break;  // not a pure integer
        ++first;
    }
    if (first >= tokens.size())
        throw ReadError("EnSight: malformed 'model:' line in case file");
    const std::string& filename = tokens[first];
    if (filename.find('*') != std::string::npos)
        throw ReadError("EnSight: transient (wildcard) geometry is not supported");

    return ensight_dirname(rCasePath) + filename;
}

// ---------------------------------------------------------------------------
// Geometry reading
// ---------------------------------------------------------------------------

// A staged cell block (added to the mesh only after AssignPoints).
struct EnsightBlock {
    std::string mType;
    NDArray mConn{DType::Int64, {}};                                       // rectangular
    std::vector<std::vector<std::int64_t>> mPolygonRows;                   // nsided
    std::vector<std::vector<std::vector<std::int64_t>>> mPolyhedronCells;  // nfaced
    int mKind = 0;  // 0 rectangular, 1 polygon, 2 polyhedron
    std::int64_t mPartId = 0;
    std::size_t mNumCells = 0;
};

// "given" and "ignore" both put id arrays in the file; only the presence
// matters — Gold connectivity is positional, so ids are always skipped.
bool ensight_ids_in_file(const std::string& rRecord, const char* pWhat) {
    // rRecord is e.g. "node id assign"; the mode is the last token.
    std::istringstream iss(rRecord);
    std::string tok, mode;
    while (iss >> tok)
        mode = tok;
    if (mode == "given" || mode == "ignore")
        return true;
    if (mode == "off" || mode == "assign")
        return false;
    throw ReadError(std::string("EnSight: malformed '") + pWhat + " id' record: " + rRecord);
}

Mesh ensight_parse_geo(EnsightCursor& rCur) {
    constexpr std::int64_t plausible_max = 100000000;  // generous id/count bound

    rCur.NextRecord();  // description line 1
    rCur.NextRecord();  // description line 2
    std::string node_id_rec = rCur.NextRecord();
    if (!ensight_starts_with(node_id_rec, "node id"))
        throw ReadError("EnSight: expected 'node id' record, got: " + node_id_rec);
    std::string elem_id_rec = rCur.NextRecord();
    if (!ensight_starts_with(elem_id_rec, "element id"))
        throw ReadError("EnSight: expected 'element id' record, got: " + elem_id_rec);
    const bool node_ids_in_file = ensight_ids_in_file(node_id_rec, "node");
    const bool elem_ids_in_file = ensight_ids_in_file(elem_id_rec, "element");

    if (ensight_starts_with(rCur.PeekRecord(), "extents")) {
        rCur.NextRecord();
        double extents[6];
        rCur.ReadFloats(6, extents);
    }

    std::vector<double> coords;  // xyz-interleaved, all parts concatenated
    std::vector<EnsightBlock> blocks;
    std::int64_t num_parts = 0;

    while (!rCur.AtEnd()) {
        std::string rec = rCur.NextRecord();
        if (!ensight_starts_with(rec, "part"))
            throw ReadError("EnSight: expected 'part' record, got: " + rec);
        rCur.CheckSwap(plausible_max, /*PreferSmaller=*/true);
        const std::int64_t part_id = rCur.NextInt();
        ++num_parts;
        rCur.NextRecord();  // part description

        rec = rCur.NextRecord();
        if (!ensight_starts_with(rec, "coordinates"))
            throw ReadError("EnSight: expected 'coordinates' record, got: " + rec);
        rCur.CheckSwap(plausible_max, /*PreferSmaller=*/false);
        const std::int64_t nn = rCur.NextInt();
        if (nn < 0)
            throw ReadError("EnSight: negative node count");
        const std::int64_t point_offset = static_cast<std::int64_t>(coords.size() / 3);

        if (node_ids_in_file)
            rCur.SkipInts(static_cast<std::size_t>(nn));
        std::vector<double> x(static_cast<std::size_t>(nn));
        std::vector<double> y(static_cast<std::size_t>(nn));
        std::vector<double> z(static_cast<std::size_t>(nn));
        rCur.ReadFloats(static_cast<std::size_t>(nn), x.data());
        rCur.ReadFloats(static_cast<std::size_t>(nn), y.data());
        rCur.ReadFloats(static_cast<std::size_t>(nn), z.data());
        coords.reserve(coords.size() + static_cast<std::size_t>(nn) * 3);
        for (std::int64_t i = 0; i < nn; ++i) {
            coords.push_back(x[static_cast<std::size_t>(i)]);
            coords.push_back(y[static_cast<std::size_t>(i)]);
            coords.push_back(z[static_cast<std::size_t>(i)]);
        }

        // 1-based positional index within this part -> global 0-based index.
        auto resolve = [&](std::int64_t v) -> std::int64_t {
            const std::int64_t local = v - 1;
            if (local < 0 || local >= nn)
                throw ReadError("EnSight: connectivity index out of range");
            return local + point_offset;
        };

        // Element sections until the next part / EOF.
        while (!rCur.AtEnd()) {
            std::string kw = rCur.PeekRecord();
            if (kw.empty() || ensight_starts_with(kw, "part"))
                break;
            rCur.NextRecord();
            // Ghost-cell sections carry the same data as their base type.
            if (ensight_starts_with(kw, "g_"))
                kw = kw.substr(2);

            rCur.CheckSwap(plausible_max, /*PreferSmaller=*/false);
            const std::int64_t ne = rCur.NextInt();
            if (ne < 0)
                throw ReadError("EnSight: negative element count");
            if (elem_ids_in_file)
                rCur.SkipInts(static_cast<std::size_t>(ne));

            if (kw == "nsided") {
                std::vector<std::int64_t> sizes(static_cast<std::size_t>(ne));
                rCur.ReadInts(sizes.size(), sizes.data());
                std::vector<std::vector<std::int64_t>> rows(static_cast<std::size_t>(ne));
                std::vector<std::int64_t> flat;
                std::size_t total = 0;
                for (auto s : sizes)
                    total += static_cast<std::size_t>(s);
                flat.resize(total);
                rCur.ReadInts(total, flat.data());
                std::size_t at = 0;
                for (std::size_t c = 0; c < rows.size(); ++c) {
                    rows[c].resize(static_cast<std::size_t>(sizes[c]));
                    for (std::size_t j = 0; j < rows[c].size(); ++j)
                        rows[c][j] = resolve(flat[at++]);
                }
                EnsightBlock b;
                b.mType = "polygon";
                b.mKind = 1;
                b.mPartId = part_id;
                b.mNumCells = rows.size();
                b.mPolygonRows = std::move(rows);
                blocks.push_back(std::move(b));
            } else if (kw == "nfaced") {
                std::vector<std::int64_t> nfaces(static_cast<std::size_t>(ne));
                rCur.ReadInts(nfaces.size(), nfaces.data());
                std::size_t total_faces = 0;
                for (auto f : nfaces)
                    total_faces += static_cast<std::size_t>(f);
                std::vector<std::int64_t> fsizes(total_faces);
                rCur.ReadInts(total_faces, fsizes.data());
                std::size_t total_nodes = 0;
                for (auto s : fsizes)
                    total_nodes += static_cast<std::size_t>(s);
                std::vector<std::int64_t> flat(total_nodes);
                rCur.ReadInts(total_nodes, flat.data());

                std::vector<std::vector<std::vector<std::int64_t>>> cells(
                    static_cast<std::size_t>(ne));
                std::size_t face_at = 0, node_at = 0;
                for (std::size_t c = 0; c < cells.size(); ++c) {
                    cells[c].resize(static_cast<std::size_t>(nfaces[c]));
                    for (auto& face : cells[c]) {
                        face.resize(static_cast<std::size_t>(fsizes[face_at++]));
                        for (auto& v : face)
                            v = resolve(flat[node_at++]);
                    }
                }
                // Group by unique node count into "polyhedron<N>" blocks (the
                // openfoam convention), preserving first-seen order.
                std::vector<std::size_t> node_counts(cells.size());
                for (std::size_t c = 0; c < cells.size(); ++c) {
                    std::vector<std::int64_t> uniq;
                    for (const auto& face : cells[c])
                        uniq.insert(uniq.end(), face.begin(), face.end());
                    std::sort(uniq.begin(), uniq.end());
                    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
                    node_counts[c] = uniq.size();
                }
                std::vector<std::size_t> group_order;
                std::map<std::size_t, std::vector<std::size_t>> groups;
                for (std::size_t c = 0; c < cells.size(); ++c) {
                    if (groups.find(node_counts[c]) == groups.end())
                        group_order.push_back(node_counts[c]);
                    groups[node_counts[c]].push_back(c);
                }
                for (std::size_t n : group_order) {
                    std::vector<std::vector<std::vector<std::int64_t>>> group_cells;
                    for (std::size_t c : groups[n])
                        group_cells.push_back(std::move(cells[c]));
                    EnsightBlock b;
                    b.mType = "polyhedron" + std::to_string(n);
                    b.mKind = 2;
                    b.mPartId = part_id;
                    b.mNumCells = group_cells.size();
                    b.mPolyhedronCells = std::move(group_cells);
                    blocks.push_back(std::move(b));
                }
            } else {
                const EnsightTypeEntry* entry = ensight_entry_from_keyword(kw);
                if (entry == nullptr)
                    throw ReadError("EnSight: unsupported element keyword: " + kw);
                const std::size_t npc = static_cast<std::size_t>(entry->mNumNodes);
                NDArray conn = NDArray::Uninit(DType::Int64, {static_cast<std::size_t>(ne), npc});
                std::int64_t* cp = conn.As<std::int64_t>();
                rCur.ReadInts(static_cast<std::size_t>(ne) * npc, cp);
                const std::vector<int>* perm = ensight_permutation(entry->mMeshioType);
                if (perm != nullptr) {
                    std::vector<std::int64_t> tmp(npc);
                    for (std::int64_t r = 0; r < ne; ++r) {
                        std::int64_t* row = cp + r * static_cast<std::int64_t>(npc);
                        for (std::size_t j = 0; j < npc; ++j)
                            tmp[j] = row[(*perm)[j]];
                        for (std::size_t j = 0; j < npc; ++j)
                            row[j] = tmp[j];
                    }
                }
                const std::size_t nvals = static_cast<std::size_t>(ne) * npc;
                for (std::size_t k = 0; k < nvals; ++k)
                    cp[k] = resolve(cp[k]);
                EnsightBlock b;
                b.mType = entry->mMeshioType;
                b.mKind = 0;
                b.mPartId = part_id;
                b.mNumCells = static_cast<std::size_t>(ne);
                b.mConn = std::move(conn);
                blocks.push_back(std::move(b));
            }
        }
    }

    Mesh mesh;
    const std::size_t npoints = coords.size() / 3;
    NDArray pts = NDArray::Uninit(DType::Float64, {npoints, 3});
    if (npoints > 0)
        std::memcpy(pts.Data(), coords.data(), npoints * 3 * sizeof(double));
    mesh.AssignPoints(std::move(pts));

    for (auto& b : blocks) {
        if (b.mKind == 0)
            mesh.AddCellBlock(b.mType, std::move(b.mConn));
        else if (b.mKind == 1)
            mesh.AddPolygonBlock(b.mType, std::move(b.mPolygonRows));
        else
            mesh.AddPolyhedronBlock(b.mType, std::move(b.mPolyhedronCells));
    }

    if (num_parts >= 2) {
        std::vector<NDArray> tags;
        tags.reserve(blocks.size());
        for (const auto& b : blocks) {
            NDArray a = NDArray::Uninit(DType::Int64, {b.mNumCells});
            std::int64_t* ap = a.As<std::int64_t>();
            for (std::size_t i = 0; i < b.mNumCells; ++i)
                ap[i] = b.mPartId;
            tags.push_back(std::move(a));
        }
        mesh.AddCellData("ensight:part", std::move(tags));
    }

    return mesh;
}

}  // namespace

Mesh read_ensight(const std::string& rPath) {
    std::string geo_path = rPath;
    if (ensight_has_suffix(rPath, ".case"))
        geo_path = ensight_parse_case(rPath);

    std::string data = ensight_read_whole_file(geo_path, "geometry file");
    if (ensight_starts_with(data, "Fortran Binary"))
        throw ReadError("EnSight: Fortran-binary geometry files are not supported");
    if (data.size() >= 80 && ensight_starts_with(data, "C Binary")) {
        EnsightBinaryCursor cur(std::move(data));
        return ensight_parse_geo(cur);
    }
    EnsightAsciiCursor cur(std::move(data));
    return ensight_parse_geo(cur);
}

namespace {

// ---------------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------------

void ensight_append_str80(std::vector<char>& rOut, const std::string& rStr) {
    char buf[80] = {};
    rStr.copy(buf, std::min<std::size_t>(rStr.size(), 79));
    rOut.insert(rOut.end(), buf, buf + 80);
}

void ensight_append_i32(std::vector<char>& rOut, std::int64_t v) {
    const std::int32_t i = static_cast<std::int32_t>(v);
    const char* p = reinterpret_cast<const char*>(&i);
    rOut.insert(rOut.end(), p, p + 4);
}

// Validate the mesh and return one keyword entry per cell block.
std::vector<const EnsightTypeEntry*> ensight_writable_blocks(const Mesh& rMesh) {
    std::vector<const EnsightTypeEntry*> entries;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsRagged())
            throw WriteError("EnSight: writing nsided/nfaced (ragged) blocks is not supported");
        const EnsightTypeEntry* entry = ensight_entry_from_meshio(cb.Type());
        if (entry == nullptr)
            throw WriteError("EnSight: cell type '" + cb.Type() + "' has no EnSight keyword");
        entries.push_back(entry);
    }
    return entries;
}

void ensight_write_geo_ascii(std::ostream& rOs, const Mesh& rMesh,
                             const std::vector<const EnsightTypeEntry*>& rEntries) {
    const NDArray& points = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    const std::size_t np = rMesh.NumPoints();

    std::string out;
    out.reserve(200 + np * 42);
    out += "EnSight Gold Geometry File\n";
    out += "Written by meshio++\n";
    out += "node id assign\n";
    out += "element id assign\n";
    out += "part\n";

    char buf[64];
    std::snprintf(buf, sizeof(buf), "%10d\n", 1);
    out += buf;
    out += "Mesh\n";
    out += "coordinates\n";
    std::snprintf(buf, sizeof(buf), "%10lld\n", static_cast<long long>(np));
    out += buf;
    for (std::size_t c = 0; c < 3; ++c) {
        for (std::size_t i = 0; i < np; ++i) {
            const double v = c < dim ? detail::read_double(points, i * dim + c) : 0.0;
            std::snprintf(buf, sizeof(buf), "%12.5e\n", v);
            out += buf;
        }
    }

    for (std::size_t bi = 0; bi < rMesh.NumCellBlocks(); ++bi) {
        const auto cb = rMesh.Cells(bi);
        const EnsightTypeEntry* entry = rEntries[bi];
        const std::size_t npc = cb.NodesPerCell();
        const std::size_t ne = cb.NumCells();
        const NDArray& conn = cb.Conn();
        const std::vector<int>* perm = ensight_permutation(cb.Type());

        out += entry->mKeyword;
        out += "\n";
        std::snprintf(buf, sizeof(buf), "%10lld\n", static_cast<long long>(ne));
        out += buf;
        for (std::size_t r = 0; r < ne; ++r) {
            for (std::size_t j = 0; j < npc; ++j) {
                const std::size_t src = perm != nullptr ? static_cast<std::size_t>((*perm)[j]) : j;
                const long long v =
                    static_cast<long long>(detail::read_int(conn, r * npc + src)) + 1;
                std::snprintf(buf, sizeof(buf), "%10lld", v);
                out += buf;
            }
            out += "\n";
        }
    }

    rOs.write(out.data(), static_cast<std::streamsize>(out.size()));
}

void ensight_write_geo_binary(std::ostream& rOs, const Mesh& rMesh,
                              const std::vector<const EnsightTypeEntry*>& rEntries) {
    const NDArray& points = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    const std::size_t np = rMesh.NumPoints();

    constexpr std::size_t i32_max =
        static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max());
    if (np > i32_max)
        throw WriteError("EnSight: mesh too large for 32-bit binary EnSight output");

    std::vector<char> out;
    out.reserve(80 * 8 + np * 12 + 64);
    ensight_append_str80(out, "C Binary");
    ensight_append_str80(out, "EnSight Gold Geometry File");
    ensight_append_str80(out, "Written by meshio++");
    ensight_append_str80(out, "node id assign");
    ensight_append_str80(out, "element id assign");
    ensight_append_str80(out, "part");
    ensight_append_i32(out, 1);
    ensight_append_str80(out, "Mesh");
    ensight_append_str80(out, "coordinates");
    ensight_append_i32(out, static_cast<std::int64_t>(np));
    {
        std::vector<float> col(np * 3);
        for (std::size_t c = 0; c < 3; ++c)
            for (std::size_t i = 0; i < np; ++i)
                col[c * np + i] =
                    c < dim ? static_cast<float>(detail::read_double(points, i * dim + c)) : 0.0f;
        const char* p = reinterpret_cast<const char*>(col.data());
        out.insert(out.end(), p, p + col.size() * sizeof(float));
    }

    for (std::size_t bi = 0; bi < rMesh.NumCellBlocks(); ++bi) {
        const auto cb = rMesh.Cells(bi);
        const EnsightTypeEntry* entry = rEntries[bi];
        const std::size_t npc = cb.NodesPerCell();
        const std::size_t ne = cb.NumCells();
        if (ne > i32_max)
            throw WriteError("EnSight: mesh too large for 32-bit binary EnSight output");
        const NDArray& conn = cb.Conn();
        const std::vector<int>* perm = ensight_permutation(cb.Type());

        ensight_append_str80(out, entry->mKeyword);
        ensight_append_i32(out, static_cast<std::int64_t>(ne));
        std::vector<std::int32_t> flat(ne * npc);
        for (std::size_t r = 0; r < ne; ++r)
            for (std::size_t j = 0; j < npc; ++j) {
                const std::size_t src = perm != nullptr ? static_cast<std::size_t>((*perm)[j]) : j;
                flat[r * npc + j] =
                    static_cast<std::int32_t>(detail::read_int(conn, r * npc + src)) + 1;
            }
        const char* p = reinterpret_cast<const char*>(flat.data());
        out.insert(out.end(), p, p + flat.size() * sizeof(std::int32_t));
    }

    rOs.write(out.data(), static_cast<std::streamsize>(out.size()));
}

}  // namespace

void write_ensight(const std::string& rPath, const Mesh& rMesh, bool binary) {
    bool ok = false;
    auto paths = ensight_case_geo_paths(rPath, ok);
    if (!ok)
        throw WriteError("EnSight: must specify a .case or .geo file");
    const std::string& case_path = paths.first;
    const std::string& geo_path = paths.second;

    if (rMesh.PointDim() > 3)
        throw WriteError("EnSight: points must have at most three components");
    const std::vector<const EnsightTypeEntry*> entries = ensight_writable_blocks(rMesh);

    {
        std::ofstream cf(case_path, std::ios::binary);
        if (!cf)
            throw WriteError("Could not open file for writing: " + case_path);
        std::string out;
        out += "FORMAT\n";
        out += "type: ensight gold\n";
        out += "\n";
        out += "GEOMETRY\n";
        out += "model: " + ensight_basename(geo_path) + "\n";
        cf.write(out.data(), static_cast<std::streamsize>(out.size()));
    }

    std::ofstream gf(geo_path, std::ios::binary);
    if (!gf)
        throw WriteError("Could not open file for writing: " + geo_path);
    if (binary)
        ensight_write_geo_binary(gf, rMesh, entries);
    else
        ensight_write_geo_ascii(gf, rMesh, entries);
}

}  // namespace meshioplusplus
