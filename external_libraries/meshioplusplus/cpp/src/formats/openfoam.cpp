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
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/openfoam.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

namespace fs = std::filesystem;

namespace meshioplusplus {

namespace {

using Face = std::vector<std::int64_t>;

struct FoamFormat {
    bool mBinary = false;
    int mLabelBytes = 8;
    int mScalarBytes = 8;
};

std::string read_whole(const std::string& rPath) {
    std::ifstream f(rPath, std::ios::binary);
    if (!f)
        throw ReadError("Could not open OpenFOAM file: " + rPath);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

std::string strip(const std::string& rS) {
    std::size_t a = rS.find_first_not_of(" \t\r\n");
    if (a == std::string::npos)
        return "";
    std::size_t b = rS.find_last_not_of(" \t\r\n");
    return rS.substr(a, b - a + 1);
}

// Parse the FoamFile header for format/arch (label/scalar byte widths).
FoamFormat detect_format(const std::string& rPath) {
    FoamFormat fmt;
    std::ifstream f(rPath, std::ios::binary);
    if (!f)
        return fmt;
    std::string line;
    while (std::getline(f, line)) {
        std::string s = strip(line);
        // format <word>;
        std::size_t p = s.find("format");
        if (p == 0) {
            std::string rest = strip(s.substr(6));
            if (!rest.empty() && rest.back() == ';')
                rest.pop_back();
            rest = strip(rest);
            if (rest == "binary")
                fmt.mBinary = true;
            else if (rest == "ascii")
                fmt.mBinary = false;
        }
        if (s.rfind("arch", 0) == 0) {
            std::size_t lp = s.find("label=");
            if (lp != std::string::npos) {
                int bits = std::atoi(s.c_str() + lp + 6);
                if (bits)
                    fmt.mLabelBytes = bits / 8;
            }
            std::size_t sp = s.find("scalar=");
            if (sp != std::string::npos) {
                int bits = std::atoi(s.c_str() + sp + 7);
                if (bits)
                    fmt.mScalarBytes = bits / 8;
            }
        }
        if (s == "}")
            break;
    }
    return fmt;
}

// Strip C-style /* */ and // comments and drop the FoamFile { ... } block.
std::string strip_comments_and_header(const std::string& rText) {
    std::string out;
    out.reserve(rText.size());
    // remove /* */ and //
    for (std::size_t i = 0; i < rText.size();) {
        if (i + 1 < rText.size() && rText[i] == '/' && rText[i + 1] == '*') {
            std::size_t e = rText.find("*/", i + 2);
            i = (e == std::string::npos) ? rText.size() : e + 2;
        } else if (i + 1 < rText.size() && rText[i] == '/' && rText[i + 1] == '/') {
            std::size_t e = rText.find('\n', i + 2);
            i = (e == std::string::npos) ? rText.size() : e;
        } else {
            out.push_back(rText[i++]);
        }
    }
    // drop FoamFile { ... }
    std::istringstream ss(out);
    std::string line, result;
    bool in_header = false;
    int depth = 0;
    while (std::getline(ss, line)) {
        std::string s = strip(line);
        if (s.find("FoamFile") != std::string::npos)
            in_header = true;
        if (in_header) {
            for (char c : s) {
                if (c == '{')
                    ++depth;
                else if (c == '}')
                    --depth;
            }
            if (depth <= 0)
                in_header = false;
            continue;
        }
        result += line;
        result.push_back('\n');
    }
    return result;
}

// ---- ASCII parsers ----

std::vector<std::array<double, 3>> parse_points_ascii(const std::string& rBody) {
    std::vector<std::array<double, 3>> pts;
    std::istringstream ss(rBody);
    std::string line;
    bool in_block = false;
    bool have_n = false;
    while (std::getline(ss, line)) {
        std::string s = strip(line);
        if (s.empty())
            continue;
        if (!have_n && s.find_first_not_of("0123456789") == std::string::npos) {
            have_n = true;
            continue;
        }
        if (s == "(" && have_n) {
            in_block = true;
            continue;
        }
        if (s == ")" && in_block)
            break;
        if (in_block) {
            // extract up to 3 numbers from within parentheses
            std::string t = s;
            for (char& c : t)
                if (c == '(' || c == ')')
                    c = ' ';
            std::istringstream ns(t);
            double a, b, c;
            if (ns >> a >> b >> c)
                pts.push_back({a, b, c});
        }
    }
    return pts;
}

std::vector<Face> parse_faces_ascii(const std::string& rBody) {
    std::vector<Face> faces;
    std::istringstream ss(rBody);
    std::string line;
    bool in_block = false, have_n = false;
    while (std::getline(ss, line)) {
        std::string s = strip(line);
        if (s.empty())
            continue;
        if (!have_n && s.find_first_not_of("0123456789") == std::string::npos) {
            have_n = true;
            continue;
        }
        if (s == "(" && have_n) {
            in_block = true;
            continue;
        }
        if (s == ")" && in_block)
            break;
        if (in_block) {
            // form: <count>(<ids...>)
            std::size_t lp = s.find('(');
            std::size_t rp = s.find(')', lp);
            if (lp == std::string::npos || rp == std::string::npos)
                continue;
            std::string inside = s.substr(lp + 1, rp - lp - 1);
            std::istringstream ns(inside);
            Face f;
            std::int64_t v;
            while (ns >> v)
                f.push_back(v);
            faces.push_back(std::move(f));
        }
    }
    return faces;
}

std::vector<std::int64_t> parse_int_list_ascii(const std::string& rBody) {
    std::vector<std::int64_t> out;
    std::istringstream ss(rBody);
    std::string line;
    bool in_block = false, have_n = false;
    while (std::getline(ss, line)) {
        std::string s = strip(line);
        if (s.empty())
            continue;
        if (!have_n && s.find_first_not_of("0123456789") == std::string::npos) {
            have_n = true;
            continue;
        }
        if (s == "(") {
            in_block = true;
            continue;
        }
        if (s == ")")
            break;
        if (in_block) {
            std::istringstream ns(s);
            std::int64_t v;
            while (ns >> v)
                out.push_back(v);
        }
    }
    return out;
}

// Boundary patch descriptor. `mNFaces`/`mStartFace` deliberately mirror
// OpenFOAM's own on-disk `boundary` field names (`nFaces`/`startFace`).
struct Patch {
    std::string mName;
    std::int64_t mNFaces = 0;
    std::int64_t mStartFace = 0;
};

std::vector<Patch> parse_boundary(const std::string& rBody) {
    // Find `name { ... }` blocks with nFaces/startFace.
    std::vector<Patch> patches;
    std::size_t i = 0, n = rBody.size();
    auto skip_ws = [&](std::size_t& p) {
        while (p < n && std::isspace(static_cast<unsigned char>(rBody[p])))
            ++p;
    };
    while (i < n) {
        skip_ws(i);
        // read a token (patch name)
        std::size_t start = i;
        while (i < n && !std::isspace(static_cast<unsigned char>(rBody[i])) && rBody[i] != '{' &&
               rBody[i] != '(' && rBody[i] != ')')
            ++i;
        std::string name = rBody.substr(start, i - start);
        skip_ws(i);
        if (i < n && rBody[i] == '{') {
            std::size_t close = rBody.find('}', i);
            if (close == std::string::npos)
                break;
            std::string block = rBody.substr(i + 1, close - i - 1);
            Patch pt;
            pt.mName = name;
            bool has_n = false, has_s = false;
            std::size_t np = block.find("nFaces");
            if (np != std::string::npos) {
                pt.mNFaces = std::atoll(block.c_str() + np + 6);
                has_n = true;
            }
            std::size_t sp = block.find("startFace");
            if (sp != std::string::npos) {
                pt.mStartFace = std::atoll(block.c_str() + sp + 9);
                has_s = true;
            }
            if (has_n && has_s && !name.empty())
                patches.push_back(pt);
            i = close + 1;
        } else if (i < n && (rBody[i] == '(' || rBody[i] == ')')) {
            ++i;  // skip list delimiters
        } else if (name.empty()) {
            ++i;
        }
    }
    return patches;
}

// ---- binary parsers ----

// Return (N, offset just after the outer '(').
std::pair<std::int64_t, std::size_t> data_start(const std::string& rRaw) {
    std::size_t end = rRaw.find('}');
    if (end == std::string::npos)
        throw ReadError("OpenFOAM: no FoamFile header");
    std::size_t lp = rRaw.find('(', end);
    if (lp == std::string::npos)
        throw ReadError("OpenFOAM: no data list '('");
    // last integer between end and lp
    std::int64_t n = 0;
    bool found = false;
    std::size_t i = end;
    while (i < lp) {
        if (std::isdigit(static_cast<unsigned char>(rRaw[i]))) {
            std::int64_t v = 0;
            while (i < lp && std::isdigit(static_cast<unsigned char>(rRaw[i])))
                v = v * 10 + (rRaw[i++] - '0');
            n = v;
            found = true;
        } else {
            ++i;
        }
    }
    if (!found)
        throw ReadError("OpenFOAM: no element count before '('");
    return {n, lp + 1};
}

template <typename T>
T read_le(const char* pP) {
    T v;
    std::memcpy(&v, pP, sizeof(T));
    return v;
}

std::vector<std::array<double, 3>> read_binary_points(const std::string& rRaw, int scalar_bytes) {
    auto [n, start] = data_start(rRaw);
    std::vector<std::array<double, 3>> pts(static_cast<std::size_t>(n));
    const char* base = rRaw.data() + start;
    for (std::int64_t i = 0; i < n; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::size_t off =
                (static_cast<std::size_t>(i) * 3 + j) * static_cast<std::size_t>(scalar_bytes);
            pts[i][j] = scalar_bytes == 4 ? static_cast<double>(read_le<float>(base + off))
                                          : read_le<double>(base + off);
        }
    }
    return pts;
}

std::vector<std::int64_t> read_binary_labels(const std::string& rRaw, int label_bytes) {
    auto [n, start] = data_start(rRaw);
    std::vector<std::int64_t> out(static_cast<std::size_t>(n));
    const char* base = rRaw.data() + start;
    for (std::int64_t i = 0; i < n; ++i) {
        std::size_t off = static_cast<std::size_t>(i) * static_cast<std::size_t>(label_bytes);
        out[i] = label_bytes == 4 ? static_cast<std::int64_t>(read_le<std::int32_t>(base + off))
                                  : read_le<std::int64_t>(base + off);
    }
    return out;
}

std::vector<Face> read_binary_faces(const std::string& rRaw, int label_bytes) {
    auto [nfaces, pos] = data_start(rRaw);
    std::vector<Face> faces(static_cast<std::size_t>(nfaces));
    std::size_t p = pos;
    for (std::int64_t i = 0; i < nfaces; ++i) {
        std::size_t lp = rRaw.find('(', p);
        if (lp == std::string::npos)
            throw ReadError("OpenFOAM: missing '(' in faces");
        std::int64_t count = std::atoll(rRaw.substr(p, lp - p).c_str());
        std::size_t blob = lp + 1;
        Face f(static_cast<std::size_t>(count));
        for (std::int64_t j = 0; j < count; ++j) {
            std::size_t off =
                blob + static_cast<std::size_t>(j) * static_cast<std::size_t>(label_bytes);
            f[j] = label_bytes == 4
                       ? static_cast<std::int64_t>(read_le<std::int32_t>(rRaw.data() + off))
                       : read_le<std::int64_t>(rRaw.data() + off);
        }
        faces[i] = std::move(f);
        p = blob + static_cast<std::size_t>(count) * static_cast<std::size_t>(label_bytes) + 1;
    }
    return faces;
}

// ---- dispatch readers ----

std::vector<std::array<double, 3>> read_points(const fs::path& rPath) {
    FoamFormat fmt = detect_format(rPath.string());
    std::string raw = read_whole(rPath.string());
    if (fmt.mBinary)
        return read_binary_points(raw, fmt.mScalarBytes);
    return parse_points_ascii(strip_comments_and_header(raw));
}

std::vector<Face> read_faces(const fs::path& rPath) {
    FoamFormat fmt = detect_format(rPath.string());
    std::string raw = read_whole(rPath.string());
    if (fmt.mBinary)
        return read_binary_faces(raw, fmt.mLabelBytes);
    return parse_faces_ascii(strip_comments_and_header(raw));
}

std::vector<std::int64_t> read_int_list(const fs::path& rPath) {
    FoamFormat fmt = detect_format(rPath.string());
    std::string raw = read_whole(rPath.string());
    if (fmt.mBinary)
        return read_binary_labels(raw, fmt.mLabelBytes);
    return parse_int_list_ascii(strip_comments_and_header(raw));
}

// ---- geometry ----

double triple(const std::array<double, 3>& rA, const std::array<double, 3>& rB,
              const std::array<double, 3>& rC) {
    // a . (b x c)
    double cx = rB[1] * rC[2] - rB[2] * rC[1];
    double cy = rB[2] * rC[0] - rB[0] * rC[2];
    double cz = rB[0] * rC[1] - rB[1] * rC[0];
    return rA[0] * cx + rA[1] * cy + rA[2] * cz;
}

std::array<double, 3> sub(const std::array<double, 3>& rA, const std::array<double, 3>& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}

std::size_t unique_node_count(const std::vector<Face>& rFaces) {
    std::unordered_set<std::int64_t> s;
    for (const auto& f : rFaces)
        for (std::int64_t v : f)
            s.insert(v);
    return s.size();
}

std::unordered_map<std::int64_t, std::unordered_set<std::int64_t>> node_adjacency(
    const std::vector<Face>& rFaces) {
    std::unordered_map<std::int64_t, std::unordered_set<std::int64_t>> adj;
    for (const auto& f : rFaces) {
        std::size_t m = f.size();
        for (std::size_t i = 0; i < m; ++i) {
            std::int64_t a = f[i], b = f[(i + 1) % m];
            adj[a].insert(b);
            adj[b].insert(a);
        }
    }
    return adj;
}

// Returns the ordered top ring, or empty if ambiguous.
std::vector<std::int64_t> match_top(const Face& rBottom, const std::vector<Face>& rOriented) {
    auto adj = node_adjacency(rOriented);
    std::unordered_set<std::int64_t> base(rBottom.begin(), rBottom.end());
    std::vector<std::int64_t> top;
    for (std::int64_t b : rBottom) {
        std::vector<std::int64_t> cand;
        for (std::int64_t x : adj[b])
            if (!base.count(x))
                cand.push_back(x);
        if (cand.size() != 1)
            return {};
        top.push_back(cand[0]);
    }
    return top;
}

using P3 = std::vector<std::array<double, 3>>;

Face build_tetra(const std::vector<Face>& rOriented, const P3& rP) {
    const Face& base = rOriented[0];
    std::unordered_set<std::int64_t> all;
    for (const auto& f : rOriented)
        for (std::int64_t v : f)
            all.insert(v);
    for (std::int64_t v : base)
        all.erase(v);
    std::int64_t apex = *all.begin();
    Face n = {base[0], base[1], base[2], apex};
    if (triple(sub(rP[n[1]], rP[n[0]]), sub(rP[n[2]], rP[n[0]]), sub(rP[n[3]], rP[n[0]])) < 0)
        n = {base[0], base[2], base[1], apex};
    return n;
}

Face build_pyramid(const std::vector<Face>& rOriented, const P3& rP) {
    Face quad;
    for (const auto& f : rOriented)
        if (f.size() == 4) {
            quad = f;
            break;
        }
    std::unordered_set<std::int64_t> all;
    for (const auto& f : rOriented)
        for (std::int64_t v : f)
            all.insert(v);
    for (std::int64_t v : quad)
        all.erase(v);
    std::int64_t apex = *all.begin();
    Face n = {quad[0], quad[1], quad[2], quad[3], apex};
    if (triple(sub(rP[n[1]], rP[n[0]]), sub(rP[n[3]], rP[n[0]]), sub(rP[n[4]], rP[n[0]])) < 0)
        n = {quad[0], quad[3], quad[2], quad[1], apex};
    return n;
}

Face build_wedge(const std::vector<Face>& rOriented, const P3& rP) {
    Face bottom;
    for (const auto& f : rOriented)
        if (f.size() == 3) {
            bottom = f;
            break;
        }
    std::vector<std::int64_t> top = match_top(bottom, rOriented);
    if (top.empty())
        return {};
    Face n = {bottom[0], bottom[1], bottom[2], top[0], top[1], top[2]};
    if (triple(sub(rP[n[1]], rP[n[0]]), sub(rP[n[2]], rP[n[0]]), sub(rP[n[3]], rP[n[0]])) < 0)
        n = {bottom[0], bottom[2], bottom[1], top[0], top[2], top[1]};
    return n;
}

Face build_hexahedron(const std::vector<Face>& rOriented, const P3& rP) {
    Face bottom;
    for (const auto& f : rOriented)
        if (f.size() == 4) {
            bottom = f;
            break;
        }
    std::vector<std::int64_t> top = match_top(bottom, rOriented);
    if (top.empty())
        return {};
    Face n = {bottom[0], bottom[1], bottom[2], bottom[3], top[0], top[1], top[2], top[3]};
    if (triple(sub(rP[n[1]], rP[n[0]]), sub(rP[n[3]], rP[n[0]]), sub(rP[n[4]], rP[n[0]])) < 0)
        n = {bottom[0], bottom[3], bottom[2], bottom[1], top[0], top[3], top[2], top[1]};
    return n;
}

// Classify a cell. Returns {meshio type, connectivity}. For "polyhedron" the
// connectivity is empty (the caller keeps the oriented faces).
std::pair<std::string, Face> reconstruct_cell(const std::vector<Face>& rOriented, const P3& rP) {
    std::size_t nf = rOriented.size();
    std::size_t np = unique_node_count(rOriented);
    if (nf == 4 && np == 4)
        return {"tetra", build_tetra(rOriented, rP)};
    if (nf == 5 && np == 5)
        return {"pyramid", build_pyramid(rOriented, rP)};
    if (nf == 5 && np == 6)
        return {"wedge", build_wedge(rOriented, rP)};
    if (nf == 6 && np == 8)
        return {"hexahedron", build_hexahedron(rOriented, rP)};
    return {"polyhedron", {}};
}

}  // namespace

Mesh read_openfoam(const std::string& rPathIn, OpenFoamInfo& rInfo) {
    // resolve polyMesh directory
    fs::path path(rPathIn);
    fs::path poly;
    if (path.extension() == ".foam") {
        fs::path c = path.parent_path() / "constant" / "polyMesh";
        if (fs::exists(c))
            poly = c;
    }
    if (poly.empty() && path.filename() == "polyMesh" && fs::is_directory(path))
        poly = path;
    if (poly.empty()) {
        for (const fs::path& c : {path / "constant" / "polyMesh", path / "polyMesh"}) {
            if (fs::exists(c)) {
                poly = c;
                break;
            }
        }
    }
    if (poly.empty())
        throw ReadError(detail::format_compat(
            "Could not locate polyMesh from '{}'. Expected <case>/constant/polyMesh/.", rPathIn));
    log::info("Reading polyMesh from {}", poly.string());

    P3 points = read_points(poly / "points");
    std::vector<Face> faces = read_faces(poly / "faces");
    std::vector<std::int64_t> owner = read_int_list(poly / "owner");
    std::vector<std::int64_t> neighbour;
    if (fs::exists(poly / "neighbour"))
        neighbour = read_int_list(poly / "neighbour");
    std::vector<Patch> boundary;
    if (fs::exists(poly / "boundary"))
        boundary =
            parse_boundary(strip_comments_and_header(read_whole((poly / "boundary").string())));

    std::int64_t owner_max = -1, neigh_max = -1;
    for (std::int64_t v : owner)
        owner_max = std::max(owner_max, v);
    for (std::int64_t v : neighbour)
        neigh_max = std::max(neigh_max, v);
    std::int64_t n_cells = owner.empty() ? 0 : std::max(owner_max, neigh_max) + 1;
    log::info("{} points, {} faces, {} cells, {} patches", points.size(), faces.size(), n_cells,
              boundary.size());

    // cell -> face ids
    std::vector<std::vector<std::int64_t>> cell_faces(static_cast<std::size_t>(n_cells));
    for (std::size_t fid = 0; fid < owner.size(); ++fid)
        cell_faces[static_cast<std::size_t>(owner[fid])].push_back(static_cast<std::int64_t>(fid));
    for (std::size_t fid = 0; fid < neighbour.size(); ++fid)
        if (neighbour[fid] >= 0)
            cell_faces[static_cast<std::size_t>(neighbour[fid])].push_back(
                static_cast<std::int64_t>(fid));

    // reconstruct volume cells
    std::vector<std::string> vol_order;
    std::map<std::string, std::vector<Face>> vol_buckets;
    // polyhedra grouped by unique node count -> "polyhedron<N>"
    std::vector<std::string> poly_order;
    std::map<std::string, std::vector<std::vector<Face>>> poly_buckets;

    // Per-cell geometric reconstruction is the expensive part and every cell
    // only reads faces/owner/points -> compute all cells in parallel into a
    // pre-sized result array, then do the (ordered) bucket grouping
    // sequentially.
    struct CellResult {
        std::string mType;         // "" = degenerate (skipped)
        Face mConn;                // named types
        std::vector<Face> mFaces;  // oriented faces, polyhedra only
    };
    std::vector<CellResult> results(static_cast<std::size_t>(n_cells));
    parallel_for(static_cast<std::size_t>(n_cells), [&](std::size_t cs) {
        const std::int64_t cid = static_cast<std::int64_t>(cs);
        std::vector<Face> oriented;
        for (std::int64_t fid : cell_faces[cs]) {
            Face f = faces[static_cast<std::size_t>(fid)];
            if (owner[static_cast<std::size_t>(fid)] != cid)
                std::reverse(f.begin(), f.end());
            oriented.push_back(std::move(f));
        }
        auto [mtype, conn] = reconstruct_cell(oriented, points);
        if (mtype == "polyhedron") {
            results[cs] = {"polyhedron", {}, std::move(oriented)};
        } else if (conn.empty()) {
            results[cs] = {};  // degenerate topology
        } else {
            results[cs] = {std::move(mtype), std::move(conn), {}};
        }
    });

    std::size_t n_skipped = 0;
    std::size_t n_polyhedra = 0;
    for (auto& res : results) {
        if (res.mType == "polyhedron") {
            std::size_t nn = unique_node_count(res.mFaces);
            std::string key = "polyhedron" + std::to_string(nn);
            if (!poly_buckets.count(key))
                poly_order.push_back(key);
            poly_buckets[key].push_back(std::move(res.mFaces));
            ++n_polyhedra;
        } else if (res.mType.empty()) {
            ++n_skipped;
        } else {
            if (!vol_buckets.count(res.mType))
                vol_order.push_back(res.mType);
            vol_buckets[res.mType].push_back(std::move(res.mConn));
        }
    }
    if (n_skipped > 0)
        log::warn("{} cell(s) skipped (degenerate topology).", n_skipped);
    if (n_polyhedra > 0)
        log::info("{} general polyhedron cell(s) found.", n_polyhedra);

    Mesh mesh;
    std::size_t npts = points.size();
    {
        NDArray pts(DType::Float64, {npts, 3});
        double* pdst = pts.As<double>();
        parallel_for(npts, [&](std::size_t i) {
            for (std::size_t j = 0; j < 3; ++j)
                pdst[i * 3 + j] = points[i][j];
        });
        mesh.AssignPoints(std::move(pts));
    }

    std::vector<NDArray> cell_tags;  // one per block, in final block order

    // rectangular volume blocks
    for (const std::string& t : vol_order) {
        const auto& rows = vol_buckets[t];
        std::size_t nc = rows.size();
        std::size_t k = nc ? rows[0].size() : 0;
        NDArray data(DType::Int64, {nc, k});
        std::int64_t* dp = data.As<std::int64_t>();
        parallel_for(nc, [&](std::size_t r) {
            for (std::size_t c = 0; c < k; ++c)
                dp[r * k + c] = rows[r][c];
        });
        mesh.AddCellBlock(t, std::move(data));
        cell_tags.emplace_back(DType::Int64, std::vector<std::size_t>{nc});  // zeros
    }
    // ragged polyhedron blocks
    for (const std::string& key : poly_order) {
        std::vector<std::vector<std::vector<std::int64_t>>> cells;
        for (const auto& cell : poly_buckets[key]) {
            std::vector<std::vector<std::int64_t>> ph;
            for (const auto& face : cell)
                ph.push_back(face);
            cells.push_back(std::move(ph));
        }
        std::size_t nc = cells.size();
        mesh.AddPolyhedronBlock(key, std::move(cells));
        cell_tags.emplace_back(DType::Int64, std::vector<std::size_t>{nc});  // zeros
    }

    // boundary cells grouped by size, with patch family tags
    std::map<int, std::vector<Face>> bysize;  // 3 -> triangles, 4 -> quads
    std::map<int, std::vector<std::int64_t>> tagsize;
    std::vector<Face> poly_faces;
    std::vector<std::int64_t> poly_tags;
    for (std::size_t pidx = 0; pidx < boundary.size(); ++pidx) {
        std::int64_t fam = -(static_cast<std::int64_t>(pidx) + 1);
        rInfo.mCellTags[fam] = {boundary[pidx].mName};
        for (std::int64_t fid = boundary[pidx].mStartFace;
             fid < boundary[pidx].mStartFace + boundary[pidx].mNFaces; ++fid) {
            if (fid < 0 || static_cast<std::size_t>(fid) >= faces.size())
                continue;
            const Face& f = faces[static_cast<std::size_t>(fid)];
            if (f.size() == 3) {
                bysize[3].push_back(f);
                tagsize[3].push_back(fam);
            } else if (f.size() == 4) {
                bysize[4].push_back(f);
                tagsize[4].push_back(fam);
            } else {
                poly_faces.push_back(f);
                poly_tags.push_back(fam);
            }
        }
    }
    auto add_boundary_block = [&](const std::string& type, const std::vector<Face>& rows,
                                  const std::vector<std::int64_t>& tags) {
        std::size_t nc = rows.size();
        std::size_t k = nc ? rows[0].size() : 0;
        NDArray data(DType::Int64, {nc, k});
        NDArray tag(DType::Int64, {nc});
        std::int64_t* dp = data.As<std::int64_t>();
        std::int64_t* tp = tag.As<std::int64_t>();
        parallel_for(nc, [&](std::size_t r) {
            for (std::size_t c = 0; c < k; ++c)
                dp[r * k + c] = rows[r][c];
            tp[r] = tags[r];
        });
        mesh.AddCellBlock(type, std::move(data));
        cell_tags.push_back(std::move(tag));
    };
    if (!bysize[3].empty())
        add_boundary_block("triangle", bysize[3], tagsize[3]);
    if (!bysize[4].empty())
        add_boundary_block("quad", bysize[4], tagsize[4]);
    if (!poly_faces.empty()) {
        // group boundary polygons by vertex count -> polygon<N>
        std::map<std::size_t, std::vector<Face>> by_n;
        std::map<std::size_t, std::vector<std::int64_t>> tag_n;
        for (std::size_t i = 0; i < poly_faces.size(); ++i) {
            by_n[poly_faces[i].size()].push_back(poly_faces[i]);
            tag_n[poly_faces[i].size()].push_back(poly_tags[i]);
        }
        for (auto& kv : by_n)
            add_boundary_block("polygon" + std::to_string(kv.first), kv.second, tag_n[kv.first]);
    }

    if (!cell_tags.empty())
        mesh.AddCellData("cell_tags", std::move(cell_tags));
    return mesh;
}

}  // namespace meshioplusplus
