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

/**
 * @file registry.cpp
 * @brief The shared format-dispatch tables (see registry.hpp). Bodies hoisted
 *        verbatim from `bindings_js/js_bindings.cpp`, extended with the
 *        HDF5/netCDF-conditional entries native (non-WASM) builds can serve.
 */

// Project includes
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/abaqus.hpp"
#include "meshioplusplus/formats/ansys.hpp"
#include "meshioplusplus/formats/ansysinp.hpp"
#include "meshioplusplus/formats/avsucd.hpp"
#include "meshioplusplus/formats/cgns.hpp"
#include "meshioplusplus/formats/dex.hpp"
#include "meshioplusplus/formats/dolfin.hpp"
#include "meshioplusplus/formats/exodus.hpp"
#include "meshioplusplus/formats/flac3d.hpp"
#include "meshioplusplus/formats/flux.hpp"
#include "meshioplusplus/formats/freefem.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/h5m.hpp"
#include "meshioplusplus/formats/hmf.hpp"
#include "meshioplusplus/formats/ip.hpp"
#include "meshioplusplus/formats/med.hpp"
#include "meshioplusplus/formats/medit.hpp"
#include "meshioplusplus/formats/mff.hpp"
#include "meshioplusplus/formats/mfm.hpp"
#include "meshioplusplus/formats/mphtxt.hpp"
#include "meshioplusplus/formats/nastran.hpp"
#include "meshioplusplus/formats/netgen.hpp"
#include "meshioplusplus/formats/obj_off.hpp"
#include "meshioplusplus/formats/openfoam.hpp"
#include "meshioplusplus/formats/permas.hpp"
#include "meshioplusplus/formats/ply.hpp"
#include "meshioplusplus/formats/stl.hpp"
#include "meshioplusplus/formats/su2.hpp"
#include "meshioplusplus/formats/tecplot.hpp"
#include "meshioplusplus/formats/tetgen.hpp"
#include "meshioplusplus/formats/ugrid.hpp"
#include "meshioplusplus/formats/unv.hpp"
#include "meshioplusplus/formats/vtk.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/formats/wkt.hpp"
#include "meshioplusplus/formats/xdmf.hpp"

namespace meshioplusplus {

const std::map<std::string, ReadFn>& registry_readers() {
    static const std::map<std::string, ReadFn> m = {
        {"abaqus", meshioplusplus::read_abaqus},
        {"ansys", meshioplusplus::read_ansys},
        {"avsucd", meshioplusplus::read_avsucd},
        {"dolfin", meshioplusplus::read_dolfin},
        {"flac3d", meshioplusplus::read_flac3d},
        {"dex", meshioplusplus::read_dex},
        {"flux", meshioplusplus::read_flux},
        {"freefem", meshioplusplus::read_freefem},
        {"gmsh", meshioplusplus::read_gmsh},
        {"ip", meshioplusplus::read_ip},
        {"medit", meshioplusplus::read_medit_ascii},
        {"mff", meshioplusplus::read_mff},
        {"mfm", meshioplusplus::read_mfm},
        {"mphtxt", meshioplusplus::read_mphtxt},
        {"nastran", meshioplusplus::read_nastran},
        {"netgen", meshioplusplus::read_netgen},
        {"obj", meshioplusplus::read_obj},
        {"off", meshioplusplus::read_off},
        {"permas", meshioplusplus::read_permas},
        {"ply", meshioplusplus::read_ply},
        {"stl", meshioplusplus::read_stl},
        {"su2", meshioplusplus::read_su2},
        {"tecplot", meshioplusplus::read_tecplot},
        {"tetgen", meshioplusplus::read_tetgen},
        {"ugrid", meshioplusplus::read_ugrid},
        {"unv", [](const std::string& path) { return meshioplusplus::read_unv(path); }},
        {"vtk", meshioplusplus::read_vtk},
        {"vtu", meshioplusplus::read_vtu},
        {"wkt", meshioplusplus::read_wkt},
        {"xdmf", meshioplusplus::read_xdmf},
        // Side-channel info (point_sets/cell_sets, cell-tag family names) is
        // not carried by the flat bindings -- v1 limitation, see doc/wasm.md
        // and doc/c_api.md.
        {"ansysinp",
         [](const std::string& path) {
             meshioplusplus::AnsysInfo info;
             return meshioplusplus::read_ansysinp(path, info);
         }},
        {"openfoam",
         [](const std::string& path) {
             meshioplusplus::OpenFoamInfo info;
             return meshioplusplus::read_openfoam(path, info);
         }},
#ifdef MESHIOPLUSPLUS_HAS_HDF5
        {"cgns", meshioplusplus::read_cgns},
        {"h5m", meshioplusplus::read_h5m},
        {"hmf", meshioplusplus::read_hmf},
        {"med",
         [](const std::string& path) {
             meshioplusplus::MedInfo info;  // families/tags side channel dropped in v1
             return meshioplusplus::read_med(path, info);
         }},
#endif
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
        {"exodus", meshioplusplus::read_exodus},
#endif
    };
    return m;
}

const std::map<std::string, WriteFn>& registry_writers() {
    static const std::map<std::string, WriteFn> m = {
        {"abaqus", meshioplusplus::write_abaqus},
        {"ansys", [](const std::string& p,
                     const Mesh& mm) { meshioplusplus::write_ansys(p, mm, /*binary=*/true); }},
        {"avsucd", meshioplusplus::write_avsucd},
        {"dolfin", meshioplusplus::write_dolfin},
        {"flac3d",
         [](const std::string& p, const Mesh& mm) {
             meshioplusplus::write_flac3d(p, mm, ".16e", /*binary=*/false);
         }},
        {"dex", meshioplusplus::write_dex},
        {"flux", meshioplusplus::write_flux},
        {"freefem", meshioplusplus::write_freefem},
        {"gmsh", [](const std::string& p,
                    const Mesh& mm) { meshioplusplus::write_gmsh41(p, mm, /*binary=*/true); }},
        {"ip", meshioplusplus::write_ip},
        {"medit", meshioplusplus::write_medit_ascii},
        {"mff", meshioplusplus::write_mff},
        {"mfm",
         [](const std::string& p, const Mesh& mm) { meshioplusplus::write_mfm(p, mm, ".16e"); }},
        {"mphtxt", meshioplusplus::write_mphtxt},
        {"nastran", meshioplusplus::write_nastran},
        {"netgen",
         [](const std::string& p, const Mesh& mm) { meshioplusplus::write_netgen(p, mm, ".16e"); }},
        {"obj", meshioplusplus::write_obj},
        {"off", meshioplusplus::write_off},
        {"permas", meshioplusplus::write_permas},
        {"ply", [](const std::string& p,
                   const Mesh& mm) { meshioplusplus::write_ply(p, mm, /*binary=*/true); }},
        {"stl", [](const std::string& p,
                   const Mesh& mm) { meshioplusplus::write_stl(p, mm, /*binary=*/false); }},
        {"su2", meshioplusplus::write_su2},
        {"tecplot", meshioplusplus::write_tecplot},
        {"tetgen", meshioplusplus::write_tetgen},
        {"ugrid", meshioplusplus::write_ugrid},
        {"unv", [](const std::string& p, const Mesh& mm) { meshioplusplus::write_unv(p, mm); }},
        {"vtk",
         [](const std::string& p, const Mesh& mm) {
             meshioplusplus::write_vtk(p, mm, /*binary=*/true, /*v51=*/true);
         }},
        {"vtu",
         [](const std::string& p, const Mesh& mm) {
             meshioplusplus::write_vtu(p, mm, /*binary=*/true, /*zlib=*/true);
         }},
        {"wkt", meshioplusplus::write_wkt},
        // XDMF's heavy-data format follows the build: HDF companion file when
        // HDF5 is available (the Python writer's default), inline XML text
        // otherwise (the only always-available option; what WASM ships).
        {"xdmf",
         [](const std::string& p, const Mesh& mm) {
#ifdef MESHIOPLUSPLUS_HAS_HDF5
             meshioplusplus::write_xdmf(p, mm, "HDF");
#else
             meshioplusplus::write_xdmf(p, mm, "XML");
#endif
         }},
        {"ansysinp",
         [](const std::string& p, const Mesh& mm) {
             meshioplusplus::AnsysInfo info;  // no point_sets/cell_sets side channel in v1
             meshioplusplus::write_ansysinp(p, mm, info);
         }},
    // openfoam is read-only in the C++ core (see openfoam.hpp) -> no writer entry.
#ifdef MESHIOPLUSPLUS_HAS_HDF5
        {"cgns", [](const std::string& p,
                    const Mesh& mm) { meshioplusplus::write_cgns(p, mm, /*gzip_level=*/4); }},
        {"h5m",
         [](const std::string& p, const Mesh& mm) {
             meshioplusplus::write_h5m(p, mm, /*add_global_ids=*/true, /*gzip_level=*/4);
         }},
        {"hmf", [](const std::string& p,
                   const Mesh& mm) { meshioplusplus::write_hmf(p, mm, /*gzip_level=*/4); }},
        {"med",
         [](const std::string& p, const Mesh& mm) {
             meshioplusplus::MedInfo info;  // families/tags side channel dropped in v1
             meshioplusplus::write_med(p, mm, info);
         }},
#endif
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
        {"exodus", meshioplusplus::write_exodus},
#endif
    };
    return m;
}

// Extension -> canonical format key for the non-ambiguous cases; `.msh`
// defaults to gmsh and `.inp` to abaqus (matching this repo's own import
// order in src/meshioplusplus/__init__.py). Pass an explicit `format` to
// select ansys/freefem (.msh) or ansysinp (.inp) instead. Optional-dependency
// extensions are mapped even in builds where the format is compiled out, so
// the resulting error names the missing dependency (registry_compiled_out())
// rather than claiming the extension is unknown.
const std::map<std::string, std::string>& registry_extension_defaults() {
    static const std::map<std::string, std::string> m = {
        {".inp", "abaqus"},  {".avs", "avsucd"},  {".xml", "dolfin"},    {".f3grid", "flac3d"},
        {".dex", "dex"},     {".ip", "ip"},       {".mff", "mff"},       {".pf3", "flux"},
        {".mesh", "medit"},  {".mfm", "mfm"},     {".mphtxt", "mphtxt"}, {".bdf", "nastran"},
        {".nas", "nastran"}, {".fem", "nastran"}, {".vol", "netgen"},    {".obj", "obj"},
        {".off", "off"},     {".post", "permas"}, {".dato", "permas"},   {".ply", "ply"},
        {".stl", "stl"},     {".su2", "su2"},     {".dat", "tecplot"},   {".tec", "tecplot"},
        {".ele", "tetgen"},  {".node", "tetgen"}, {".ugrid", "ugrid"},   {".unv", "unv"},
        {".vtk", "vtk"},     {".vtu", "vtu"},     {".wkt", "wkt"},       {".xdmf", "xdmf"},
        {".xmf", "xdmf"},    {".msh", "gmsh"},    {".cgns", "cgns"},     {".h5m", "h5m"},
        {".hmf", "hmf"},     {".med", "med"},     {".e", "exodus"},      {".exo", "exodus"},
        {".ex2", "exodus"},
    };
    return m;
}

namespace {

std::string extension_of(const std::string& rPath) {
    auto pos = rPath.find_last_of('.');
    return pos == std::string::npos ? "" : rPath.substr(pos);
}

}  // namespace

std::string resolve_format(const std::string& rPath, const std::string& rFormat) {
    if (!rFormat.empty())
        return rFormat;
    auto it = registry_extension_defaults().find(extension_of(rPath));
    if (it == registry_extension_defaults().end())
        throw meshioplusplus::ReadError("meshio++: cannot infer format from '" + rPath +
                                        "' -- pass an explicit format argument");
    return it->second;
}

const char* registry_compiled_out(const std::string& rFormat) {
#ifndef MESHIOPLUSPLUS_HAS_HDF5
    if (rFormat == "cgns" || rFormat == "h5m" || rFormat == "hmf" || rFormat == "med")
        return "HDF5";
#endif
#ifndef MESHIOPLUSPLUS_HAS_NETCDF
    if (rFormat == "exodus")
        return "netCDF";
#endif
    return nullptr;
}

}  // namespace meshioplusplus
