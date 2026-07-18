//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/io.h"

namespace Kratos
{
class IntegrationValuesExtrapolationToNodesProcess; // forward declaration (gauss point extrapolation)

///@name Kratos Classes
///@{

/**
 * @brief Transient (temporal collection) XDMF writer.
 * @details Writes an XDMF v3 file whose static mesh grid is written once and
 * whose temporal collection grid is extended with one grid per solution step.
 * Appending re-parses the existing file, so a writer constructed on an
 * already-written time series continues it instead of overwriting it.
 * The produced layout matches the meshio `TimeSeriesWriter`/`TimeSeriesReader`
 * convention and is readable by ParaView.
 * @note This class is an implementation detail of @ref MeshioPlusPlusIO
 * (exposed here so it can also be used standalone from C++).
 */
class KRATOS_API(KRATOS_CORE) XdmfTimeSeriesWriter
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of XdmfTimeSeriesWriter
    KRATOS_CLASS_POINTER_DEFINITION(XdmfTimeSeriesWriter);

    /// A named data array (row-major, NumberOfRows x NumberOfComponents),
    /// used for point data (one row per node) and cell data (one row per cell)
    struct DataArray
    {
        std::string Name;
        std::size_t NumberOfComponents = 1;
        std::vector<double> Values;
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor.
     * @param rFileName The .xdmf file to create or extend.
     * @param rDataFormat How heavy data is stored: "XML" (inline text),
     *                    "Binary" (sibling .bin files) or "HDF" (sibling .h5
     *                    file, requires an HDF5-enabled build).
     * @param Precision Digits used for inline "XML" floating point data.
     */
    XdmfTimeSeriesWriter(
        const std::filesystem::path& rFileName,
        const std::string& rDataFormat,
        const std::size_t Precision = 16
        );

    /// Destructor.
    ~XdmfTimeSeriesWriter();

    /// Copying is disabled (the writer owns file state).
    XdmfTimeSeriesWriter(const XdmfTimeSeriesWriter&) = delete;
    XdmfTimeSeriesWriter& operator=(const XdmfTimeSeriesWriter&) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Whether the target file already holds a valid time-series mesh.
     * @details True when the file exists and contains the static mesh grid and
     * the temporal collection ("the buffer already exists"): appending steps
     * will extend it and @ref WriteMesh is not needed.
     */
    bool MeshWritten() const;

    /**
     * @brief Writes (or overwrites) the file with the static mesh and an empty
     * temporal collection.
     * @param rModelPart The model part providing nodes/elements/conditions.
     * @param WriteElements If false, elements are not written.
     * @param WriteConditions If false, conditions are not written.
     * @param WriteDeformedConfiguration If true the current node coordinates
     *                                   are written; otherwise the initial ones.
     */
    void WriteMesh(
        const ModelPart& rModelPart,
        const bool WriteElements = true,
        const bool WriteConditions = true,
        const bool WriteDeformedConfiguration = false
        );

    /**
     * @brief Appends one step grid (time value + point/cell data) to the
     * temporal collection of the existing file.
     * @param TimeValue The value written as <Time Value="...">.
     * @param rPointData Nodal arrays, in the node order of the mesh written by
     *                   @ref WriteMesh.
     * @param rCellData Cell arrays, in the cell order of the mesh written by
     *                  @ref WriteMesh (elements first, then conditions).
     */
    void AppendStep(
        const double TimeValue,
        const std::vector<DataArray>& rPointData,
        const std::vector<DataArray>& rCellData = {}
        );

    ///@}

private:
    ///@name Member Variables
    ///@{

    class Impl;
    std::unique_ptr<Impl> mpImpl; /// Hides the XML/meshio++ machinery from this header

    ///@}

}; // Class XdmfTimeSeriesWriter

/**
 * @brief Multi-format mesh input/output based on the meshio++ library.
 * @details Reads and writes ~40 mesh file formats (vtu, vtk, gmsh, med, xdmf,
 * abaqus, ...; see @ref Format) converting to/from Kratos model parts through
 * the meshio++ Kratos bridge. Availability of the HDF5-backed formats (med,
 * cgns, h5m, hmf, the HDF data path of xdmf) and the netCDF-backed ones
 * (exodus) depends on the build; query @ref GetSupportedFormats /
 * @ref IsFormatAvailable.
 *
 * When writing a volume mesh to a surface-only format (stl, ply) the boundary
 * skin is extracted and written by default; set the "skin" setting to false
 * for the legacy behavior (volume cells dropped, only existing surface cells
 * written).
 *
 * Reading: elements are created from the highest-dimension cell blocks and
 * conditions from the lower-dimension ones; integer tag arrays (gmsh physical
 * groups etc.) become sub model parts.
 *
 * Writing supports transient output: repeated calls to @ref WriteModelPart on
 * the same instance extend the current output instead of overwriting it. For
 * XDMF the steps are appended to the temporal collection of a single file
 * (checking whether the file "buffer" already exists); for the other formats a
 * file series <stem>_<label>.<ext> is produced. See the "time_series" setting.
 */
class KRATOS_API(KRATOS_CORE) MeshioPlusPlusIO
    : public IO
{
public:
    ///@name Type Definitions
    ///@{

    /// The index type definition
    using IndexType = std::size_t;

    /// Pointer definition of MeshioPlusPlusIO
    KRATOS_CLASS_POINTER_DEFINITION(MeshioPlusPlusIO);

    ///@}
    ///@name  Enum's
    ///@{

    /**
     * @brief The file formats this IO is compatible with.
     * @details AUTOMATIC resolves the format from the file extension. Formats
     * backed by an optional dependency (HDF5: CGNS, H5M, HMF, MED and the HDF
     * data path of XDMF; netCDF: EXODUS) are only available when the build
     * enables them - check with @ref IsFormatAvailable. OPENFOAM is read only;
     * SVG and TIKZ are write only.
     */
    enum class Format
    {
        AUTOMATIC, /// Resolve from the file extension
        ABAQUS,    /// Abaqus .inp
        ANSYS,     /// ANSYS .msh
        ANSYSINP,  /// ANSYS Mechanical APDL .inp (cdb)
        AVSUCD,    /// AVS-UCD .avs
        CGNS,      /// CGNS .cgns (requires HDF5)
        DEX,       /// Dexelas .dex
        DOLFIN,    /// DOLFIN XML .xml
        ENSIGHT,   /// EnSight Gold .case/.geo
        EXODUS,    /// Exodus II .e/.exo (requires netCDF)
        FLAC3D,    /// FLAC3D .f3grid
        FLUX,      /// Flux .pf3
        FREEFEM,   /// FreeFEM .msh
        GMSH,      /// Gmsh .msh
        H5M,       /// MOAB .h5m (requires HDF5)
        HMF,       /// HMF .hmf (requires HDF5)
        IP,        /// Ipreo .ip
        MED,       /// salome MED .med (requires HDF5)
        MEDIT,     /// Medit .mesh/.meshb
        MFF,       /// MFF .mff
        MFM,       /// MFM .mfm
        MPHTXT,    /// COMSOL .mphtxt
        NASTRAN,   /// Nastran .bdf/.nas
        NETGEN,    /// Netgen .vol
        OBJ,       /// Wavefront .obj
        OFF,       /// Object File Format .off
        OPENFOAM,  /// OpenFOAM polyMesh (read only)
        PERMAS,    /// PERMAS .post/.dato
        PLY,       /// Polygon File Format .ply
        STL,       /// Stereolithography .stl
        SU2,       /// SU2 .su2
        SVG,       /// Scalable Vector Graphics .svg (write only)
        TECPLOT,   /// Tecplot .dat
        TETGEN,    /// TetGen .node/.ele
        TIKZ,      /// LaTeX TikZ/PGF .tikz (write only)
        TRIANGLE,  /// Shewchuk Triangle .node/.ele/.poly (.node/.ele resolve to TETGEN by extension - select via the "format" setting)
        UGRID,     /// AFLR3 .ugrid
        UNV,       /// I-deas universal .unv
        VTK,       /// VTK legacy .vtk
        VTP,       /// VTK PolyData XML .vtp
        VTU,       /// VTK unstructured XML .vtu
        WKT,       /// Well-known text .wkt
        XDMF       /// XDMF v3 .xdmf/.xmf (HDF data path requires HDF5)
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructs a MeshioPlusPlusIO object for a file.
     * @param rFileName The path of the file to read or write.
     * @param ThisParameters Optional configuration, validated against
     *                       @ref GetDefaultParameters.
     */
    MeshioPlusPlusIO(
        const std::filesystem::path& rFileName,
        Parameters ThisParameters = Parameters()
        );

    /// Destructor.
    ~MeshioPlusPlusIO() override;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Retrieves the default parameters.
     * @return The default parameters.
     */
    static Parameters GetDefaultParameters();

    /**
     * @brief The names of every format this build can read or write.
     * @details Reflects the compiled-in optional dependencies (HDF5/netCDF).
     */
    static std::vector<std::string> GetSupportedFormats();

    /**
     * @brief The names of every format this build can read.
     */
    static std::vector<std::string> GetSupportedReadFormats();

    /**
     * @brief The names of every format this build can write.
     * @details Read-only formats (e.g. "openfoam") are not listed.
     */
    static std::vector<std::string> GetSupportedWriteFormats();

    /**
     * @brief Converts a format name to the @ref Format enum ("gmsh" -> GMSH).
     * @details "auto"/"automatic"/"" map to AUTOMATIC. Throws on unknown names.
     */
    static Format FormatFromString(const std::string& rFormatName);

    /**
     * @brief The canonical name of a @ref Format value (GMSH -> "gmsh").
     */
    static std::string FormatName(const Format Value);

    /**
     * @brief Resolves the format from a file path extension (".vtu" -> VTU).
     * @details Throws if the extension is not associated with any format.
     */
    static Format ResolveFormat(const std::filesystem::path& rPath);

    /**
     * @brief Whether the format is usable in this build.
     * @details False when the format was compiled out for lack of an optional
     * dependency (HDF5/netCDF); true otherwise. AUTOMATIC is always true.
     */
    static bool IsFormatAvailable(const Format Value);

    /**
     * @brief Reads a model part from the file.
     * @details The target model part must be empty. Elements are created from
     * the highest-dimension cell blocks, conditions from lower-dimension ones,
     * and integer tag arrays (e.g. gmsh physical groups) become sub model parts.
     * @param rThisModelPart Reference to the model part to read into.
     */
    void ReadModelPart(ModelPart& rThisModelPart) override;

    /**
     * @brief Writes the model part to the file.
     * @details Repeated calls extend the current output: for XDMF (with
     * "time_series" set to "automatic") the steps are appended to the temporal
     * collection of the file - if the file already exists with a valid time
     * series it is extended, otherwise it is created; for other formats one
     * file per call is written as <stem>_<label>.<ext>. With "time_series" set
     * to "single_file" every call overwrites the file.
     * @param rThisModelPart Const reference to the model part to write from.
     */
    void WriteModelPart(const ModelPart& rThisModelPart) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::filesystem::path mFileName; /// The file to read from / write to
    Parameters mParameters;          /// The configuration parameters

    // Transient write state ("does the output buffer already exist?")
    IndexType mOutputStep = 0; /// Number of WriteModelPart calls performed

    /// Time-series writers for XDMF transient output, one per output target
    /// ("" = the root model part; sub model part suffixes when
    /// "output_sub_model_parts" is enabled)
    std::map<std::string, std::unique_ptr<XdmfTimeSeriesWriter>> mXdmfWriters;

    /// Extrapolates the "gauss_point_variables_extrapolated_to_nodes" to nodal
    /// data before every write (created lazily at the first WriteModelPart)
    std::unique_ptr<IntegrationValuesExtrapolationToNodesProcess> mpGaussToNodesProcess;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief The effective format name: the "format" setting, or the extension
     * default when it is "auto". Throws a descriptive error for formats
     * compiled out of this build.
     * @param CheckWritable If true additionally require a writer for the format.
     */
    std::string ResolveEffectiveFormat(const bool CheckWritable) const;

    /**
     * @brief The step/time label used for file-series names and XDMF time values.
     */
    std::string GetOutputLabel(const ModelPart& rModelPart) const;

    /**
     * @brief The <Time Value> written for a transient XDMF step.
     */
    double GetOutputTimeValue(const ModelPart& rModelPart) const;

    /**
     * @brief Whether elements are written (from the "entity_type" setting).
     */
    bool WritesElements() const;

    /**
     * @brief Whether conditions are written (from the "entity_type" setting).
     */
    bool WritesConditions() const;

    /**
     * @brief Composes the output path from the configured directory/stem,
     * "custom_name_prefix"/"custom_name_postfix", an optional target suffix
     * (sub model part output) and an optional step/time label:
     * <dir>/<prefix><stem><TargetSuffix><postfix>[_<Label>]<ext>.
     */
    std::filesystem::path ComposeOutputPath(
        const std::string& rTargetSuffix,
        const std::string& rLabel
        ) const;

    /**
     * @brief Writes one output target (the root model part or one sub model
     * part) applying the "time_series" logic.
     */
    void WriteTarget(
        const ModelPart& rModelPart,
        const std::string& rTargetSuffix
        );

    /**
     * @brief Writes the model part to a concrete path in the given format
     * (single-shot, no transient logic).
     */
    void WriteStatic(
        const std::filesystem::path& rPath,
        const std::string& rFormatName,
        const ModelPart& rThisModelPart
        ) const;

    /**
     * @brief Transient XDMF write: creates or extends the temporal collection
     * of the target's file.
     */
    void WriteXdmfStep(
        const ModelPart& rThisModelPart,
        const std::string& rTargetSuffix
        );

    /**
     * @brief Collects the configured nodal variables, flags, node ids and
     * extrapolated gauss point results as flat point-data arrays (node
     * container order).
     */
    std::vector<XdmfTimeSeriesWriter::DataArray> CollectPointData(const ModelPart& rThisModelPart) const;

    /**
     * @brief Collects the configured elemental/conditional variables, flags,
     * ids and gauss point results as flat cell-data arrays (element rows
     * first, then condition rows, zero-filled for the entity kind a variable
     * does not apply to).
     */
    std::vector<XdmfTimeSeriesWriter::DataArray> CollectCellData(const ModelPart& rThisModelPart) const;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    MeshioPlusPlusIO& operator=(MeshioPlusPlusIO const& rOther) = delete;

    /// Copy constructor.
    MeshioPlusPlusIO(MeshioPlusPlusIO const& rOther) = delete;

    ///@}

}; // Class MeshioPlusPlusIO

///@}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const MeshioPlusPlusIO& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.
