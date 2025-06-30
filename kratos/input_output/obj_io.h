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
#include <string>
#include <iostream>
#include <filesystem>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/io.h"
#include "includes/variables.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class ObjIO
 * @ingroup KratosCore
 * @brief This class reads from OBJ file format and creates elements in given model_part
 * @details The current version reads vertices, vertex normals, and faces from the OBJ file.
 * For the OBJ file format, please check http://fegemo.github.io/cefet-cg/attachments/obj-spec.pdf
 * A sample file format with vertices, vertex normals, and faces:
    v 0.123 0.234 0.345
    v 0.234 0.345 0.456
    v 0.345 0.456 0.567
    vn 0.707 0.000 0.707
    vn 0.000 1.000 0.000
    vn -0.707 0.000 0.707
    f 1 2 3
 * The vertices are defined by the "v" keyword followed by the x, y, and z coordinates.
 * The vertex normals are defined by the "vn" keyword followed by the x, y, and z coordinates.
 * The faces are defined by the "f" keyword followed by the vertex index, vertex normal index, and texture index.
 * The vertex normal index is optional and can be omitted. (e.g., "f 1 2 3"). We will assume that the vertex normal index is omitted.
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) ObjIO
    : public IO
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of supported geometries
    static constexpr std::array<const char*, 2> SupportedGeometries = {"Triangle3D3", "Quadrilateral3D4"};

    /// The index type definition
    using IndexType = std::size_t;

    /// Pointer definition of ObjIO
    KRATOS_CLASS_POINTER_DEFINITION(ObjIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructs an ObjIO object using a filename.
     * @details This constructor will create an ObjIO object and open a file with the provided filename and open options. The default open option is read mode.
     * @param rFilename The path of the file to open.
     * @param ThisParameters Optional. Additional parameters for the ObjIO object.
     *                       Defaults to an empty Parameters object.
     */
    ObjIO(
        const std::filesystem::path& rFilename,
        Parameters ThisParameters = Parameters()
        );

    /**
     * @brief Constructs an ObjIO object using an input/output stream.
     * @details This constructor will create an ObjIO object using a provided input/output stream.
     * @param pInputStream A shared pointer to the input/output stream to use.
     * @param ThisParameters Optional. Additional parameters for the ObjIO object.
     *                       Defaults to an empty Parameters object.
     */
    ObjIO(
        Kratos::shared_ptr<std::iostream> pInputStream,
        Parameters ThisParameters = Parameters()
        );

    /// Destructor.
    virtual ~ObjIO(){}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Retrieves the default parameters.
     * @return The default parameters.
     */
    static Parameters GetDefaultParameters();

    /**
     * @brief Reads a model part.
     * @details Reads a model part from the source and stores it in the provided model part object.
     * @param rThisModelPart Reference to the model part to read into.
     */
    void ReadModelPart(ModelPart& rThisModelPart) override;

    /**
     * @brief Writes a model part.
     * @details Writes the provided model part to the destination.
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
protected:
    ///@name Protected member Variables
    ///@{

    Parameters mParameters;          /// The configuration parameters

    IndexType mFirstNodeId = 0;      /// The first node ID
    IndexType mNextNodeId = 0;       /// The next node ID

    IndexType mFirstElementId = 0;   /// The first element ID
    IndexType mNextElementId = 0;    /// The next element ID

    IndexType mFirstConditionId = 0; /// The first condition ID
    IndexType mNextConditionId = 0;  /// The next condition ID

    IndexType mNormalCounter = 0;    /// The normal counter

    ///@}
private:
    ///@name Member Variables
    ///@{

    Kratos::shared_ptr<std::iostream> mpInputStream; /// The input stream
    Flags mOptions; /// The flags of the IO

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Reads vertices, vertex normals, and faces from the OBJ file.
     * @details This function reads the OBJ file line by line, creating nodes and elements in the model part.
     * @param rThisModelPart Reference to the model part to read into.
     * @param rEntityType The entity type to create in the model part. Can be "element" or "condition".
     * @param NormalAsHistoricalVariable If true, the normals are stored as historical variables in the nodes.
     * @param DecomposeQuadrilateral If true, the quadrilateral faces are decomposed into triangles.
     */
    void ReadVerticesAndFaces(
        ModelPart& rThisModelPart,
        const std::string& rEntityType = "element",
        const bool NormalAsHistoricalVariable = false,
        const bool DecomposeQuadrilateral = false
        );

    /**
     * @brief Parses a vertex line from the OBJ file.
     * @details Creates a node in the model part for the vertex.
     * @param rThisModelPart Reference to the model part.
     * @param rLine The line containing the vertex data.
     */
    void ParseVertexLine(
        ModelPart& rThisModelPart,
        const std::string& rLine
        );

    /**
     * @brief Parses a vertex normal line from the OBJ file.
     * @details Stores the vertex normal in the normals list.
     * @param rThisModelPart Reference to the model part.
     * @param rLine The line containing the vertex normal data.
     * @param NormalAsHistoricalVariable If true, the normals are stored as historical variables in the nodes.
     */
    void ParseNormalLine(
        ModelPart& rThisModelPart,
        const std::string& rLine,
        const bool NormalAsHistoricalVariable = false
        );

    /**
     * @brief Parses a face line from the OBJ file.
     * @details Creates an element or condition in the model part for the face.
     * @param rThisModelPart Reference to the model part.
     * @param rLine The line containing the face data.
     * @param rEntityType The entity type to create in the model part. Can be "element" or "condition".
     * @param DecomposeQuadrilateral If true, the quadrilateral faces are decomposed into triangles.
     */
    void ParseFaceLine(
        ModelPart& rThisModelPart,
        const std::string& rLine,
        const std::string& rEntityType = "element",
        const bool DecomposeQuadrilateral = false
        );

    /**
     * @brief Splits a string into tokens based on whitespace.
     * @param rLine The string to split.
     * @return A vector of tokens.
     */
    std::vector<std::string> Tokenize(const std::string& rLine);

    ///@}

}; // Class ObjIO

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const ObjIO& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.