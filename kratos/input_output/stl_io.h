//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include <filesystem>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/io.h"

namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// This class reads from STL file format and creates triangular elements in given model_part
/** The current version only reads triangles from the STL and not higher order polygons
 * The nodes corresponging to given vertices are not collapsed
 * A SubModelPart for each additional solid block will be created
 * For definition STL format please check https://en.wikipedia.org/wiki/STL_(file_format)
 * A sample file format with 3 triangles:
   solid 3 triangles
   facet normal  1.000000 0.000000 0.000000 
       outer loop 
          vertex 0.1 -2.56114e-08 0.1
          vertex 0.1 -0.499156 -0.0352136
          vertex 0.1 -0.473406 -0.0446259
       endloop 
   endfacet 
   facet normal  1.000000 -0.000000 0.000000 
       outer loop 
          vertex 0.1 -0.473406 -0.0446259
          vertex 0.1 -0.447464 -0.0534931
          vertex 0.1 -2.56114e-08 0.1
       endloop 
   endfacet 
   facet normal  1.000000 0.000000 0.000000 
       outer loop 
          vertex 0.1 -0.6 0.1
          vertex 0.1 -0.524702 -0.0252604
          vertex 0.1 -0.499156 -0.0352136
       endloop 
   endfacet 
   endsolid 3 triangles

*/
class KRATOS_API(KRATOS_CORE) StlIO 
    : public IO
{
public:
    ///@name Type Definitions
    ///@{

    /// Geometries map type definition
    using GeometriesMapType = ModelPart::GeometriesMapType;

    /// The nodes array type definition
    using NodesArrayType = Element::NodesArrayType;

    /// The index type definition
    using IndexType = std::size_t;

    /// Pointer definition of StlIO
    KRATOS_CLASS_POINTER_DEFINITION(StlIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /** 
     * @brief Constructs a StlIO object using a filename.
     * @details This constructor will create a StlIO object and open a file with the provided filename and open options. The default open option is read mode.
     * @param rFilename The path of the file to open.
     * @param ThisParameters Optional. Additional parameters for the StlIO object.
     *                       Defaults to an empty Parameters object.
     */
    StlIO(
        const std::filesystem::path& rFilename,
        Parameters ThisParameters = Parameters()
        );

    /** 
     * @brief Constructs a StlIO object using an input/output stream.
     * @details This constructor will create a StlIO object using a provided input/output stream.
     * @param pInputStream A shared pointer to the input/output stream to use.
     * @param ThisParameters Optional. Additional parameters for the StlIO object.
     *                       Defaults to an empty Parameters object.
     */
    StlIO(
        Kratos::shared_ptr<std::iostream> pInputStream,
        Parameters ThisParameters = Parameters()
        );

    /// Destructor.
    virtual ~StlIO(){}

    ///@}
    ///@name Operators
    ///@{

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
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

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
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    Parameters mParameters;         /// The configuration parameters

    IndexType mNextNodeId = 0;      /// The next node ID
    IndexType mNextElementId = 0;   /// The next element ID
    IndexType mNextConditionId = 0; /// The next condition ID

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    Kratos::shared_ptr<std::iostream> mpInputStream; /// The input stream
    Flags mOptions; /// The flags of the IO

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Read a solid object from a model part.
     * @details This function reads a solid object from the given model part using the provided function to create the entity.
     * @param rThisModelPart Reference to the model part to read the solid from.
     * @param rCreateEntityFunctor A functor to create entities.
     */
    void ReadSolid(
        ModelPart& rThisModelPart,
        const std::function<void(ModelPart&, NodesArrayType&)>& rCreateEntityFunctor 
        );

    /**
     * @brief Read a facet from a model part.
     * @details This function reads a facet from the given model part using the provided function to create the entity.
     * @param rThisModelPart Reference to the model part to read the facet from.
     * @param rCreateEntityFunctor A functor to create entities.
     */
    void ReadFacet(
        ModelPart& rThisModelPart,
        const std::function<void(ModelPart&, NodesArrayType&)>& rCreateEntityFunctor
        );

    /**
     * @brief Read a loop from a model part.
     * @details This function reads a loop from the given model part using the provided function to create the entity.
     * @param rThisModelPart Reference to the model part to read the loop from.
     * @param rCreateEntityFunctor A functor to create entities.
     */
    void ReadLoop(
        ModelPart& rThisModelPart,
        const std::function<void(ModelPart&, NodesArrayType&)>& rCreateEntityFunctor
        );

    /**
     * @brief Reads a point.
     * @return A point read from the source.
    */
    Point ReadPoint();

    /**
     * @brief Reads a keyword.
     * @param Keyword A string reference containing the keyword to read.
     */
    void ReadKeyword(const std::string& Keyword);

    /**
     * @brief Writes an entity block.
     * @details Writes an entity block of a given type to the destination.
     * @tparam TContainerType The type of the container for entities.
     * @param rThisEntities Reference to the entities to write.
     */
    template<class TContainerType>
    void WriteEntityBlock(const TContainerType& rThisEntities);

    /**
     * @brief Writes a geometry block.
     * @details Writes a geometry block to the destination.
     * @param rThisGeometries Reference to the geometries to write.
     */
    void WriteGeometryBlock(const GeometriesMapType& rThisGeometries);

    /**
     * @brief Writes an entity block (MPI version).
     * @details Writes an entity block of a given type to the destination.
     * @tparam TContainerType The type of the container for entities.
     * @param rThisEntities Reference to the entities to write.
     * @param rDataCommunicator The data communicator considered for MPI.
     */
    template<class TContainerType>
    void WriteEntityBlockMPI(
        const TContainerType& rThisEntities,
        const DataCommunicator& rDataCommunicator
        );

    /**
     * @brief Writes a geometry block (MPI version).
     * @details Writes a geometry block to the destination.
     * @param rThisGeometries Reference to the geometries to write.
     * @param rDataCommunicator The data communicator considered for MPI.
     */
    void WriteGeometryBlockMPI(
        const GeometriesMapType& rThisGeometries,
        const DataCommunicator& rDataCommunicator
        );

    /**
     * @brief Writes a facet.
     * @details Writes a facet to the destination.
     * @param rGeom Reference to the geometry to write.
     * @param rStream The stream considered
     * @tparam TStreamType The stream type considered.
     */
    template<class TStreamType>
    void WriteFacet(
        const GeometryType& rGeom,
        TStreamType& rStream
        );

    /**
     * @brief Checks the validity of a geometry.
     * @details Checks if a given geometry is valid or not.
     * @param rGeometry Reference to the geometry to check.
     * @param rNumDegenerateGeos Reference to store the number of degenerate geometries.
     * @return True if the geometry is valid, false otherwise.
     */
    bool IsValidGeometry(
        const Geometry<Node>& rGeometry,
        IndexType& rNumDegenerateGeos
        ) const;

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    StlIO& operator=(StlIO const& rOther);

    /// Copy constructor.
    StlIO(StlIO const& rOther);

    ///@}

}; // Class StlIO

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                StlIO& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const StlIO& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.