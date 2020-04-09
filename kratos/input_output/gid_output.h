//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Janosch Stascheit
//                   Pooyan Dadvand
//

#if !defined(KRATOS_GID_OUTPUT_H_INCLUDED)
#define  KRATOS_GID_OUTPUT_H_INCLUDED

// System includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstddef>
#include <iomanip>

// External includes
#define USE_CONST
#include "gidpost/source/gidpost.h"

// Project includes
#include "includes/define.h"
#include "includes/io.h"
#include "includes/gid_gauss_point_container.h"
#include "includes/gid_mesh_container.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///Flags for mesh writing
enum WriteDeformedMeshFlag {WriteDeformed, WriteUndeformed};
enum WriteConditionsFlag {WriteConditions, WriteElementsOnly, WriteConditionsOnly};
enum MultiFileFlag {SingleFile, MultipleFiles};

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class GidIOBase
 * @ingroup KratosCore
 * @brief Base class for GidIO
 * @author Riccardo Rossi
 * @author Janosch Stascheit
 * @author Pooyan Dadvand
 */
class KRATOS_API(KRATOS_CORE) GidIOBase
    : public IO
{
protected:
    /**
     * Counter of live GidIO instances
     * (to ensure GiD_PostInit and GiD_PostDone are properly called)
     */
    int data;

    // Private constructor so that no objects can be created.
    GidIOBase() {
        data = 0;
    }

public:
    static GidIOBase& GetInstance();

    int GetData();
    void SetData(int data);

private:
    static void Create();

    static GidIOBase* mpInstance;
};

/**
 * @class GidIO
 * @ingroup KratosCore
 * @brief This class defines an interface to the GiDPost library in order to provide GiD compliant I/O functionality
 * @tparam TGaussPointContainer The gauss point container considered
 * @tparam TMeshContainer The mesh container considered
 * @author Riccardo Rossi
 * @author Janosch Stascheit
 * @author Pooyan Dadvand
 */
template<class TGaussPointContainer = GidGaussPointsContainer, class TMeshContainer = GidMeshContainer>
class KRATOS_API(KRATOS_CORE) GidIO
    : public GidIOBase
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GidIO
    KRATOS_CLASS_POINTER_DEFINITION(GidIO);

    /// Base class definition
    typedef IO BaseType;

    /// Containers definition
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    /// Integration method definition
    typedef GeometryData::IntegrationMethod IntegrationMethodType;

    /// Geometry family definition
    typedef GeometryData::KratosGeometryFamily KratosGeometryFamily;

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor. Single stream IO constructor
     */
    GidIO(
        const std::string& rDatafilename,
        const GiD_PostMode Mode,
        const MultiFileFlag UseMultipleFilesFlag,
        const WriteDeformedMeshFlag WriteDeformedFlag,
        const WriteConditionsFlag WriteConditions
         );

    ///Destructor.
    ~GidIO() override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates the mesh containers for all different element types.
     * @note The containers are not filled yet in here!
     */
    void SetUpMeshContainers();

    /**
     * @brief Creates the gauss point containers for all different element types.
     * @note The containers are not filled yet in here!
     */
    virtual void SetUpGaussPointContainers();


    /// General GidIO related functions ///

    /**
     * @brief
     * @todo To be removed
     */
    void ChangeOutputName(const std::string& rDatafilename );

    /**
     * sets up the file names and opens the result file in case there
     * is ASCII mode and only one file written
     */
    void InitializeResultFile( std::string const& rResultFileName );

    /**
     * TODO: check whether this is still necessary!
     */
    void  CloseResultFile();

    /**
     * TODO: check whether this is still necessary!
     */
    void Flush();

    /**
     * Turn back information as a string.
     */
    std::string Info() const override
    {
        return "gid io";
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///result functions
    /**
     * This has to be called for each solution step BEFORE any results
     * (on nodes and on gauss points) is written
     * @param SolutionTag the current solution step (i.e. time)
     * @param conditions_flag states whether results should also be written on conditions
     */
    virtual void InitializeResults(
        const double name,
        const MeshType& rThisMesh
        );
    /**
     * This has to be called for each solution step AFTER all the results
     * have been written
     */
    void FinalizeResults();

    ///functions for writing nodal results

    ///////////////////////////////////////////////////////////////////////
    //////                  HISTORICAL DATABASE BLOCK                 /////
    ///////////////////////////////////////////////////////////////////////
     /**
     * writes nodal results for variables of type bool
     */
    void WriteNodalResults(
        Variable<bool> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag,
        const std::size_t SolutionStepNumber
        );

    ///functions for writing nodal results
    /**
     * writes nodal results for variables of type double
     */
    void WriteNodalResults(
        Variable<double> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag,
        const std::size_t SolutionStepNumber
        );

    /**
     * writes nodal results for variables of type int
     */
    void WriteNodalResults(
        Variable<int> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag,
        const std::size_t SolutionStepNumber
        );


    /**
     * writes nodal results for variables of type array_1d<double, 3>
     * (e.g. DISPLACEMENT)
     */
    void WriteNodalResults(
        Variable<array_1d<double, 3> > const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag,
        const std::size_t SolutionStepNumber
        );

    /**
     * writes nodal results for variables of type Vector
     * (note that only vectors with 3 components can be printed)
     */
    void WriteNodalResults(
        Variable<Vector> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag,
        const std::size_t SolutionStepNumber
        );

    /**
     * writes nodal results for variables of type Matrix
     */
    void WriteNodalResults(
        Variable<Matrix> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag,
        const std::size_t SolutionStepNumber
        );

    void WriteLocalAxesOnNodes(
        Variable<array_1d<double, 3> > const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag,
        const std::size_t SolutionStepNumber
        );

    ///////////////////////////////////////////////////////////////////////
    //////                 NON- HISTORICAL DATABASE BLOCK             /////
    ///////////////////////////////////////////////////////////////////////

   /**
    * Writes nodal flags
    */
    void WriteNodalFlags(
        const Kratos::Flags& rFlag,
        const std::string& rFlagName,
        const NodesContainerType& rNodes,
        const double SolutionTag
        );

    /**
     * writes nodal results for variables of type bool
     */
    void WriteNodalResultsNonHistorical(
        Variable<bool> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag
        );

    ///functions for writing nodal results
    /**
     * writes nodal results for variables of type double
     */
    void WriteNodalResultsNonHistorical(
        Variable<double> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag
        );

    /**
     * writes nodal results for variables of type int
     */
    void WriteNodalResultsNonHistorical(
        Variable<int> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag
        );

    /**
     * writes nodal results for variables of type array_1d<double, 3>
     * (e.g. DISPLACEMENT)
     */
    void WriteNodalResultsNonHistorical(
        Variable<array_1d<double, 3> > const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag
        );

    /**
     * writes nodal results for variables of type Vector
     * (note that only vectors with 3 components can be printed)
     */
    void WriteNodalResultsNonHistorical(
        Variable<Vector> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag
        );

    /**
     * writes nodal results for variables of type Matrix
     */
    void WriteNodalResultsNonHistorical(
        Variable<Matrix> const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag
        );

    void WriteLocalAxesOnNodesNonHistorical(
        Variable<array_1d<double, 3> > const& rVariable,
        const NodesContainerType& rNodes,
        const double SolutionTag
        );

    ///mesh writing functions
    /**
     * opens a new mesh group
     */
    void InitializeMesh( const double name );

    /**
     * closes a mesh group
     */
    void FinalizeMesh();

    /**
     * Writes a node mesh.
     * @param rThisMesh the given mesh to be written to the output file
     * @param solution_step the current solution step
     * @param deformed_flag indicates whether the mesh shall be written in deformed
     * or undeformed state
     * @param Mode either GiD_PostAscii (default) or GiD_PostBinary
     */
    void WriteNodeMesh(MeshType& rThisMesh) override;

    void WriteSphereMesh(const MeshType& rThisMesh);


    void WriteCircleMesh(const MeshType& rThisMesh);

    void WriteClusterMesh(const MeshType& rThisMesh);


    /**
     * This is a multi-purpose function that writes arbitrary meshes of elements
     * and conditions in either deformed or undeformed state
     * @param rThisMesh the current mesh to be written
     * @param deformed_flag states whether the mesh should be written in deformed configuration
     * @param conditions_flag states whether conditions should also be written
     */
    void WriteMesh(MeshType& rThisMesh) override;


    ///functions for printing results on gauss points

    /**
    * @brief Writes elemental and conditional flags
    * @param rFlag the flag
    * @param rFlagName the given flag name
    * @param rModelPart the current model part
    */
    void PrintFlagsOnGaussPoints(
        const Kratos::Flags& rFlag,
        const std::string& rFlagName,
        const ModelPart& rModelPart,
        const double SolutionTag
        );

    /**
     * Prints variables of type int on gauss points of the complete mesh
     * @param rVariable the given variable name
     * @param rModelPart the current model part
     */
    virtual void PrintOnGaussPoints(
        const Variable<bool>& rVariable,
        const ModelPart& rModelPart,
        const double SolutionTag, const int ValueIndex = 0 );

    /**
     * Prints variables of type int on gauss points of the complete mesh
     * @param rVariable the given variable name
     * @param rModelPart the current model part
     */
    virtual void PrintOnGaussPoints(
        const Variable<int>& rVariable,
        const ModelPart& rModelPart,
        const double SolutionTag,
        const int ValueIndex = 0
        );

    /**
     * Prints variables of type double on gauss points of the complete mesh
     * @param rVariable the given variable name
     * @param rModelPart the current model part
     */
    virtual void PrintOnGaussPoints(
        const Variable<double>& rVariable,
        const ModelPart& rModelPart,
        const double SolutionTag,
        const int ValueIndex = 0
        );

    /**
     * Prints variables of type double on gauss points of the complete mesh
     * @param rVariable the given variable name
     * @param rModelPart the current model part
     */
    virtual void PrintOnGaussPoints(
        const Variable<array_1d<double,3> >& rVariable,
        const ModelPart& rModelPart,
        const double SolutionTag,
        const int ValueIndex = 0
        );

    /**
     * Prints variables of type double on gauss points of the complete mesh
     * @param rVariable the given variable name
     * @param rModelPart the current model part
     */
    virtual void PrintOnGaussPoints(
        const Variable<Vector>& rVariable,
        const ModelPart& rModelPart,
        const double SolutionTag,
        const int ValueIndex = 0
        );

    /**
     * Prints variables of type double on gauss points of the complete mesh
     * @param rVariable the given variable name
     * @param rModelPart the current model part
     */
    virtual void PrintOnGaussPoints(
        const Variable<Matrix>& rVariable,
        const ModelPart& rModelPart,
        const double SolutionTag,
        const int ValueIndex = 0
        );

protected:
    /**
     * File names
     */
    std::string mResultFileName;
    std::string mMeshFileName;

    GiD_FILE mMeshFile;
    GiD_FILE mResultFile;

    /**
     * Flags
     */
    WriteDeformedMeshFlag mWriteDeformed;
    WriteConditionsFlag mWriteConditions;
    MultiFileFlag mUseMultiFile;
    GiD_PostMode mMode;

    /**
     * member variables
     */
    std::vector<TMeshContainer> mGidMeshContainers;
    std::vector<TGaussPointContainer> mGidGaussPointContainers;
    bool mMeshFileOpen;
    bool mResultFileOpen;

private:
    /**
     * assignment operator
     */
    GidIO& operator=(GidIO const& rOther);

    /**
     * Copy constructor
     */
    GidIO(GidIO const& rOther);
}; // Class GidIO

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) GidIO<GidGaussPointsContainer,GidMeshContainer>;

///@}
///@name Input and output
///@{

inline std::ostream& operator << (std::ostream& rOStream, const GidIO<>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

///@}

#endif // KRATOS_GID_OUTPUT_H_INCLUDED  defined
