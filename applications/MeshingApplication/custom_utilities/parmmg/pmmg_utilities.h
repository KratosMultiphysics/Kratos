// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:         BSD License
//                   license: MeshingApplication/license.txt
//
//  Main authors:    Marc Nunez
//                   Carlos Roig
//                   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_PMMG_UTILITIES)
#define KRATOS_PMMG_UTILITIES

// System includes

// External includes

// Project includes
#include "custom_utilities/mmg/mmg_utilities.h"

// NOTE: The following contains the license of the PMMG library
/* =============================================================================
**  Copyright (c) Bx INP/Inria/UBordeaux, 2017- .
**
**  parmmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  parmmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// Index definition
    typedef std::size_t                                               IndexType;

    /// Size definition
    typedef std::size_t                                                SizeType;

    /// Index vector
    typedef std::vector<IndexType>                              IndexVectorType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @struct PMMGMeshInfo
 * @ingroup MeshingApplication
 * @brief Stores the ParMmg mesh information
 * @author Vicente Mataix Ferrandiz
 */
template<PMMGLibrary TPMMGLibrary>
struct PMMGMeshInfo : public MMGMeshInfo<MMGLibrary::MMG3D>
{
};

/**
 * @class ParMmgUtilities
 * @ingroup MeshingApplication
 * @brief Provides the Kratos interface to the PMMG library API
 * @author Marc Nunez
 * @author Carlos Roig
 * @author Vicente Mataix Ferrandiz
 */
template<PMMGLibrary TPMMGLibrary>
class KRATOS_API(MESHING_APPLICATION) ParMmgUtilities
    : public MmgUtilities<MMGLibrary::MMG3D>
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of ParMmgUtilities
    KRATOS_CLASS_POINTER_DEFINITION(ParMmgUtilities);

    typedef MmgUtilities<MMGLibrary::MMG3D> BaseType;

    /// Node definition
    typedef Node                                                   NodeType;
    // Geometry definition
    typedef Geometry<NodeType>                                     GeometryType;

    /// Spatial dimension
    static constexpr SizeType Dimension = 3;

    /// The type of array considered for the tensor
    typedef typename std::conditional<Dimension == 2, array_1d<double, 3>, array_1d<double, 6>>::type TensorArrayType;

    /// Double vector
    typedef std::vector<double> DoubleVectorType;

    /// Double vector map
    typedef std::unordered_map<DoubleVectorType, IndexType, KeyHasherRange<DoubleVectorType>, KeyComparorRange<DoubleVectorType> > DoubleVectorMapType;

    /// Index vector map
    typedef std::unordered_map<IndexVectorType, IndexType, KeyHasherRange<IndexVectorType>, KeyComparorRange<IndexVectorType> > IndexVectorMapType;

    /// Colors map
    typedef std::unordered_map<IndexType,IndexType> ColorsMapType;

    /// Index pair
    typedef std::pair<IndexType,IndexType> IndexPairType;

    /// Index and string vector pair
    typedef std::pair<IndexType, std::vector<std::string>> IndexStringVectorPairType;

    /// Definition of the zero tolerance
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

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
     * @brief It prints info about the current mesh
     * @param[in,out] rPMMGMeshInfo The number of nodes, conditions and elements
     */
    void PrintAndGetParMmgMeshInfo(PMMGMeshInfo<TPMMGLibrary>& rPMMGMeshInfo);

    /**
     * @brief Returns a vector of ids of spatially repeated nodes
     * @param[in,out] rModelPart The model part whose nodes are checked
     */
    IndexVectorType FindDuplicateNodeIds(const ModelPart& rModelPart) override;

    /**
     * @brief Returns a vector of ids of repeated conditions
     * @details For 2D and surface meshes it returns the ids of repeated edges found in the PMMG mesh data structure. In 3D it returns ids of triangles.
     */
    IndexVectorType CheckFirstTypeConditions() override;

    /**
     * @brief Returns a vector of ids of repeated elements
     * @details For 2D and surface meshes it returns the ids of repeated triangles found in the PMMG mesh data structure. In 3D it returns ids of tetrahedra.
     */
    IndexVectorType CheckFirstTypeElements() override;

    /**
     * @brief It blocks certain nodes before remesh the model
     * @param[in] iNode The index of the node
     */
    void BlockNode(const IndexType iNode) override;

    /**
     * @brief It blocks certain conditions before remesh the model
     * @param[in] iCondition The index of the condition
     */
    void BlockCondition(const IndexType iCondition) override;

    /**
     * @brief It blocks certain elements before remesh the model
     * @param[in] iElement The index of the element
     */
    void BlockElement(const IndexType iElement) override;

    /**
     * @brief It creates the new node
     * @details Each call to this function increments the internal counter of the PMMG mesh data structure, extracts the coordinates of the current vertex and creates the new node with the extracted coordinates.
     * @param[in,out] rModelPart The model part which owns the new node
     * @param[in] iNode The index of the new node
     * @param[out] Ref PMMG point reference
     * @param[out] IsRequired PMMG required entity (0 or 1)
     * @return pNode The pointer to the new node created
     */
    NodeType::Pointer CreateNode(
        ModelPart& rModelPart,
        IndexType iNode,
        int& Ref,
        int& IsRequired
        ) override;

    /**
     * @brief It creates the new condition (first type, depends if the library work in 2D/3D/Surfaces)
     * @details Each call to this function increments the internal counter of the PMMG mesh data structure and extracts the vertices of the current  edge.
     * @param[in,out] rModelPart The model part of interest to study
     * @param[in,out] rMapPointersRefCondition The pointer to the condition of reference
     * @param[in] CondId The id of the new condition
     * @param[out] Ref PMMG edge reference
     * @param[out] IsRequired PMMG required entity (0 or 1)
     * @param[in] SkipCreation Skips the creation of the new condition
     * @return pCondition The pointer to the new condition created
     */
    Condition::Pointer CreateFirstTypeCondition(
        ModelPart& rModelPart,
        std::unordered_map<IndexType,Condition::Pointer>& rMapPointersRefCondition,
        const IndexType CondId,
        int& Ref,
        int& IsRequired,
        bool SkipCreation
        ) override;

    /**
     * @brief It creates the new element (first type, depends if the library work in 2D/3D/Surfaces)
     * @param[in,out] rModelPart The model part of interest to study
     * @param[in,out] rMapPointersRefElement The pointer to the element of reference
     * @param[in] ElemId The id of the element
     * @param[out] Ref PMMG edge reference
     * @param[out] IsRequired PMMG required entity (0 or 1)
     * @return pElement The pointer to the new condition created
     */
    Element::Pointer CreateFirstTypeElement(
        ModelPart& rModelPart,
        std::unordered_map<IndexType,Element::Pointer>& rMapPointersRefElement,
        const IndexType ElemId,
        int& Ref,
        int& IsRequired,
        bool SkipCreation
        ) override;

    /**
     * @brief Initialisation of mesh and sol structures
     * @details Initialisation of mesh and sol structures args of InitMesh:
     * -# MMG5_ARG_start we start to give the args of a variadic func
     * -# MMG5_ARG_ppMesh next arg will be a pointer over a MMG5_pMesh
     * -# &mmgMesh pointer toward your MMG5_pMesh (that store your mesh)
     * -# MMG5_ARG_ppMet next arg will be a pointer over a MMG5_pSol storing a metric
     * -# &mmgSol pointer toward your MMG5_pSol (that store your metric)
     * @param[in] rDataCommunicator Data communicator to be used to initialize ParMmg
     *
     */
    void InitMesh(const DataCommunicator& rDataCommunicator);

    /**
     * @brief Here the verbosity is set
     */
    void InitVerbosity() override;

    /**
     * @brief Here the API mode is set using the API
     * @param[in] API Mode sets the mode in which the parallel communicator works
     */
    void InitAPIModeParameter(const IndexType VerbosityPMMG);

    /**
     * @brief Ask for output node global numbering
     * @param[in] 1 to ask it, 0 to skip it
     */
    void InitNodeGloNumParameter(const IndexType nodeGloNum);

    /**
     * @brief Here the verbosity is set using the API
     * @param[in] VerbosityPMMG The equivalent verbosity level in the PMMG API
     */
    void InitVerbosityParameter(const IndexType VerbosityPMMG) override;

    /**
     * @brief This sets the size of the mesh
     * @param[in,out] rPMMGMeshInfo The number of nodes, conditions and elements
     */
    void SetMeshSize(PMMGMeshInfo<TPMMGLibrary>& rPMMGMeshInfo);

    /**
     * @brief This sets the size of the solution for the scalar case
     * @param[in] NumNodes Number of nodes
     */
    void SetSolSizeScalar(const SizeType NumNodes) override;

    /**
     * @brief This sets the size of the solution for the vector case
     * @param[in] NumNodes Number of nodes
     */
    void SetSolSizeVector(const SizeType NumNodes) override;

    /**
     * @brief This sets the size of the solution for the tensor case
     * @param[in] NumNodes Number of nodes
     */
    void SetSolSizeTensor(const SizeType NumNodes) override;

    /**
     * @brief This sets the size of the displacement for lagrangian movement
     * @param[in] NumNodes Number of nodes
     */
    void SetDispSizeVector(const SizeType NumNodes) override;

    /**
     * @brief This checks the mesh data and prints if it is OK
     */
    void CheckMeshData() override;

    /**
     * @brief This sets the output mesh
     * @param[in] rInputName The input name
     */
    void InputMesh(const std::string& rInputName) override;

    /**
     * @brief This sets the output sol
     * @param[in] rInputName The input name
     */
    void InputSol(const std::string& rInputName) override;

    /**
     * @brief This sets the output mesh
     * @param[in] rOutputName The output name
     */
    void OutputMesh(const std::string& rOutputName) override;

    /**
     * @brief This sets the output sol
     * @param[in] rOutputName The output name
     */
    void OutputSol(const std::string& rOutputName) override;

    /**
     * @brief This sets the output displacement
     * @param[in] rOutputName The output name
     */
    void OutputDisplacement(const std::string& rOutputName) override;

    /**
     * @brief This method generates the maps of reference for conditions and elements from an existing json
     * @param[in] rOutputName The name of the files
     * @param[in] rRefCondition The conditions of reference
     * @param[in] rRefElement The elements of reference
     */
    void OutputReferenceEntitities(
        const std::string& rOutputName,
        const std::unordered_map<IndexType,Condition::Pointer>& rRefCondition,
        const std::unordered_map<IndexType,Element::Pointer>& rRefElement
        ) override;

    /**
     * @brief This frees the PMMG structures
     */
    void FreeAll() override;

    /**
     * @brief This loads the solution
     */
    void PMMGLibCallMetric(Parameters ConfigurationParameters);

    /**
     * @brief This loads the solution
     */
    void PMMGLibCallIsoSurface(Parameters ConfigurationParameters);

    /**
     * @brief This sets the nodes of the mesh
     * @param[in] X Coordinate X
     * @param[in] Y Coordinate Y
     * @param[in] Z Coordinate Z
     * @param[in] Color Reference of the node(submodelpart)
     * @param[in] Index The index number of the node
     */
    void SetNodes(
        const double X,
        const double Y,
        const double Z,
        const IndexType Color,
        const IndexType Index
        ) override;

    /**
     * @brief This sets the conditions of the mesh
     * @param[in,out] rGeometry The geometry of the condition
     * @param[in] Color Reference of the node(submodelpart)
     * @param[in] Index The index number of the node
     */
    void SetConditions(
        GeometryType& rGeometry,
        const IndexType Color,
        const IndexType Index
        ) override;

    /**
     * @brief This sets elements of the mesh
     * @param[in,out] rGeometry The geometry of the element
     * @param[in] Color Reference of the node(submodelpart)
     * @param[in] Index The index number of the node
     */
    void SetElements(
        GeometryType& rGeometry,
        const IndexType Color,
        const IndexType Index
        ) override;

    /**
     * @brief This function is used to set the metric scalar
     * @param[in] Metric The inverse of the size node
     * @param[in] NodeId The id of the node
     */
    void SetMetricScalar(
        const double Metric,
        const IndexType NodeId
        ) override;

    /**
     * @brief This function is used to set the metric vector (x, y, z)
     * @param[in] rMetric This array contains the components of the metric vector
     * @param[in] NodeId The id of the node
     */
    void SetMetricVector(
        const array_1d<double, Dimension>& rMetric,
        const IndexType NodeId
        ) override;

    /**
     * @brief This function is used to set the Hessian metric tensor, note that when using the Hessian, more than one metric can be defined simultaneously, so in consecuence we need to define the elipsoid which defines the volume of maximal intersection
     * @param[in] rMetric This array contains the components of the metric tensor in the PMMG defined order
     * @param[in] NodeId The id of the node
     */
    void SetMetricTensor(
        const TensorArrayType& rMetric,
        const IndexType NodeId
        ) override;

    /**
     * @brief This function is used to set the displacement vector (x, y, z)
     * @param[in] rMetric This array contains the components of the displacement vector
     * @param[in] NodeId The id of the node
     */
    void SetDisplacementVector(
        const array_1d<double, 3>& rDisplacement,
        const IndexType NodeId
        ) override;

    /**
     * @brief This function is used to retrieve the metric scalar
     * @param[in,out] rMetric The inverse of the size node
     */
    void GetMetricScalar(double& rMetric) override;

    /**
     * @brief This function is used to retrieve the metric vector (x, y, z)
     * @param[in,out] rMetric This array contains the components of the metric vector
     */
    void GetMetricVector(array_1d<double, Dimension>& rMetric) override;

    /**
     * @brief This function is used to retrieve the Hessian metric tensor, note that when using the Hessian, more than one metric can be defined simultaneously, so in consecuence we need to define the elipsoid which defines the volume of maximal intersection
     * @param[in,out] rMetric This array contains the components of the metric tensor in the PMMG defined order
     */
    void GetMetricTensor(TensorArrayType& rMetric) override;

    /**
     * @brief This function is used to retrieve the displacement vector (x, y, z)
     * @param[in,out] rDisplacement This array contains the components of the displacement vector
     */
    void GetDisplacementVector(array_1d<double, 3>& rDisplacement) override;

    /**
     * @brief This method generates mesh data from an existing model part
     * @param[in,out] rModelPart The model part of interest to study
     * @param[in,out] rColors Where the sub model parts IDs are stored
     * @param[in,out] rColorMapCondition Auxiliar color map for conditions
     * @param[in,out] rColorMapElement Auxiliar color map for elements
     * @param[in] Framework The framework considered
     * @param[in] CollapsePrismElements If the prisms elements are going to be collapsed
     */
    void GenerateMeshDataFromModelPart(
        ModelPart& rModelPart,
        std::unordered_map<IndexType,std::vector<std::string>>& rColors,
        ColorsMapType& rColorMapCondition,
        ColorsMapType& rColorMapElement,
        const FrameworkEulerLagrange Framework = FrameworkEulerLagrange::EULERIAN,
        const bool CollapsePrismElements = false
        ) override;

    /**
     * @brief This method generates the interface data for the parallel communicator
     * @param[in,out] rModelPart The model part with the kratos communicator.
     */
    void GenerateParallelInterfaces(
        ModelPart& rModelPart
    );

    /**
     * @brief This method prints the interface data for the parallel communicator
     * @param[in,out] rModelPart The model part with the kratos communicator.
     */
    void PrintParallelInterfaces(
        ModelPart& rModelPart
    );

    /**
     * @brief This method generates the maps of reference for conditions and elements
     * @param[in] rModelPart The model part of interest to study
     * @param[in] rColorMapCondition Auxiliar color map for conditions
     * @param[in] rColorMapElement Auxiliar color map for elements
     * @param[in,out] rRefCondition The conditions of reference
     * @param[in,out] rRefElement The elements of reference
     */
    void GenerateReferenceMaps(
        ModelPart& rModelPart,
        const ColorsMapType& rColorMapCondition,
        const ColorsMapType& rColorMapElement,
        std::unordered_map<IndexType,Condition::Pointer>& rRefCondition,
        std::unordered_map<IndexType,Element::Pointer>& rRefElement
        ) override;

    /**
     * @brief This method generates solution (metric) data from an existing model part
     * @param[in,out] rModelPart The model part of interest to study
     */
    void GenerateSolDataFromModelPart(ModelPart& rModelPart) override;

    /**
     * @brief This method generates displacement data from an existing model part
     * @param[in,out] rModelPart The model part of interest to study
     */
    void GenerateDisplacementDataFromModelPart(ModelPart& rModelPart) override;

    /**
     * @brief This method writes mesh data to an existing model part
     * @param[in,out] rModelPart The model part of interest to study
     * @param[in] rColors Where the sub model parts IDs are stored
     * @param[in] rDofs Storage for the dof of the node
     * @param[in] rPMMGMeshInfo The resulting mesh data
     * @param[in] rMapPointersRefCondition The map of the conditions by reference (color)
     * @param[in] rMapPointersRefElement The map of the elements by reference (color)
     */
    void WriteMeshDataToModelPart(
        ModelPart& rModelPart,
        const std::unordered_map<IndexType,std::vector<std::string>>& rColors,
        const NodeType::DofsContainerType& rDofs,
        const PMMGMeshInfo<TPMMGLibrary>& rPMMGMeshInfo,
        std::unordered_map<IndexType,Condition::Pointer>& rMapPointersRefCondition,
        std::unordered_map<IndexType,Element::Pointer>& rMapPointersRefElement
        );

    /**
     * @brief This method writes the maps of reference for conditions and elements from an existing json
     * @param[in] rModelPart The model part of interest to study
     * @param[in] rFilename The name of the files
     * @param[in,out] rRefCondition The conditions of reference
     * @param[in,out] rRefElement The elements of reference
     */
    void WriteReferenceEntitities(
        ModelPart& rModelPart,
        const std::string& rFilename,
        std::unordered_map<IndexType,Condition::Pointer>& rRefCondition,
        std::unordered_map<IndexType,Element::Pointer>& rRefElement
        ) override;

    /**
     * @brief This function generates a list of submodelparts to be able to reassign flags after remesh
     * @param[in,out] rModelPart The model part of interest to study
     */
    void CreateAuxiliarSubModelPartForFlags(ModelPart& rModelPart) override;

    /**
     * @brief This function assigns the flags and clears the auxiliar sub model part for flags
     * @param[in,out] rModelPart The model part of interest to study
     */
    void AssignAndClearAuxiliarSubModelPartForFlags(ModelPart& rModelPart) override;

    std::unordered_map<int, int> GetNodalLocalToGlobalMap();

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
    std::string Info() const // override
    {
        return "ParMmgUtilities";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const // override
    {
        rOStream << "ParMmgUtilities";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const // override
    {
    }

protected:

    ///@name Protected Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

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

    std::unordered_map<int, int> mGlobalToLocalNodePreMap;                           /// Map of nodal global ids to local ids before remeshing
    std::unordered_map<int, int> mGlobalToLocalElemPreMap;                           /// Map of elemental global ids to local ids before remeshing
    std::unordered_map<int, int> mGlobalToLocalCondPreMap;                           /// Map of condition global ids to local ids before remeshing
    std::unordered_map<int, int> mLocalToGlobalNodePostMap;                          /// Map of nodal local ids to global ids after remeshing

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

//     /// Assignment operator.
//     ParMmgUtilities& operator=(ParMmgUtilities const& rOther);

//     /// Copy constructor.
//     ParMmgUtilities(ParMmgUtilities const& rOther);

    ///@}

};// class ParMmgUtilities
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<PMMGLibrary TPMMGLibrary>
inline std::istream& operator >> (std::istream& rIStream,
                                  ParMmgUtilities<TPMMGLibrary>& rThis);

/// output stream function
template<PMMGLibrary TPMMGLibrary>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ParMmgUtilities<TPMMGLibrary>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}// namespace Kratos.
#endif /* KRATOS_PMMG_UTILITIES defined */
