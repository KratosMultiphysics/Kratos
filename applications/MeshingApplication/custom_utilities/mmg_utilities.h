// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_MMG_UTILITIES)
#define KRATOS_MMG_UTILITIES

// System includes

// External includes

// Project includes
#include "meshing_application.h"
#include "includes/key_hash.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/variable_utils.h"
#include "utilities/assign_unique_model_part_collection_tag_utility.h"
#include "processes/fast_transfer_between_model_parts_process.h"

// NOTE: The following contains the license of the MMG library
/* =============================================================================
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
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
 * @struct MMGMeshInfo
 * @ingroup MeshingApplication
 * @brief Stores the Mmg mesh information
 * @author Vicente Mataix Ferrandiz
 */
template<MMGLibrary TMMGLibrary>
struct MMGMeshInfo
{
    // Info stored
    SizeType NumberOfNodes;
    SizeType NumberOfLines;
    SizeType NumberOfTriangles;
    SizeType NumberOfQuadrilaterals;
    SizeType NumberOfPrism;
    SizeType NumberOfTetrahedra;

    /**
     * @brief It returns the number of the first type of conditions
     */
    SizeType NumberFirstTypeConditions() const;

    /**
     * @brief It returns the number of the second type of conditions
     */
    SizeType NumberSecondTypeConditions() const;

    /**
     * @brief It returns the number of the first type of elements
     */
    SizeType NumberFirstTypeElements() const;

    /**
     * @brief It returns the number of the second type of elements
     */
    SizeType NumberSecondTypeElements() const;
};

/**
 * @class MmgUtilities
 * @ingroup MeshingApplication
 * @brief Provides the Kratos interface to the MMG library API
 * @author Vicente Mataix Ferrandiz
 */
template<MMGLibrary TMMGLibrary>
class KRATOS_API(MESHING_APPLICATION) MmgUtilities
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of MmgUtilities
    KRATOS_CLASS_POINTER_DEFINITION(MmgUtilities);

    /// Node definition
    typedef Node <3>                                                   NodeType;
    // Geometry definition
    typedef Geometry<NodeType>                                     GeometryType;

    /// Spatial dimension
    static constexpr SizeType Dimension = (TMMGLibrary == MMGLibrary::MMG2D) ? 2 : 3;

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
     * @brief This method sets the echo level
     * @param[in] EchoLevel Sets the echo level
     */
    void SetEchoLevel(const SizeType EchoLevel);

    /**
     * @brief This method gets the echo level
     * @return mEchoLevel Gets the echo level
     */
    SizeType GetEchoLevel();

    /**
     * @brief This method sets the discretization method
     * @param[in] Discretization Sets the discretization method
     */
    void SetDiscretization(const DiscretizationOption Discretization);

    /**
     * @brief This method gets the discretization method
     * @return mDiscretization Gets the discretization method
     */
    DiscretizationOption GetDiscretization();

    /**
     * @brief This method sets if the regions must be removed
     * @param[in] RemoveRegions Sets if the regions must be removed
     */
    void SetRemoveRegions(const bool RemoveRegions);

    /**
     * @brief This method gets if the regions must be removed
     * @return mRemoveRegions Gets if the regions must be removed
     */
    bool GetRemoveRegions();

    /**
     * @brief It prints info about the current mesh
     * @param[in,out] rMMGMeshInfo The number of nodes, conditions and elements
     */
    void PrintAndGetMmgMeshInfo(MMGMeshInfo<TMMGLibrary>& rMMGMeshInfo);

    /**
     * @brief Returns a vector of ids of spatially repeated nodes
     * @param[in,out] rModelPart The model part whose nodes are checked
     */
    IndexVectorType FindDuplicateNodeIds(const ModelPart& rModelPart);

    /**
     * @brief Returns a vector of ids of repeated conditions
     * @details For 2D and surface meshes it returns the ids of repeated edges found in the MMG mesh data structure. In 3D it returns ids of triangles.
     */
    IndexVectorType CheckFirstTypeConditions();

    /**
     * @brief Returns a vector of ids of repeated conditions
     * @details For 3D meshes it returns the ids of repeated quadrilaterals found in the MMG mesh data structure. Otherwise, it returns an empty vector.
     */
    IndexVectorType CheckSecondTypeConditions();

    /**
     * @brief Returns a vector of ids of repeated elements
     * @details For 2D and surface meshes it returns the ids of repeated triangles found in the MMG mesh data structure. In 3D it returns ids of tetrahedra.
     */
    IndexVectorType CheckFirstTypeElements();

    /**
     * @brief Returns a vector of ids of repeated elements
     * @details For 3D meshes it returns the ids of repeated prisms found in the MMG mesh data structure. Otherwise, it returns an empty vector.
     */
    IndexVectorType CheckSecondTypeElements();

    /**
     * @brief It blocks certain nodes before remesh the model
     * @param[in] iNode The index of the node
     */
    void BlockNode(const IndexType iNode);

    /**
     * @brief It blocks certain conditions before remesh the model
     * @param[in] iCondition The index of the condition
     */
    void BlockCondition(const IndexType iCondition);

    /**
     * @brief It blocks certain elements before remesh the model
     * @param[in] iElement The index of the element
     */
    void BlockElement(const IndexType iElement);

    /**
     * @brief It creates the new node
     * @details Each call to this function increments the internal counter of the MMG mesh data structure, extracts the coordinates of the current vertex and creates the new node with the extracted coordinates.
     * @param[in,out] rModelPart The model part which owns the new node
     * @param[in] iNode The index of the new node
     * @param[out] Ref MMG point reference
     * @param[out] IsRequired MMG required entity (0 or 1)
     * @return pNode The pointer to the new node created
     */
    NodeType::Pointer CreateNode(
        ModelPart& rModelPart,
        IndexType iNode,
        int& Ref,
        int& IsRequired
        );

    /**
     * @brief It creates the new condition (first type, depends if the library work in 2D/3D/Surfaces)
     * @details Each call to this function increments the internal counter of the MMG mesh data structure and extracts the vertices of the current  edge.
     * @param[in,out] rModelPart The model part of interest to study
     * @param[in,out] rMapPointersRefCondition The pointer to the condition of reference
     * @param[in] CondId The id of the new condition
     * @param[out] Ref MMG edge reference
     * @param[out] IsRequired MMG required entity (0 or 1)
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
        );

    /**
     * @brief It creates the new condition (second type, depends if the library work in 2D/3D/Surfaces)
     * @param[in,out] rModelPart The model part of interest to study
     * @param[in,out] rMapPointersRefCondition The pointer to the condition of reference
     * @param[in] CondId The id of the condition
     * @param[out] Ref MMG edge reference
     * @param[out] IsRequired MMG required entity (0 or 1)
     * @return pCondition The pointer to the new condition created
     */
    Condition::Pointer CreateSecondTypeCondition(
        ModelPart& rModelPart,
        std::unordered_map<IndexType,Condition::Pointer>& rMapPointersRefCondition,
        const IndexType CondId,
        int& Ref,
        int& IsRequired,
        bool SkipCreation
        );

    /**
     * @brief It creates the new element (first type, depends if the library work in 2D/3D/Surfaces)
     * @param[in,out] rModelPart The model part of interest to study
     * @param[in,out] rMapPointersRefElement The pointer to the element of reference
     * @param[in] ElemId The id of the element
     * @param[out] Ref MMG edge reference
     * @param[out] IsRequired MMG required entity (0 or 1)
     * @return pElement The pointer to the new condition created
     */
    Element::Pointer CreateFirstTypeElement(
        ModelPart& rModelPart,
        std::unordered_map<IndexType,Element::Pointer>& rMapPointersRefElement,
        const IndexType ElemId,
        int& Ref,
        int& IsRequired,
        bool SkipCreation
        );

    /**
     * @brief It creates the new element (second type, depends if the library work in 2D/3D/Surfaces)
     * @param[in,out] rModelPart The model part of interest to study
     * @param[in,out] rMapPointersRefElement The pointer to the element of reference
     * @param[in] ElemId The id of the element
     * @param[out] Ref MMG edge reference
     * @param[out] IsRequired MMG required entity (0 or 1)
     * @return pElement The pointer to the new condition created
     */
    Element::Pointer CreateSecondTypeElement(
        ModelPart& rModelPart,
        std::unordered_map<IndexType,Element::Pointer>& rMapPointersRefElement,
        const IndexType ElemId,
        int& Ref,
        int& IsRequired,
        bool SkipCreation
        );

    /**
     * @brief Initialisation of mesh and sol structures
     * @details Initialisation of mesh and sol structures args of InitMesh:
     * -# MMG5_ARG_start we start to give the args of a variadic func
     * -# MMG5_ARG_ppMesh next arg will be a pointer over a MMG5_pMesh
     * -# &mmgMesh pointer toward your MMG5_pMesh (that store your mesh)
     * -# MMG5_ARG_ppMet next arg will be a pointer over a MMG5_pSol storing a metric
     * -# &mmgSol pointer toward your MMG5_pSol (that store your metric)
     */
    void InitMesh();

    /**
     * @brief Here the verbosity is set
     */
    void InitVerbosity();

    /**
     * @brief Here the verbosity is set using the API
     * @param[in] VerbosityMMG The equivalent verbosity level in the MMG API
     */
    void InitVerbosityParameter(const IndexType VerbosityMMG);

    /**
     * @brief This sets the size of the mesh
     * @param[in,out] rMMGMeshInfo The number of nodes, conditions and elements
     */
    void SetMeshSize(MMGMeshInfo<TMMGLibrary>& rMMGMeshInfo);

    /**
     * @brief This sets the size of the solution for the scalar case
     * @param[in] NumNodes Number of nodes
     */
    void SetSolSizeScalar(const SizeType NumNodes);

    /**
     * @brief This sets the size of the solution for the vector case
     * @param[in] NumNodes Number of nodes
     */
    void SetSolSizeVector(const SizeType NumNodes);

    /**
     * @brief This sets the size of the solution for the tensor case
     * @param[in] NumNodes Number of nodes
     */
    void SetSolSizeTensor(const SizeType NumNodes);

    /**
     * @brief This sets the size of the displacement for lagrangian movement
     * @param[in] NumNodes Number of nodes
     */
    void SetDispSizeVector(const SizeType NumNodes);

    /**
     * @brief This checks the mesh data and prints if it is OK
     */
    void CheckMeshData();

    /**
     * @brief This sets the output mesh
     * @param[in] rInputName The input name
     */
    void InputMesh(const std::string& rInputName);

    /**
     * @brief This sets the output sol
     * @param[in] rInputName The input name
     */
    void InputSol(const std::string& rInputName);

    /**
     * @brief This sets the output mesh
     * @param[in] rOutputName The output name
     */
    void OutputMesh(const std::string& rOutputName);

    /**
     * @brief This sets the output sol
     * @param[in] rOutputName The output name
     */
    void OutputSol(const std::string& rOutputName);

    /**
     * @brief This sets the output displacement
     * @param[in] rOutputName The output name
     */
    void OutputDisplacement(const std::string& rOutputName);

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
        );

    /**
     * @brief This frees the MMG structures
     */
    void FreeAll();

    /**
     * @brief This loads the solution
     */
    void MMGLibCallMetric(Parameters ConfigurationParameters);

    /**
     * @brief This loads the solution
     */
    void MMGLibCallIsoSurface(Parameters ConfigurationParameters);

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
        );

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
        );

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
        );

    /**
     * @brief This function is used to set the metric scalar
     * @param[in] Metric The inverse of the size node
     * @param[in] NodeId The id of the node
     */
    void SetMetricScalar(
        const double Metric,
        const IndexType NodeId
        );

    /**
     * @brief This function is used to set the metric vector (x, y, z)
     * @param[in] rMetric This array contains the components of the metric vector
     * @param[in] NodeId The id of the node
     */
    void SetMetricVector(
        const array_1d<double, Dimension>& rMetric,
        const IndexType NodeId
        );

    /**
     * @brief This function is used to set the Hessian metric tensor, note that when using the Hessian, more than one metric can be defined simultaneously, so in consecuence we need to define the elipsoid which defines the volume of maximal intersection
     * @param[in] rMetric This array contains the components of the metric tensor in the MMG defined order
     * @param[in] NodeId The id of the node
     */
    void SetMetricTensor(
        const TensorArrayType& rMetric,
        const IndexType NodeId
        );

    /**
     * @brief This function is used to set the displacement vector (x, y, z)
     * @param[in] rMetric This array contains the components of the displacement vector
     * @param[in] NodeId The id of the node
     */
    void SetDisplacementVector(
        const array_1d<double, 3>& rDisplacement,
        const IndexType NodeId
        );

    /**
     * @brief This function is used to retrieve the metric scalar
     * @param[in,out] rMetric The inverse of the size node
     */
    void GetMetricScalar(double& rMetric);

    /**
     * @brief This function is used to retrieve the metric vector (x, y, z)
     * @param[in,out] rMetric This array contains the components of the metric vector
     */
    void GetMetricVector(array_1d<double, Dimension>& rMetric);

    /**
     * @brief This function is used to retrieve the Hessian metric tensor, note that when using the Hessian, more than one metric can be defined simultaneously, so in consecuence we need to define the elipsoid which defines the volume of maximal intersection
     * @param[in,out] rMetric This array contains the components of the metric tensor in the MMG defined order
     */
    void GetMetricTensor(TensorArrayType& rMetric);

    /**
     * @brief This function is used to retrieve the displacement vector (x, y, z)
     * @param[in,out] rDisplacement This array contains the components of the displacement vector
     */
    void GetDisplacementVector(array_1d<double, 3>& rDisplacement);

    /**
     * @brief This function reorder the nodes, conditions and elements to avoid problems with non-consecutive ids
     * @param[in,out] rModelPart The model part of interest to study
     */
    void ReorderAllIds(ModelPart& rModelPart);

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
        );

    /**
     * @brief This method generates solution (metric) data from an existing model part
     * @param[in,out] rModelPart The model part of interest to study
     */
    void GenerateSolDataFromModelPart(ModelPart& rModelPart);

    /**
     * @brief This method generates displacement data from an existing model part
     * @param[in,out] rModelPart The model part of interest to study
     */
    void GenerateDisplacementDataFromModelPart(ModelPart& rModelPart);

    /**
     * @brief This method writes mesh data to an existing model part
     * @param[in,out] rModelPart The model part of interest to study
     * @param[in] rColors Where the sub model parts IDs are stored
     * @param[in] rDofs Storage for the dof of the node
     * @param[in] rMMGMeshInfo The resulting mesh data
     * @param[in] rMapPointersRefCondition The map of the conditions by reference (color)
     * @param[in] rMapPointersRefElement The map of the elements by reference (color)
     */
    void WriteMeshDataToModelPart(
        ModelPart& rModelPart,
        const std::unordered_map<IndexType,std::vector<std::string>>& rColors,
        const NodeType::DofsContainerType& rDofs,
        const MMGMeshInfo<TMMGLibrary>& rMMGMeshInfo,
        std::unordered_map<IndexType,Condition::Pointer>& rMapPointersRefCondition,
        std::unordered_map<IndexType,Element::Pointer>& rMapPointersRefElement
        );

    /**
     * @brief This method writes sol data to an existing model part
     * @param[in,out] rModelPart The model part of interest to study
     */
    void WriteSolDataToModelPart(ModelPart& rModelPart);

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
        );

    /**
     * @brief This function generates a list of submodelparts to be able to reassign flags after remesh
     * @param[in,out] rModelPart The model part of interest to study
     */
    void CreateAuxiliarSubModelPartForFlags(ModelPart& rModelPart);

    /**
     * @brief This function assigns the flags and clears the auxiliar sub model part for flags
     * @param[in,out] rModelPart The model part of interest to study
     */
    void AssignAndClearAuxiliarSubModelPartForFlags(ModelPart& rModelPart);

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
        return "MmgUtilities";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const // override
    {
        rOStream << "MmgUtilities";
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

    SizeType mEchoLevel = 0;                                               /// The echo level of the utilities
    bool mRemoveRegions = false;                                           /// Cuttig-out specified regions during surface remeshing
    DiscretizationOption mDiscretization = DiscretizationOption::STANDARD; /// Discretization The discretization type

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Sets a flag according to a given status over all submodelparts
     * @param rFlag flag to be set
     * @param FlagValue flag value to be set
     */
    void ResursivelyAssignFlagEntities(
        ModelPart& rModelPart,
        const Flags& rFlag,
        const bool FlagValue
        )
    {
        // We call it recursively
        for (auto& r_sub_model_part : rModelPart.SubModelParts()) {
            VariableUtils().SetFlag(rFlag, FlagValue, r_sub_model_part.Conditions());
            VariableUtils().SetFlag(rFlag, FlagValue, r_sub_model_part.Elements());
            ResursivelyAssignFlagEntities(r_sub_model_part, rFlag, FlagValue);
        }
    }

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
//     MmgUtilities& operator=(MmgUtilities const& rOther);

//     /// Copy constructor.
//     MmgUtilities(MmgUtilities const& rOther);

    ///@}

};// class MmgUtilities
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<MMGLibrary TMMGLibrary>
inline std::istream& operator >> (std::istream& rIStream,
                                  MmgUtilities<TMMGLibrary>& rThis);

/// output stream function
template<MMGLibrary TMMGLibrary>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MmgUtilities<TMMGLibrary>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}// namespace Kratos.
#endif /* KRATOS_MMG_UTILITIES defined */
