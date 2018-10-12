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

#if !defined(KRATOS_MMG_PROCESS)
#define KRATOS_MMG_PROCESS

// System includes
#include <unordered_map>

// External includes
// The includes related with the MMG library
// #include "mmg/libmmg.h"

// Project includes
#include "processes/process.h"
#include "includes/key_hash.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "containers/variables_list.h"

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

    /**
     * @brief This enum defines the type of MMG libray used
     */
    enum class MMGLibray
    {
        MMG2D = 0,
        MMG3D = 1,
        MMGS  = 2
    };

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class MmgProcess
 * @ingroup MeshingApplication
 * @brief This class is a remesher which uses the MMG library
 * @details This class is a remesher which uses the MMG library. The class uses a class for the 2D and 3D cases.
 * The remesher keeps the previous submodelparts and interpolates the nodal values between the old and new mesh
 * @author Vicente Mataix Ferrandiz
 */
template<MMGLibray TMMGLibray>
class MmgProcess
    : public Process
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of MmgProcess
    KRATOS_CLASS_POINTER_DEFINITION(MmgProcess);

    /// Node containers definition
    typedef ModelPart::NodesContainerType                        NodesArrayType;
    /// Elements containers definition
    typedef ModelPart::ElementsContainerType                  ElementsArrayType;
    /// Conditions containers definition
    typedef ModelPart::ConditionsContainerType              ConditionsArrayType;

    /// Node definition
    typedef Node <3>                                                   NodeType;
    // Geometry definition
    typedef Geometry<NodeType>                                     GeometryType;
    /// Properties definition
    typedef Properties                                           PropertiesType;
    /// Element definition
    typedef Element                                                 ElementType;
    /// Condition definition
    typedef Condition                                             ConditionType;

    /// DoF definition
    typedef Dof<double>                                                 DofType;

    /// Mesh definition
    typedef Mesh<NodeType, PropertiesType, ElementType, ConditionType> MeshType;
    /// Properties container definition
    typedef MeshType::PropertiesContainerType           PropertiesContainerType;
    /// Nodes container definition
    typedef MeshType::NodeConstantIterator                 NodeConstantIterator;
    /// Conditions container definition
    typedef MeshType::ConditionConstantIterator       ConditionConstantIterator;
    /// Elements container definition
    typedef MeshType::ElementConstantIterator           ElementConstantIterator;

    /// Conditions array size
    static constexpr SizeType Dimension = (TMMGLibray == MMGLibray::MMG2D) ? 2 : 3;

    /// Conditions array size
    static constexpr SizeType ConditionsArraySize = (Dimension == 2) ? 1 : 2;

    /// Elements array size
    static constexpr SizeType ElementsArraySize = (Dimension == 2) ? 1 : 2;

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

    ///@}
    ///@name  Enum's
    ///@{

    /**
     * @brief This enums allows to differentiate the working framework
     */
    enum class FrameworkEulerLagrange {EULERIAN = 0, LAGRANGIAN = 1, ALE = 2};

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor

    /**
     * @brief This is the default constructor, which is used to read the input files
     * @param rThisModelPart The model part
     * @param ThisParameters The parameters
     */
    MmgProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    ~MmgProcess() override = default;

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

    void operator()();

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Execute method is used to execute the Process algorithms.
     */
    void Execute() override;

    /**
     * @brief This function is designed for being execute once before the solution loop but after all of the solvers where built
     */
    void ExecuteInitialize() override;

    /**
     * @brief This function is designed for being execute once before the solution loop but after all of the solvers where built
     */
    void ExecuteBeforeSolutionLoop() override;

    /**
     * @brief This function will be executed at every time step BEFORE performing the solve phase
     */
    void ExecuteInitializeSolutionStep() override;

    /**
     * @brief This function will be executed at every time step AFTER performing the solve phase
     */
    void ExecuteFinalizeSolutionStep() override;

    /**
     * @brief This function will be executed at every time step BEFORE  writing the output
     */
    void ExecuteBeforeOutputStep() override;

    /**
     * @brief This function will be executed at every time step AFTER writing the output
     */
    void ExecuteAfterOutputStep() override;

    /**
     * @brief This function is designed for being called at the end of the computations right after reading the model and the groups
     */
    void ExecuteFinalize() override;

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
    std::string Info() const override
    {
        return "MmgProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MmgProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

protected:

    ///@name Protected static Member Variables
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

    ModelPart& mrThisModelPart;                                      /// The model part to compute
    Parameters mThisParameters;                                      /// The parameters (can be used for general pourposes)
    NodeType::DofsContainerType  mDofs;                              /// Storage for the dof of the node

    char* mFilename;                                                 /// I/O file name
    std::string mStdStringFilename;                                  /// I/O file name (string)
    IndexType mEchoLevel;                                            /// The echo level

    FrameworkEulerLagrange mFramework;                               /// The framework

    std::unordered_map<IndexType,std::vector<std::string>> mColors;  /// Where the sub model parts IDs are stored

    std::unordered_map<IndexType,Element::Pointer>   mpRefElement;   /// Reference condition
    std::unordered_map<IndexType,Condition::Pointer> mpRefCondition; /// Reference element

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This converts the framework string to an enum
     * @param Str The string
     * @return FrameworkEulerLagrange: The equivalent enum
     */

    static inline FrameworkEulerLagrange ConvertFramework(const std::string& Str)
    {
        if(Str == "Lagrangian" || Str == "LAGRANGIAN")
            return FrameworkEulerLagrange::LAGRANGIAN;
        else if(Str == "Eulerian" || Str == "EULERIAN")
            return FrameworkEulerLagrange::EULERIAN;
        else if(Str == "ALE")
            return FrameworkEulerLagrange::ALE;
        else
            return FrameworkEulerLagrange::EULERIAN;
    }

    /**
     * @brief This function generates the mesh MMG5 structure from a Kratos Model Part
     */

    void InitializeMeshData();

    /**
     *@brief This function generates the metric MMG5 structure from a Kratos Model Part
     */

    void InitializeSolData();

    /**
     * @brief We execute the MMg library and build the new model part from the old model part
     */

    void ExecuteRemeshing();

    /**
     * @brief This function reorder the nodes, conditions and elements to avoid problems with non-consecutive ids
     */

    void ReorderAllIds();

    /**
     * @brief After we have transfer the information from the previous modelpart we initilize the elements and conditions
     */

    void InitializeElementsAndConditions();

    /**
     * @brief It checks if the nodes are repeated and remove the repeated ones
     */

    IndexVectorType CheckNodes();

    /**
     * @brief It checks if the conditions are repeated and remove the repeated ones
     */

    IndexVectorType CheckConditions0();

    /**
     * @brief It checks if the conditions are repeated and remove the repeated ones
     */

    IndexVectorType CheckConditions1();

    /**
     * @brief It checks if the elemenst are removed and remove the repeated ones
     */

    IndexVectorType CheckElements0();

    /**
     * @brief It checks if the elemenst are removed and remove the repeated ones
     */

    IndexVectorType CheckElements1();

    /**
     * @brief It blocks certain nodes before remesh the model
     * @param iNode The index of the noode
     */

    void BlockNode(IndexType iNode);

    /**
     * @brief It creates the new node
     * @param iNode The index of the new noode
     * @param Ref The submodelpart id
     * @param IsRequired MMG value (I don't know that it does)
     * @return pNode The pointer to the new node created
     */

    NodeType::Pointer CreateNode(
        IndexType iNode,
        int& Ref,
        int& IsRequired
        );

    /**
     * @brief It creates the new condition
     * @param CondId The id of the condition
     * @param PropId The submodelpart id
     * @param IsRequired MMG value (I don't know that it does)
     * @return pCondition The pointer to the new condition created
     */

    ConditionType::Pointer CreateCondition0(
        const IndexType CondId,
        int& PropId,
        int& IsRequired,
        bool SkipCreation
        );

    /**
     * @brief It creates the new condition
     * @param CondId The id of the condition
     * @param PropId The submodelpart id
     * @param IsRequired MMG value (I don't know that it does)
     * @return pCondition The pointer to the new condition created
     */

    ConditionType::Pointer CreateCondition1(
        const IndexType CondId,
        int& PropId,
        int& IsRequired,
        bool SkipCreation
        );

    /**
     * @brief It creates the new element
     * @param ElemId The id of the element
     * @param PropId The submodelpart id
     * @param IsRequired MMG value (I don't know that it does)
     * @return pElement The pointer to the new condition created
     */

    ElementType::Pointer CreateElement0(
        const IndexType ElemId,
        int& PropId,
        int& IsRequired,
        bool SkipCreation
        );

    /**
     * @brief It creates the new element
     * @param ElemId The id of the element
     * @param PropId The submodelpart id
     * @param IsRequired MMG value (I don't know that it does)
     * @return pElement The pointer to the new condition created
     */

    ElementType::Pointer CreateElement1(
        const IndexType ElemId,
        int& PropId,
        int& IsRequired,
        bool SkipCreation
        );

    /**
     * @brief It saves the solution and mesh to files (for debugging pourpose g.e)
     * @param PostOutput If the file to save is after or before remeshing
     */

    void SaveSolutionToFile(const bool PostOutput);

    /**
     * @brief It frees the memory used during all the process
     */

    void FreeMemory();

    /**
     * @brief Initialisation of mesh and sol structures
     * @detauils Initialisation of mesh and sol structures args of InitMesh:
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
     * @param VerbosityMMG The equivalent verbosity level in the MMG API
     */

    void InitVerbosityParameter(const IndexType VerbosityMMG);

    /**
     * @brief This sets the size of the mesh
     * @param NumNodes Number of nodes
     * @param NumArrayElements Number of Elements
     * @param NumArrayConditions Number of Conditions
     */

    void SetMeshSize(
        const SizeType NumNodes,
        const array_1d<SizeType, ElementsArraySize>& NumArrayElements,
        const array_1d<SizeType, ConditionsArraySize>& NumArrayConditions
        );

    /**
     * @brief This sets the size of the solution for the scalar case
     * @param NumNodes Number of nodes
     */

    void SetSolSizeScalar(const SizeType NumNodes);

    /**
     * @brief This sets the size of the solution for the vector case
     * @param NumNodes Number of nodes
     */

    void SetSolSizeVector(const SizeType NumNodes);

    /**
     * @brief This sets the size of the solution for the tensor case
     * @param NumNodes Number of nodes
     */

    void SetSolSizeTensor(const SizeType NumNodes);

    /**
     * @brief This checks the mesh data and prints if it is OK
     */

    void CheckMeshData();

    /**
     * @brief This sets the output mesh
     * @param PostOutput If the ouput file is the solution after take into account the metric or not
     * @param Step The step to postprocess
     */

    void OutputMesh(
        const bool PostOutput,
        const IndexType Step
        );

    /**
     * @brief This sets the output mesh in a .mdpa format
     */
    void OutputMdpa();

    /**
     * @brief This sets the output sol
     * @param PostOutput If the ouput file is the solution after take into account the metric or not
     * @param Step The step to postprocess
     */
    void OutputSol(
        const bool PostOutput,
        const IndexType Step
        );

    /**
     * @brief This loads the solution
     */
    void MMGLibCall();

    /**
     * @brief This frees the MMG structures
     */
    void FreeAll();

    /**
     * @brief This sets the nodes of the mesh
     * @param X Coordinate X
     * @param Y Coordinate Y
     * @param Z Coordinate Z
     * @param Color Reference of the node(submodelpart)
     * @param Index The index number of the node
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
     * @param Geom The geometry of the condition
     * @param Color Reference of the node(submodelpart)
     * @param Index The index number of the node
     */
    void SetConditions(
        GeometryType& Geom,
        const IndexType Color,
        const IndexType Index
        );

    /**
     * @brief This sets elements of the mesh
     * @param Geom The geometry of the element
     * @param Color Reference of the node(submodelpart)
     * @param Index The index number of the node
     */
    void SetElements(
        GeometryType& Geom,
        const IndexType Color,
        const IndexType Index
        );

    /**
     * @brief This function is used to compute the metric scalar
     * @param Metric The inverse of the size node
     */
    void SetMetricScalar(
        const double Metric,
        const IndexType NodeId
        );

    /**
     * @brief This function is used to compute the metric vector (x, y, z)
     * @param Metric This array contains the components of the metric vector
     */
    void SetMetricVector(
        const array_1d<double, Dimension>& Metric,
        const IndexType NodeId
        );

    /**
     * @brief This function is used to compute the Hessian metric tensor, note that when using the Hessian, more than one metric can be defined simultaneously, so in consecuence we need to define the elipsoid which defines the volume of maximal intersection
     * @param Metric This array contains the components of the metric tensor in the MMG defined order
     */
    void SetMetricTensor(
        const TensorArrayType& Metric,
        const IndexType NodeId
        );

    /**
     * @brief This function generates a list of submodelparts to be able to reassign flags after remesh
     */
    void CreateAuxiliarSubModelPartForFlags();

    /**
     * @brief This function assigns the flags and clears the auxiliar sub model part for flags
     */
    void AssignAndClearAuxiliarSubModelPartForFlags();

    /**
     * @brief This function creates an before/after remesh output file
     * @param rOldModelPart The old model part before remesh
     */
    void CreateDebugPrePostRemeshOutput(ModelPart& rOldModelPart);

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
//     MmgProcess& operator=(MmgProcess const& rOther);

//     /// Copy constructor.
//     MmgProcess(MmgProcess const& rOther);

    ///@}

};// class MmgProcess
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<MMGLibray TMMGLibray>
inline std::istream& operator >> (std::istream& rIStream,
                                  MmgProcess<TMMGLibray>& rThis);

/// output stream function
template<MMGLibray TMMGLibray>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MmgProcess<TMMGLibray>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}// namespace Kratos.
#endif /* KRATOS_MMG_PROCESS defined */
