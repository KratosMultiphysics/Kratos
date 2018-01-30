// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
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

    // Containers definition
    typedef ModelPart::NodesContainerType                        NodesArrayType;
    typedef ModelPart::ElementsContainerType                  ElementsArrayType;
    typedef ModelPart::ConditionsContainerType              ConditionsArrayType;
    
    // Components definition
    typedef Node <3>                                                   NodeType;
    typedef Properties                                           PropertiesType;
    typedef Element                                                 ElementType;
    typedef Condition                                             ConditionType;
    
    // Index defintion
    typedef std::size_t                                               IndexType;
    typedef std::size_t                                                SizeType;
    
    // DoF definition
    typedef Dof<double>                                                 DofType;
    
    // Mesh definition
    typedef Mesh<NodeType, PropertiesType, ElementType, ConditionType> MeshType;
    typedef MeshType::PropertiesContainerType           PropertiesContainerType;
    typedef MeshType::NodeConstantIterator                 NodeConstantIterator;
    typedef MeshType::ConditionConstantIterator       ConditionConstantIterator;
    typedef MeshType::ElementConstantIterator           ElementConstantIterator;
    
    #if !defined(HASH_COMBINE)
    #define HASH_COMBINE
    template <class TClassType>
    inline void HashCombine(std::size_t& Seed, const TClassType& Value)
    {
        std::hash<TClassType> hasher;
        Seed ^= hasher(Value) + 0x9e3779b9 + (Seed<<6) + (Seed>>2);
    }
    #endif
    
    #if !defined(HASH_RANGE)
    #define HASH_RANGE
    template <class TClassType>
    inline std::size_t HashRange(TClassType First, TClassType Last)
    {
        std::size_t seed = 0;

        while (First!=Last)
        {
            HashCombine(seed, *First);
            ++First;
        }
        
        return seed;
    }
    #endif
    
    #if !defined(KEY_COMPAROR_RANGE)
    #define KEY_COMPAROR_RANGE
    template<class TClassType>
    struct KeyComparorRange
    {
        bool operator()(const TClassType& lhs, const TClassType& rhs) const
        {
            if(lhs.size() != rhs.size())
            {
                return false;
            }

            auto it_lhs = lhs.begin();
            auto it_rhs = rhs.begin();

            while(it_lhs != lhs.end()) // NOTE: We already checked that are same size
            {
                if(*it_lhs != *it_rhs) 
                {
                    return false;
                }
                if(it_lhs != lhs.end())
                {
                    ++it_lhs;
                    ++it_rhs;
                }
            }

            return true;
        }
    };
    #endif
    
    #if !defined(KEY_HASHER_RANGE)
    #define KEY_HASHER_RANGE
    template<class TClassType>
    struct KeyHasherRange
    {
        std::size_t operator()(const TClassType& rRange) const
        {
            return HashRange(rRange.begin(), rRange.end());
        }
    };
    #endif

///@}
///@name  Enum's
///@{

    /**
     * This enums are used to simplify the computation of the std::vector containing the conditions and elements
     */
    #if !defined(MMG_GEOMETRY)
    #define MMG_GEOMETRY
        enum CondGeometries2D {Line = 0};
        
        enum ElemGeometries2D {Triangle2D = 0};
        
        enum CondGeometries3D {Triangle3D = 0, Quadrilateral3D = 1};
        
        enum ElemGeometries3D {Tetrahedra = 0, Prism = 1};
    #endif
    
    #if !defined(FRAMEWORK_EULER_LAGRANGE)
    #define FRAMEWORK_EULER_LAGRANGE
        enum FrameworkEulerLagrange {Eulerian = 0, Lagrangian = 1};
    #endif
    
///@}
///@name  Functions
///@{
    
///@}
///@name Kratos Classes
///@{

//// This class is a remesher which uses the MMG library 
// The class uses a class for the 2D and 3D cases 

template<unsigned int TDim>  
class MmgProcess 
    : public Process
{
public:

    ///@name Type Definitions
    ///@{
    
    ///@}
    ///@name Life Cycle
    ///@{
     
    // Constructor
    
    /**
     * This is the default constructor, which is used to read the input files 
     * @param rThisModelPart: The model part
     * @param ThisParameters: The parameters
     */
    
    MmgProcess(ModelPart& rThisModelPart, Parameters ThisParameters = Parameters(R"({})"));

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
     * Instead of using an files already created we read an existing model part
     */
    
    void Execute() override;

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
    
    // The model part to compute
    ModelPart& mrThisModelPart;                      
    
    // The parameters (can be used for general pourposes)
    Parameters mThisParameters;
    
    // Storage for the dof of the node
    Node<3>::DofsContainerType  mDofs;
    
    // I/O information
    char* mFilename;
    std::string mStdStringFilename;
    unsigned int mEchoLevel;
    
    // The framework
    FrameworkEulerLagrange mFramework;
    
    // Where the sub model parts IDs are stored
    std::unordered_map<int,std::vector<std::string>> mColors;
    
    // Reference element and condition
    std::unordered_map<int,Element::Pointer>   mpRefElement; 
    std::unordered_map<int,Condition::Pointer> mpRefCondition;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * This function generates the mesh MMG5 structure from a Kratos Model Part
     */
    
    void InitializeMeshData();
    
    /**
     * This function generates the metric MMG5 structure from a Kratos Model Part
     */
    
    void InitializeSolData();
    
    /**
     * We execute the MMg library and build the new model part from the old model part
     */
    
    void ExecuteRemeshing();
    
    /**
     * This function reorder the nodes, conditions and elements to avoid problems with non-consecutive ids
     */
    
    void ReorderAllIds();
    
    /**
     * After we have transfer the information from the previous modelpart we initilize the elements and conditions
     */
    
    void InitializeElementsAndConditions();
    
    /**
     * It checks if the nodes are repeated and remove the repeated ones
     */
    
    std::vector<unsigned int> CheckNodes();
    
    /**
     * It checks if the conditions are repeated and remove the repeated ones
     */
    
    std::vector<unsigned int> CheckConditions0();
    std::vector<unsigned int> CheckConditions1();
    
    /**
     * It checks if the elemenst are removed and remove the repeated ones
     */
    
    std::vector<unsigned int> CheckElements0();
    std::vector<unsigned int> CheckElements1();
    
    /**
     * It blocks certain nodes before remesh the model
     * @param iNode: The index of the noode
     */
    
    void BlockNode(unsigned int iNode);
    
    /**
     * It creates the new node
     * @param iNode: The index of the new noode
     * @param Ref: The submodelpart id
     * @param IsRequired: MMG value (I don't know that it does)
     * @return pNode: The pointer to the new node created
     */
    
    NodeType::Pointer CreateNode(
        unsigned int iNode,
        int& Ref, 
        int& IsRequired
        );
    
    /**
     * It creates the new condition
     * @param CondId: The id of the condition
     * @param PropId: The submodelpart id
     * @param IsRequired: MMG value (I don't know that it does)
     * @return pCondition: The pointer to the new condition created
     */
    
    ConditionType::Pointer CreateCondition0(
        const unsigned int CondId,
        int& PropId, 
        int& IsRequired,
        bool SkipCreation
        );
    
    ConditionType::Pointer CreateCondition1(
        const unsigned int CondId,
        int& PropId, 
        int& IsRequired,
        bool SkipCreation
        );
    
    /**
     * It creates the new element
     * @param CondId: The id of the element
     * @param PropId: The submodelpart id
     * @param IsRequired: MMG value (I don't know that it does)
     * @return pElement: The pointer to the new condition created
     */
    
    ElementType::Pointer CreateElement0(
        const unsigned int ElemId,
        int& PropId, 
        int& IsRequired,
        bool SkipCreation
        );
    
    ElementType::Pointer CreateElement1(
        const unsigned int ElemId,
        int& PropId, 
        int& IsRequired,
        bool SkipCreation
        );
    
    /**
     * It saves the solution and mesh to files (for debugging pourpose g.e)
     * @param PostOutput: If the file to save is after or before remeshing
     */
    
    void SaveSolutionToFile(const bool PostOutput);
    
    /**
     * It frees the memory used during all the process
     */
    
    void FreeMemory();
    
    /** 
     * Initialisation of mesh and sol structures args of InitMesh:
     * @param MMG5_ARG_start: we start to give the args of a variadic func
     * @param MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
     * @param &mmgMesh: pointer toward your MMG5_pMesh (that store your mesh)
     * @param MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
     * @param &mmgSol: pointer toward your MMG5_pSol (that store your metric) 
     */
    
    void InitMesh();
    
    /** 
     * Here the verbosity is set 
     */
    
    void InitVerbosity();
    
    /** 
     * Here the verbosity is set using the API
     * @param verbosityMMG: The equivalent verbosity level in the MMG API
     */
        
    void InitVerbosityParameter(const int& VerbosityMMG);
    
    /**
     * This sets the size of the mesh
     * @param NumNodes: Number of nodes
     * @param NumElements: Number of Elements
     * @param NumConditions: Number of Conditions
     */
    
    void SetMeshSize(
        const SizeType NumNodes,
        const array_1d<SizeType, TDim - 1> NumArrayElements,  // NOTE: We do this tricky thing to take into account the prisms
        const array_1d<SizeType, TDim - 1> NumArrayConditions // NOTE: We do this tricky thing to take into account the quadrilaterals
        );
    
    /**
     * This sets the size of the solution for the scalar case
     * @param NumNodes: Number of nodes
     */
    
    void SetSolSizeScalar(const int NumNodes);
    
    /**
     * This sets the size of the solution for the vector case
     * @param NumNodes: Number of nodes
     */
    
    void SetSolSizeVector(const int NumNodes);
    
    /**
     * This sets the size of the solution for the tensor case
     * @param NumNodes: Number of nodes
     */
    
    void SetSolSizeTensor(const int NumNodes);
    
    /**
     * This checks the mesh data and prints if it is OK
     */
    
    void CheckMeshData();
    
    /**
     * This sets the output mesh
     * @param PostOutput: If the ouput file is the solution after take into account the metric or not
     * @param step: The step to postprocess
     */
    
    void OutputMesh(
        const bool PostOutput, 
        const unsigned int Step
        );
    
    /**
     * This sets the output mesh in a .mdpa format
     */
    void OutputMdpa();

    /**
     * This sets the output sol
     * @param PostOutput: If the ouput file is the solution after take into account the metric or not
     * @param step: The step to postprocess
     */
    
    void OutputSol(
        const bool PostOutput, 
        const unsigned int Step
        );
    
    /**
     * This loads the solution
     */
    
    void MMGLibCall();
    
    /**
     * This frees the MMG structures
     */
    
    void FreeAll();
    
    /**
     * This sets the nodes of the mesh
     * @param X: Coordinate X
     * @param Y: Coordinate Y
     * @param Z: Coordinate Z
     * @param Color: Reference of the node(submodelpart)
     * @param Index: The index number of the node 
     */
    
    void SetNodes(
        const double X,
        const double Y,
        const double Z,
        const int Color,
        const int Index
        );
    
    /**
     * This sets the conditions of the mesh
     * @param Geom: The geometry of the condition
     * @param Color: Reference of the node(submodelpart)
     * @param Index: The index number of the node 
     */
    
    void SetConditions(
        Geometry<Node<3> > & Geom,
        const int Color,
        const int Index
        );
    
    /**
     * This sets elements of the mesh
     * @param Geom: The geometry of the element
     * @param Color: Reference of the node(submodelpart)
     * @param Index: The index number of the node 
     */
    
    void SetElements(
        Geometry<Node<3> > & Geom,
        const int Color,
        const int Index
        );
    
    /**
     * This functions gets the "colors", parts of a model part to process
     * @param NodeColors: Map where the submodelparts and nodes are stored
     * @param CondColors: Map where the submodelparts and conditions are stored
     * @param ElemColors: Map where the submodelparts and elements are stored
     */
    
    void ComputeColors(
        std::unordered_map<int,int>& NodeColors,
        std::unordered_map<int,int>& CondColors,
        std::unordered_map<int,int>& ElemColors
        );

    /**
     * This function is used to compute the metric scalar
     * @param Metric: The inverse of the size node
     */

    void SetMetricScalar(
        const double& Metric,
        const int NodeId 
        );
    
    /**
     * This function is used to compute the metric vector (x, y, z)
     * @param Metric: This array contains the components of the metric vector
     */

    void SetMetricVector(
        const array_1d<double, 3>& Metric,
        const int NodeId 
        );
    
    /**
     * This function is used to compute the Hessian metric tensor, note that when using the Hessian, more than one metric can be defined simultaneously, so in consecuence we need to define the elipsoid which defines the volume of maximal intersection
     * @param Metric: This array contains the components of the metric tensor in the MMG defined order
     */

    void SetMetricTensor(
        const Vector& Metric,
        const int NodeId 
        );
    
    /**
     * This converts the framework string to an enum
     * @param Str: The string
     * @return FrameworkEulerLagrange: The equivalent enum
     */
        
    FrameworkEulerLagrange ConvertFramework(const std::string& Str);

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
// inline std::istream& operator >> (std::istream& rIStream,
//                                   MmgProcess& rThis);
//
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const MmgProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
//
//     return rOStream;
// }

}// namespace Kratos.
#endif /* KRATOS_MMG_PROCESS defined */
