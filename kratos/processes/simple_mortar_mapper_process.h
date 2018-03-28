//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_SIMPLE_MORTAR_MAPPER_PROCESS)
#define KRATOS_SIMPLE_MORTAR_MAPPER_PROCESS

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "processes/process.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "geometries/point.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

/* Custom includes */
#include "includes/mortar_classes.h"

/* Custom utilities */
#include "utilities/exact_mortar_segmentation_utility.h"

/* Tree structures */
// #include "spatial_containers/bounding_volume_tree.h" // k-DOP
#include "spatial_containers/spatial_containers.h" // kd-tree 

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
    
///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
    
/** 
 * @ingroup KratosCore
 * @class PointMapper
 * @brief Custom Point container to be used by the mapper
 * @details The main difference with this point and the base one is that it contains the pointer to condition where the center of the points belongs
 * @author Vicente Mataix Ferrandiz 
 */
class PointMapper
    : public Point
{
public:
    ///@name Type Definitions
    ///@{
    
    typedef Point BaseType; 
    
    /// Counted pointer of PointMapper
    KRATOS_CLASS_POINTER_DEFINITION( PointMapper );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    PointMapper():
        BaseType(),
        mpOriginCond(nullptr)
    {}

    PointMapper(const array_1d<double, 3>& Coords)
        :BaseType(Coords),
         mpOriginCond(nullptr)
    {}
    
    PointMapper(Condition::Pointer pCond):
        mpOriginCond(pCond)
    {
        UpdatePoint();
    }
    
    PointMapper(
        const array_1d<double, 3>& Coords,
        Condition::Pointer pCond
    ):
        BaseType(Coords),
        mpOriginCond(pCond)
    {}

    ///Copy constructor  (not really required)
    PointMapper(const PointMapper& rhs):
        BaseType(rhs),
        mpOriginCond(rhs.mpOriginCond)
    {
    }

    /// Destructor.
    ~PointMapper() override= default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Returns the point
     * @return The point
     */
    BaseType GetPoint()
    {
        BaseType Point(this->Coordinates());
        return Point;
    }
    
    /**
     * @brief Set the point
     * @param Point The point
     */
    void SetPoint(const BaseType Point)
    {
        this->Coordinates() = Point.Coordinates();
    }

    /**
     * @brief Sets the condition associated to the point
     * @param pCond The pointer to the condition
     */
    void SetCondition(Condition::Pointer pCond)
    {
        mpOriginCond = pCond;
    }
    
    /**
     * @brief Returns the condition associated to the point
     * @return mpOriginCond The pointer to the condition associated to the point
     */
    Condition::Pointer GetCondition()
    {
        KRATOS_DEBUG_ERROR_IF(mpOriginCond == nullptr) << "Condition no initialized in the PointMapper class" << std::endl;
        return mpOriginCond;
    }
    
    /**
     * @brief This method checks everything is right
     */
    void Check()
    {
        KRATOS_TRY;
        
        auto aux_coord = Kratos::make_shared<array_1d<double, 3>>(this->Coordinates());
        KRATOS_ERROR_IF(!aux_coord) << "Coordinates no initialized in the PointMapper class" << std::endl;
        KRATOS_ERROR_IF(mpOriginCond == nullptr) << "Condition no initialized in the PointMapper class" << std::endl;
        
        KRATOS_CATCH("");
    }
    
    /**
     * @brief This function updates the database, using as base for the coordinates the condition center
     */
    void UpdatePoint()
    {
        noalias(this->Coordinates()) = mpOriginCond->GetGeometry().Center().Coordinates();
    }

private:
    ///@name Member Variables
    ///@{
    Condition::Pointer mpOriginCond; /// Condition pointer  
    ///@}

}; // Class PointMapper 
    
/**
 * @ingroup KratosCore
 * @class SimpleMortarMapperProcess
 * @brief This is basic mapper of values between domains using mortar formulation
 * @details Using the dual mortar formulation the resolution of the system of equations is not needed. 
 * Several types of constructors are avaible depending of the needs.
 * If the pairs sets are not provided a serach will be performed using a KDTree
 * @author Vicente Mataix Ferrandiz
 */
template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination = THistOrigin> 
class KRATOS_API(KRATOS_CORE) SimpleMortarMapperProcess
        : public Process
{
public:
    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of SimpleMortarMapperProcess
    KRATOS_CLASS_POINTER_DEFINITION(SimpleMortarMapperProcess);
    
    typedef Point                                        PointType;
    typedef Node<3>                                       NodeType;
    typedef Geometry<NodeType>                        GeometryType;
    typedef Geometry<PointType>                  GeometryPointType;
    typedef ModelPart::NodesContainerType           NodesArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    
    /// Type definition for integration methods
    typedef GeometryData::IntegrationMethod      IntegrationMethod;
    
    /// Auxiliar geometries
    typedef Line2D2<PointType>                            LineType;
    typedef Triangle3D3<PointType>                    TriangleType;
    typedef typename std::conditional<TDim == 2, LineType, TriangleType >::type DecompType;
    
    /// Component type
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;

    /// Linear solver
    typedef UblasSpace<double, CompressedMatrix, Vector>    SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector>               LocalSpaceType;
    typedef typename SparseSpaceType::MatrixType                 MatrixType;
    typedef typename SparseSpaceType::VectorType                 VectorType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    
    /// Component type
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentType;  
    
    /// An integer map
    typedef std::unordered_map<int, int>                             IntMap;
    
    /// BoundedMatrix
    typedef bounded_matrix<double, TNumNodes, TNumNodes>  BoundedMatrixType;

    // Type definitions for the tree
    typedef PointMapper                                     PointMapperType;
    typedef PointMapperType::Pointer                       PointTypePointer;
    typedef std::vector<PointTypePointer>                       PointVector;
    typedef PointVector::iterator                             PointIterator;
    typedef std::vector<double>                              DistanceVector;
    typedef DistanceVector::iterator                       DistanceIterator;
    
    // KDtree definitions
    typedef Bucket< 3ul, PointMapperType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTreeType;
    
    ///@}
    ///@name Life Cycle
    ///@{
    
    /**
     * @brief Default constructor
     * @param rOriginModelPart The origin model part to compute 
     * @param rDestinationModelPart The destination model part to compute 
     * @param ThisVariable The variable to transfer and be transfered
     * @param ThisParameters The configuration parameters
     * @param pThisLinearSolver The pointer to the linear to be used (in case of implicit resolution)
     */
    SimpleMortarMapperProcess( 
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        TVarType& ThisVariable, 
        Parameters ThisParameters = Parameters(R"({})" ),
        LinearSolverType::Pointer pThisLinearSolver = nullptr
        );
    
    /**
     * @brief A constructor where two different variables can be considered for each subdomain
     * @param rOriginModelPart The origin model part to compute 
     * @param rDestinationModelPart The destination model part to compute 
     * @param OriginVariable The variable to transfer
     * @param DestinationVariable The variable to be transfered
     * @param ThisParameters The configuration parameters
     * @param pThisLinearSolver The pointer to the linear to be used (in case of implicit resolution)
     */
    SimpleMortarMapperProcess( 
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        TVarType& OriginVariable,
        TVarType& DestinationVariable,
        Parameters ThisParameters = Parameters(R"({})" ),
        LinearSolverType::Pointer pThisLinearSolver = nullptr
        );

    /// Destructor.
    ~SimpleMortarMapperProcess() override = default;

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

    void operator()()
    {
        Execute();
    }
    
    ///@}
    ///@name Operations
    ///@{
    
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
        return "SimpleMortarMapperProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SimpleMortarMapperProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

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
    
    ModelPart& mOriginModelPart;                  /// The origin model part to compute
    ModelPart& mDestinationModelPart;             /// The destination model part to compute
    TVarType mOriginVariable;                     /// The origin variable to map
    TVarType mDestinationVariable;                /// The destiny variable to map
    
    unsigned int mEchoLevel;                      /// The verbosity level    
    Parameters mThisParameters;                   /// The configuration parameters
    
    LinearSolverType::Pointer mpThisLinearSolver; // The linear solver used to compute the solution

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    /**
     * @brief Check if the pairs has been created
     */
    void CheckAndPerformSearch();
    
    /**
     * @brief This method resets the nodal area
     */
    void ResetNodalArea();
    
    /**
     * @brief This method gets the max area of the conditions from the modelpart
     */
    double GetReferenceArea();
    
    /**
     * @brief This method assemble locally the mortar operators
     * @param ConditionsPointSlave The list of points that form the triangle decomposition
     * @param SlaveGeometry The slave geometry
     * @param MasterGeometry The master geometry
     * @param MasterNormal The normal vector of the master geometry
     * @param ThisKinematicVariables The kinematic variables of the geometries, needed to integrate the mortar operators
     * @param ThisMortarOperators The mortar operators
     * @param ThisIntegrationMethod The integration method used, determines the integration order
     * @param Ae The dual lagrange multiplier operator
     */
    void AssemblyMortarOperators(
        const std::vector<array_1d<PointType,TDim>>& ConditionsPointSlave,
        GeometryType& SlaveGeometry,
        GeometryType& MasterGeometry,
        const array_1d<double, 3>& MasterNormal,
        MortarKinematicVariables<TNumNodes>& ThisKinematicVariables,
        MortarOperator<TNumNodes>& ThisMortarOperators,
        const IntegrationMethod& ThisIntegrationMethod,
        const BoundedMatrixType Ae = IdentityMatrix(TNumNodes)
        );
    
    /**
     * @brief This method computes the Ae matrix
     * @param SlaveGeometry The slave geometry
     * @param ThisKinematicVariables The kinematic variables
     * @param ConditionsPointsSlave The list of decomposed triangles
     * @param ThisIntegrationMethod The integration method considered
     * @return Ae: The matrix of dual LM
     */
    static inline BoundedMatrixType CalculateAe(
        GeometryType& SlaveGeometry,
        MortarKinematicVariables<TNumNodes>& ThisKinematicVariables,
        std::vector<array_1d<PointType,TDim>>& ConditionsPointsSlave,
        const IntegrationMethod& ThisIntegrationMethod
        );
        
    /**
     * @brief This method inverts a diagonal matrix
     * @param InputMatrix The matrix to invert
     * @return The matrix inverted
     */
    static inline BoundedMatrixType InvertDiagonalMatrix(const BoundedMatrixType& InputMatrix);

    /**
     * @brief This method inverts a diagonal matrix
     * @param InputMatrix The matrix to invert
     * @param InvertedMatrix The matrix inverted
     */
    static inline void InvertDiagonalMatrix(
        const BoundedMatrixType& InputMatrix,
        BoundedMatrixType& InvertedMatrix
        );
    
    /**
     * @brief This method lumps a matrix
     * @param InputMatrix The matrix to lump
     */
    void LumpMatrix(BoundedMatrixType& InputMatrix);
    
    /**
     * @brief This method computes the size of the system
     * @param SizeSystem The size of the system
     */
        
    void GetSystemSize(std::size_t& SizeSystem);

    /**
     * @brief This method creates a slave database needed to assemble the system
     * @param SizeSystem The size of the system
     * @param ConectivityDatabase The database that will be used to assemble the system
     * @param InverseConectivityDatabase The inverse database that will be used to assemble the system
     */
        
    void CreateSlaveConectivityDatabase(
        std::size_t& SizeSystem,
        IntMap& ConectivityDatabase,
        IntMap& InverseConectivityDatabase
        );
    
    /**
     * @brief This method returns the corresponding integration order considered
     * @return The integration order considered
     */
    IntegrationMethod GetIntegrationMethod();
    
    /**
     * @brief This method checks if all components of a vector are true
     * @param VectorToCheck The vector to check
     * @return result True if all componets are true
     */
    bool CheckWholeVector(std::vector<bool> VectorToCheck);
    
    /**
     * @brief This method computes the residual matrix of the mapping
     * @param ResidualMatrix The matrix containing the residual of the mappping
     * @param SlaveGeometry The slave geometry
     * @param MasterGeometry The master geometry
     * @param ThisMortarOperators The mortar operators
     */
    void ComputeResidualMatrix(       
        Matrix& ResidualMatrix,
        GeometryType& SlaveGeometry,
        GeometryType& MasterGeometry,
        const MortarOperator<TNumNodes>& ThisMortarOperators
        );
    
    /**
     * @brief This method assembles the LHS and the RHS
     * @param A The LHS of the system
     * @param b The RHS of the system
     * @param VariableSize The size of the variable
     * @param ResidualMatrix The matrix containing the residual of the mappping
     * @param SlaveGeometry The slave geometry
     * @param InverseConectivityDatabase The inverse database that will be used to assemble the system
     * @param ThisMortarOperators The mortar operators
     */
    void AssembleRHSAndLHS(
        MatrixType& A,
        std::vector<VectorType>& b,
        const unsigned int& VariableSize,
        const Matrix& ResidualMatrix,
        GeometryType& SlaveGeometry,
        IntMap& InverseConectivityDatabase,
        const MortarOperator<TNumNodes>& ThisMortarOperators
        );
    
    /**
     * @brief This method assembles the RHS
     * @param b The RHS of the system
     * @param VariableSize The size of the variable
     * @param ResidualMatrix The matrix containing the residual of the mappping
     * @param SlaveGeometry The slave geometry
     * @param InverseConectivityDatabase The inverse database that will be used to assemble the system
     */
    void AssembleRHS(
        std::vector<VectorType>& b,
        const unsigned int& VariableSize,
        const Matrix& ResidualMatrix,
        GeometryType& SlaveGeometry,
        IntMap& InverseConectivityDatabase
        );
    
    /**
     * @brief This method executes the explicit mapping (when no linear solver is avalaible)
     */
    void ExecuteExplicitMapping();
    
    /**
     * @brief This method executes the mapping when a linear solver is avalaible and a system of equations can be solved
     */
    void ExecuteImplicitMapping();
    
    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    Parameters GetDefaultParameters();
        
    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    
    /// Assignment operator.
    SimpleMortarMapperProcess& operator=(SimpleMortarMapperProcess const& rOther) = delete;

    /// Copy constructor.
    //SimpleMortarMapperProcess(SimpleMortarMapperProcess const& rOther);
    
    ///@}
};// class SimpleMortarMapperProcess

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
#endif /* KRATOS_SIMPLE_MORTAR_MAPPER_PROCESS defined */
