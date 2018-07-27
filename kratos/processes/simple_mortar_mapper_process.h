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
template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
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
    
    /// Index type definition
    typedef std::size_t                                           IndexType;

    /// Size type definition
    typedef std::size_t                                            SizeType;

    /// A map for integers
    typedef std::unordered_map<IndexType, IndexType>                 IntMap;
    
    /// BoundedMatrix
    typedef BoundedMatrix<double, TNumNodes, TNumNodes>  BoundedMatrixType;

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
    
    bool mOriginHistorical;                       /// A bool that defines if the origin variables is historical
    bool mDestinationHistorical;                  /// A bool that defines if the destination variables is historical

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
        const SizeType& VariableSize,
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
        const SizeType& VariableSize,
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
     * @brief This method computes common methods between the implicit and explicit formulation
     */
    template<class TClassType, bool TImplicit = false>
    void PerformMortarOperations(
        MatrixType& A,
        std::vector<VectorType>& b,
        IntMap& InverseConectivityDatabase,
        typename TClassType::Pointer pIndexesPairs,
        ConditionsArrayType::iterator itCond,
        ExactMortarIntegrationUtility<TDim, TNumNodes>& IntegrationUtility,
        MortarKinematicVariables<TNumNodes>& rThisKineticVariables,
        MortarOperator<TNumNodes>& rThisMortarOperators,
        const IndexType Iteration
        )
    {
        // The root model part
        ModelPart& root_model_part = mOriginModelPart.GetRootModelPart();

        // Getting the auxiliar variable
        TVarType aux_variable = MortarUtilities::GetAuxiliarVariable<TVarType>();

        // Indexes of the pair to be removed
        std::vector<IndexType> indexes_to_remove, conditions_to_erase;

        // Geometrical values
        const array_1d<double, 3>& slave_normal = itCond->GetValue(NORMAL);
        GeometryType& slave_geometry = itCond->GetGeometry();

        for (auto it_pair = pIndexesPairs->begin(); it_pair != pIndexesPairs->end(); ++it_pair ) {
            const IndexType master_id = pIndexesPairs->GetId(it_pair);
            Condition::Pointer p_cond_master = mOriginModelPart.pGetCondition(master_id); // MASTER
            const array_1d<double, 3>& master_normal = p_cond_master->GetValue(NORMAL);
            GeometryType& master_geometry = p_cond_master->GetGeometry();

            const IntegrationMethod& this_integration_method = GetIntegrationMethod();

            // Reading integration points
            std::vector<array_1d<PointType,TDim>> conditions_points_slave; // These are the segmentation points, with this points it is possible to create the lines or triangles used on the mapping
            const bool is_inside = IntegrationUtility.GetExactIntegration(slave_geometry, slave_normal, master_geometry, master_normal, conditions_points_slave);

            if (is_inside == true) {
                // Initialize general variables for the current master element
                rThisKineticVariables.Initialize();

                // Initialize the mortar operators
                rThisMortarOperators.Initialize();

                const BoundedMatrixType Ae = CalculateAe(slave_geometry, rThisKineticVariables, conditions_points_slave, this_integration_method);

                AssemblyMortarOperators( conditions_points_slave, slave_geometry, master_geometry,master_normal, rThisKineticVariables, rThisMortarOperators, this_integration_method, Ae);

                /* We compute the residual */
                const IndexType size_to_compute = MortarUtilities::SizeToCompute<TDim, TVarType>();
                Matrix residual_matrix(TNumNodes, size_to_compute);
                ComputeResidualMatrix(residual_matrix, slave_geometry, master_geometry, rThisMortarOperators);

                MortarUtilities::AddValue<TVarType, NonHistorical>(slave_geometry, aux_variable, residual_matrix);

                // We check if DOperator is diagonal
                if (mEchoLevel > 1) {
                    BoundedMatrixType aux_copy_D = rThisMortarOperators.DOperator;
                    LumpMatrix(aux_copy_D);
                    const BoundedMatrixType aux_diff = aux_copy_D - rThisMortarOperators.DOperator;
                    const double norm_diff = norm_frobenius(aux_diff);
                    if (norm_diff > 1.0e-4)
                        KRATOS_WARNING("D OPERATOR") << " THE MORTAR OPERATOR D IS NOT DIAGONAL" << std::endl;
                    if (mEchoLevel == 3) {
                        KRATOS_WATCH(norm_diff);
                        KRATOS_WATCH(rThisMortarOperators.DOperator);
                    }
                }

                if (Iteration == 0) { // Just assembled the first iteration
                    if (TImplicit) {
                        /* We compute the residual and assemble */
                        const SizeType variable_size = MortarUtilities::SizeToCompute<TDim, TVarType>();
                        AssembleRHSAndLHS(A, b, variable_size, residual_matrix, slave_geometry, InverseConectivityDatabase, rThisMortarOperators);
                    } else {
                        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                            slave_geometry[i_node].GetValue(NODAL_AREA) += rThisMortarOperators.DOperator(i_node, i_node);
                        }
                    }
                } else if (TImplicit) {
                    const SizeType variable_size = MortarUtilities::SizeToCompute<TDim, TVarType>();
                    AssembleRHS(b, variable_size, residual_matrix, slave_geometry, InverseConectivityDatabase);
                }
            } else { // NOTE: The condition considered maybe is to tight
                indexes_to_remove.push_back(master_id);
                const IndexType other_id = pIndexesPairs->GetOtherId(it_pair);
                if (std::is_same<TClassType, IndexMap>::value && other_id != 0) {
                    conditions_to_erase.push_back(other_id);
                }
            }
        }

        // Clear indexes
        for (IndexType i_to_remove = 0; i_to_remove < indexes_to_remove.size(); ++i_to_remove) {
            for (auto& id : conditions_to_erase ) {
                Condition::Pointer p_cond = root_model_part.pGetCondition(id);
                p_cond->Set(TO_ERASE, true);
            }
            pIndexesPairs->RemoveId(indexes_to_remove[i_to_remove]);
        }
    }

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
