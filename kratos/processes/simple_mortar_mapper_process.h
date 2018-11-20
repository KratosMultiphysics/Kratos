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

    /// The definition of the size type
    typedef std::size_t SizeType;

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
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it 
        this->Coordinates() = mpOriginCond->GetGeometry().Center().Coordinates();
#else
        noalias(this->Coordinates()) = mpOriginCond->GetGeometry().Center().Coordinates();
#endif // ifdef KRATOS_USE_AMATRIX
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
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template< const SizeType TDim, const SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster = TNumNodes>
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
    typedef std::size_t                                          IndexType;

    /// A map for integers
    typedef std::unordered_map<IndexType, IndexType>                IntMap;

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

    /// Mortar definition
    typedef MortarKinematicVariables<TNumNodes, TNumNodesMaster>                        MortarKinematicVariablesType;
    typedef MortarOperator<TNumNodes, TNumNodesMaster>                                            MortarOperatorType;
    typedef DualLagrangeMultiplierOperators<TNumNodes, TNumNodesMaster>          DualLagrangeMultiplierOperatorsType;
    typedef ExactMortarIntegrationUtility<TDim, TNumNodes, false, TNumNodesMaster> ExactMortarIntegrationUtilityType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rOriginModelPart The origin model part to compute
     * @param rDestinationModelPart The destination model part to compute
     * @param ThisParameters The configuration parameters
     * @param pThisLinearSolver The pointer to the linear to be used (in case of implicit resolution)
     */
    SimpleMortarMapperProcess(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        Parameters ThisParameters = Parameters(R"({})" ),
        LinearSolverType::Pointer pThisLinearSolver = nullptr
        );

    /**
     * @brief Default constructor
     * @param rOriginModelPart The origin model part to compute
     * @param rDestinationModelPart The destination model part to compute
     * @param rThisVariable The variable to transfer and be transfered
     * @param ThisParameters The configuration parameters
     * @param pThisLinearSolver The pointer to the linear to be used (in case of implicit resolution)
     */
    SimpleMortarMapperProcess(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        TVarType& rThisVariable,
        Parameters ThisParameters = Parameters(R"({})" ),
        LinearSolverType::Pointer pThisLinearSolver = nullptr
        );

    /**
     * @brief A constructor where two different variables can be considered for each subdomain
     * @param rOriginModelPart The origin model part to compute
     * @param rDestinationModelPart The destination model part to compute
     * @param rOriginVariable The variable to transfer
     * @param rDestinationVariable The variable to be transfered
     * @param ThisParameters The configuration parameters
     * @param pThisLinearSolver The pointer to the linear to be used (in case of implicit resolution)
     */
    SimpleMortarMapperProcess(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        TVarType& rOriginVariable,
        TVarType& rDestinationVariable,
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

    double mMappingCoefficient = 1.0;             /// The mapping coefficient

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
     * @param rConditionsPointSlave The list of points that form the triangle decomposition
     * @param rSlaveGeometry The slave geometry
     * @param rMasterGeometry The master geometry
     * @param rMasterNormal The normal vector of the master geometry
     * @param rThisKinematicVariables The kinematic variables of the geometries, needed to integrate the mortar operators
     * @param rThisMortarOperators The mortar operators
     * @param rThisIntegrationMethod The integration method used, determines the integration order
     * @param Ae The dual lagrange multiplier operator
     */
    void AssemblyMortarOperators(
        const std::vector<array_1d<PointType,TDim>>& rConditionsPointSlave,
        GeometryType& rSlaveGeometry,
        GeometryType& rMasterGeometry,
        const array_1d<double, 3>& rMasterNormal,
        MortarKinematicVariablesType& rThisKinematicVariables,
        MortarOperatorType& rThisMortarOperators,
        const IntegrationMethod& rThisIntegrationMethod,
        const BoundedMatrixType Ae = IdentityMatrix(TNumNodes)
        );

    /**
     * @brief This method computes the Ae matrix
     * @param rSlaveGeometry The slave geometry
     * @param rThisKinematicVariables The kinematic variables
     * @param rConditionsPointsSlave The list of decomposed triangles
     * @param rThisIntegrationMethod The integration method considered
     * @return Ae: The matrix of dual LM
     */
    static inline BoundedMatrixType CalculateAe(
        GeometryType& rSlaveGeometry,
        MortarKinematicVariablesType& rThisKinematicVariables,
        std::vector<array_1d<PointType,TDim>>& rConditionsPointsSlave,
        const IntegrationMethod& rThisIntegrationMethod
        );

    /**
     * @brief This method inverts a diagonal matrix
     * @param rInputMatrix The matrix to invert
     * @return The matrix inverted
     */
    static inline BoundedMatrixType InvertDiagonalMatrix(const BoundedMatrixType& rInputMatrix);

    /**
     * @brief This method inverts a diagonal matrix
     * @param rInputMatrix The matrix to invert
     * @param rInvertedMatrix The matrix inverted
     */
    static inline void InvertDiagonalMatrix(
        const BoundedMatrixType& rInputMatrix,
        BoundedMatrixType& rInvertedMatrix
        );

    /**
     * @brief This method lumps a matrix
     * @param rInputMatrix The matrix to lump
     */
    void LumpMatrix(BoundedMatrixType& rInputMatrix);

    /**
     * @brief This method computes the size of the system
     * @param rSizeSystem The size of the system
     */
    void GetSystemSize(SizeType& rSizeSystem);

    /**
     * @brief This method creates a slave database needed to assemble the system
     * @param rSizeSystem The size of the system
     * @param rConectivityDatabase The database that will be used to assemble the system
     * @param rInverseConectivityDatabase The inverse database that will be used to assemble the system
     */
    void CreateSlaveConectivityDatabase(
        SizeType& rSizeSystem,
        IntMap& rConectivityDatabase,
        IntMap& rInverseConectivityDatabase
        );

    /**
     * @brief This method returns the corresponding integration order considered
     * @return The integration order considered
     */
    IntegrationMethod GetIntegrationMethod();

    /**
     * @brief This method checks if all components of a vector are true
     * @param rVectorToCheck The vector to check
     * @return result True if all componets are true
     */
    bool CheckWholeVector(std::vector<bool>& rVectorToCheck);

    /**
     * @brief This method computes the residual matrix of the mapping
     * @param rResidualMatrix The matrix containing the residual of the mappping
     * @param rSlaveGeometry The slave geometry
     * @param rMasterGeometry The master geometry
     * @param rThisMortarOperators The mortar operators
     */
    void ComputeResidualMatrix(
        Matrix& rResidualMatrix,
        GeometryType& rSlaveGeometry,
        GeometryType& rMasterGeometry,
        const MortarOperatorType& rThisMortarOperators
        );

    /**
     * @brief This method assembles the LHS and the RHS
     * @param rA The LHS of the system
     * @param rb The RHS of the system
     * @param VariableSize The size of the variable
     * @param rResidualMatrix The matrix containing the residual of the mappping
     * @param rSlaveGeometry The slave geometry
     * @param rInverseConectivityDatabase The inverse database that will be used to assemble the system
     * @param rThisMortarOperators The mortar operators
     */
    void AssembleRHSAndLHS(
        MatrixType& rA,
        std::vector<VectorType>& rb,
        const SizeType VariableSize,
        const Matrix& rResidualMatrix,
        GeometryType& rSlaveGeometry,
        IntMap& rInverseConectivityDatabase,
        const MortarOperatorType& rThisMortarOperators
        );

    /**
     * @brief This method assembles the RHS
     * @param rb The RHS of the system
     * @param VariableSize The size of the variable
     * @param rResidualMatrix The matrix containing the residual of the mappping
     * @param rSlaveGeometry The slave geometry
     * @param rInverseConectivityDatabase The inverse database that will be used to assemble the system
     */
    void AssembleRHS(
        std::vector<VectorType>& rb,
        const SizeType VariableSize,
        const Matrix& rResidualMatrix,
        GeometryType& rSlaveGeometry,
        IntMap& rInverseConectivityDatabase
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
     * @param rA The LHS of the system
     * @param rb The RHS of the system
     * @param rInverseConectivityDatabase The inverse database that will be used to assemble the system
     * @param pIndexesPairs The pointer to indexed objects
     * @param itCond Iterator of a condition
     * @param rIntegrationUtility An integration utility for mortar
     * @param rThisKineticVariables Kinematic variables (shape functions)
     * @param rThisMortarOperators The mortar operators
     * @param Iteration The current non-linear iteration
     * @tparam TClassType The class of index pairs considered
     * @tparam TImplicit If we solve with lamping or we use a linear solver
     */
    template<class TClassType, bool TImplicit = false>
    void PerformMortarOperations(
        MatrixType& rA,
        std::vector<VectorType>& rb,
        IntMap& rInverseConectivityDatabase,
        typename TClassType::Pointer pIndexesPairs,
        ConditionsArrayType::iterator itCond,
        ExactMortarIntegrationUtilityType& rIntegrationUtility,
        MortarKinematicVariablesType& rThisKineticVariables,
        MortarOperatorType& rThisMortarOperators,
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
            const bool is_inside = rIntegrationUtility.GetExactIntegration(slave_geometry, slave_normal, master_geometry, master_normal, conditions_points_slave);

            if (is_inside) {
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

                if (!TImplicit) {
                    MortarUtilities::AddValue<TVarType, NonHistorical>(slave_geometry, aux_variable, residual_matrix);
                }

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
                        AssembleRHSAndLHS(rA, rb, variable_size, residual_matrix, slave_geometry, rInverseConectivityDatabase, rThisMortarOperators);
                    } else {
                        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                            slave_geometry[i_node].GetValue(NODAL_AREA) += rThisMortarOperators.DOperator(i_node, i_node);
                        }
                    }
                } else if (TImplicit) {
                    const SizeType variable_size = MortarUtilities::SizeToCompute<TDim, TVarType>();
                    AssembleRHS(rb, variable_size, residual_matrix, slave_geometry, rInverseConectivityDatabase);
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

/**
 * @ingroup KratosCore
 * @class SimpleMortarMapperProcessWrapper
 * @brief This class wraps automatically the different types mof mappers
 * @author Vicente Mataix Ferrandiz
 */
class SimpleMortarMapperProcessWrapper
        : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SimpleMortarMapperProcessWrapper
    KRATOS_CLASS_POINTER_DEFINITION(SimpleMortarMapperProcessWrapper);

    /// Linear solver
    typedef UblasSpace<double, CompressedMatrix, Vector>    SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector>               LocalSpaceType;
    typedef typename SparseSpaceType::MatrixType                 MatrixType;
    typedef typename SparseSpaceType::VectorType                 VectorType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    /// Index type definition
    typedef std::size_t                                          IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rOriginModelPart The origin model part to compute
     * @param rDestinationModelPart The destination model part to compute
     * @param rThisVariable The variable to transfer and be transfered
     * @param ThisParameters The configuration parameters
     * @param pThisLinearSolver The pointer to the linear to be used (in case of implicit resolution)
     */
    SimpleMortarMapperProcessWrapper(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        Parameters ThisParameters = Parameters(R"({})" ),
        LinearSolverType::Pointer pThisLinearSolver = nullptr
        )
    {
        // The default parameters
        Parameters default_parameters = GetDefaultParameters();
        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // The condition iterators
        auto it_cond_origin_begin = rOriginModelPart.Conditions().begin();
        auto it_cond_destination_begin = rOriginModelPart.Conditions().begin();

        // The dimensions
        const SizeType dimension = it_cond_origin_begin->GetGeometry().WorkingSpaceDimension();
        const SizeType size_1 = it_cond_origin_begin->GetGeometry().size();
        const SizeType size_2 = it_cond_destination_begin->GetGeometry().size();

        // The variable names
        const std::string& origin_variable_name = ThisParameters["origin_variable"].GetString();
        const std::string& destination_variable_name = ThisParameters["destination_variable"].GetString();

        bool double_variable = true;
        if(KratosComponents<Variable<double>>::Has(origin_variable_name)) {
            if (!(KratosComponents<Variable<double>>::Has(destination_variable_name)))
                KRATOS_ERROR << "The destination variable is not the same type (double) as the origin" << std::endl;
        } else if (KratosComponents< Variable< array_1d< double, 3> > >::Has(origin_variable_name)) {
            double_variable = false;
            if (!(KratosComponents<Variable<array_1d< double, 3>>>::Has(destination_variable_name)))
                KRATOS_ERROR << "The destination variable is not the same type (array_1d< double, 3>) as the origin" << std::endl;
        } else {
            KRATOS_ERROR << "The types of the variables are not supported array_1d< double, 3> or double" << std::endl;
        }

        // Creating the mapper
        if (dimension == 2) {
            // 2D
            if (double_variable) {
                mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<2, 2, Variable<double>>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
            } else {
                mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<2, 2, Variable<array_1d< double, 3>>>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
            }
        } else {
            // 3D
            if (size_1 == 3 && size_2 == 3) {
                if (double_variable) {
                    mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<3, 3, Variable<double>>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
                } else {
                    mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<3, 3, Variable<array_1d< double, 3>>>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
                }
            } else if (size_1 == 4 && size_2 == 4) {
                if (double_variable) {
                    mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<3, 4, Variable<double>>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
                } else {
                    mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<3, 4, Variable<array_1d< double, 3>>>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
                }
            } else if (size_1 == 3 && size_2 == 4) {
                if (double_variable) {
                    mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<3, 3, Variable<double>, 4>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
                } else {
                    mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<3, 3, Variable<array_1d< double, 3>>, 4>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
                }
            } else if (size_1 == 4 && size_2 == 3) {
                if (double_variable) {
                    mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<3, 4, Variable<double>, 3>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
                } else {
                    mpMapperProcess = Kratos::make_shared<SimpleMortarMapperProcess<3, 4, Variable<array_1d< double, 3>>, 3>>(rOriginModelPart, rDestinationModelPart, ThisParameters, pThisLinearSolver);
                }
            }
        }
    }

    /// Destructor.
    ~SimpleMortarMapperProcessWrapper() override = default;

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

    void Execute() override
    {
        mpMapperProcess->Execute();
    }

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
        return "SimpleMortarMapperProcessWrapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SimpleMortarMapperProcessWrapper";
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
    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    Parameters GetDefaultParameters()
    {
        Parameters default_parameters = Parameters(R"(
        {
            "echo_level"                       : 0,
            "absolute_convergence_tolerance"   : 1.0e-9,
            "relative_convergence_tolerance"   : 1.0e-4,
            "max_number_iterations"            : 10,
            "integration_order"                : 2,
            "distance_threshold"               : 1.0e24,
            "origin_variable"                  : "TEMPERATURE",
            "destination_variable"             : "",
            "origin_variable_historical"       : true,
            "destination_variable_historical"  : true,
            "mapping_coefficient "             : 1.0,
            "search_parameters"                : {
                "allocation_size"                  : 1000,
                "bucket_size"                      : 4,
                "search_factor"                    : 3.5
            }
        })" );

        return default_parameters;
    }

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

    Process::Pointer mpMapperProcess = nullptr; /// The real mapper process

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

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{

    ///@}
};// class SimpleMortarMapperProcessWrapper


///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
#endif /* KRATOS_SIMPLE_MORTAR_MAPPER_PROCESS defined */
