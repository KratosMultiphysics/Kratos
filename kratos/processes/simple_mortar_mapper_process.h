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
#include <unordered_map>

// External includes

// Project includes
#include "processes/process.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "includes/mortar_classes.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/atomic_utilities.h"
#include "utilities/exact_mortar_segmentation_utility.h"

/* Tree structures */
#include "spatial_containers/point_object.h"
#include "spatial_containers/spatial_containers.h" // kd-tree

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The definition of the size type
    using SizeType = std::size_t;

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
 * @class SimpleMortarMapperProcess
 * @brief This is basic mapper of values between domains using mortar formulation
 * @details Using the dual mortar formulation the resolution of the system of equations is not needed.
 * Several types of constructors are avaible depending of the needs.
 * If the pairs sets are not provided a search will be performed using a KDTree
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

    // DEFINITION OF FLAGS TO CONTROL THE BEHAVIOUR
    KRATOS_DEFINE_LOCAL_FLAG(AVERAGE_NORMAL);                      /// If using average normal
    KRATOS_DEFINE_LOCAL_FLAG(DISCONTINOUS_INTERFACE);              /// If interface is discontinous
    KRATOS_DEFINE_LOCAL_FLAG(ORIGIN_IS_HISTORICAL);                /// If the origin variables is historical
    KRATOS_DEFINE_LOCAL_FLAG(DESTINATION_IS_HISTORICAL);           /// If the destination variables is historical
    KRATOS_DEFINE_LOCAL_FLAG(ORIGIN_SKIN_IS_CONDITION_BASED);      /// If the entities to take into account on the origin model part are the conditions, otherwise we will take elements into consideration
    KRATOS_DEFINE_LOCAL_FLAG(DESTINATION_SKIN_IS_CONDITION_BASED); /// If the entities to take into account on the destination model part are the conditions, otherwise we will take elements into consideration
    KRATOS_DEFINE_LOCAL_FLAG(CONSIDER_TESELLATION);                /// If we consider the tesellation in the mortar integration

    /// Pointer definition of SimpleMortarMapperProcess
    KRATOS_CLASS_POINTER_DEFINITION(SimpleMortarMapperProcess);

    //// Type definition for Point
    using PointType = Point;

    /// Type definition for node geometry
    using GeometryType = Geometry<Node>;

    /// Type definition for point geometry
    using GeometryPointType = Geometry<PointType>;

    //// Type definition for integration methods
    using IntegrationMethod = GeometryData::IntegrationMethod;

    /// Type definition for LineType
    using LineType = Line2D2<PointType>;

    /// Type definition for TriangleType
    using TriangleType = Triangle3D3<PointType>;

    /// Type definition for DecompositionType based on TDim value
    using DecompositionType = typename std::conditional<TDim == 2, LineType, TriangleType>::type;

    /// Type definition for sparse space type
    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;

    /// Type definition for local space type
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

    /// Type definition for matrix
    using MatrixType = typename SparseSpaceType::MatrixType;

    /// Type definition for vector
    using VectorType = typename SparseSpaceType::VectorType;

    /// Type definition for linear solver
    using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;

    /// Index type definition
    using IndexType = std::size_t;

    /// A map for integers
    using IntMap = std::unordered_map<IndexType, IndexType>;

    /// BoundedMatrix
    using BoundedMatrixType = BoundedMatrix<double, TNumNodes, TNumNodes>;

    /// Type definitions for the tree
    using PointMapperType = PointObject<GeometricalObject>;
    using PointTypePointer = typename PointMapperType::Pointer;
    using PointVector = std::vector<PointTypePointer>;

    // KDtree definitions
    using BucketType = Bucket<3ul, PointMapperType, PointVector>;
    using KDTreeType = Tree<KDTreePartition<BucketType>>;

    /// Mortar definition
    using MortarKinematicVariablesType = MortarKinematicVariables<TNumNodes, TNumNodesMaster>;
    using MortarOperatorType = MortarOperator<TNumNodes, TNumNodesMaster>;
    using DualLagrangeMultiplierOperatorsType = DualLagrangeMultiplierOperators<TNumNodes, TNumNodesMaster>;
    using ExactMortarIntegrationUtilityType = ExactMortarIntegrationUtility<TDim, TNumNodes, false, TNumNodesMaster>;

    /// Auxiliar struct for mapping
    struct TLS {
        MortarKinematicVariablesType this_kinematic_variables;    // Create and initialize condition variables:
        MortarOperatorType this_mortar_operators;                 // Create the mortar operators
        ExactMortarIntegrationUtilityType integration_utility;    // We call the exact integration utility
    };

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

    /**
     * @brief Execute method is used to execute the Process algorithms.
     */
    void Execute() override;

    /**
     * @details This function will be executed at every time step BEFORE performing the solve phase
     */
    void ExecuteInitializeSolutionStep() override;

    /**
     * @brief This method is a direct map between the origin and destination modelpart with custom variables
     * @param rOriginVariable The origin model part
     * @param rDestinationVariable The destination model part
     * @param Flag The flags to special settings. Right now does nothing
     */
    void Map(
        TVarType& rOriginVariable,
        TVarType& rDestinationVariable,
        const Flags Flag = Flags()
        )
    {
        // Reassign the variables
        mpOriginVariable = &rOriginVariable;
        mpDestinationVariable = &rDestinationVariable;

        // Execute the process
        Execute();
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

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
private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mOriginModelPart;                  /// The origin model part to compute
    ModelPart& mDestinationModelPart;             /// The destination model part to compute
    const TVarType* mpOriginVariable;             /// The origin variable to map
    const TVarType* mpDestinationVariable;        /// The destiny variable to map

    double mMappingCoefficient = 1.0;             /// The mapping coefficient
    Flags mOptions;                               /// Local flags

    unsigned int mEchoLevel = 0;                  /// The verbosity level
    Parameters mThisParameters;                   /// The configuration parameters

    LinearSolverType::Pointer mpThisLinearSolver; /// The linear solver used to compute the solution

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
     * @param rGeometricalObjectsPointSlave The list of points that form the triangle decomposition
     * @param rSlaveGeometry The slave geometry
     * @param rMasterGeometry The master geometry
     * @param rMasterNormal The normal vector of the master geometry
     * @param rThisKinematicVariables The kinematic variables of the geometries, needed to integrate the mortar operators
     * @param rThisMortarOperators The mortar operators
     * @param rThisIntegrationMethod The integration method used, determines the integration order
     * @param Ae The dual lagrange multiplier operator
     */
    void AssemblyMortarOperators(
        const std::vector<array_1d<PointType,TDim>>& rGeometricalObjectsPointSlave,
        const GeometryType& rSlaveGeometry,
        const GeometryType& rMasterGeometry,
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
        const GeometryType& rSlaveGeometry,
        const GeometryType& rMasterGeometry,
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
        const GeometryType& rSlaveGeometry,
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
        const GeometryType& rSlaveGeometry,
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
     * @param rGeometricalObject Reference of a geometrical object
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
        GeometricalObject& rGeometricalObject,
        ExactMortarIntegrationUtilityType& rIntegrationUtility,
        MortarKinematicVariablesType& rThisKineticVariables,
        MortarOperatorType& rThisMortarOperators,
        const IndexType Iteration
        )
    {
        // The root model part
        ModelPart& r_root_model_part = mOriginModelPart.GetRootModelPart();

        // Getting the auxiliar variable
        const TVarType& r_aux_variable = KratosComponents<TVarType>::Get(MortarUtilities::GetAuxiliarVariable<TVarType>());

        // Indexes of the pair to be removed
        std::vector<IndexType> indexes_to_remove, geometrical_objects_to_erase;

        // Getting discontinous factor
        const double discontinous_interface_factor = mOptions.Is(DISCONTINOUS_INTERFACE) ? mThisParameters["discontinous_interface_factor"].GetDouble() : 1.0;

        // Declare auxiliar coordinates
        GeometryType::CoordinatesArrayType aux_coords;

        // Geometrical values
        auto& r_slave_geometry = rGeometricalObject.GetGeometry();
        r_slave_geometry.PointLocalCoordinates(aux_coords, r_slave_geometry.Center());
        const array_1d<double, 3> slave_normal = r_slave_geometry.UnitNormal(aux_coords);

        // The model part as const to avoid race conditions
        const auto& r_const_origin_model_part = mOriginModelPart;
        const auto& r_const_destination_model_part = mDestinationModelPart;

        for (auto it_pair = pIndexesPairs->begin(); it_pair != pIndexesPairs->end(); ++it_pair ) {
            const IndexType master_id = pIndexesPairs->GetId(it_pair); // MASTER

            const auto& r_master_geometry = mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED) ? r_const_origin_model_part.GetCondition(master_id).GetGeometry() : r_const_origin_model_part.GetElement(master_id).GetGeometry();
            r_master_geometry.PointLocalCoordinates(aux_coords, r_master_geometry.Center());
            const array_1d<double, 3> master_normal = r_master_geometry.UnitNormal(aux_coords);

            const IntegrationMethod& r_integration_method = GetIntegrationMethod();

            // Reading integration points
            std::vector<array_1d<PointType,TDim>> geometrical_objects_points_slave; // These are the segmentation points, with this points it is possible to create the lines or triangles used on the mapping
            const bool is_inside = rIntegrationUtility.GetExactIntegration(r_slave_geometry, slave_normal, r_master_geometry, master_normal, geometrical_objects_points_slave);

            if (is_inside) {
                // Initialize general variables for the current master element
                rThisKineticVariables.Initialize();

                // Initialize the mortar operators
                rThisMortarOperators.Initialize();

                const BoundedMatrixType Ae = CalculateAe(r_slave_geometry, rThisKineticVariables, geometrical_objects_points_slave, r_integration_method);

                AssemblyMortarOperators( geometrical_objects_points_slave, r_slave_geometry, r_master_geometry,master_normal, rThisKineticVariables, rThisMortarOperators, r_integration_method, Ae);

                /* We compute the residual */
                const IndexType size_to_compute = MortarUtilities::SizeToCompute<TDim, TVarType>();
                Matrix residual_matrix(TNumNodes, size_to_compute);
                ComputeResidualMatrix(residual_matrix, r_slave_geometry, r_master_geometry, rThisMortarOperators);

                if (!TImplicit) {
                    MortarUtilities::AddValue<TVarType, MortarUtilitiesSettings::SaveAsNonHistoricalVariable>(r_slave_geometry, r_aux_variable, residual_matrix);
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
                    if constexpr (TImplicit) {
                        /* We compute the residual and assemble */
                        const SizeType variable_size = MortarUtilities::SizeToCompute<TDim, TVarType>();
                        AssembleRHSAndLHS(rA, rb, variable_size, residual_matrix, r_slave_geometry, rInverseConectivityDatabase, rThisMortarOperators);
                    } else {
                        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                            double& r_nodal_area = r_slave_geometry[i_node].GetValue(NODAL_AREA);
                            AtomicAdd(r_nodal_area, rThisMortarOperators.DOperator(i_node, i_node));
                        }
                        // In case of discontinous interface we add contribution to near nodes
                        if (mOptions.Is(DISCONTINOUS_INTERFACE)) {
                            const double element_length = r_slave_geometry.Length();

                            // Iterating over nodes
                            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                                const double nodal_area_contribution = rThisMortarOperators.DOperator(i_node, i_node);

                                // The original node coordinates
                                const auto& r_slave_node = r_slave_geometry[i_node];

                                // Iterating over other paired geometrical objects
                                const auto& r_index_masp_master = mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED) ? r_const_origin_model_part.GetCondition(master_id).GetValue(INDEX_SET) : r_const_origin_model_part.GetElement(master_id).GetValue(INDEX_SET);
                                for (auto it_master_pair = r_index_masp_master->begin(); it_master_pair != r_index_masp_master->end(); ++it_master_pair ) {

                                    const IndexType auxiliar_slave_id = r_index_masp_master->GetId(it_master_pair);
                                    if (rGeometricalObject.Id() != auxiliar_slave_id) {
                                        GeometryType& r_auxiliar_slave_geometry =  const_cast<GeometryType&>(mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED) ? r_const_destination_model_part.GetCondition(auxiliar_slave_id).GetGeometry() : r_const_destination_model_part.GetElement(auxiliar_slave_id).GetGeometry());

                                        for (IndexType j_node = 0; j_node < TNumNodes; ++j_node) {
                                            // The auxiliary node distance
                                            auto& r_auxiliary_slave_node = r_auxiliar_slave_geometry[j_node];
                                            const double distance = r_auxiliary_slave_node.Distance(r_slave_node);
                                            const double contribution_coeff = 1.0/std::pow((1.0 + distance/(discontinous_interface_factor * element_length)), 2);

                                            double& r_nodal_area = r_auxiliary_slave_node.GetValue(NODAL_AREA);
                                            AtomicAdd(r_nodal_area, contribution_coeff * nodal_area_contribution);
                                        }
                                    }
                                }
                            }
                        }
                    }
                } else if constexpr (TImplicit) {
                    const SizeType variable_size = MortarUtilities::SizeToCompute<TDim, TVarType>();
                    AssembleRHS(rb, variable_size, residual_matrix, r_slave_geometry, rInverseConectivityDatabase);
                }
            } else { // NOTE: The geometrical object considered maybe is to tight
                indexes_to_remove.push_back(master_id);
                const IndexType other_id = pIndexesPairs->GetOtherId(it_pair);
                if (std::is_same<TClassType, IndexMap>::value && other_id != 0) {
                    geometrical_objects_to_erase.push_back(other_id);
                }
            }
        }

        // Clear indexes
        for (IndexType i_to_remove = 0; i_to_remove < indexes_to_remove.size(); ++i_to_remove) {
            if (mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED)) {
                for (auto& id : geometrical_objects_to_erase ) {
                    auto& r_cond = r_root_model_part.GetCondition(id);
                    r_cond.Set(TO_ERASE, true);
                }
            } else {
                for (auto& id : geometrical_objects_to_erase ) {
                    auto& r_elem = r_root_model_part.GetElement(id);
                    r_elem.Set(TO_ERASE, true);
                }
            }
            pIndexesPairs->RemoveId(indexes_to_remove[i_to_remove]);
        }
    }

    /**
     * @brief This method can be used to clear the unused indexes
     * @param pIndexesPairs The pointer to indexed objects
     * @param rGeometricalObject Reference of a geometrical object
     * @param rIntegrationUtility An integration utility for mortar
     * @tparam TClassType The class of index pairs considered
     */
    template<class TClassType>
    void ClearIndexes(
        typename TClassType::Pointer pIndexesPairs,
        GeometricalObject& rGeometricalObject,
        ExactMortarIntegrationUtilityType& rIntegrationUtility
        )
    {
        // The root model part
        ModelPart& r_root_model_part = mOriginModelPart.GetRootModelPart();

        // Indexes of the pair to be removed
        std::vector<IndexType> indexes_to_remove, geometrical_objects_to_erase;

        // Declare auxiliar coordinates
        GeometryType::CoordinatesArrayType aux_coords;

        // Geometrical values
        auto& r_slave_geometry = rGeometricalObject.GetGeometry();
        r_slave_geometry.PointLocalCoordinates(aux_coords, r_slave_geometry.Center());
        const array_1d<double, 3> slave_normal = r_slave_geometry.UnitNormal(aux_coords);

        // The model part as const to avoid race conditions
        const auto& r_const_origin_model_part = mOriginModelPart;

        for (auto it_pair = pIndexesPairs->begin(); it_pair != pIndexesPairs->end(); ++it_pair ) {
            const IndexType master_id = pIndexesPairs->GetId(it_pair); // MASTER

            const auto& r_master_geometry = mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED) ? r_const_origin_model_part.GetCondition(master_id).GetGeometry() : r_const_origin_model_part.GetElement(master_id).GetGeometry();
            r_master_geometry.PointLocalCoordinates(aux_coords, r_master_geometry.Center());
            const array_1d<double, 3> master_normal = r_master_geometry.UnitNormal(aux_coords);

            // Reading integration points
            std::vector<array_1d<PointType,TDim>> geometrical_objects_points_slave; // These are the segmentation points, with this points it is possible to create the lines or triangles used on the mapping
            const bool is_inside = rIntegrationUtility.GetExactIntegration(r_slave_geometry, slave_normal, r_master_geometry, master_normal, geometrical_objects_points_slave);

            if (!is_inside) {
                indexes_to_remove.push_back(master_id);
                const IndexType other_id = pIndexesPairs->GetOtherId(it_pair);
                if (std::is_same<TClassType, IndexMap>::value && other_id != 0) {
                    geometrical_objects_to_erase.push_back(other_id);
                }
            }
        }

        // Clear indexes
        for (IndexType i_to_remove = 0; i_to_remove < indexes_to_remove.size(); ++i_to_remove) {
            if (mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED)) {
                for (auto& id : geometrical_objects_to_erase ) {
                    auto& r_cond = r_root_model_part.GetCondition(id);
                    r_cond.Set(TO_ERASE, true);
                }
            } else {
                for (auto& id : geometrical_objects_to_erase ) {
                    auto& r_elem = r_root_model_part.GetElement(id);
                    r_elem.Set(TO_ERASE, true);
                }
            }
            pIndexesPairs->RemoveId(indexes_to_remove[i_to_remove]);
        }
    }

    /**
     * @brief This method fills the database
     * @param rGeometricalObject Reference of a geometrical object
     * @param rTreePoints The search tree
     * @param AllocationSize The allocation size of the tree
     * @param SearchFactor The search factor of the tree
     */
    template<class TEntity>
    void FillDatabase(
        TEntity& rGeometricalObject,
        KDTreeType& rTreePoints,
        const SizeType AllocationSize,
        const double SearchFactor
        )
    {
        // Initialize values
        PointVector points_found(AllocationSize);

        GeometryType& r_geometry = rGeometricalObject.GetGeometry();
        const Point center = r_geometry.Center();

        double radius = 0.0;
        for(IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)  {
            const array_1d<double, 3> aux_vector = center.Coordinates() - r_geometry[i_node].Coordinates();
            const double aux_value = inner_prod(aux_vector, aux_vector);
            if(aux_value > radius) radius = aux_value;
        }

        const double search_radius = SearchFactor * std::sqrt(radius);

        const SizeType number_points_found = rTreePoints.SearchInRadius(center, search_radius, points_found.begin(), AllocationSize);

        if (number_points_found > 0) {
            // In case of missing is created
            if (!rGeometricalObject.Has(INDEX_SET))
                rGeometricalObject.SetValue(INDEX_SET, Kratos::make_shared<IndexSet>());

            // Accessing to the index set
            IndexSet::Pointer indexes_set = rGeometricalObject.GetValue(INDEX_SET);

            for (IndexType i_point = 0; i_point < number_points_found; ++i_point ) {
                auto p_geometrical_object_master = points_found[i_point]->pGetObject();
                indexes_set->AddId(p_geometrical_object_master->Id());
            }
        }
    }

    /**
     * @brief This method creates an inverse database
     */
    void CreateInverseDatabase();

    /**
     * @brief Reset the interface database
     * @details This method resets the mapping database saved in the destination database.
     * @note Note that this needs to be done if such modelpart has changed its number of nodes or geometrical objects. This needs to be done even though the mapping instance is deleted since such information is saved in the destination nodes and geometrical objects.
     */
    void UpdateInterface();

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
