//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Ruben Zorrilla
//
//

#if !defined(KRATOS_FIXED_MESH_ALE_UTILITIES_H_INCLUDED )
#define  KRATOS_FIXED_MESH_ALE_UTILITIES_H_INCLUDED

// System includes


// External includes


// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "factories/linear_solver_factory.h"
#include "linear_solvers/linear_solver.h"
#include "processes/calculate_embedded_nodal_variable_from_skin_process.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "spaces/ublas_space.h"

// Application includes


namespace Kratos
{
///@addtogroup MeshMovingApplication
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

/**
 * @brief Utility to perform the FM-ALE algorithm operations
 * This utility implements the Fixed Mesh - Arbitrary Lagrangian Eulerian (FM-ALE)
 * algorithm operations. After setting a virtual mesh, which is a copy of the
 * problem background element mesh, it is moved in accordante to the immersed
 * object movement. The virtual mesh movement is solved using a common ALE mesh
 * solver (in this case the Laplacian mesh solver is used). Once the mesh movement,
 * including the mesh velocity, have been computed, the origin mesh historical
 * values (velocity and pressure), as well as the mesh velocity, are computed
 * as a projection from the virtual mesh. This is required to consistently
 * initialize the historical values when nodes change its topological status.
 */
class FixedMeshALEUtilities
{
public:
    ///@name Type Definitions
    ///@{

    typedef Element::NodesArrayType NodesArrayType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
    typedef LinearSolverFactory<SparseSpaceType, LocalSpaceType> LinearSolverFactoryType;
    typedef ResidualBasedIncrementalUpdateStaticScheme<SparseSpaceType, LocalSpaceType> SchemeType;
    typedef ResidualBasedLinearStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> StrategyType;
    typedef ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> BuilderAndSolverType;
    typedef CalculateEmbeddedNodalVariableFromSkinProcess<array_1d<double, 3>, SparseSpaceType, LocalSpaceType, LinearSolverType> EmbeddedNodalVariableProcessArrayType;

    typedef typename SchemeType::Pointer SchemePointerType;
    typedef typename StrategyType::Pointer StrategyPointerType;
    typedef typename BuilderAndSolverType::Pointer BuilderAndSolverPointerType;

    /// Pointer definition of FixedMeshALEUtilities
    KRATOS_CLASS_POINTER_DEFINITION(FixedMeshALEUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    FixedMeshALEUtilities(
        ModelPart &rVirtualModelPart,
        ModelPart &rStructureModelPart);

    /// Constructor with model and parameters
    FixedMeshALEUtilities(
        Model &rModel,
        Parameters &rParameters);

    /// Destructor.
    virtual ~FixedMeshALEUtilities() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Initializes the FM-ALE utility
     * This method fills the virtual model part as a copy of the origin model part.
     * In case it is used, it also creates and initializes the mesh moving strategy
     * @param rOriginModelPart model part from where the nodes and elements are copied
     */
    virtual void Initialize(ModelPart &rOriginModelPart);

    /**
     * @brief Set the Virtual Mesh Values From Origin Mesh object
     * This method sets the VELOCITY and PRESSURE historical values in the virtual mesh.
     * The values are retrieved from the origin mesh. This needs to be called in the
     * InitializeSolutionStep of the solver. Note that if sub-iteration is done (e.g.
     * FSI) this must be called before and just once.
     */
    virtual void SetVirtualMeshValuesFromOriginMesh();

    /**
     * @brief Compute the virtual mesh movement
     * This method computes the virtual mesh movement in accordance to the immersed structure
     * DISPLACEMENT values. To that purpose it sets the fixity and creates & solves the
     * mesh moving strategy.
     * @param DeltaTime time step value (required for the computation of the MESH_VELOCITY)
     */
    virtual void ComputeMeshMovement(const double DeltaTime);

    /**
     * @brief Revert the virtual mesh movement
     * This method reverts the virtual mesh movement to recover its original configuration,
     * which coincides with the background mesh one. It has to be called once the values
     * projection from the virtual mesh to the origin has been performed.
     */
    virtual void UndoMeshMovement();

    /**
     * @brief This method projects the virtual model part mesh values to the origin mesh
     * Once the FM-ALE operations have been performed, this method projects the nodal values
     * from the virtual mesh to the origin mesh. The projected variables are PRESSURE,
     * VELOCITY and MESH_VELOCITY.
     * @tparam TDim Template parameter containing the problem domain size
     * @param rOriginModelPart Reference to the model part to which the values are projected
     * @param BufferSize Buffer values that are projected
     */
    template <unsigned int TDim>
    void ProjectVirtualValues(
        ModelPart &rOriginModelPart,
        const unsigned int BufferSize = 3);

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
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;

    ///@}
    ///@name Friends
    ///@{


    ///@}
protected:
    ///@}
    ///@name Life Cycle
    ///@{


    ///@}
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart &mrVirtualModelPart;
    ModelPart &mrStructureModelPart;
    ModelPart *mpOriginModelPart = nullptr;

    Parameters mEmbeddedNodalVariableSettings;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    /**
    * This method fills the mrVirtualModelPart with the nodes and elmens of a given model part
    * It has to be performed once since after each values projection the virtual mesh configuration
    * is reverted to its original status.
    * @param rOriginModelPart model part from where the nodes and elements are copied
    */
    virtual void FillVirtualModelPart(ModelPart &rOriginModelPart);

    /**
     * @brief Create the virtual model part elements
     * This method creates the elements in the virtual model part
     * @param rOriginModelPart Origin model part to mimic the elements from
     */
    virtual void CreateVirtualModelPartElements(const ModelPart &rOriginModelPart);

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    LinearSolverType::Pointer mpLinearSolver = nullptr;
    StrategyPointerType mpMeshMovingStrategy = nullptr;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Get the Default Settings object
     * Return the default FM-ALE settings
     * @return Parameters json string encapsulation the default settings
     */
    Parameters GetDefaultParameters();

    /**
     * @brief Set the Linear Solver Pointer object
     * This methods sets a pointer to the mesh moving strategy linear solver
     * @param rLinearSolverSettings Settings of the linear solver
     */
    void SetLinearSolverPointer(const Parameters &rLinearSolverSettings);

    /**
     * @brief Set the Mesh Moving Strategy object
     * This methods sets the mesh moving linear solver strategy
     */
    void SetMeshMovingStrategy();

    /**
     * @brief Initialize the MESH_DISPLACEMENT and MESH_VELOCITY values
     * Initialize the MESH_DISPLACEMENT and MESH_VELOCITY values
     * Note that both positions of the buffer are initialized to zero. This is important
     * in case the CloneTimeStep() is done in the virtual model part, since the method
     * assumes that the mesh is in the origin position when computing the MESH_VELOCITY.
     */
    void InitializeVirtualMeshValues();

    /**
     * @brief Initialize the mesh displacement fixity
     * This methods frees all the mesh displacement DOFs. It has
     * to be called before the current problem iteration fixity check.
     */
    void InitializeMeshDisplacementFixity();

    /**
     * @brief Set the virtual mesh origin model part based fixity
     * This method checks the mesh displacement fixity in the origin
     * model part and applies it to the virtual mesh before solving.
     */
    void SetMeshDisplacementFixityFromOriginModelPart();

    /**
     * @brief Set the embedded nodal mesh displacement
     * This method calls the utility that computes the nodal mesh
     * displacement from the immersed structure displacement. Then
     * it fixes such values before the ALE strategy solve call.
     */
    void SetEmbeddedNodalMeshDisplacement();

    /**
     * @brief Set the and solve the mesh movement strategy
     * After all the mesh BCs have been set, this method is called to
     * set and solve the ALE mesh movement strategy. The mesh velocity
     * calculation and virtual mesh update are performe din here.
     * @param DeltaTime time step value (required for the computation of the MESH_VELOCITY)
     */
    void SolveMeshMovementStrategy(const double DeltaTime);

    /**
     * @brief Revert the MESH_DISPLACEMENT fixity
     * Once the mesh moving strategy has been solved, this method
     * has to be called to free the MESH_DISPLACEMENT DOFs. It will
     * be set again in accordance to the immersed body displacement
     * before the next mesh movement computation
     */
    void RevertMeshDisplacementFixity();

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
    FixedMeshALEUtilities& operator=(FixedMeshALEUtilities const& rOther) = delete;

    /// Copy constructor.
    FixedMeshALEUtilities(FixedMeshALEUtilities const& rOther) = delete;

    ///@}
}; // Class FixedMeshALEUtilities

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const FixedMeshALEUtilities& rThis);

///@}
///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_FIXED_MESH_ALE_UTILITIES_H_INCLUDED  defined
