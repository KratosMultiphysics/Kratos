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
#include "includes/define.h"
#include "includes/model_part.h"
#include "linear_solvers/linear_solver.h"
#include "processes/calculate_embedded_nodal_variable_from_skin_process.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/strategies/solving_strategy.h"

// Application includes
#include "custom_strategies/strategies/laplacian_meshmoving_strategy.h"

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

  /// Utility to initialize the historical data in moving boundary problems
  /** This utility is based on the Fixed Mesh - Arbitrary Lagrangian Eulerian
   * (FM-ALE) method but solving the mesh problem in an explicit manner. Thus,
   * a virtual mesh is set. This virtual mesh is moved according to the embedded
   * object movement. The virtual mesh movement, is computed in an explicit
   * manner as a weighted average. Such weights are computed by means of a
   * kernel function. Once the mesh movement (and velocity) have been computed,
   * the origin mesh historical values (velocity and pressure) are computed as
   * an interpolation in the virtualmodel part.
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
    typedef LaplacianMeshMovingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> LaplacianMeshMovingStrategyType;
    typedef CalculateEmbeddedNodalVariableFromSkinProcess<array_1d<double, 3>, SparseSpaceType, LocalSpaceType, LinearSolverType> EmbeddedNodalVariableProcessArrayType;

    /// Pointer definition of FixedMeshALEUtilities
    KRATOS_CLASS_POINTER_DEFINITION(FixedMeshALEUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    FixedMeshALEUtilities(
        ModelPart &rVirtualModelPart,
        ModelPart &rStructureModelPart,
        const std::string LevelSetType);

    /// Destructor.
    ~FixedMeshALEUtilities() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // /**
    // * This method performs the explicit mesh movement (computes the MESH_DISPLACEMENT value and moves
    // * the mesh accordingly) and computes the MESH_VELOCITY values.
    // * @param DeltaTime time step value (used in the computation of the MESH_VELOCITY values)
    // */
    void ComputeMeshMovement(const double DeltaTime);

    /**
    * This method fills the mrVirtualModelPart with the nodes and elmens of a given model part
    * It is supposed to be performed once.
    * @param rOriginModelPart model part from where the nodes and elements are copied
    */
    void FillVirtualModelPart(ModelPart& rOriginModelPart);

    /**
    * This method undoes the performed mesh movement to recover the original mesh in
    */
    void UndoMeshMovement();

    /**
    * This method projects the virtual model part mesh values to the origin mesh
    * @param rOriginModelPart model part to where the values are projected
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

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    const std::string mLevelSetType;

    ModelPart &mrVirtualModelPart;
    ModelPart &mrStructureModelPart;
    ModelPart *mpOriginModelPart = nullptr;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    const Vector SetDistancesVector(ModelPart::ElementIterator ItElem);

    inline bool IsSplit(const Vector &rDistances);

    void GetOriginModelPartMeshDisplacementFixity();

    void SetEmbeddedMeshDisplacement();

    void SetAndSolveMeshMovementStrategy(const double DeltaTime);

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
    FixedMeshALEUtilities& operator=(FixedMeshALEUtilities const& rOther);

    /// Copy constructor.
    FixedMeshALEUtilities(FixedMeshALEUtilities const& rOther);

    ///@}

}; // Class FixedMeshALEUtilities

///@}

///@name Type Definitions
///@{


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
