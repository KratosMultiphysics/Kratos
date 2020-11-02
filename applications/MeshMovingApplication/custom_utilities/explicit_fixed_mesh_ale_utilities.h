//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Ruben Zorrilla
//

#if !defined(KRATOS_EXPLICIT_FIXED_MESH_ALE_UTILITIES_H_INCLUDED )
#define  KRATOS_EXPLICIT_FIXED_MESH_ALE_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "spatial_containers/spatial_search.h"
#include "spatial_containers/bins_dynamic_objects.h"
#include "spatial_containers/configures/node_configure.h"

// Application includes
#include "custom_utilities/fixed_mesh_ale_utilities.h"

namespace Kratos
{
  ///@addtogroup FluidDynamicsApplication
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
   * an interpolation in the virtual model part.
   */
  class KRATOS_API(MESH_MOVING_APPLICATION) ExplicitFixedMeshALEUtilities : public FixedMeshALEUtilities
  {
  public:

    ///@name Type Definitions
    ///@{

    typedef Element::GeometryType::PointsArrayType                             PointsArrayType;
    typedef BinsObjectDynamic<NodeConfigure>                                      NodeBinsType;
    typedef std::vector<double>                                             DistanceVectorType;
    typedef SpatialSearch::ResultNodesContainerType                   ResultNodesContainerType;
    typedef std::vector<std::vector<double>>                       DistanceVectorContainerType;
    typedef SpatialSearch::VectorResultNodesContainerType       VectorResultNodesContainerType;

    /// Pointer definition of ExplicitFixedMeshALEUtilities
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitFixedMeshALEUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ExplicitFixedMeshALEUtilities(
        ModelPart &rVirtualModelPart,
        ModelPart &rStructureModelPart,
        const double SearchRadius);

    /// Constructor with model and parameters
    ExplicitFixedMeshALEUtilities(
        Model &rModel,
        Parameters &rParameters);

    /// Destructor.
    ~ExplicitFixedMeshALEUtilities() = default;

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
    void Initialize(ModelPart &rOriginModelPart) override;

    /**
    * This method performs the explicit mesh movement (computes the MESH_DISPLACEMENT value and moves
    * the mesh accordingly) and computes the MESH_VELOCITY values.
    * @param DeltaTime time step value (used in the computation of the MESH_VELOCITY values)
    */
    void ComputeMeshMovement(const double DeltaTime) override;

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
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Create the virtual model part elements
     * This method creates the elements in the virtual model part
     * @param rOriginModelPart Origin model part to mimic the elements from
     */
    void CreateVirtualModelPartElements(const ModelPart &rOriginModelPart) override;

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

    const double mSearchRadius;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
    * According to mSearchRadius, performs the bins search of the close structure nodes for each fluid node
    * @return rSearchResults vector containing the the rStructureModelPart nodes inside
    * @return rSearchDistanceResults vector containing the the rStructureModelPart nodes
    * inside SearchRadius distance values for each rModelPart nodes
    */
    void SearchStructureNodes(
        VectorResultNodesContainerType &rSearchResults,
        DistanceVectorContainerType &rSearchDistanceResults);

    /**
    * Computes the MESH_DISPLACEMENT value for each fluid node. This operation is explicitly computed
    * as a weighted average of the structure nodes DISPLACEMENT values within the mSearchRadius. The
    * weights are computed using a kernel function.
    * @param rSearchResults vector containing the the rStructureModelPart nodes inside
    * @param rSearchDistanceResults vector containing the the rStructureModelPart nodes
    * inside SearchRadius distance values for each rModelPart nodes
    */
    void ComputeExplicitMeshDisplacement(
        const VectorResultNodesContainerType &rSearchResults,
        const DistanceVectorContainerType &rSearchDistanceResults);

    /**
    * Computes the kernel function value for a given normalised distance value
    * @param NormalisedDistance structure node distance value normalised with the mSearchRadius
    */
    inline double ComputeKernelValue(const double NormalisedDistance);

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
    ExplicitFixedMeshALEUtilities& operator=(ExplicitFixedMeshALEUtilities const& rOther) = delete;

    /// Copy constructor.
    ExplicitFixedMeshALEUtilities(ExplicitFixedMeshALEUtilities const& rOther) = delete;

    ///@}
}; // Class ExplicitFixedMeshALEUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const ExplicitFixedMeshALEUtilities& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_EXPLICIT_FIXED_MESH_ALE_UTILITIES_H_INCLUDED  defined
