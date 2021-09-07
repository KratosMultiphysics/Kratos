// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:		 BSD License
//					 license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_DERIVATIVES_UTILITIES)
#define KRATOS_DERIVATIVES_UTILITIES

// System includes

// External includes

// Project includes
#include "contact_structural_mechanics_application_variables.h"

/* Includes */
#include "includes/model_part.h"
#include "includes/mortar_classes.h"

/* Geometries */
#include "geometries/line_2d_2.h"
#include "geometries/triangle_3d_3.h"

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
/**
 * @class DerivativesUtilities
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This utilities are used in order to compute the directional derivatives during mortar contact
 * @details The derivatives take the same argument templates than the contact conditions
 * @author Vicente Mataix Ferrandiz
 * @author Gabriel Valdes Alonzo
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TNormalVariation If the normal variation is considered
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template< const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster = TNumNodes>
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) DerivativesUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// The index type
    typedef std::size_t                                                                                                  IndexType;

    /// The geometry of nodes
    typedef Geometry<NodeType>                                                                                        GeometryType;

    /// The array of nodes contained in a geometry
    typedef Geometry<NodeType>::PointsArrayType                                                                     NodesArrayType;

    /// The belong type (for derivatives definition)
    typedef typename std::conditional<TNumNodes == 2, PointBelongsLine2D2N, typename std::conditional<TNumNodes == 3, typename std::conditional<TNumNodesMaster == 3, PointBelongsTriangle3D3N, PointBelongsTriangle3D3NQuadrilateral3D4N>::type, typename std::conditional<TNumNodesMaster == 3, PointBelongsQuadrilateral3D4NTriangle3D3N, PointBelongsQuadrilateral3D4N>::type>::type>::type BelongType;

    /// The points used for derivatives definition
    typedef PointBelong<TNumNodes, TNumNodesMaster>                                                                PointBelongType;

    /// A geometry defined by the point belongs (the points used for derivatives definition)
    typedef Geometry<PointBelongType>                                                                      GeometryPointBelongType;

    /// An array of belong point to define the geometries of belong points
    typedef array_1d<PointBelongType,TDim>                                                                      ConditionArrayType;

    /// The definition of an array pf conditions
    typedef typename std::vector<ConditionArrayType>                                                        ConditionArrayListType;

    /// The line definition
    typedef Line2D2<PointType>                                                                                            LineType;

    /// The triangle definition
    typedef Triangle3D3<PointType>                                                                                    TriangleType;

    /// The geometry for decomposition (line in 2D and triangle for 3D)
    typedef typename std::conditional<TDim == 2, LineType, TriangleType >::type                                  DecompositionType;

    /// The derivative data type
    typedef typename std::conditional<TFrictional, DerivativeDataFrictional<TDim, TNumNodes, TNumNodesMaster>, DerivativeData<TDim, TNumNodes, TNumNodesMaster> >::type DerivativeDataType;

    /// The kinematic variables
    typedef MortarKinematicVariablesWithDerivatives<TDim, TNumNodes, TNumNodesMaster>                             GeneralVariables;

    /// The dual LM operators
    typedef DualLagrangeMultiplierOperatorsWithDerivatives<TDim, TNumNodes, TFrictional, TNumNodesMaster>                   AeData;

    /// The mortar operators
    typedef MortarOperatorWithDerivatives<TDim, TNumNodes, TFrictional, TNumNodesMaster>                   MortarConditionMatrices;

    /// Definition of epsilon
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    /// Pointer definition of DerivativesUtilities
    KRATOS_CLASS_POINTER_DEFINITION( DerivativesUtilities );

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
    ///@name Operations
    ///@{

    /**
     * @brief This method is used to compute the directional derivatives of the jacobian determinant
     * @param DecompGeom The triangle used to decompose the geometry
     * @param rVariables The kinematic variables
     * @param rDerivativeData The derivatives container
     */
    static void CalculateDeltaDetjSlave(
        const DecompositionType& DecompGeom,
        const GeneralVariables& rVariables,
        DerivativeDataType& rDerivativeData
        );

    /**
     * @brief This method is used to compute the local increment of the normal
     * @param rJacobian The jacobian on the GP
     * @param rDNDe The local gradient
     * @return The matrix containing the delta normals
     * @note Not the mean, look in the contact utilities
     */
    static inline array_1d<array_1d<double, 3>, TDim * TNumNodes> GPDeltaNormalSlave(
        const Matrix& rJacobian,
        const Matrix& rDNDe
        );

    /**
     * @brief This method is used to compute the local increment of the normal
     * @param rJacobian The jacobian on the GP
     * @param rDNDe The local gradient
     * @return The matrix containing the delta normals
     * @note Not the mean, look in the contact utilities
     * @note Hardcopied for performance
     */
    static inline array_1d<array_1d<double, 3>, TDim * TNumNodesMaster> GPDeltaNormalMaster(
        const Matrix& rJacobian,
        const Matrix& rDNDe
        );

    /**
     * @brief It computes the delta normal of the center of the geometry
     * @param rThisGeometry The geometry where the delta normal is computed
     */
    static array_1d<array_1d<double, 3>, TDim * TNumNodes> DeltaNormalCenter(const GeometryType& rThisGeometry);

    /**
     * @brief Calculates the increment of the normal and in the slave condition
     * @param rDeltaNormal The derivative of the normal
     * @param rThisGeometry The geometry of the salve side
     */
    static void CalculateDeltaNormalSlave(
        array_1d<BoundedMatrix<double, TNumNodes, TDim>, TNumNodes * TDim>& rDeltaNormal,
        GeometryType& rThisGeometry
        );

    /**
     * @brief Calculates the increment of the normal and in the master condition
     * @param rDeltaNormal The derivative of the normal
     * @param rThisGeometry The geometry of the master side
     * @note Hardcopied for performance
     */
    static void CalculateDeltaNormalMaster(
        array_1d<BoundedMatrix<double, TNumNodesMaster, TDim>, TNumNodesMaster * TDim>& rDeltaNormal,
        const GeometryType& rThisGeometry
        );

    /**
     * @brief This method is used to compute the directional derivatives of the cell vertex
     * @param rVariables The kinematic variables
     * @param rDerivativeData The derivatives container
     * @param rTheseBelongs The belongs list used in the derivatives
     * @param ConsiderNormalVariation If consider the normal derivative
     * @param rSlaveGeometry The slave geometry
     * @param rMasterGeometry The master geometry
     * @param rNormal The normal vector of the slave geometry
     * @details The  procedure will be the following in order to compute the derivative of the clipping
     * The expression of the clipping is the following:
     *      xclipp = xs1 - num/denom * diff3
     * Being:
     *      diff1 = xs1 - xs2;
     *      diff2 = xe2 - xs2;
     *      diff3 = xe1 - xs1;
     *      num = (diff1 x diff2) · n0
     *      denom = (diff3 x diff2) · n0
     * The derivative can be defined then as:
     *     delta_xclipp = delta_xs1 - (delta_num * denom - num * delta_denom)/delta_denom**2 * diff3 - num/ denom * delta_diff3
     * And here:
     *     delta_num = num · delta_n0 + n0 · (delta_diff1 x diff2 + diff1 x delta_diff2)
     *     delta_denom = denom · delta_n0 + n0 · (delta_diff3 x diff2 + diff3 x delta_diff2)
     */
    static void CalculateDeltaCellVertex(
        const GeneralVariables& rVariables,
        DerivativeDataType& rDerivativeData,
        const array_1d<BelongType, TDim>& rTheseBelongs,
        const NormalDerivativesComputation ConsiderNormalVariation,
        const GeometryType& rSlaveGeometry,
        const GeometryType& rMasterGeometry,
        const array_1d<double, 3>& rNormal
        );

    /**
     * @brief Calculates the increment of the shape functions
     * @param rVariables The kinematic variables
     * @param rDerivativeData The derivatives container
     * @param rSlaveGeometry The geometry of the slave side
     * @param rMasterGeometry The geometry of the master side
     * @param rSlaveNormal The normal of the slave side
     * @param rDecompGeom The triangle used to decompose the geometry
     * @param rLocalPointDecomp The local coordinates in the decomposed geometry
     * @param rLocalPointParent The local coordinates in the slave geometry
     * @param ConsiderNormalVariation If consider the normal derivative
     */
    static inline void CalculateDeltaN1(
        const GeneralVariables& rVariables,
        DerivativeDataType& rDerivativeData,
        const GeometryType& rSlaveGeometry,
        const GeometryType& rMasterGeometry,
        const array_1d<double, 3> rSlaveNormal,
        const DecompositionType& rDecompGeom,
        const PointType& rLocalPointDecomp,
        const PointType& rLocalPointParent,
        const NormalDerivativesComputation ConsiderNormalVariation = NO_DERIVATIVES_COMPUTATION
        );

    /**
     * @brief Calculates the increment of the shape functions
     * @param rVariables The kinematic variables
     * @param rDerivativeData The derivatives container
     * @param rSlaveGeometry The geometry of the slave side
     * @param rMasterGeometry The geometry of the master side
     * @param rSlaveNormal The normal of the slave side
     * @param rMasterNormal The normal of the master side
     * @param rDecompGeom The triangle used to decompose the geometry
     * @param rLocalPointDecomp The local coordinates in the decomposed geometry
     * @param rLocalPointParent The local coordinates in the slave geometry
     * @param ConsiderNormalVariation If consider the normal derivative
     * @param DualLM If the dual Lm formulation is considered
     */
    static void CalculateDeltaN(
        const GeneralVariables& rVariables,
        DerivativeDataType& rDerivativeData,
        const GeometryType& rSlaveGeometry,
        const GeometryType& rMasterGeometry,
        const array_1d<double, 3>& rSlaveNormal,
        const array_1d<double, 3>& rMasterNormal,
        const DecompositionType& rDecompGeom,
        const PointType& rLocalPointDecomp,
        const PointType& rLocalPointParent,
        const NormalDerivativesComputation ConsiderNormalVariation = NO_DERIVATIVES_COMPUTATION,
        const bool DualLM = false
        );

    /**
     * @brief Returns a matrix with the increment of displacements, that can be used for compute the Jacobian reference (current) configuration
     * @param rDeltaPosition The matrix with the increment of displacements
     * @param rThisGeometry The geometry considered
     * @param rLocalCoordinates The array containing the local coordinates of the exact integration segment
     */
    Matrix& CalculateDeltaPosition(
        Matrix& rDeltaPosition,
        const GeometryType& rThisGeometry,
        const ConditionArrayType& rLocalCoordinates
        );

    /**
     * @brief Returns a matrix with the increment of displacements
     * @param rDeltaPosition The matrix with the increment of displacements
     * @param ThisGeometry The geometry considered
     */
    static inline Matrix& CalculateDeltaPosition(
        Matrix& rDeltaPosition,
        const GeometryType& ThisGeometry
        );

    /**
     * @brief Returns a vector with the increment of displacements
     * @param rDeltaPosition The resulting vector with the increment of position
     * @param rSlaveGeometry The slave geometry
     * @param rMasterGeometry The master geometry
     * @param IndexNode The node index
     */
    static inline void CalculateDeltaPosition(
        Vector& rDeltaPosition,
        const GeometryType& rSlaveGeometry,
        const GeometryType& rMasterGeometry,
        const IndexType IndexNode
        );

    /**
     * @brief Returns a vector with the increment of displacements
     * @param rDeltaPosition The resulting vector with the increment of position
     * @param rSlaveGeometry The slave geometry
     * @param rMasterGeometry The master geometry
     * @param IndexNode The node index
     * @param iDoF The degree of freedom index
     */
    static inline void CalculateDeltaPosition(
        Vector& rDeltaPosition,
        const GeometryType& rSlaveGeometry,
        const GeometryType& rMasterGeometry,
        const IndexType IndexNode,
        const IndexType iDoF
        );

    /**
     * @brief Returns a double with the increment of displacements
     * @param rDeltaPosition The resulting double with the increment of position
     * @param rSlaveGeometry The slave geometry
     * @param rMasterGeometry The master geometry
     * @param IndexNode The node index
     * @param iDoF The degree of freedom index
     */
    static inline void CalculateDeltaPosition(
        double& rDeltaPosition,
        const GeometryType& rSlaveGeometry,
        const GeometryType& rMasterGeometry,
        const IndexType IndexNode,
        const IndexType iDoF
        );

    /**
     * @brief Calculate Ae and DeltaAe matrices
     * @param rSlaveGeometry The geometry of the slave side
     * @param rSlaveNormal The normal of the slave side
     * @param rMasterGeometry The master side geometry
     * @param rDerivativeData The derivatives container
     * @param rVariables The kinematic variables
     * @param ConsiderNormalVariation If consider the normal derivative
     * @param rConditionsPointsSlave The points that configure the exact decomposition of the geometry
     * @param ThisIntegrationMethod The integration method considered
     * @param AxiSymCoeff The axisymmetric coefficient
     */
    static bool CalculateAeAndDeltaAe(
        const GeometryType& rSlaveGeometry,
        const array_1d<double, 3>& rSlaveNormal,
        const GeometryType& rMasterGeometry,
        DerivativeDataType& rDerivativeData,
        GeneralVariables& rVariables,
        const NormalDerivativesComputation ConsiderNormalVariation,
        ConditionArrayListType& rConditionsPointsSlave,
        GeometryData::IntegrationMethod ThisIntegrationMethod,
        const double AxiSymCoeff = 1.0
        );

private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method is used to compute the directional derivatives of the cell vertex (locally)
     * @param rNormal The normal of the slave surface
     * @param rDeltaNormal The derivative of the normal vector
     * @param iDoF The DoF computed index
     * @param iBelong The belong (intersection, node, etc..) index
     * @param ConsiderNormalVariation If the normal variation is considered
     * @param rSlaveGeometry The geometry of the slave side
     * @param rMasterGeometry The geometry of the master side
     * @param Coeff The coefficient considered in proportion
     * @return The local vertex derivative
     */
    static inline array_1d<double, 3> LocalDeltaVertex(
        const array_1d<double, 3>& rNormal,
        const array_1d<double, 3>& rDeltaNormal,
        const IndexType iDoF,
        const IndexType iBelong,
        const NormalDerivativesComputation ConsiderNormalVariation,
        const GeometryType& rSlaveGeometry,
        const GeometryType& rMasterGeometry,
        double Coeff = 1.0
        );

    /**
     * @brief This method computes the auxiliar matrix used to keep unitary the normal
     * @param rDiffVector The auxiliar vector of difference of two normal vectors
     * @param rDeltaNormal The vector containing the delta normal
     * @return The auxiliar matrix computed
     */
    static inline BoundedMatrix<double, 3, 3> ComputeRenormalizerMatrix(
        const array_1d<double, 3>& rDiffVector,
        const array_1d<double, 3>& rDeltaNormal
        );

    /**
     * @brief This method computes the auxiliar matrix used to keep unitary the normal
     * @param rDiffMatrix The auxiliar matrix of difference of two normal matrices
     * @param rDeltaNormal The matrix containing the delta normal
     * @param iGeometry The index of the node of the geometry computed
     * @return The auxiliar matrix computed
     */
    static inline BoundedMatrix<double, 3, 3> ComputeRenormalizerMatrix(
        const BoundedMatrix<double, TNumNodes, TDim>& rDiffMatrix,
        const BoundedMatrix<double, TNumNodes, TDim>& rDeltaNormal,
        const IndexType iGeometry
        );

    /**
     * @brief This method computes the normal in the previous configuration
     * @param rThisGeometry The geometry where compute
     * @param rPointLocal The local coordinates of the point
     * @return The normal in the previous configuration
     */
    static inline array_1d<double, 3> PreviousNormalGeometry(
        const GeometryType& rThisGeometry,
        const GeometryType::CoordinatesArrayType& rPointLocal
        );

    /**
     * @brief This method computes the equivalent indexes to the auxiliar hash
     * @param AuxIndex The auxiliar index to decompose
     * @param riBelongSlaveStart The index of the first/slave segment and first node
     * @param riBelongSlaveEnd The index of the first/slave segment and end node
     * @param riBelongMasterStart The index of the second/master segment and first node
     * @param riBelongMasterEnd The index of the second/master segment and end node
     */
    static inline void ConvertAuxHashIndex(
        const IndexType AuxIndex,
        IndexType& riBelongSlaveStart,
        IndexType& riBelongSlaveEnd,
        IndexType& riBelongMasterStart,
        IndexType& riBelongMasterEnd
        );

    /**
     * @brief This method computes the increment of local coordinates
     * @param rResult The solution obtained
     * @param rDeltaPoint The increment of position in the points
     * @param rThisGeometry The geometry considered
     * @param rThisNormal The normal of the geometry
     * @note Hardcopied for performance
     */
    static inline void DeltaPointLocalCoordinatesSlave(
        array_1d<double, 2>& rResult,
        const array_1d<double, 3>& rDeltaPoint,
        const Matrix& rDNDe,
        const GeometryType& rThisGeometry,
        const array_1d<double, 3>& rThisNormal
        );
    /**
     * @brief This method computes the increment of local coordinates
     * @param rResult The solution obtained
     * @param rDeltaPoint The increment of position in the points
     * @param rThisGeometry The geometry considered (master)
     * @param rThisNormal The normal of the geometry
     * @note Hardcopied for performance
     */
    static inline void DeltaPointLocalCoordinatesMaster(
        array_1d<double, 2>& rResult,
        const array_1d<double, 3>& rDeltaPoint,
        const Matrix& rDNDe,
        const GeometryType& rThisGeometry,
        const array_1d<double, 3>& rThisNormal
        );

    /**
     * @brief This method is used to compute the directional derivatives of the mortar segment on the slave side
     * @param rDeltaNormal All derivatives of normals in points of both geometries
     * @param rSlaveNormal Normal in the center of slave geometry
     * @param rSlaveGeometry Slave geometry where to compute the value
     * @param rMasterGeometry Master node projected to the Slave Geometry
     * @param rN1 Shape function at SlaveGeometry
     * @param rDNDe1 Gradient of shape function N1
     * @param MortarNode Index of mortar node where computation is being carried
     * @param iNode Actual node where computation is being carried
     * @param iDoF Direction in which computation is carried
     * @param ConsiderNormalVariation Flag to determine if consider delta_normal
     * @return The mortar node derivative
     */
    static inline double LocalDeltaSegmentN1(
        const array_1d<array_1d<double, 3>, TDim * TNumNodes>& rDeltaNormal,
        const array_1d<double, 3>& rSlaveNormal,
        const GeometryType& rSlaveGeometry,
        const GeometryType& rMasterGeometry,
        const Vector& rN1,
        const Matrix& rDNDe1,
        const IndexType MortarNode,
        const IndexType iNode,
        const IndexType iDoF,
        const NormalDerivativesComputation ConsiderNormalVariation
        );

    /**
     * @brief This method is used to compute the directional derivatives of the mortar segment on the master side
     * @param rDeltaNormal All derivatives of normals in points of both geometries
     * @param rSlaveNormal Normal in the center of slave geometry
     * @param rSlaveGeometry Slave geometry where to compute the value
     * @param rMasterGeometry Master node projected to the Slave Geometry
     * @param rN2 Shape function at MasterGeometry
     * @param rDNDe2 Gradient of shape function N2
     * @param MortarNode Index of mortar node where computation is being carried
     * @param iNode Actual node where computation is being carried
     * @param iDoF Direction in which computation is carried
     * @param ConsiderNormalVariation flag to determine if consider delta_normal
     * @return The mortar node derivative
     */
     static inline double LocalDeltaSegmentN2(
         const array_1d<array_1d<double, 3>, TDim * TNumNodes>& rDeltaNormal,
         const array_1d<double, 3>& rSlaveNormal,
         const GeometryType& rSlaveGeometry,
         const GeometryType& rMasterGeometry,
         const Vector& rN2,
         const Matrix& rDNDe2,
         const IndexType MortarNode,
         const IndexType iNode,
         const IndexType iDoF,
         const NormalDerivativesComputation ConsiderNormalVariation
         );

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
};// class DerivativesUtilities

}
#endif /* KRATOS_DERIVATIVES_UTILITIES defined */
