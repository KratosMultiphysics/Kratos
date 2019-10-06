//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#if !defined(KRATOS_GENERIC_SMALL_STRAIN_FEMDEM_ELEMENT_H_INCLUDED)
#define KRATOS_GENERIC_SMALL_STRAIN_FEMDEM_ELEMENT_H_INCLUDED


// System includes


// External include

// Project includes

#include "custom_elements/generic_total_lagrangian_femdem_element.h"
#include "custom_utilities/constitutive_law_utilities.h"

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
 * @class GenericSmallStrainFemDemElement
 * @ingroup FemToDemApplication
 * @brief Small Displacement element for the 2D and 3D cases
 * @author Alejandro Cornejo
 */
template<unsigned int TDim, unsigned int TyieldSurf>
class GenericSmallStrainFemDemElement 
    : public GenericTotalLagrangianFemDemElement<TDim,TyieldSurf> // Derived Element from SolidMechanics
{
public:
    ///@name Type Definitions
    ///@{

    ///definition of element type
    typedef Element ElementType;

    ///base type: an GeometricalObject that automatically has a unique number
    typedef GenericTotalLagrangianFemDemElement<TDim,TyieldSurf> BaseType;

    ///definition of node type (default is: Node<3>)
    typedef Node < 3 > NodeType;

    /**
     * Properties are used to store any parameters
     * related to the constitutive law
     */
    typedef Properties PropertiesType;

    ///definition of the geometry type with given NodeType
    typedef Geometry<NodeType> GeometryType;

    ///definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    typedef GeometryData GeometryDataType;

    /// We define the dimension
    static constexpr SizeType VoigtSize = (TDim == 3) ? 6 : 3;

    /// We define the number of edges
    static constexpr SizeType NumberOfEdges = (TDim == 3) ? 6 : 3;

    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    /// The zero tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    /// Counted pointer of GenericSmallStrainFemDemElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GenericSmallStrainFemDemElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    GenericSmallStrainFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry);
    GenericSmallStrainFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    GenericSmallStrainFemDemElement(GenericSmallStrainFemDemElement const &rOther);

    /// Destructor.
    virtual ~GenericSmallStrainFemDemElement();

    /// Assignment operator.
    GenericSmallStrainFemDemElement &operator=(GenericSmallStrainFemDemElement const &rOther);
    Element::Pointer Create(IndexType NewId, NodesArrayType const &ThisNodes, PropertiesType::Pointer pProperties) const override;
    Element::Pointer Clone(IndexType NewId, NodesArrayType const &ThisNodes) const override;

    GenericSmallStrainFemDemElement()
    {
    }

    /**
     * this is called at the end of each solution step
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;


    // ************** Methods to compute the tangent constitutive tensor via numerical derivation ************** 
    /**
     * this computes the Tangent tensor via numerical derivation (perturbations)
     */
    void CalculateTangentTensor(Matrix& rTangentTensor,const Vector& rStrainVectorGP,const Vector& rStressVectorGP,const Matrix& rElasticMatrix, ConstitutiveLaw::Parameters& rValues);

protected:

    /**
     * @brief This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix The LHS
     * @param rRightHandSideVector The RHS
     * @param rCurrentProcessInfo The current process info instance
     * @param CalculateStiffnessMatrixFlag The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag The flag to set if compute the RHS
     */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        ) override;

    /**
     * @brief This functions updates the kinematics variables
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param PointNumber The integration point considered
     * @param rIntegrationMethod The integration method considered
     */
    void CalculateKinematicVariables(
        KinematicVariables& rThisKinematicVariables,
        const IndexType PointNumber,
        const GeometryType::IntegrationMethod& rIntegrationMethod
        ) override;

    /**
     * @brief This method computes the deformation matrix B
     * @param rB The deformation matrix
     * @param rF The deformation gradient
     * @param rDN_DX The gradient derivative of the shape function
     */
    void CalculateB(Matrix& rB, const Matrix& rDN_DX);
    void Calculate2DB(Matrix& rB, const Matrix& rDN_DX);
    void Calculate3DB(Matrix& rB, const Matrix& rDN_DX);
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

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
}; // Class GenericSmallStrainFemDemElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif