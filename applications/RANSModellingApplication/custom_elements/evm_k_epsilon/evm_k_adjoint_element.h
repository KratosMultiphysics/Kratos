//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_EVM_K_ADJOINT_ELEMENT_H_INCLUDED)
#define KRATOS_EVM_K_ADJOINT_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/stabilized_convection_diffusion_reaction_adjoint_element.h"
#include "includes/checks.h"
#include "includes/element.h"
#include "includes/properties.h"

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

struct EvmKAdjointElementData
{
    double KinematicViscosity;
    double EffectiveKinematicViscosity;
    double WallDistance;
    double Gamma;
    double Fmu;
    double y_plus;
    double VelocityDivergence;

    double TurbulentKineticEnergy;
    double TurbulentKinematicViscosity;

    Matrix ShapeFunctionDerivatives;
    Vector ShapeFunctions;

    Vector NodalTurbulentKineticEnergy;
    Vector NodalTurbulentEnergyDissipationRate;
    Vector NodalFmu;
    Vector NodalYPlus;
    Matrix NodalVelocity;
};

template <unsigned int TDim, unsigned int TNumNodes>
class EvmKAdjointElement
    : public StabilizedConvectionDiffusionReactionAdjointElement<TDim, TNumNodes, EvmKAdjointElementData>
{
public:
    ///@name Type Definitions
    ///@{

    typedef StabilizedConvectionDiffusionReactionAdjointElement<TDim, TNumNodes, EvmKAdjointElementData> BaseType;

    /// Node type (default is: Node<3>)
    typedef Node<3> NodeType;

    /**
     * Properties are used to store any parameters
     * related to the constitutive law
     */
    typedef Properties PropertiesType;

    /// Geometry type (using with given NodeType)
    typedef Geometry<NodeType> GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    /// Vector type for local contributions to the linear system
    typedef Vector VectorType;

    /// Matrix type for local contributions to the linear system
    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector<Dof<double>::Pointer> DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    /// Type for shape function values container
    typedef MatrixRow<Matrix> ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    typedef Kratos::Matrix ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of EvmKAdjointElement
    KRATOS_CLASS_POINTER_DEFINITION(EvmKAdjointElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    EvmKAdjointElement(IndexType NewId = 0);

    /**
     * Constructor using an array of nodes
     */
    EvmKAdjointElement(IndexType NewId, const NodesArrayType& ThisNodes);

    /**
     * Constructor using Geometry
     */
    EvmKAdjointElement(IndexType NewId, GeometryType::Pointer pGeometry);

    /**
     * Constructor using Properties
     */
    EvmKAdjointElement(IndexType NewId,
                       GeometryType::Pointer pGeometry,
                       PropertiesType::Pointer pProperties);

    /**
     * Copy Constructor
     */
    EvmKAdjointElement(EvmKAdjointElement const& rOther);

    /**
     * Destructor
     */
    ~EvmKAdjointElement() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    EvmKAdjointElement& operator=(EvmKAdjointElement const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /**
     * ELEMENTS inherited from this class have to implement next
     * Create and Clone methods: MANDATORY
     */

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override;

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

    void Calculate(const Variable<Matrix>& rVariable,
                   Matrix& Output,
                   const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief GetIntegrationMethod Return the integration order to be used.
     * @return Gauss Order
     */
    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the contitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     * this method is: MANDATORY
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

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
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    const Variable<double>& GetPrimalVariable() const override;

    const Variable<double>& GetPrimalRelaxedRateVariable() const override;

    const Variable<double>& GetAdjointVariable() const override;

    const Variable<double>& GetAdjointSecondVariable() const override;

    void CalculateElementData(EvmKAdjointElementData& rData,
                              const Vector& rShapeFunctions,
                              const Matrix& rShapeFunctionDerivatives,
                              const ProcessInfo& rCurrentProcessInfo) const override;

    double CalculateEffectiveKinematicViscosity(const EvmKAdjointElementData& rCurrentData,
                                                const ProcessInfo& rCurrentProcessInfo) const override;

    double CalculateReactionTerm(const EvmKAdjointElementData& rCurrentData,
                                 const ProcessInfo& rCurrentProcessInfo) const override;

    double CalculateSourceTerm(const EvmKAdjointElementData& rCurrentData,
                               const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateEffectiveKinematicViscosityScalarDerivatives(
        Vector& rOutput,
        const Variable<double>& rDerivativeVariable,
        const EvmKAdjointElementData& rCurrentData,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateReactionTermScalarDerivatives(Vector& rOutput,
                                                const Variable<double>& rDerivativeVariable,
                                                const EvmKAdjointElementData& rCurrentData,
                                                const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateSourceTermScalarDerivatives(Vector& rOutput,
                                              const Variable<double>& rDerivativeVariable,
                                              const EvmKAdjointElementData& rCurrentData,
                                              const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateEffectiveKinematicViscosityVelocityDerivatives(
        Matrix& rOutput,
        const EvmKAdjointElementData& rCurrentData,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateReactionTermVelocityDerivatives(Matrix& rOutput,
                                                  const EvmKAdjointElementData& rCurrentData,
                                                  const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateSourceTermVelocityDerivatives(Matrix& rOutput,
                                                const EvmKAdjointElementData& rCurrentData,
                                                const ProcessInfo& rCurrentProcessInfo) const override;

    double CalculateEffectiveKinematicViscosityShapeSensitivity(
        const EvmKAdjointElementData& rCurrentData,
        const ShapeParameter& rShapeDerivative,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
        const ProcessInfo& rCurrentProcessInfo) const override;

    double CalculateReactionTermShapeSensitivity(
        const EvmKAdjointElementData& rCurrentData,
        const ShapeParameter& rShapeDerivative,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
        const ProcessInfo& rCurrentProcessInfo) const override;

    double CalculateSourceTermShapeSensitivity(
        const EvmKAdjointElementData& rCurrentData,
        const ShapeParameter& rShapeDerivative,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
        const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class EvmTurbulentKineticEnergyElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_EVM_K_ADJOINT_ELEMENT_H_INCLUDED  defined
