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

#if !defined(KRATOS_RANS_EVM_K_EPSILON_LOW_RE_K_H_INCLUDED)
#define KRATOS_RANS_EVM_K_EPSILON_LOW_RE_K_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/stabilized_convection_diffusion_reaction.h"

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

struct RansEvmKEpsilonLowReKElementData
{
    double KinematicViscosity;
    double WallDistance;
    double Gamma;
    double TurbulentKineticEnergy;
    double TurbulentKinematicViscosity;
};

template <unsigned int TDim, unsigned int TNumNodes>
class RansEvmKEpsilonLowReKElement
    : public StabilizedConvectionDiffusionReaction<TDim, TNumNodes, RansEvmKEpsilonLowReKElementData>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType =
        StabilizedConvectionDiffusionReaction<TDim, TNumNodes, RansEvmKEpsilonLowReKElementData>;

    /// Node type (default is: Node<3>)
    using NodeType = Node<3>;

    /**
     * Properties are used to store any parameters
     * related to the constitutive law
     */
    using PropertiesType = Properties;

    /// Geometry type (using with given NodeType)
    using GeometryType = Geometry<NodeType>;

    /// Definition of nodes container type, redefined from GeometryType
    using NodesArrayType = Geometry<NodeType>::PointsArrayType;

    /// Vector type for local contributions to the linear system
    using VectorType = Vector;

    using IndexType = std::size_t;

    using EquationIdVectorType = std::vector<IndexType>;

    using DofsVectorType = std::vector<Dof<double>::Pointer>;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of RansEvmKEpsilonLowReKElement
    KRATOS_CLASS_POINTER_DEFINITION(RansEvmKEpsilonLowReKElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit RansEvmKEpsilonLowReKElement(IndexType NewId = 0) : BaseType(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    RansEvmKEpsilonLowReKElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : BaseType(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    RansEvmKEpsilonLowReKElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    RansEvmKEpsilonLowReKElement(IndexType NewId,
                                 GeometryType::Pointer pGeometry,
                                 PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    RansEvmKEpsilonLowReKElement(RansEvmKEpsilonLowReKElement const& rOther)
        : BaseType(rOther)
    {
    }

    /**
     * Destructor
     */
    ~RansEvmKEpsilonLowReKElement() override = default;

    ///@}
    ///@name Operators
    ///@{

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

    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult: the elemental equation ID vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) override;

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo) override;

    void GetValuesVector(VectorType& rValues, int Step = 0) override;

    void GetFirstDerivativesVector(VectorType& values, int Step = 0) override;

    void GetSecondDerivativesVector(VectorType& values, int Step = 0) override;

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

    void CalculateElementData(RansEvmKEpsilonLowReKElementData& rData,
                              const Vector& rShapeFunctions,
                              const Matrix& rShapeFunctionDerivatives,
                              const ProcessInfo& rCurrentProcessInfo,
                              const int Step = 0) const override;

    double CalculateEffectiveKinematicViscosity(const RansEvmKEpsilonLowReKElementData& rData,
                                                const Vector& rShapeFunctions,
                                                const Matrix& rShapeFunctionDerivatives,
                                                const ProcessInfo& rCurrentProcessInfo,
                                                const int Step = 0) const override;

    double CalculateReactionTerm(const RansEvmKEpsilonLowReKElementData& rData,
                                 const Vector& rShapeFunctions,
                                 const Matrix& rShapeFunctionDerivatives,
                                 const ProcessInfo& rCurrentProcessInfo,
                                 const int Step = 0) const override;

    double CalculateSourceTerm(const RansEvmKEpsilonLowReKElementData& rData,
                               const Vector& rShapeFunctions,
                               const Matrix& rShapeFunctionDerivatives,
                               const ProcessInfo& rCurrentProcessInfo,
                               const int Step = 0) const override;

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

#endif // KRATOS_RANS_EVM_K_EPSILON_LOW_RE_K_H_INCLUDED  defined
