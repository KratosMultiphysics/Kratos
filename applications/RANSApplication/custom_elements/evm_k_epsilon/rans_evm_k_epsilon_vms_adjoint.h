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

#if !defined(KRATOS_RANS_EVM_K_EPSILON_VMS_ADJOINT_ELEMENT_H_INCLUDED)
#define KRATOS_RANS_EVM_K_EPSILON_VMS_ADJOINT_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/rans_evm_vms_adjoint.h"

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

template <std::size_t TNumNodes>
struct RANSEvmVMSAdjointData
{
    Matrix ShapeFunctionDerivatives;
    Vector ShapeFunctions;

    BoundedVector<double, TNumNodes> TurbulentKinematicViscositySensitivitiesK;
    BoundedVector<double, TNumNodes> TurbulentKinematicViscositySensitivitiesEpsilon;
};

template <unsigned int TDim, unsigned int TNumNodes = TDim + 1>
class KRATOS_API(RANS_APPLICATION) RansEvmKEpsilonVMSAdjoint
    : public RANSEvmVMSAdjoint<TDim, RANSEvmVMSAdjointData<TNumNodes>, TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    constexpr static unsigned int TFluidLocalSize = (TDim + 1) * TNumNodes;

    typedef RANSEvmVMSAdjoint<TDim, RANSEvmVMSAdjointData<TNumNodes>, TNumNodes> BaseType;

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

    /// Type for shape function values container
    typedef MatrixRow<Matrix> ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    typedef Kratos::Matrix ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of RansEvmKEpsilonVMSAdjoint
    KRATOS_CLASS_POINTER_DEFINITION(RansEvmKEpsilonVMSAdjoint);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit RansEvmKEpsilonVMSAdjoint(IndexType NewId = 0) : BaseType(NewId)
    {
    }

    /**
     * Constructor using Geometry
     */
    RansEvmKEpsilonVMSAdjoint(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    RansEvmKEpsilonVMSAdjoint(IndexType NewId,
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Destructor
     */
    ~RansEvmKEpsilonVMSAdjoint() override = default;

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

    void Calculate(const Variable<Matrix>& rVariable,
                   Matrix& Output,
                   const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief GetIntegrationMethod Return the integration order to be used.
     * @return Gauss Order
     */

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

    void CalculateElementData(RANSEvmVMSAdjointData<TNumNodes>& rData,
                              const Vector& rShapeFunctions,
                              const Matrix& rShapeFunctionDerivatives,
                              const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateTurbulentKinematicViscosityScalarDerivatives(
        BoundedVector<double, TNumNodes>& rOutput,
        const Variable<double>& rDerivativeVariable,
        const RANSEvmVMSAdjointData<TNumNodes>& rCurrentData,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateTurbulentKinematicViscosityVelocityDerivatives(
        BoundedMatrix<double, TNumNodes, TDim>& rOutput,
        const RANSEvmVMSAdjointData<TNumNodes>& rCurrentData,
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

#endif // KRATOS_RANS_EVM_K_EPSILON_VMS_ADJOINT_ELEMENT_H_INCLUDED  defined
