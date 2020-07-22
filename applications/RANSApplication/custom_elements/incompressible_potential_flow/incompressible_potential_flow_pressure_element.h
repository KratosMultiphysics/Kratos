//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_PRESSURE_ELEMENT_H_INCLUDED)
#define KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_PRESSURE_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/incompressible_potential_flow/laplace_element.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/variables.h"

// Application includes
#include "rans_application_variables.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes>
class IncompressiblePotentialFlowPressureElement
: public LaplaceElement<TDim, TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = LaplaceElement<TDim, TNumNodes>;

    using NodeType = Node<3>;

    /// Geometry type (using with given NodeType)
    using GeometryType = Geometry<NodeType>;

    using PropertiesType = Properties;

    using MatrixType = Element::MatrixType;

    using VectorType = Element::VectorType;

    /// Definition of nodes container type, redefined from GeometryType
    using NodesArrayType = GeometryType::PointsArrayType;

    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    using IndexType = std::size_t;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of IncompressiblePotentialFlowPressureElement
    KRATOS_CLASS_POINTER_DEFINITION(IncompressiblePotentialFlowPressureElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit IncompressiblePotentialFlowPressureElement(
        IndexType NewId = 0)
    : BaseType(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    IncompressiblePotentialFlowPressureElement(
        IndexType NewId,
        const NodesArrayType& ThisNodes)
    : BaseType(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    IncompressiblePotentialFlowPressureElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    IncompressiblePotentialFlowPressureElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    IncompressiblePotentialFlowPressureElement(
        IncompressiblePotentialFlowPressureElement const& rOther)
    : BaseType(rOther)
    {
    }

    /**
     * Destructor
     */
    ~IncompressiblePotentialFlowPressureElement() override = default;

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
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<IncompressiblePotentialFlowPressureElement>(
            NewId, Element::GetGeometry().Create(ThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<IncompressiblePotentialFlowPressureElement>(
            NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& ThisNodes) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<IncompressiblePotentialFlowPressureElement>(
            NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
        KRATOS_CATCH("");
    }

    const Variable<double>& GetVariable() const override
    {
        return PRESSURE_POTENTIAL;
    }

    /**
     * ELEMENTS inherited from this class have to implement next
     * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide
     * methods they can be managed internally with a private method to do the
     * same calculations only once: MANDATORY
     */

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        // Calculate RHS
        this->CalculateRightHandSideVelocityContribution(rRightHandSideVector);

        // Calculate LHS
        this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

        VectorType values;
        this->GetValuesVector(values, 0);

        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, values);
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        this->CalculateRightHandSideVelocityContribution(rRightHandSideVector);

        Matrix lhs;
        this->CalculateLeftHandSide(lhs, rCurrentProcessInfo);

        Vector values;
        this->GetValuesVector(values, 0);

        noalias(rRightHandSideVector) -= prod(lhs, values);

        KRATOS_CATCH("");
    }

    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        int check = BaseType::Check(rCurrentProcessInfo);

        for (IndexType i_node = 0; i_node < this->GetGeometry().size(); ++i_node) {
            const NodeType& r_node = this->GetGeometry()[i_node];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY_POTENTIAL, r_node);
        }

        return check;

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "IncompressiblePotentialFlowPressureElement #" << this->Id();
        return buffer.str();
    }
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IncompressiblePotentialFlowPressureElement #" << this->Id();
    }

    ///@}

protected:
    ///@name Private Operations
    ///@{

    void CalculateRightHandSideVelocityContribution(Vector& rRHS) const
    {
        if (rRHS.size() != TNumNodes) {
            rRHS.resize(TNumNodes, false);
        }
        noalias(rRHS) = ZeroVector(TNumNodes);

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const IndexType num_gauss_points = gauss_weights.size();

        const auto& r_geometry = this->GetGeometry();

        for (IndexType g = 0; g < num_gauss_points; ++g) {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector gauss_shape_functions = row(shape_functions, g);

            const double density = RansCalculationUtilities::EvaluateInPoint(
                r_geometry, DENSITY, gauss_shape_functions);

            array_1d<double, 3> kinetic_energy_gradient;
            // VELOCITY_POTENTIAL contains velocity magnitude square
            RansCalculationUtilities::CalculateGradient(
                kinetic_energy_gradient, r_geometry, VELOCITY_POTENTIAL, r_shape_derivatives);
            noalias(kinetic_energy_gradient) =
                kinetic_energy_gradient * (gauss_weights[g] * 0.5 * density);

            for (IndexType a = 0; a < TNumNodes; ++a) {
                double value = 0.0;
                for (IndexType d = 0; d < TDim; ++d) {
                    value += r_shape_derivatives(a, d) * kinetic_energy_gradient[d];
                }

                rRHS[a] -= value;
            }
        }
    }

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_TRY

        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);

        KRATOS_CATCH("");
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_TRY

        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);

        KRATOS_CATCH("");
    }

    ///@}
}; // Class IncompressiblePotentialFlowPressureElement

///@}

///@name Input and output
///@{

template <unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator>>(
    std::istream& rIStream,
    IncompressiblePotentialFlowPressureElement<TDim, TNumNodes>& rThis);

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const IncompressiblePotentialFlowPressureElement<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_PRESSURE_ELEMENT_H_INCLUDED defined
