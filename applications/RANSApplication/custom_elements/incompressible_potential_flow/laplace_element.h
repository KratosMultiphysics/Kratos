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

#if !defined(KRATOS_LAPLACE_ELEMENT_H_INCLUDED)
#define KRATOS_LAPLACE_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/element.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes>
class LaplaceElement : public Element
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = Element;

    /// Node type (default is: Node<3>)
    using NodeType = BaseType::NodeType;

    /// Geometry type (using with given NodeType)
    using GeometryType = Geometry<NodeType>;

    using PropertiesType = Properties;

    /// Definition of nodes container type, redefined from GeometryType
    using NodesArrayType = GeometryType::PointsArrayType;

    /// Vector type for local contributions to the linear system
    using VectorType = BaseType::VectorType;

    /// Matrix type for local contributions to the linear system
    using MatrixType = BaseType::MatrixType;

    using EquationIdVectorType = BaseType::EquationIdVectorType;

    using DofsVectorType = BaseType::DofsVectorType;

    using IndexType = std::size_t;

    /// Type for an array of shape function gradient matrices
    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of LaplaceElement
    KRATOS_CLASS_POINTER_DEFINITION(LaplaceElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit LaplaceElement(
        IndexType NewId = 0)
    : Element(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    LaplaceElement(
        IndexType NewId,
        const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    LaplaceElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    LaplaceElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    LaplaceElement(LaplaceElement const& rOther)
    : Element(rOther)
    {
    }

    /**
     * Destructor
     */
    ~LaplaceElement() override = default;

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
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to Create base "
                        "LaplaceElement instances."
                     << std::endl;
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
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to Create base "
                        "LaplaceElement instances."
                     << std::endl;
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
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to Clone base "
                        "LaplaceElement instances."
                     << std::endl;
        KRATOS_CATCH("");
    }

    virtual const Variable<double>& GetVariable() const
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to call base "
                        "LaplaceElement "
                        "GetVariable method. Please implement it in the "
                        "derrived class."
                     << std::endl;
        KRATOS_CATCH("");
    }

    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult: the elemental equation ID vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& CurrentProcessInfo) const override
    {
        if (rResult.size() != TNumNodes) {
            rResult.resize(TNumNodes, false);
        }

        const Variable<double>& r_variable = this->GetVariable();

        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rResult[i] = Element::GetGeometry()[i].GetDof(r_variable).EquationId();
        }
    }

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& CurrentProcessInfo) const override
    {
        if (rElementalDofList.size() != TNumNodes) {
            rElementalDofList.resize(TNumNodes);
        }

        const Variable<double>& r_variable = this->GetVariable();

        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rElementalDofList[i] = Element::GetGeometry()[i].pGetDof(r_variable);
        }
    }

    void GetValuesVector(
        VectorType& rValues,
        int Step) const override
    {
        if (rValues.size() != TNumNodes) {
            rValues.resize(TNumNodes, false);
        }

        const Variable<double>& r_variable = this->GetVariable();

        const auto& r_geometry = this->GetGeometry();
        IndexType LocalIndex = 0;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            rValues[LocalIndex++] =
                r_geometry[i_node].FastGetSolutionStepValue(r_variable, Step);
        }
    }

    GeometryData::IntegrationMethod GetIntegrationMethod() const override
    {
        return GeometryData::GI_GAUSS_2;
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
    virtual void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        // Calculate RHS
        if (rRightHandSideVector.size() != TNumNodes) {
            rRightHandSideVector.resize(TNumNodes, false);
        }

        noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

        // Calculate LHS
        this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

        VectorType values;
        this->GetValuesVector(values, 0);

        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, values);
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rLeftHandSideMatrix.size1() != TNumNodes ||
            rLeftHandSideMatrix.size2() != TNumNodes) {
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
        }

        noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const IndexType num_gauss_points = gauss_weights.size();

        for (IndexType g = 0; g < num_gauss_points; ++g) {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector gauss_shape_functions = row(shape_functions, g);

            for (IndexType a = 0; a < TNumNodes; ++a) {
                for (IndexType b = 0; b < TNumNodes; ++b) {
                    double dNa_dNb = 0.0;
                    for (IndexType i = 0; i < TDim; ++i) {
                        dNa_dNb += r_shape_derivatives(a, i) * r_shape_derivatives(b, i);
                    }

                    rLeftHandSideMatrix(a, b) += gauss_weights[g] * dNa_dNb;
                }
            }
        }

        KRATOS_CATCH("");
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    virtual void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rRightHandSideVector.size() != TNumNodes) {
            rRightHandSideVector.resize(TNumNodes, false);
        }

        Matrix lhs;
        this->CalculateLeftHandSide(lhs, rCurrentProcessInfo);

        Vector values;
        this->GetValuesVector(values, 0);

        noalias(rRightHandSideVector) = prod(lhs, values) * -1.0;

        KRATOS_CATCH("");
    }

    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the contitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     * this method is: MANDATORY
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        int check = BaseType::Check(rCurrentProcessInfo);

        const Variable<double>& r_variable = this->GetVariable();

        for (IndexType i_node = 0; i_node < this->GetGeometry().size(); ++i_node) {
            const NodeType& r_node = this->GetGeometry()[i_node];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_variable, r_node);
            KRATOS_CHECK_DOF_IN_NODE(r_variable, r_node);
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
        buffer << "LaplaceElement #" << this->Id();
        return buffer.str();
    }
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "LaplaceElement #" << this->Id();
    }

    ///@}

protected:
    ///@name Protected Operations
    ///@{
    /**
     * @brief Calculates shape function data for this element
     *
     * @param rGaussWeights Gauss point weights list
     * @param rNContainer   Shape function values. Each row contains shape functions for respective gauss point
     * @param rDN_DX        List of matrices containing shape function derivatives for each gauss point
     */
    virtual void CalculateGeometryData(
        Vector& rGaussWeights,
        Matrix& rNContainer,
        ShapeFunctionDerivativesArrayType& rDN_DX) const
    {
        const auto& r_geometry = this->GetGeometry();

        RansCalculationUtilities::CalculateGeometryData(
            r_geometry, this->GetIntegrationMethod(), rGaussWeights, rNContainer, rDN_DX);
    }

    ///@}

private:
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
}; // Class LaplaceElement

///@}

///@name Input and output
///@{

template <unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator>>(
    std::istream& rIStream,
    LaplaceElement<TDim, TNumNodes>& rThis);

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const LaplaceElement<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_LAPLACE_ELEMENT_H_INCLUDED defined
