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

#if !defined(KRATOS_CONVECTION_DIFFUSION_REACTION_ELEMENT_H_INCLUDED)
#define KRATOS_CONVECTION_DIFFUSION_REACTION_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/element.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
class ConvectionDiffusionReactionElement : public Element
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = Element;

    /// Node type (default is: Node<3>)
    using NodeType = Node<3>;

    /// Geometry type (using with given NodeType)
    using GeometryType = Geometry<NodeType>;

    /// Definition of nodes container type, redefined from GeometryType
    using NodesArrayType = Geometry<NodeType>::PointsArrayType;

    /// Vector type for local contributions to the linear system
    using VectorType = Vector;

    /// Matrix type for local contributions to the linear system
    using MatrixType = Matrix;

    using IndexType = std::size_t;

    using EquationIdVectorType = std::vector<IndexType>;

    using DofsVectorType = std::vector<Dof<double>::Pointer>;

    /// Type for an array of shape function gradient matrices
    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    using ConvectionDiffusionReactionDataType = TConvectionDiffusionReactionData;

    using CurrentElementType =
        ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ConvectionDiffusionReactionElement
    KRATOS_CLASS_POINTER_DEFINITION(ConvectionDiffusionReactionElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit ConvectionDiffusionReactionElement(IndexType NewId = 0)
        : Element(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    ConvectionDiffusionReactionElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    ConvectionDiffusionReactionElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    ConvectionDiffusionReactionElement(IndexType NewId,
                                       GeometryType::Pointer pGeometry,
                                       PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    ConvectionDiffusionReactionElement(ConvectionDiffusionReactionElement const& rOther)
        : Element(rOther)
    {
    }

    /**
     * Destructor
     */
    ~ConvectionDiffusionReactionElement() override = default;

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
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<CurrentElementType>(
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
    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<CurrentElementType>(NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<CurrentElementType>(
            NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
        KRATOS_CATCH("");
    }

    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& CurrentProcessInfo) const override
    {
        if (rResult.size() != TNumNodes)
            rResult.resize(TNumNodes, false);

        const Variable<double>& r_variable =
            TConvectionDiffusionReactionData::GetScalarVariable();

        for (unsigned int i = 0; i < TNumNodes; ++i)
            rResult[i] = Element::GetGeometry()[i].GetDof(r_variable).EquationId();
    }

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList,
                    const ProcessInfo& CurrentProcessInfo) const override
    {
        if (rElementalDofList.size() != TNumNodes)
            rElementalDofList.resize(TNumNodes);

        const Variable<double>& r_variable =
            TConvectionDiffusionReactionData::GetScalarVariable();

        for (unsigned int i = 0; i < TNumNodes; ++i)
            rElementalDofList[i] = Element::GetGeometry()[i].pGetDof(r_variable);
    }

    void GetValuesVector(Vector& rValues, int Step = 0) const override
    {
        this->GetFirstDerivativesVector(rValues, Step);
    }

    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override
    {
        if (rValues.size() != TNumNodes)
            rValues.resize(TNumNodes, false);

        const GeometryType& rGeom = this->GetGeometry();
        const Variable<double>& r_variable =
            TConvectionDiffusionReactionData::GetScalarVariable();

        IndexType LocalIndex = 0;
        for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
        {
            rValues[LocalIndex++] =
                rGeom[iNode].FastGetSolutionStepValue(r_variable, Step);
        }
    }

    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override
    {
        if (rValues.size() != TNumNodes)
            rValues.resize(TNumNodes, false);

        const GeometryType& rGeom = this->GetGeometry();
        const Variable<double>& r_variable =
            TConvectionDiffusionReactionData::GetScalarRateVariable();

        IndexType LocalIndex = 0;
        for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
        {
            rValues[LocalIndex++] =
                rGeom[iNode].FastGetSolutionStepValue(r_variable, Step);
        }
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
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override
    {
        // Check sizes and initialize matrix
        if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);

        // Calculate RHS
        this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes, false);

        noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const IndexType num_gauss_points = gauss_weights.size();

        const GeometryType& r_geometry = this->GetGeometry();
        TConvectionDiffusionReactionData r_current_data(r_geometry);

        r_current_data.CalculateConstants(rCurrentProcessInfo);

        for (IndexType g = 0; g < num_gauss_points; ++g)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector gauss_shape_functions = row(shape_functions, g);

            r_current_data.CalculateGaussPointData(gauss_shape_functions, r_shape_derivatives);
            const double source = r_current_data.CalculateSourceTerm(
                gauss_shape_functions, r_shape_derivatives);

            for (IndexType a = 0; a < TNumNodes; ++a)
            {
                double value = 0.0;

                value += gauss_shape_functions[a] * source;

                rRightHandSideVector[a] += gauss_weights[g] * value;
            }
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief CalculateLocalVelocityContribution Calculate the local contribution in terms of velocity and pressure.
     * @param rDampMatrix Local finite element system matrix (output)
     * @param rRightHandSideVector Local finite element residual vector (output)
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    void CalculateLocalVelocityContribution(MatrixType& rDampingMatrix,
                                            VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateDampingMatrix(rDampingMatrix, rCurrentProcessInfo);
    }

    /**
     * ELEMENTS inherited from this class must implement this methods
     * if they need to add dynamic element contributions
     * note: second derivatives means the accelerations if the displacements are the dof of the analysis
     * note: time integration parameters must be set in the rCurrentProcessInfo before calling these methods
     * CalculateSecondDerivativesContributions,
     * CalculateSecondDerivativesLHS, CalculateSecondDerivativesRHS methods are : OPTIONAL
     */

    /**
     * this is called during the assembling process in order
     * to calculate the elemental mass matrix
     * @param rMassMatrix: the elemental mass matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override
    {
        BoundedMatrix<double, TNumNodes, TNumNodes> local_matrix;
        this->CalculatePrimalMassMatrix(local_matrix, rCurrentProcessInfo);

        if (rMassMatrix.size1() != TNumNodes || rMassMatrix.size2() != TNumNodes)
            rMassMatrix.resize(TNumNodes, TNumNodes, false);

        noalias(rMassMatrix) = local_matrix;
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental damping matrix
     * @param rDampingMatrix: the elemental damping matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateDampingMatrix(MatrixType& rDampingMatrix,
                                const ProcessInfo& rCurrentProcessInfo) override
    {
        BoundedMatrix<double, TNumNodes, TNumNodes> local_matrix;
        this->CalculatePrimalDampingMatrix(local_matrix, rCurrentProcessInfo);

        if (rDampingMatrix.size1() != TNumNodes || rDampingMatrix.size2() != TNumNodes)
            rDampingMatrix.resize(TNumNodes, TNumNodes, false);

        noalias(rDampingMatrix) = local_matrix;
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
        TConvectionDiffusionReactionData::Check(this->GetGeometry(), rCurrentProcessInfo);

        return check;

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates scalar value for given gauss point
     *
     * @param rVariable      Scalar variable
     * @param rShapeFunction Gauss point shape functions
     * @param Step           Step
     * @return double        Gauss point scalar value
     */
    double EvaluateInPoint(const Variable<double>& rVariable,
                           const Vector& rShapeFunction,
                           const int Step = 0) const
    {
        return RansCalculationUtilities::EvaluateInPoint(
            this->GetGeometry(), rVariable, rShapeFunction, Step);
    }

    /**
     * @brief Calculates vector value for given gauss point
     *
     * @param rVariable            Vector variable
     * @param rShapeFunction       Gauss point shape functions
     * @param Step                 Step
     * @return array_1d<double, 3> Gauss point vector value
     */
    array_1d<double, 3> EvaluateInPoint(const Variable<array_1d<double, 3>>& rVariable,
                                        const Vector& rShapeFunction,
                                        const int Step = 0) const
    {
        return RansCalculationUtilities::EvaluateInPoint(
            this->GetGeometry(), rVariable, rShapeFunction, Step);
    }

    /**
     * @brief Get the Divergence Operator object
     *
     * Calculates divergence of a vector at a gauss point
     *
     * @param rVariable          Vector variable
     * @param rShapeDerivatives  Shape derivatives at gauss point
     * @param Step               time step
     * @return double            Divergence of the variable
     */
    double GetDivergenceOperator(const Variable<array_1d<double, 3>>& rVariable,
                                 const Matrix& rShapeDerivatives,
                                 const int Step = 0) const
    {
        double value = 0.0;
        const GeometryType& r_geometry = this->GetGeometry();

        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            const array_1d<double, 3>& r_value =
                r_geometry[i].FastGetSolutionStepValue(rVariable, Step);
            for (IndexType j = 0; j < TDim; ++j)
            {
                value += r_value[j] * rShapeDerivatives(i, j);
            }
        }

        return value;
    }

    void CalculatePrimalDampingMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& rDampingMatrix,
                                      const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        rDampingMatrix.clear();

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const IndexType num_gauss_points = gauss_weights.size();

        const GeometryType& r_geometry = this->GetGeometry();
        TConvectionDiffusionReactionData r_current_data(r_geometry);

        r_current_data.CalculateConstants(rCurrentProcessInfo);

        for (IndexType g = 0; g < num_gauss_points; ++g)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector gauss_shape_functions = row(shape_functions, g);

            r_current_data.CalculateGaussPointData(gauss_shape_functions, r_shape_derivatives);
            const array_1d<double, 3>& velocity = r_current_data.CalculateEffectiveVelocity(
                gauss_shape_functions, r_shape_derivatives);
            BoundedVector<double, TNumNodes> velocity_convective_terms;
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            const double effective_kinematic_viscosity =
                r_current_data.CalculateEffectiveKinematicViscosity(
                    gauss_shape_functions, r_shape_derivatives);

            const double reaction = r_current_data.CalculateReactionTerm(
                gauss_shape_functions, r_shape_derivatives);

            for (IndexType a = 0; a < TNumNodes; ++a)
            {
                for (IndexType b = 0; b < TNumNodes; ++b)
                {
                    double dNa_dNb = 0.0;
                    for (IndexType i = 0; i < TDim; ++i)
                        dNa_dNb += r_shape_derivatives(a, i) * r_shape_derivatives(b, i);

                    double value = 0.0;

                    value += gauss_shape_functions[a] * velocity_convective_terms[b];
                    value += gauss_shape_functions[a] * reaction *
                             gauss_shape_functions[b];
                    value += effective_kinematic_viscosity * dNa_dNb;

                    rDampingMatrix(a, b) += gauss_weights[g] * value;
                }
            }
        }

        r_current_data.UpdateElementDataValueContainer(*this);

        KRATOS_CATCH("");
    }

    void CalculatePrimalMassMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& rMassMatrix,
                                   const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        rMassMatrix.clear();

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const IndexType num_gauss_points = gauss_weights.size();

        for (IndexType g = 0; g < num_gauss_points; ++g)
        {
            const double mass = gauss_weights[g] * (1.0 / TNumNodes);
            this->AddLumpedMassMatrix(rMassMatrix, mass);
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Get the Convection Operator object
     *
     * Calculates convection operator given by following equation
     *
     * \[
     *  w_i\frac{\partial N^a}{\partial x_i}
     * \]
     *
     * $w_i$ being the $i^{th}$ dimension of $\underline{w}$ vector, $N^a$ being the
     * shape function of $a^{th}$ node, $x_i$ being the $i^{th}$ dimension
     * of local coordinates
     *
     * @param rOutput           Vector of results
     * @param rVector           Input vector (i.e. $\underline{w}$)
     * @param rShapeDerivatives Shape function derivatives w.r.t. physical coordinates
     */
    void GetConvectionOperator(BoundedVector<double, TNumNodes>& rOutput,
                               const array_1d<double, 3>& rVector,
                               const Matrix& rShapeDerivatives) const
    {
        rOutput.clear();
        for (IndexType i = 0; i < TNumNodes; ++i)
            for (IndexType j = 0; j < TDim; ++j)
            {
                rOutput[i] += rVector[j] * rShapeDerivatives(i, j);
            }
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ConvectionDiffusionReactionElement #" << Id();
        return buffer.str();
    }
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CDR" << TConvectionDiffusionReactionData::GetName();
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
    virtual void CalculateGeometryData(Vector& rGaussWeights,
                                       Matrix& rNContainer,
                                       ShapeFunctionDerivativesArrayType& rDN_DX) const
    {
        const GeometryType& r_geometry = this->GetGeometry();

        RansCalculationUtilities::CalculateGeometryData(
            r_geometry, TConvectionDiffusionReactionData::GetIntegrationMethod(),
            rGaussWeights, rNContainer, rDN_DX);
    }

    void AddLumpedMassMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& rMassMatrix,
                             const double Mass) const
    {
        for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
            rMassMatrix(iNode, iNode) += Mass;
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

}; // Class ConvectionDiffusionReactionElement

///@}
///@name Input and output
///@{

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
inline std::istream& operator>>(
    std::istream& rIStream,
    ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>& rThis);

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_CONVECTION_DIFFUSION_REACTION_ELEMENT_H_INCLUDED defined
