// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//


#pragma once

#include "adjoint_finite_difference_base_element.h"

namespace Kratos
{

/** \brief AdjointFiniteDifferencingBaseElement
 *
 * This element is a wrapper for a primal linear truss element. It is responsible to deliver local stresses and
 * the stress displacement derivative. It is designed to be used in adjoint
 * sensitivity analysis.
 */
template <typename TPrimalElement>
class AdjointFiniteDifferenceCableElement
    : public AdjointFiniteDifferencingBaseElement<TPrimalElement>
{
public:

    // redefine the typedefs because of templated base class
    typedef AdjointFiniteDifferencingBaseElement<TPrimalElement> BaseType;
    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::GeometryType GeometryType;
    typedef typename BaseType::PropertiesType PropertiesType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::VectorType VectorType;
    typedef typename BaseType::MatrixType MatrixType;
    typedef typename BaseType::EquationIdVectorType EquationIdVectorType;
    typedef typename BaseType::DofsVectorType DofsVectorType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::IntegrationMethod IntegrationMethod;
    typedef typename BaseType::GeometryDataType GeometryDataType;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AdjointFiniteDifferenceCableElement);

    AdjointFiniteDifferenceCableElement(IndexType NewId = 0)
    : BaseType(NewId)
    {
    }

    AdjointFiniteDifferenceCableElement(IndexType NewId, typename GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
    {
    }

    AdjointFiniteDifferenceCableElement(IndexType NewId,
                        typename GeometryType::Pointer pGeometry,
                        typename PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
    {
    }

    Element::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointFiniteDifferenceCableElement<TPrimalElement>>(
            NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId,
                              typename GeometryType::Pointer pGeometry,
                              typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointFiniteDifferenceCableElement<TPrimalElement>>(
            NewId, pGeometry, pProperties);
    }

    // void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
    //              std::vector< array_1d<double, 3 > >& rOutput,
    //             const ProcessInfo& rCurrentProcessInfo) override;

    // void CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
    //                                 Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       const ProcessInfo& rCurrentProcessInfo) override;

    // void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rDesignVariable,
    //                                 Matrix& rOutput,
    //                                 const ProcessInfo& rCurrentProcessInfo) override;
    // Sensitivity functions

    /**
     * Calculates the pseudo-load contribution of the element w.r.t. all properties which are available at the element.
     * This is done by finite differencing of the RHS of the primal element when perturbing a property value.
     * This operation is thread-save, because the property pointer is exchanged by a local one before pertubation.
     */
    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculates the pseudo-load of the design variable SHAPE_SENSITIVITY (coordinates of nodes) contribution of the element.
     * This is done by finite differencing of the RHS of the primal element when perturbing a nodal coordinate.
     * This operation is currently NOT thread-save!
     */
    void CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;


private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};


}
