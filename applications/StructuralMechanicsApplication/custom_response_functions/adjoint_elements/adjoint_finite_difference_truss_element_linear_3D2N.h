// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//


#if !defined(ADJOINT_FINITE_DIFFERENCE_TRUSS_ELEMENT_LINEAR_H_INCLUDED )
#define  ADJOINT_FINITE_DIFFERENCE_TRUSS_ELEMENT_LINEAR_H_INCLUDED

#include "adjoint_finite_difference_truss_element_3D2N.h"

namespace Kratos
{

/** \brief AdjointFiniteDifferencingBaseElement
 *
 * This element is a wrapper for a primal linear truss element. It is responsible to deliver local stresses and
 * the stress displacement derivative. It is designed to be used in adjoint
 * sensitivity analysis.
 */
template <typename TPrimalElement>
class AdjointFiniteDifferenceTrussElementLinear
    : public AdjointFiniteDifferenceTrussElement<TPrimalElement>
{
public:

    // redefine the typedefs because of templated base class
    typedef AdjointFiniteDifferenceTrussElement<TPrimalElement> BaseType;
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

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AdjointFiniteDifferenceTrussElementLinear);

    AdjointFiniteDifferenceTrussElementLinear(IndexType NewId = 0)
    : BaseType(NewId)
    {
    }

    AdjointFiniteDifferenceTrussElementLinear(IndexType NewId, typename GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
    {
    }

    AdjointFiniteDifferenceTrussElementLinear(IndexType NewId,
                        typename GeometryType::Pointer pGeometry,
                        typename PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
    {
    }

    Element::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>>(
            NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId,
                              typename GeometryType::Pointer pGeometry,
                              typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>>(
            NewId, pGeometry, pProperties);
    }

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                 std::vector< array_1d<double, 3 > >& rOutput,
                const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo) override;


private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};


}

#endif
