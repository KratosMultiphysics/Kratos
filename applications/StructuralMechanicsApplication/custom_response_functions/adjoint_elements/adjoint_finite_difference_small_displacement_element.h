// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//                   Daniel Baumgaertner, https://github.com/dbaumgaertner
//


#if !defined(ADJOINT_FINITE_DIFFERENCE_SMALL_DISPLACEMENT_ELEMENT_H_INCLUDED )
#define  ADJOINT_FINITE_DIFFERENCE_SMALL_DISPLACEMENT_ELEMENT_H_INCLUDED

#include "adjoint_finite_difference_base_element.h"

namespace Kratos
{


/** \brief AdjointFiniteDifferencingSmallDisplacementElement
 *
 * This element is inherited from AdjointFiniteDifferencingBaseElement.
 * It overwrites some functions necessary to do proper finite differencing with small displacement elements
 */
template <typename TPrimalElement>
class AdjointFiniteDifferencingSmallDisplacementElement : public AdjointFiniteDifferencingBaseElement<TPrimalElement>
{
public:

    ///@name Type Definitions
    ///@{

    // redefine the typedefs because of templated base class
    typedef AdjointFiniteDifferencingBaseElement<TPrimalElement> BaseType;
    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::GeometryType GeometryType;
    typedef typename BaseType::NodeType NodeType;
    typedef typename BaseType::PropertiesType PropertiesType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::VectorType VectorType;
    typedef typename BaseType::MatrixType MatrixType;
    typedef typename BaseType::EquationIdVectorType EquationIdVectorType;
    typedef typename BaseType::DofsVectorType DofsVectorType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::IntegrationMethod IntegrationMethod;
    typedef typename BaseType::GeometryDataType GeometryDataType;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AdjointFiniteDifferencingSmallDisplacementElement);

    AdjointFiniteDifferencingSmallDisplacementElement(IndexType NewId = 0)
    : BaseType(NewId, false)
    {
    }

    AdjointFiniteDifferencingSmallDisplacementElement(IndexType NewId, typename GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry, false)
    {
    }

    AdjointFiniteDifferencingSmallDisplacementElement(IndexType NewId,
                                                    typename GeometryType::Pointer pGeometry,
                                                    typename PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties, false)
    {
    }

    /**
     * Calculates the stress-displacement derivative of the given rStressVariable.
     * For the linear case this happens by calculating the stress using a unit-displacement for each dof.
     */
    void CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                               Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointFiniteDifferencingSmallDisplacementElement<TPrimalElement>>(
            NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId,
                            typename GeometryType::Pointer pGeometry,
                            typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointFiniteDifferencingSmallDisplacementElement<TPrimalElement>>(
            NewId, pGeometry, pProperties);
    }

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

private:

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};


}

#endif