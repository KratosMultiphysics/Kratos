/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Anna Bauer
//                  Thomas Oberbichler
//                  Tobias Teschemacher
*/

#if !defined(KRATOS_TRUSS_DISCRETE_ELEMENT_H_INCLUDED )
#define  KRATOS_TRUSS_DISCRETE_ELEMENT_H_INCLUDED

// System includes

// External includes
#include "includes/define.h"
#include "includes/element.h"

// Project includes
#include "custom_elements/curve_base_discrete_element.h"

namespace Kratos
{

/**
 * @class TrussDiscreteElement
 *
 * @brief This is a 3D-X-node isogeometric truss element with 3 translational
 * dofs per node
 */
class TrussDiscreteElement
    : public CurveBaseDiscreteElement
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TrussDiscreteElement);

    TrussDiscreteElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : CurveBaseDiscreteElement(NewId, pGeometry)
    {};

    TrussDiscreteElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : CurveBaseDiscreteElement(NewId, pGeometry, pProperties)
    {};

    TrussDiscreteElement()
        : CurveBaseDiscreteElement()
    {};

    virtual ~TrussDiscreteElement() override
    {};

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared<TrussDiscreteElement>(NewId,
            GetGeometry().Create(ThisNodes), pProperties);
    };

    static constexpr std::size_t DofsPerNode();

    std::size_t NumberOfNodes() const;

    std::size_t NumberOfDofs() const;

    void Initialize() override;

    /**
    * @brief This functions calculates both the RHS and the LHS
    * @param rLeftHandSideMatrix: The LHS
    * @param rRightHandSideVector: The RHS
    * @param rCurrentProcessInfo: The current process info instance
    * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
    * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
    */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag) override;

protected:

private:
    friend class Serializer;

    void save(
        Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
    }

    void load(
        Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
    }
}; // class TrussDiscreteElement

} // namespace Kratos

#endif // !defined(KRATOS_TRUSS_DISCRETE_ELEMENT_H_INCLUDED)
