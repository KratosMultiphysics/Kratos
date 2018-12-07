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

#if !defined(KRATOS_IGA_TRUSS_ELEMENT_H_INCLUDED)
#define KRATOS_IGA_TRUSS_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "iga_base_element.h"


namespace Kratos
{

class IgaTrussElement
    : public IgaBaseElement<3>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION( IgaTrussElement );

    using IgaBaseElementType::IgaBaseElementType;

    ~IgaTrussElement() override
    {
    };

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override;

    void GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo) override;

    void Initialize() override;

    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool ComputeLeftHandSide,
        const bool ComputeRightHandSide) override;

    void PrintInfo(std::ostream& rOStream) const override;

private:
    Vector3 mReferenceBaseVector;

    Vector3 GetActualBaseVector() const;
};

} // namespace Kratos

#endif // !defined(KRATOS_IGA_TRUSS_ELEMENT_H_INCLUDED)
