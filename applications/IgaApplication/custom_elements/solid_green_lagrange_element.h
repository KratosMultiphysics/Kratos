//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:       Thomas Lillieskold
//  Co-author:          Andrea Gorgi
//
#pragma once

#include "includes/element.h"

namespace Kratos
{

/// IGA solid in the total Lagrangian formulation.
/// The 2D formulation uses the XY plane. The 3D formulation is pending.
class KRATOS_API(IGA_APPLICATION) SolidGreenLagrangeElement : public Element
{
public:
    //--------------------------------------------[1. Creación del elemento]--------------------------------------------

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SolidGreenLagrangeElement);

    SolidGreenLagrangeElement() = default;

    SolidGreenLagrangeElement(IndexType NewId, GeometryType::Pointer pGeometry);

    SolidGreenLagrangeElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties);

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& rThisNodes,
        PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) const override;

    //--------------------------------------------[2. Inicialización, cuadratura y grados de libertad]--------------------------------------------

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(
        DofsVectorType& rDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    IntegrationMethod GetIntegrationMethod() const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    std::string Info() const override;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;
    
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

//--------------------------------------------[3. Cinemática de Green–Lagrange]--------------------------------------------

protected:
    struct KinematicVariables
    {
        Matrix F = IdentityMatrix(2);
        double DeterminantF = 1.0;
        Vector Strain = ZeroVector(3);
    };

    KinematicVariables CalculateKinematics() const;

private:
    void CheckGeometryDimension() const;

    // One geometry (and therefore one quadrature point) per IGA element.
    Matrix mDN_DX;
    //--------------------------------------------[4. Respuesta constitutiva (pendiente)]--------------------------------------------

    void InitializeMaterial();

    struct ConstitutiveVariables
    
    {
        ConstitutiveLaw::StressVectorType StressVector;
        ConstitutiveLaw::VoigtSizeMatrixType D;

        explicit ConstitutiveVariables(SizeType StrainSize)
            : StressVector(ZeroVector(StrainSize)),
              D(ZeroMatrix(StrainSize, StrainSize))
        {
        }
    };

    void ApplyConstitutiveLaw(
        KinematicVariables& rKinematics,
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveVariables& rConstitutiveVariables);

    ///@}
    ///@name Protected member Variables
    ///@{

    ConstitutiveLaw::Pointer mpConstitutiveLaw; /// The pointer containing the constitutive laws


    ///@}
    ///@name Protected member Variables
    ///@{


    //--------------------------------------------[5. Ensamblado de volumen (pendiente)]--------------------------------------------

};

} // namespace Kratos
