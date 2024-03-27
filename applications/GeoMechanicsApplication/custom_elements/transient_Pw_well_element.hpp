// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   John van Esch
//

#pragma once

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) TransientPwWellElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TransientPwWellElement);

    explicit TransientPwWellElement(IndexType NewId = 0) : Element(NewId) {}

    TransientPwWellElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes)
    {
    }

    TransientPwWellElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    TransientPwWellElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    /// Destructor
    ~TransientPwWellElement() override {}

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;


    struct ElementVariables {
        /// Properties variables
        double dynamicViscosity = 0.0;
        double waterDensity = 0.0;
        double wellLength = 0.0;
        double wellDiameter = 0.0;
        double compressibilityFactor = 0.0;
        double waterCompressibility  = 0.0;

        double BiotCoefficient;
        double BiotModulusInverse;
        BoundedMatrix<double, TNumNodes, TNumNodes> compressibilityMatrix;
        BoundedMatrix<double, TNumNodes, TNumNodes> permeabilityMatrix;
        array_1d<double, TNumNodes> compressibilityVector;
        array_1d<double, TNumNodes> permeabilityVector;
        BoundedMatrix<double, 1, 1> permeabilityTensor;
        array_1d<double, TNumNodes> bodyForceVector;

        double DtPressureCoefficient;

        array_1d<double, TNumNodes> pressureVector;
        array_1d<double, TNumNodes> DtPressureVector;
        array_1d<double, TNumNodes * 1> VolumeAcceleration;

        Vector N;
        Matrix GradNT;
        Matrix GradNpTInitialConfiguration;

        Matrix F;
        double detF;
        Vector detJContainer;
        Matrix NContainer;
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer;

        double detJ;
        double detJInitialConfiguration;
        double IntegrationCoefficient;
        double IntegrationCoefficientInitialConfiguration;
    };

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateAll(MatrixType& rLeftHandSideMatrix,
                      VectorType& rRightHandSideVector,
                      const ProcessInfo& CurrentProcessInfo,
                      const bool CalculateStiffnessMatrixFlag,
                      const bool CalculateResidualVectorFlag);

    void InitializeElementVariables(ElementVariables& rVariables, const ProcessInfo& CurrentProcessInfo);

    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateKinematics(ElementVariables& rVariables, unsigned int PointNumber);

    void CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);
    void CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);
    void CalculateAndAddCompressibilityVector(VectorType& rRightHandSideVector, ElementVariables& rVariables);
    void CalculateAndAddPermeabilityVector(VectorType& rRightHandSideVector, ElementVariables& rVariables);
    //void CalculateAndAddFluidBodyVector(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void InitializeNodalPorePressureVariables(ElementVariables& rVariables);
    void InitializeNodalVolumeAccelerationVariables(ElementVariables& rVariables);
    double CalculateIntegrationCoefficient(const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
                                           unsigned int PointNumber,
                                           double detJ);

    void CalculateCompressibilityMatrix(ElementVariables& rVariables);
    void CalculatePermeabilityMatrix(ElementVariables& rVariables);
    void CalculateCompressibilityVector(ElementVariables& rVariables);
    void CalculatePermeabilityVector(ElementVariables& rVariables);

    void InitializeProperties(ElementVariables& rVariables);

    unsigned int GetNumberOfDOF() const;

    void CalculateCompressibilityFactor(ElementVariables& rVariables);
    void CalculatePermeabilityTensor(ElementVariables& rVariables);
    void CalculateBodyForceVector(ElementVariables& rVariables);
    void CalculateAndAddBodyForceVector(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo);

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;
    double CalculateIntegrationCoefficient(const Matrix& Jacobian, const double& Weight);

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    /// Member Variables
    bool mIsThermalCoupled       = false;
    bool mUpdateDensityViscosity = true;
    bool mIsInitialised          = false;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
    }

    //// Assignment operator.
    //TransientPwWellElement&
    //operator=(TransientPwWellElement const& rOther);
    //
    //// Copy constructor.
    //TransientPwWellElement(TransientPwWellElement const& rOther);

}; // Class TransientPwWellElement

} // namespace Kratos