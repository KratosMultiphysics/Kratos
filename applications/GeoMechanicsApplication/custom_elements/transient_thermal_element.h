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
//                   Gennady Markelov
//

#pragma once

#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/serializer.h"

namespace Kratos {

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) TransientThermalElement : public Element {
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TransientThermalElement);

    struct ElementVariables {
        double WaterDensity;
        double SolidDensity;
        double WaterHeatCapacity;
        double SolidHeatCapacity;
        double Porosity;
        double Saturation;

        double DtTemperatureCoefficient;
        array_1d<double, TNumNodes> TemperatureVector;
        array_1d<double, TNumNodes> DtTemperatureVector;
        Matrix ConstitutiveMatrix;
        Vector N;
        Matrix GradNT;
        Matrix GradNTInitialConfiguration;
        Vector detJContainer;
        Matrix NContainer;
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer;
        double detJ;
        double IntegrationCoefficient;
        BoundedMatrix<double, TNumNodes, TNumNodes> ConductivityMatrix;
        BoundedMatrix<double, TNumNodes, TNumNodes> CapacityMatrix;
        array_1d<double, TNumNodes> ConductivityVector;
        array_1d<double, TNumNodes> CapacityVector;
    };

    explicit TransientThermalElement(IndexType NewId = 0);

    TransientThermalElement(IndexType NewId, GeometryType::Pointer pGeometry);

    TransientThermalElement(IndexType NewId,
                            GeometryType::Pointer pGeometry,
                            PropertiesType::Pointer pProperties);

    TransientThermalElement& operator=(TransientThermalElement const& rOther) = delete;

    TransientThermalElement(TransientThermalElement const& rOther) = delete;

    ~TransientThermalElement() override;

    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& rThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& rElementalDofList,
                    const ProcessInfo& rCurrentProcessInfo) const override;

    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& rCurrentProcessInfo) const override;

protected:
    void CalculateAll(MatrixType& rLeftHandSideMatrix,
                      VectorType& rRightHandSideVector,
                      const ProcessInfo& CurrentProcessInfo,
                      const bool CalculateStiffnessMatrixFlag,
                      const bool CalculateResidualVectorFlag);

    void InitializeElementVariables(ElementVariables& rVariables,
                                    const ProcessInfo& CurrentProcessInfo);

    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateKinematics(ElementVariables& rVariables, unsigned int PointNumber);

    void CalculateAndAddConductivityMatrix(MatrixType& rLeftHandSideMatrix,
                                           ElementVariables& rVariables);

    void CalculateAndAddCapacityMatrix(MatrixType& rLeftHandSideMatrix,
                                       ElementVariables& rVariables);

    void CalculateAndAddConductivityVector(VectorType& rRightHandSideVector,
                                           ElementVariables& rVariables);

    void CalculateAndAddCapacityVector(VectorType& rRightHandSideVector,
                                       ElementVariables& rVariables);

    void InitializeNodalTemperatureVariables(ElementVariables& rVariables);

    unsigned int GetNumberOfDOF() const;

    virtual double CalculateIntegrationCoefficient(const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
                                                   unsigned int PointNumber,
                                                   double detJ);

    void InitializeProperties(ElementVariables& rVariables);

    virtual void CalculateConductivityMatrix(ElementVariables& rVariables);

    virtual void CalculateCapacityMatrix(ElementVariables& rVariables) const;

    virtual void CalculateCapacityVector(ElementVariables& rVariables) const;

    virtual void CalculateConductivityVector(ElementVariables& rVariables);

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

private:
    void VerifyProperty(Kratos::Variable<double>& rVariable) const;
    void CheckDomainSize() const;
    void CheckSolutionStepsData(int rId, Kratos::Variable<double>& rVariable) const;

    bool mIsInitialised = false;
    const ProcessInfo * mpCurrentProcessInfo;

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;
};

} // namespace Kratos
