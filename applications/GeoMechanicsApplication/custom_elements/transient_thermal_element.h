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
        double DtTemperatureCoefficient;
        array_1d<double, TNumNodes> TemperatureVector;
        array_1d<double, TNumNodes> DtTemperatureVector;
        Vector detJContainer;
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer;
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

    void GetDofList(DofsVectorType& rElementalDofList,
                    const ProcessInfo& rCurrentProcessInfo) const override;

    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;
private:
    void CalculateAll(MatrixType& rLeftHandSideMatrix,
                      VectorType& rRightHandSideVector,
                      const ProcessInfo& CurrentProcessInfo);

    void InitializeElementVariables(ElementVariables& rVariables,
                                    const ProcessInfo& CurrentProcessInfo);

    void InitializeNodalTemperatureVariables(ElementVariables& rVariables);

    Vector CalculateIntegrationCoefficients(const Vector& detJContainer) const;

    BoundedMatrix<double, TNumNodes, TNumNodes> CalculateConductivityMatrix(const GeometryType::ShapeFunctionsGradientsType& rShapeFunctionGradients,
                                                                            const Vector& rIntegrationCoefficients,
                                                                            const ProcessInfo& rCurrentProcessInfo) const;

    BoundedMatrix<double, TNumNodes, TNumNodes> CalculateCapacityMatrix(const Vector& rIntegrationCoefficients) const;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    void VerifyProperty(Kratos::Variable<double>& rVariable) const;
    void CheckDomainSize() const;
    void CheckSolutionStepsData(int rId, Kratos::Variable<double>& rVariable) const;

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;
};

} // namespace Kratos
