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

#include "custom_constitutive/thermal_law.h"
#include "includes/element.h"
#include "integration_coefficients_calculator.hpp"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) TransientThermalElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TransientThermalElement);

    explicit TransientThermalElement(IndexType NewId = 0);

    TransientThermalElement(IndexType             NewId,
                            GeometryType::Pointer pGeometry,
                            std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier = nullptr);

    TransientThermalElement(IndexType               NewId,
                            GeometryType::Pointer   pGeometry,
                            PropertiesType::Pointer pProperties,
                            std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier = nullptr);

    Element::Pointer Create(IndexType               NewId,
                            const NodesArrayType&   rThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override;

    void CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                              VectorType&        rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    std::unique_ptr<IntegrationCoefficientModifier> CloneIntegrationCoefficientModifier() const;

private:
    IntegrationCoefficientsCalculator mIntegrationCoefficientsCalculator;

    static void AddContributionsToLhsMatrix(MatrixType& rLeftHandSideMatrix,
                                            const BoundedMatrix<double, TNumNodes, TNumNodes>& rConductivityMatrix,
                                            const BoundedMatrix<double, TNumNodes, TNumNodes>& rCapacityMatrix,
                                            double DtTemperatureCoefficient);

    void AddContributionsToRhsVector(VectorType& rRightHandSideVector,
                                     const BoundedMatrix<double, TNumNodes, TNumNodes>& rConductivityMatrix,
                                     const BoundedMatrix<double, TNumNodes, TNumNodes>& rCapacityMatrix) const;

    BoundedMatrix<double, TNumNodes, TNumNodes> CalculateConductivityMatrix(
        const GeometryType::ShapeFunctionsGradientsType& rShapeFunctionGradients,
        const Vector&                                    rIntegrationCoefficients) const;

    BoundedMatrix<double, TNumNodes, TNumNodes> CalculateCapacityMatrix(const Vector& rIntegrationCoefficients) const;

    std::unique_ptr<GeoThermalLaw> CreateThermalLaw() const;

    [[nodiscard]] DofsVectorType GetDofs() const;

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;
};

} // namespace Kratos
