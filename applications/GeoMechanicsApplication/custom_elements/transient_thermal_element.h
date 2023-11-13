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
#include "custom_constitutive/thermal_dispersion_law.h"
#include "includes/serializer.h"

namespace Kratos {

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) TransientThermalElement : public Element {
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TransientThermalElement);

    explicit TransientThermalElement(IndexType NewId = 0) : Element(NewId) {}

    TransientThermalElement(IndexType             NewId,
                            GeometryType::Pointer pGeometry)
            : Element(NewId, pGeometry)
    {
    }

    TransientThermalElement(IndexType               NewId,
                            GeometryType::Pointer   pGeometry,
                            PropertiesType::Pointer pProperties)
            : Element(NewId, pGeometry, pProperties)
    {
    }

    ~TransientThermalElement() override = default;
    TransientThermalElement(const TransientThermalElement&) = delete;
    TransientThermalElement& operator=(const TransientThermalElement&) = delete;

    Element::Pointer Create(IndexType               NewId,
                            const NodesArrayType&   rThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        return Element::Pointer{new TransientThermalElement{NewId, GetGeometry().Create(rThisNodes), pProperties}};
    }

    Element::Pointer Create(IndexType               NewId,
                            GeometryType::Pointer   pGeom,
                            PropertiesType::Pointer pProperties) const override
    {
        return Element::Pointer{new TransientThermalElement{NewId, pGeom, pProperties}};
    }

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& rElementalDofList,
                    const ProcessInfo& ) const override
    {
        KRATOS_TRY

        rElementalDofList.resize(TNumNodes);
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rElementalDofList[i] = GetGeometry()[i].pGetDof(TEMPERATURE);
        }

        KRATOS_CATCH("")
    }

    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& ) const override
    {
        KRATOS_TRY

        rResult.resize(TNumNodes, false);
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rResult[i] = GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
        }

        KRATOS_CATCH("")
    }

    void CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                              VectorType&        rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        GeometryType::ShapeFunctionsGradientsType DN_DXContainer;
        Vector detJContainer;
        GetGeometry().ShapeFunctionsIntegrationPointsGradients(DN_DXContainer, detJContainer, GetIntegrationMethod());
        const auto integration_coefficients = CalculateIntegrationCoefficients(detJContainer);
        const auto conductivity_matrix =
                CalculateConductivityMatrix(DN_DXContainer, integration_coefficients, rCurrentProcessInfo);
        const auto capacity_matrix = CalculateCapacityMatrix(integration_coefficients);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);
        AddContributionsToLhsMatrix(rLeftHandSideMatrix, conductivity_matrix, capacity_matrix,
                                    rCurrentProcessInfo[DT_TEMPERATURE_COEFFICIENT]);

        noalias(rRightHandSideVector) = ZeroVector(TNumNodes);
        AddContributionsToRhsVector(rRightHandSideVector, conductivity_matrix, capacity_matrix);

        KRATOS_CATCH("")
    }

    GeometryData::IntegrationMethod GetIntegrationMethod() const override
    {
        switch (TNumNodes) {
            case 3:
            case 6:
                return GeometryData::IntegrationMethod::GI_GAUSS_2;
            case 10:
                return GeometryData::IntegrationMethod::GI_GAUSS_4;
            case 15:
                return GeometryData::IntegrationMethod::GI_GAUSS_5;
            default:
                return GeometryData::IntegrationMethod::GI_GAUSS_2;
        }
    }

private:
    static void AddContributionsToLhsMatrix(MatrixType&                                        rLeftHandSideMatrix,
                                            const BoundedMatrix<double, TNumNodes, TNumNodes>& rConductivityMatrix,
                                            const BoundedMatrix<double, TNumNodes, TNumNodes>& rCapacityMatrix,
                                            double                                             DtTemperatureCoefficient)
    {
        GeoElementUtilities::AssemblePBlockMatrix<0, TNumNodes>(rLeftHandSideMatrix, rConductivityMatrix);
        GeoElementUtilities::AssemblePBlockMatrix<0, TNumNodes>(rLeftHandSideMatrix,
                                                                DtTemperatureCoefficient * rCapacityMatrix);
    }

    void AddContributionsToRhsVector(VectorType&                                        rRightHandSideVector,
                                     const BoundedMatrix<double, TNumNodes, TNumNodes>& rConductivityMatrix,
                                     const BoundedMatrix<double, TNumNodes, TNumNodes>& rCapacityMatrix) const
    {
        const auto capacity_vector =
                array_1d<double, TNumNodes>{-prod(rCapacityMatrix, GetNodalValuesOf(DT_TEMPERATURE))};
        GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector, capacity_vector);

        const auto conductivity_vector =
                array_1d<double, TNumNodes>{-prod(rConductivityMatrix, GetNodalValuesOf(TEMPERATURE))};
        GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector, conductivity_vector);
    }

    Vector CalculateIntegrationCoefficients(const Vector& rDetJContainer) const
    {
        const auto& r_integration_points = GetGeometry().IntegrationPoints(GetIntegrationMethod());

        auto result = Vector{r_integration_points.size()};
        for (unsigned int GPoint = 0; GPoint < r_integration_points.size(); ++GPoint) {
            result[GPoint] = r_integration_points[GPoint].Weight() * rDetJContainer[GPoint];
        }

        return result;
    }

    BoundedMatrix<double, TNumNodes, TNumNodes> CalculateConductivityMatrix(const GeometryType::ShapeFunctionsGradientsType& rShapeFunctionGradients,
                                                                            const Vector& rIntegrationCoefficients,
                                                                            const ProcessInfo& rCurrentProcessInfo) const
    {
        GeoThermalDispersionLaw law{TDim};
        const auto constitutive_matrix =
                law.CalculateThermalDispersionMatrix(GetProperties(), rCurrentProcessInfo, GetGeometry());

        auto result = BoundedMatrix<double, TNumNodes, TNumNodes>{ZeroMatrix{TNumNodes, TNumNodes}};
        for (unsigned int GPoint = 0; GPoint < GetGeometry().IntegrationPointsNumber(GetIntegrationMethod()); ++GPoint) {
            BoundedMatrix<double, TDim, TNumNodes> Temp = prod(constitutive_matrix, trans(rShapeFunctionGradients[GPoint]));
            result += prod(rShapeFunctionGradients[GPoint], Temp) * rIntegrationCoefficients[GPoint];
        }

        return result;
    }

    BoundedMatrix<double, TNumNodes, TNumNodes> CalculateCapacityMatrix(const Vector& rIntegrationCoefficients) const
    {
        const auto& r_properties = GetProperties();
        const auto  cWater = r_properties[POROSITY] * r_properties[SATURATION] *
                             r_properties[DENSITY_WATER] * r_properties[SPECIFIC_HEAT_CAPACITY_WATER];
        const auto  cSolid = (1.0 - r_properties[POROSITY]) *
                             r_properties[DENSITY_SOLID] * r_properties[SPECIFIC_HEAT_CAPACITY_SOLID];

        auto result = BoundedMatrix<double, TNumNodes, TNumNodes>{ZeroMatrix{TNumNodes, TNumNodes}};
        const auto& NContainer = GetGeometry().ShapeFunctionsValues(GetIntegrationMethod());
        for (unsigned int GPoint = 0; GPoint < GetGeometry().IntegrationPointsNumber(GetIntegrationMethod()); ++GPoint) {
            const auto N = Vector{row(NContainer, GPoint)};
            result += (cWater + cSolid) * outer_prod(N, N) * rIntegrationCoefficients[GPoint];
        }

        return result;
    }

    array_1d<double, TNumNodes> GetNodalValuesOf(const Variable<double>& rNodalVariable) const
    {
        auto result = array_1d<double, TNumNodes>{};
        const auto& r_geometry = GetGeometry();
        std::transform(r_geometry.begin(), r_geometry.end(), result.begin(),
                       [&rNodalVariable](const auto& node) {
                           return node.FastGetSolutionStepValue(rNodalVariable);
                       });
        return result;
    }

    void VerifyProperty(const Kratos::Variable<double>& rVariable) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(rVariable)) << rVariable.Name()
                                                            << " does not exist in the material properties" << std::endl;
        KRATOS_ERROR_IF(GetProperties()[rVariable] < 0.0) << rVariable.Name() << " has an invalid value at element "
                                                          << Id() << std::endl;
    }

    void CheckDomainSize() const
    {
        constexpr auto min_domain_size = 1.0e-15;
        KRATOS_ERROR_IF(GetGeometry().DomainSize() < min_domain_size)
            << "DomainSize smaller than " << min_domain_size << " for element " << Id() << std::endl;
    }

    void CheckHasSolutionStepsDataFor(const Variable<double>& rVariable) const
    {
        for (const auto& node : GetGeometry()) {
            KRATOS_ERROR_IF_NOT(node.SolutionStepsDataHas(rVariable))
                << "Missing variable " << rVariable.Name() << " on node " << node.Id() << std::endl;
        }
    }

    void CheckHasDofsFor(const Variable<double>& rVariable) const
    {
        for (const auto& node : GetGeometry()) {
            KRATOS_ERROR_IF_NOT(node.HasDofFor(rVariable))
                << "Missing degree of freedom for " << rVariable.Name() << " on node " << node.Id() << std::endl;
        }
    }

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
    }
};

} // namespace Kratos
