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

#include "custom_retention/retention_law_factory.h"
#include "custom_utilities/dof_utilities.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/element.h"
#include "includes/serializer.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) TransientPwLineElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TransientPwLineElement);

    explicit TransientPwLineElement(IndexType NewId = 0) : Element(NewId) {}

    TransientPwLineElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    }

    TransientPwLineElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    Element::Pointer Create(IndexType NewId, const NodesArrayType& rThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return make_intrusive<TransientPwLineElement>(NewId, GetGeometry().Create(rThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return make_intrusive<TransientPwLineElement>(NewId, pGeom, pProperties);
    }

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const override
    {
        rElementalDofList = GetDofs();
    }

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override
    {
        rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
    }

    void CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                              VectorType&        rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        Vector det_J_container;
        GetGeometry().DeterminantOfJacobian(det_J_container, this->GetIntegrationMethod());
        GeometryType::ShapeFunctionsGradientsType dN_dX_container =
            GetGeometry().ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        std::transform(dN_dX_container.begin(), dN_dX_container.end(), det_J_container.begin(),
                       dN_dX_container.begin(), std::divides<>());
        const Matrix& r_N_container = GetGeometry().ShapeFunctionsValues(GetIntegrationMethod());

        const auto integration_coefficients = CalculateIntegrationCoefficients(det_J_container);
        const auto permeability_matrix = CalculatePermeabilityMatrix(dN_dX_container, integration_coefficients);
        const auto compressibility_matrix =
            CalculateCompressibilityMatrix(r_N_container, integration_coefficients);

        const auto fluid_body_vector =
            CalculateFluidBodyVector(r_N_container, dN_dX_container, integration_coefficients);

        AddContributionsToLhsMatrix(rLeftHandSideMatrix, permeability_matrix, compressibility_matrix,
                                    rCurrentProcessInfo[DT_PRESSURE_COEFFICIENT]);
        AddContributionsToRhsVector(rRightHandSideVector, permeability_matrix,
                                    compressibility_matrix, fluid_body_vector);

        KRATOS_CATCH("")
    }

    GeometryData::IntegrationMethod GetIntegrationMethod() const override
    {
        switch (TNumNodes) {
        case 2:
        case 3:
            return GeometryData::IntegrationMethod::GI_GAUSS_2;
        case 4:
            return GeometryData::IntegrationMethod::GI_GAUSS_3;
        case 5:
            return GeometryData::IntegrationMethod::GI_GAUSS_5;
        default:
            KRATOS_ERROR << "Can't return integration method: unexpected number of nodes: " << TNumNodes
                         << std::endl;
        }
    }

    int Check(const ProcessInfo&) const override
    {
        KRATOS_TRY

        CheckDomainSize();
        CheckHasSolutionStepsDataFor(WATER_PRESSURE);
        CheckHasSolutionStepsDataFor(DT_WATER_PRESSURE);
        CheckHasDofsFor(WATER_PRESSURE);
        CheckProperties();
        CheckForNonZeroZCoordinateIn2D();

        KRATOS_CATCH("")

        return 0;
    }

private:
    std::vector<RetentionLaw::Pointer> mRetentionLawVector;

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
                << "Missing degree of freedom for " << rVariable.Name() << " on node " << node.Id()
                << std::endl;
        }
    }

    void CheckProperties() const
    {
        CheckProperty(DENSITY_WATER);
        CheckProperty(DENSITY_SOLID);
        CheckProperty(POROSITY);
        CheckProperty(BULK_MODULUS_SOLID);
        CheckProperty(BULK_MODULUS_FLUID);
        CheckProperty(DYNAMIC_VISCOSITY);
        CheckProperty(BIOT_COEFFICIENT);
        CheckProperty(PERMEABILITY_XX);
    }

    void CheckProperty(const Kratos::Variable<double>& rVariable) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(rVariable))
            << rVariable.Name() << " does not exist in the pressure element's properties" << std::endl;
        KRATOS_ERROR_IF(GetProperties()[rVariable] < 0.0)
            << rVariable.Name() << " has an invalid value at element " << Id() << std::endl;
    }

    void CheckProperty(const Kratos::Variable<std::string>& rVariable, const std::string& rName) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(rVariable))
            << rVariable.Name() << " does not exist in the pressure element's properties" << std::endl;
        KRATOS_ERROR_IF_NOT(GetProperties()[rVariable] == rName)
            << rVariable.Name() << " has a value of (" << GetProperties()[rVariable]
            << ") instead of (" << rName << ") at element " << Id() << std::endl;
    }

    void CheckForNonZeroZCoordinateIn2D() const
    {
        if constexpr (TDim == 2) {
            const auto& r_geometry = GetGeometry();
            auto        pos        = std::find_if(r_geometry.begin(), r_geometry.end(),
                                                  [](const auto& node) { return node.Z() != 0.0; });
            KRATOS_ERROR_IF_NOT(pos == r_geometry.end())
                << " Node with non-zero Z coordinate found. Id: " << pos->Id() << std::endl;
        }
    }

    static void AddContributionsToLhsMatrix(MatrixType& rLeftHandSideMatrix,
                                            const BoundedMatrix<double, TNumNodes, TNumNodes>& rPermeabilityMatrix,
                                            const BoundedMatrix<double, TNumNodes, TNumNodes>& rCompressibilityMatrix,
                                            double DtPressureCoefficient)
    {
        rLeftHandSideMatrix = rPermeabilityMatrix + DtPressureCoefficient * rCompressibilityMatrix;
    }

    void AddContributionsToRhsVector(VectorType& rRightHandSideVector,
                                     const BoundedMatrix<double, TNumNodes, TNumNodes>& rPermeabilityMatrix,
                                     const BoundedMatrix<double, TNumNodes, TNumNodes>& rCompressibilityMatrix,
                                     const array_1d<double, TNumNodes>& rFluidBodyVector) const
    {
        const auto compressibility_vector =
            array_1d<double, TNumNodes>{-prod(rCompressibilityMatrix, GetNodalValuesOf(DT_WATER_PRESSURE))};
        const auto permeability_vector =
            array_1d<double, TNumNodes>{-prod(rPermeabilityMatrix, GetNodalValuesOf(WATER_PRESSURE))};
        rRightHandSideVector = compressibility_vector + permeability_vector + rFluidBodyVector;
    }

    Vector CalculateIntegrationCoefficients(const Vector& rDetJContainer) const
    {
        const auto& r_properties         = GetProperties();
        const auto& r_integration_points = GetGeometry().IntegrationPoints(GetIntegrationMethod());

        auto result = Vector{r_integration_points.size()};
        std::transform(r_integration_points.begin(), r_integration_points.end(), rDetJContainer.begin(),
                       result.begin(), [&r_properties](const auto& rIntegrationPoint, const auto& rDetJ) {
            return rIntegrationPoint.Weight() * rDetJ * r_properties[CROSS_AREA];
        });
        return result;
    }

    BoundedMatrix<double, TNumNodes, TNumNodes> CalculatePermeabilityMatrix(
        const GeometryType::ShapeFunctionsGradientsType& rShapeFunctionGradients,
        const Vector&                                    rIntegrationCoefficients) const
    {
        RetentionLaw::Parameters    RetentionParameters(GetProperties());
        BoundedMatrix<double, 1, 1> constitutive_matrix;
        const auto&                 r_properties = GetProperties();
        GeoElementUtilities::FillPermeabilityMatrix(constitutive_matrix, r_properties);

        auto result = BoundedMatrix<double, TNumNodes, TNumNodes>{ZeroMatrix{TNumNodes, TNumNodes}};
        for (unsigned int integration_point_index = 0;
             integration_point_index < GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
             ++integration_point_index) {
            const double RelativePermeability =
                mRetentionLawVector[integration_point_index]->CalculateRelativePermeability(RetentionParameters);
            double dynamic_viscosity_inverse = 1.0 / r_properties[DYNAMIC_VISCOSITY];
            result += GeoTransportEquationUtilities::CalculatePermeabilityMatrix<TDim, TNumNodes>(
                rShapeFunctionGradients[integration_point_index], dynamic_viscosity_inverse, constitutive_matrix,
                RelativePermeability, rIntegrationCoefficients[integration_point_index]);
        }
        return result;
    }

    BoundedMatrix<double, TNumNodes, TNumNodes> CalculateCompressibilityMatrix(const Matrix& rNContainer,
                                                                               const Vector& rIntegrationCoefficients) const
    {
        const auto&              r_properties = GetProperties();
        RetentionLaw::Parameters parameters(r_properties);
        auto                     retention_law = RetentionLawFactory::Clone(r_properties);

        auto result = BoundedMatrix<double, TNumNodes, TNumNodes>{ZeroMatrix{TNumNodes, TNumNodes}};
        for (unsigned int integration_point_index = 0;
             integration_point_index < GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
             ++integration_point_index) {
            const auto   N                  = Vector{row(rNContainer, integration_point_index)};
            const double BiotModulusInverse = CalculateBiotModulusInverse(integration_point_index);
            result += GeoTransportEquationUtilities::CalculateCompressibilityMatrix<TNumNodes>(
                N, BiotModulusInverse, rIntegrationCoefficients[integration_point_index]);
        }
        return result;
    }

    array_1d<double, TNumNodes> GetNodalValuesOf(const Variable<double>& rNodalVariable) const
    {
        auto        result     = array_1d<double, TNumNodes>{};
        const auto& r_geometry = GetGeometry();
        std::transform(r_geometry.begin(), r_geometry.end(), result.begin(), [&rNodalVariable](const auto& node) {
            return node.FastGetSolutionStepValue(rNodalVariable);
        });
        return result;
    }

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override
    {
        if (const std::size_t number_integration_points =
                GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
            mRetentionLawVector.size() != number_integration_points) {
            mRetentionLawVector.resize(number_integration_points);
        }
        for (unsigned int i = 0; i < mRetentionLawVector.size(); ++i) {
            mRetentionLawVector[i] = RetentionLawFactory::Clone(GetProperties());
            mRetentionLawVector[i]->InitializeMaterial(
                GetProperties(), GetGeometry(),
                row(GetGeometry().ShapeFunctionsValues(GetIntegrationMethod()), i));
        }
    }

    double CalculateBiotModulusInverse(const unsigned int integrationPointIndex) const
    {
        const auto&  r_properties     = GetProperties();
        const double biot_coefficient = r_properties[BIOT_COEFFICIENT];

        double bulk_fluid = TINY;
        if (!r_properties[IGNORE_UNDRAINED]) {
            bulk_fluid = r_properties[BULK_MODULUS_FLUID];
        }
        double result = (biot_coefficient - r_properties[POROSITY]) / r_properties[BULK_MODULUS_SOLID] +
                        r_properties[POROSITY] / bulk_fluid;

        RetentionLaw::Parameters RetentionParameters(GetProperties());
        const double             degree_of_saturation =
            mRetentionLawVector[integrationPointIndex]->CalculateSaturation(RetentionParameters);
        const double derivative_of_saturation =
            mRetentionLawVector[integrationPointIndex]->CalculateDerivativeOfSaturation(RetentionParameters);

        result *= degree_of_saturation;
        result -= derivative_of_saturation * r_properties[POROSITY];
        return result;
    }

    void InitializeSolutionStep(const ProcessInfo&) override
    {
        RetentionLaw::Parameters RetentionParameters(this->GetProperties());
        for (auto retention_law : mRetentionLawVector) {
            retention_law->InitializeSolutionStep(RetentionParameters);
        }
    }

    void FinalizeSolutionStep(const ProcessInfo&) override
    {
        RetentionLaw::Parameters RetentionParameters(this->GetProperties());
        for (auto retention_law : mRetentionLawVector) {
            retention_law->FinalizeSolutionStep(RetentionParameters);
        }
    }

    array_1d<double, TNumNodes> CalculateFluidBodyVector(const Matrix& rNContainer,
                                                         const GeometryType::ShapeFunctionsGradientsType& rShapeFunctionGradients,
                                                         const Vector& rIntegrationCoefficients) const
    {
        const std::size_t number_integration_points =
            GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
        GeometryType::JacobiansType J_container;
        J_container.resize(number_integration_points, false);
        for (std::size_t i = 0; i < number_integration_points; ++i) {
            J_container[i].resize(GetGeometry().WorkingSpaceDimension(),
                                  GetGeometry().LocalSpaceDimension(), false);
        }
        GetGeometry().Jacobian(J_container, this->GetIntegrationMethod());

        const auto&                 r_properties = GetProperties();
        BoundedMatrix<double, 1, 1> constitutive_matrix;
        GeoElementUtilities::FillPermeabilityMatrix(constitutive_matrix, r_properties);

        RetentionLaw::Parameters RetentionParameters(GetProperties());

        array_1d<double, TNumNodes * TDim> volume_acceleration;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(
            volume_acceleration, GetGeometry(), VOLUME_ACCELERATION);
        array_1d<double, TDim> body_acceleration;

        array_1d<double, TNumNodes> fluid_body_vector = ZeroVector(TNumNodes);
        for (unsigned int integration_point_index = 0;
             integration_point_index < GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
             ++integration_point_index) {
            GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
                body_acceleration, rNContainer, volume_acceleration, integration_point_index);

            array_1d<double, TDim> tangent_vector = column(J_container[integration_point_index], 0);
            tangent_vector /= norm_2(tangent_vector);

            array_1d<double, 1> projected_gravity = ZeroVector(1);
            projected_gravity(0) = MathUtils<double>::Dot(tangent_vector, body_acceleration);
            const auto N         = Vector{row(rNContainer, integration_point_index)};
            double     RelativePermeability =
                mRetentionLawVector[integration_point_index]->CalculateRelativePermeability(RetentionParameters);
            fluid_body_vector +=
                r_properties[DENSITY_WATER] * RelativePermeability *
                prod(prod(rShapeFunctionGradients[integration_point_index], constitutive_matrix), projected_gravity) *
                rIntegrationCoefficients[integration_point_index] / r_properties[DYNAMIC_VISCOSITY];
        }
        return fluid_body_vector;
    }

    [[nodiscard]] DofsVectorType GetDofs() const
    {
        return Geo::DofUtilities::ExtractDofsFromNodes(GetGeometry(), WATER_PRESSURE);
    }

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
        rSerializer.save("RetentionlawVector", mRetentionLawVector);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
    }
};

} // namespace Kratos
