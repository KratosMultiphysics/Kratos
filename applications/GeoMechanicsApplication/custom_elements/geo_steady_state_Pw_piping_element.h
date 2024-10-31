// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//

#pragma once

#include "custom_utilities/dof_utilities.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/cfd_variables.h"
#include "includes/element.h"
#include "includes/serializer.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoSteadyStatePwPipingElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoSteadyStatePwPipingElement);

    explicit GeoSteadyStatePwPipingElement(IndexType NewId = 0) : Element(NewId) {}

    GeoSteadyStatePwPipingElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    }

    GeoSteadyStatePwPipingElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    Element::Pointer Create(IndexType NewId, const NodesArrayType& rThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return this->Create(NewId, GetGeometry().Create(rThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) const override
    {
        return make_intrusive<GeoSteadyStatePwPipingElement>(NewId, pGeometry, pProperties);
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
        AddContributionsToLhsMatrix(rLeftHandSideMatrix, permeability_matrix);

        const auto fluid_body_vector =
            CalculateFluidBodyVector(r_N_container, dN_dX_container, integration_coefficients);
        AddContributionsToRhsVector(rRightHandSideVector, permeability_matrix, fluid_body_vector);

        KRATOS_CATCH("")
    }

    GeometryData::IntegrationMethod GetIntegrationMethod() const override
    {
        KRATOS_ERROR_IF(TNumNodes != 2)
            << "Can't return integration method: unexpected number of nodes: " << TNumNodes << std::endl;
        return GeometryData::IntegrationMethod::GI_GAUSS_2;
    }

    int Check(const ProcessInfo&) const override
    {
        KRATOS_TRY

        CheckDomainSize();
        CheckHasSolutionStepsDataFor(WATER_PRESSURE);
        CheckHasDofsFor(WATER_PRESSURE);
        CheckProperties();
        // conditional on model dimension
        CheckForNonZeroZCoordinateIn2D();

        KRATOS_CATCH("")

        return 0;
    }

private:
    void CheckDomainSize() const
    {
        constexpr auto min_domain_size = 1.0e-15;
        KRATOS_ERROR_IF(GetGeometry().DomainSize() < min_domain_size)
            << "DomainSize (" << GetGeometry().DomainSize() << ") is smaller than "
            << min_domain_size << " for element " << Id() << std::endl;
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
        // typical material parameters check, this should be in the check of the constitutive law.
        // possibly check PIPE_HEIGHT, CROSS_SECTION == 1.0
        CheckProperty(DENSITY_WATER);
        CheckProperty(DYNAMIC_VISCOSITY);
        CheckProperty(PIPE_HEIGHT);
    }

    void CheckProperty(const Kratos::Variable<double>& rVariable) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(rVariable))
            << rVariable.Name() << " does not exist in the properties of element " << Id() << std::endl;
        KRATOS_ERROR_IF(GetProperties()[rVariable] < 0.0)
            << rVariable.Name() << " (" << GetProperties()[rVariable]
            << ") is not in the range [0,-> at element " << Id() << std::endl;
    }

    void CheckForNonZeroZCoordinateIn2D() const
    {
        if constexpr (TDim == 2) {
            const auto& r_geometry = GetGeometry();
            auto        pos        = std::find_if(r_geometry.begin(), r_geometry.end(),
                                                  [](const auto& node) { return node.Z() != 0.0; });
            KRATOS_ERROR_IF_NOT(pos == r_geometry.end())
                << "Node with non-zero Z coordinate found. Id: " << pos->Id() << std::endl;
        }
    }

    static void AddContributionsToLhsMatrix(MatrixType& rLeftHandSideMatrix,
                                            const BoundedMatrix<double, TNumNodes, TNumNodes>& rPermeabilityMatrix)
    {
        rLeftHandSideMatrix = rPermeabilityMatrix;
    }

    void AddContributionsToRhsVector(VectorType& rRightHandSideVector,
                                     const BoundedMatrix<double, TNumNodes, TNumNodes>& rPermeabilityMatrix,
                                     const array_1d<double, TNumNodes>& rFluidBodyVector) const
    {
        const auto permeability_vector =
            array_1d<double, TNumNodes>{-prod(rPermeabilityMatrix, GetNodalValuesOf(WATER_PRESSURE))};
        rRightHandSideVector = permeability_vector + rFluidBodyVector;
    }

    Vector CalculateIntegrationCoefficients(const Vector& rDetJContainer) const
    {
        const auto& r_integration_points = GetGeometry().IntegrationPoints(GetIntegrationMethod());

        auto result = Vector{r_integration_points.size()};
        // all governed by PIPE_HEIGHT and element length so without CROSS_AREA
        std::transform(r_integration_points.begin(), r_integration_points.end(), rDetJContainer.begin(),
                       result.begin(), [](const auto& rIntegrationPoint, const auto& rDetJ) {
            return rIntegrationPoint.Weight() * rDetJ;
        });
        return result;
    }

    Matrix FillPermeabilityMatrix(double pipe_height) const
    {
        return ScalarMatrix{1, 1, std::pow(pipe_height, 3) / 12.0};
    }

    BoundedMatrix<double, TNumNodes, TNumNodes> CalculatePermeabilityMatrix(
        const GeometryType::ShapeFunctionsGradientsType& rShapeFunctionGradients,
        const Vector&                                    rIntegrationCoefficients) const
    {
        const auto&  r_properties              = GetProperties();
        const double dynamic_viscosity_inverse = 1.0 / r_properties[DYNAMIC_VISCOSITY];

        auto constitutive_matrix = FillPermeabilityMatrix(r_properties[PIPE_HEIGHT]);

        auto result = BoundedMatrix<double, TNumNodes, TNumNodes>{ZeroMatrix{TNumNodes, TNumNodes}};
        for (unsigned int integration_point_index = 0;
             integration_point_index < GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
             ++integration_point_index) {
            result += GeoTransportEquationUtilities::CalculatePermeabilityMatrix<TDim, TNumNodes>(
                rShapeFunctionGradients[integration_point_index], dynamic_viscosity_inverse,
                constitutive_matrix, 1.0, rIntegrationCoefficients[integration_point_index]);
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

        const auto& r_properties        = GetProperties();
        auto        constitutive_matrix = FillPermeabilityMatrix(r_properties[PIPE_HEIGHT]);

        array_1d<double, TNumNodes * TDim> volume_acceleration;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(
            volume_acceleration, GetGeometry(), VOLUME_ACCELERATION);

        array_1d<double, TNumNodes> fluid_body_vector = ZeroVector(TNumNodes);
        for (unsigned int integration_point_index = 0;
             integration_point_index < GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
             ++integration_point_index) {
            array_1d<double, TDim> body_acceleration;
            GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
                body_acceleration, rNContainer, volume_acceleration, integration_point_index);

            array_1d<double, TDim> tangent_vector = column(J_container[integration_point_index], 0);
            tangent_vector /= norm_2(tangent_vector);

            array_1d<double, 1> projected_gravity = ZeroVector(1);
            projected_gravity(0) = MathUtils<double>::Dot(tangent_vector, body_acceleration);
            const auto N         = Vector{row(rNContainer, integration_point_index)};
            fluid_body_vector +=
                r_properties[DENSITY_WATER] *
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
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
    }
};

} // namespace Kratos
