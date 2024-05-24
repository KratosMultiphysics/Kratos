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

#include "custom_constitutive/thermal_dispersion_law.h"
#include "custom_utilities/thermal_utilities.hpp"
#include "custom_constitutive/thermal_filter_law.h"
#include "custom_retention/retention_law_factory.h"
#include "custom_utilities/dof_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/element.h"
#include "includes/serializer.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) TransientThermalElement : public Element
{
private:
    bool mIsPressureCoupled = false;
    bool mUpdateDensityViscosity = true;
    Vector mWaterDensity;
    Vector mWaterViscosity;

public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TransientThermalElement);

    explicit TransientThermalElement(IndexType NewId = 0) : Element(NewId) {}

    TransientThermalElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    }

    TransientThermalElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    Element::Pointer Create(IndexType NewId, const NodesArrayType& rThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return make_intrusive<TransientThermalElement>(NewId, GetGeometry().Create(rThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return make_intrusive<TransientThermalElement>(NewId, pGeom, pProperties);
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

        GeometryType::ShapeFunctionsGradientsType dN_dX_container;
        Vector                                    det_J_container;

        // ShapreFunctionsIntegrationsPointsGradients does not allow for the line element in 2D/3D configuration
        // and will produce errors. To circumvent this, the dN_dX_container is separately computed with correct
        // dimensions for the line element.
        if (GetGeometry().LocalSpaceDimension() == 1) {
            GetGeometry().DeterminantOfJacobian(det_J_container, this->GetIntegrationMethod());
            dN_dX_container = GetGeometry().ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
            std::transform(dN_dX_container.begin(), dN_dX_container.end(), det_J_container.begin(),
                           dN_dX_container.begin(), std::divides<>());
        } 
        else {
            GetGeometry().ShapeFunctionsIntegrationPointsGradients(dN_dX_container, det_J_container,
                                                                   GetIntegrationMethod());
        }
        
        const auto integration_coefficients = CalculateIntegrationCoefficients(det_J_container);

        BoundedMatrix<double, TNumNodes, TNumNodes> convection_matrix = ZeroMatrix(TNumNodes, TNumNodes);
        if (mIsPressureCoupled) {
            if (mUpdateDensityViscosity) {
                this->UpdateMaterialProperties();
            }
            convection_matrix = CalculateConvectionMatrix(dN_dX_container, integration_coefficients);
        }

        const auto conductivity_matrix =
            CalculateConductivityMatrix(dN_dX_container, integration_coefficients, rCurrentProcessInfo);
        const auto capacity_matrix = CalculateCapacityMatrix(integration_coefficients, rCurrentProcessInfo);

        AddContributionsToLhsMatrix(rLeftHandSideMatrix, conductivity_matrix, capacity_matrix, convection_matrix,
                                    rCurrentProcessInfo[DT_TEMPERATURE_COEFFICIENT]);
        AddContributionsToRhsVector(rRightHandSideVector, conductivity_matrix, capacity_matrix, convection_matrix);

        KRATOS_CATCH("")
    }

    GeometryData::IntegrationMethod GetIntegrationMethod() const override
    {
        if (GetGeometry().LocalSpaceDimension() == 1)
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

        switch (TNumNodes) {
        case 3:
        case 4:
        case 6:
        case 8:
        case 9:
        case 20:
        case 27:
            return GeometryData::IntegrationMethod::GI_GAUSS_2;
        case 10:
            return GeometryData::IntegrationMethod::GI_GAUSS_4;
        case 15:
            return GeometryData::IntegrationMethod::GI_GAUSS_5;
        default:
            KRATOS_ERROR << "Can't return integration method: unexpected number of nodes: "
                         << TNumNodes << std::endl;
        }
    }

    int Check(const ProcessInfo&) const override
    {
        KRATOS_TRY

        CheckDomainSize();
        CheckHasSolutionStepsDataFor(TEMPERATURE);
        CheckHasSolutionStepsDataFor(DT_TEMPERATURE);
        CheckHasDofsFor(TEMPERATURE);
        CheckProperties();
        CheckForNonZeroZCoordinateIn2D();

        KRATOS_CATCH("")

        return 0;
    }

private:
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
        CheckProperty(POROSITY);
        CheckProperty(RETENTION_LAW, "SaturatedLaw");
        CheckProperty(SATURATED_SATURATION);
        CheckProperty(DENSITY_SOLID);
        CheckProperty(SPECIFIC_HEAT_CAPACITY_WATER);
        CheckProperty(SPECIFIC_HEAT_CAPACITY_SOLID);
        CheckProperty(THERMAL_CONDUCTIVITY_WATER);
        CheckProperty(THERMAL_CONDUCTIVITY_SOLID_XX);
        CheckProperty(THERMAL_CONDUCTIVITY_SOLID_YY);
        CheckProperty(THERMAL_CONDUCTIVITY_SOLID_XY);

        if constexpr (TDim == 3) {
            CheckProperty(THERMAL_CONDUCTIVITY_SOLID_ZZ);
            CheckProperty(THERMAL_CONDUCTIVITY_SOLID_YZ);
            CheckProperty(THERMAL_CONDUCTIVITY_SOLID_XZ);
        }
    }

    void CheckProperty(const Kratos::Variable<double>& rVariable) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(rVariable))
            << rVariable.Name() << " does not exist in the thermal element's properties" << std::endl;
        KRATOS_ERROR_IF(GetProperties()[rVariable] < 0.0)
            << rVariable.Name() << " has an invalid value at element " << Id() << std::endl;
    }

    void CheckProperty(const Kratos::Variable<std::string>& rVariable, const std::string& rName) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(rVariable))
            << rVariable.Name() << " does not exist in the thermal element's properties" << std::endl;
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
                                            const BoundedMatrix<double, TNumNodes, TNumNodes>& rConductivityMatrix,
                                            const BoundedMatrix<double, TNumNodes, TNumNodes>& rCapacityMatrix,
                                            const BoundedMatrix<double, TNumNodes, TNumNodes>& rConvectionMatrix,
                                            double DtTemperatureCoefficient)
    {
        rLeftHandSideMatrix = rConductivityMatrix;
        rLeftHandSideMatrix += (DtTemperatureCoefficient * rCapacityMatrix);
        rLeftHandSideMatrix += rConvectionMatrix;
    }

    void AddContributionsToRhsVector(VectorType& rRightHandSideVector,
                                     const BoundedMatrix<double, TNumNodes, TNumNodes>& rConductivityMatrix,
                                     const BoundedMatrix<double, TNumNodes, TNumNodes>& rCapacityMatrix,
                                     const BoundedMatrix<double, TNumNodes, TNumNodes>& rConvectionMatrix) const
    {
        const auto capacity_vector =
            array_1d<double, TNumNodes>{-prod(rCapacityMatrix, GetNodalValuesOf(DT_TEMPERATURE))};
        rRightHandSideVector = capacity_vector;
        const auto conductivity_vector =
            array_1d<double, TNumNodes>{-prod(rConductivityMatrix, GetNodalValuesOf(TEMPERATURE))};
        rRightHandSideVector += conductivity_vector;

        const auto convection_vector =
            array_1d<double, TNumNodes>{-prod(rConvectionMatrix, GetNodalValuesOf(TEMPERATURE))};
        rRightHandSideVector += convection_vector;
    }

    Vector CalculateIntegrationCoefficients(const Vector& rDetJContainer) const
    {
        const auto& r_properties = GetProperties();
        const auto& r_integration_points = GetGeometry().IntegrationPoints(GetIntegrationMethod());

        auto result = Vector{r_integration_points.size()};
        for (unsigned int integration_point_index = 0;
             integration_point_index < r_integration_points.size(); ++integration_point_index) {
            result[integration_point_index] = r_integration_points[integration_point_index].Weight() *
                                              rDetJContainer[integration_point_index];
            if (GetGeometry().LocalSpaceDimension() == 1) {
                result[integration_point_index] *= r_properties[CROSS_AREA];
            }
        }

        return result;
    }

    BoundedMatrix<double, TNumNodes, TNumNodes> CalculateConductivityMatrix(
        const GeometryType::ShapeFunctionsGradientsType& rShapeFunctionGradients,
        const Vector&                                    rIntegrationCoefficients,
        const ProcessInfo&                               rCurrentProcessInfo) const
    {
        Vector discharge_vector = ZeroVector(TDim);

        const auto law = CreateThermalLaw();
        auto constitutive_matrix = law->CalculateThermalDispersionMatrix(GetProperties(), rCurrentProcessInfo);

        auto result = BoundedMatrix<double, TNumNodes, TNumNodes>{ZeroMatrix{TNumNodes, TNumNodes}};
        for (unsigned int integration_point_index = 0;
             integration_point_index < GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
             ++integration_point_index) {

             if (mIsPressureCoupled) {
                 discharge_vector = this->CalculateDischargeVector(rShapeFunctionGradients, integration_point_index);
                 constitutive_matrix = law->CalculateThermalDispersionMatrix(GetProperties(), rCurrentProcessInfo, discharge_vector);
             }
 
            BoundedMatrix<double, TDim, TNumNodes> Temp =
                prod(constitutive_matrix, trans(rShapeFunctionGradients[integration_point_index]));
            result += prod(rShapeFunctionGradients[integration_point_index], Temp) *
                      rIntegrationCoefficients[integration_point_index];
        }

        return result;
    }

    BoundedMatrix<double, TNumNodes, TNumNodes> CalculateCapacityMatrix(const Vector& rIntegrationCoefficients,
                                                                        const ProcessInfo& rCurrentProcessInfo) const
    {
        const auto&              r_properties = GetProperties();
        RetentionLaw::Parameters parameters(r_properties, rCurrentProcessInfo);
        auto                     retention_law = RetentionLawFactory::Clone(r_properties);
        const double             saturation    = retention_law->CalculateSaturation(parameters);

        const auto c_solid = (1.0 - r_properties[POROSITY]) * r_properties[DENSITY_SOLID] *
                             r_properties[SPECIFIC_HEAT_CAPACITY_SOLID];

        auto result = BoundedMatrix<double, TNumNodes, TNumNodes>{ZeroMatrix{TNumNodes, TNumNodes}};
        const auto& r_N_container = GetGeometry().ShapeFunctionsValues(GetIntegrationMethod());
        for (unsigned int integration_point_index = 0;
             integration_point_index < GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
             ++integration_point_index) {
            const auto N = Vector{row(r_N_container, integration_point_index)};

            const auto c_water = r_properties[POROSITY] * saturation * mWaterDensity[integration_point_index] *
                                 r_properties[SPECIFIC_HEAT_CAPACITY_WATER];

            result += (c_water + c_solid) * outer_prod(N, N) * rIntegrationCoefficients[integration_point_index];
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

    // ============================================================================================
    // ============================================================================================
    BoundedMatrix<double, TNumNodes, TNumNodes> CalculateConvectionMatrix(
        const GeometryType::ShapeFunctionsGradientsType& rShapeFunctionGradients,
        const Vector& rIntegrationCoefficients) const
    {
        const auto& r_properties = GetProperties();
        auto result = BoundedMatrix<double, TNumNodes, TNumNodes>{ZeroMatrix{TNumNodes, TNumNodes}};
        const auto& r_N_container = GetGeometry().ShapeFunctionsValues(GetIntegrationMethod());
        for (unsigned int integration_point_index = 0;
             integration_point_index < GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
             ++integration_point_index) {

            auto discharge_vector = this->CalculateDischargeVector(rShapeFunctionGradients, integration_point_index);

            array_1d<double, TNumNodes> Temp =
                prod(rShapeFunctionGradients[integration_point_index], discharge_vector);
            const auto N = Vector{row(r_N_container, integration_point_index)};
            result += outer_prod(N, Temp) * rIntegrationCoefficients[integration_point_index] *
                      mWaterDensity[integration_point_index] * r_properties[SPECIFIC_HEAT_CAPACITY_WATER];
        }
        return result;
    }

    // ============================================================================================
    // ============================================================================================
    Vector CalculateDischargeVector(const GeometryType::ShapeFunctionsGradientsType& rShapeFunctionGradients,
                                                    unsigned int integration_point_index) const
    {
        const GeometryType& rGeom = this->GetGeometry();
        array_1d<double, TNumNodes> pressure_vector;
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            pressure_vector[i] = -rGeom[i].FastGetSolutionStepValue(WATER_PRESSURE);
        }
        array_1d<double, TDim> pressure_gradient = prod(trans(rShapeFunctionGradients[integration_point_index]), pressure_vector);

        array_1d<double, TDim> weight_vector = rGeom[0].FastGetSolutionStepValue(VOLUME_ACCELERATION) *
                                               mWaterDensity[integration_point_index];
        pressure_gradient = pressure_gradient - weight_vector;

        auto permeability_matrix = this->CalculatePermeabilityMatrix();
        array_1d<double, TDim> result = -prod(permeability_matrix, pressure_gradient) / mWaterViscosity[integration_point_index];

        return result;
    }

    // ============================================================================================
    // ============================================================================================
    BoundedMatrix<double, TDim, TDim> CalculatePermeabilityMatrix() const
    {
        const PropertiesType&             rProp = this->GetProperties();
        BoundedMatrix<double, TDim, TDim> C;

        C(0, 0) = rProp[PERMEABILITY_XX];

        if (TDim == 2) {
            C(0, 1) = rProp[PERMEABILITY_XY];
            C(1, 0) = rProp[PERMEABILITY_XY];
            C(1, 1) = rProp[PERMEABILITY_YY];
        }
        
        if (TDim == 3) {
            C(2, 2) = rProp[PERMEABILITY_ZZ];
            C(1, 2) = rProp[PERMEABILITY_YZ];
            C(2, 1) = rProp[PERMEABILITY_YZ];
            C(0, 2) = rProp[PERMEABILITY_ZX];
            C(2, 0) = rProp[PERMEABILITY_ZX];
        }
        return C;
    }

    // ============================================================================================
    // ============================================================================================
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const GeometryType& rGeom     = this->GetGeometry();
        const unsigned int NumGPoints = rGeom.IntegrationPointsNumber(this->GetIntegrationMethod());
        const PropertiesType& rProp   = this->GetProperties();

        mIsPressureCoupled      = rGeom[0].SolutionStepsDataHas(WATER_PRESSURE);
        mUpdateDensityViscosity = rCurrentProcessInfo[UPDATE_DENSITY_VISCOSITY];
        mWaterDensity           = ZeroVector(NumGPoints);
        mWaterViscosity         = ZeroVector(NumGPoints);

        if (!mIsPressureCoupled || !mUpdateDensityViscosity) {
            for (unsigned int integration_point_index = 0;
                 integration_point_index < GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
                 ++integration_point_index) {
                 mWaterDensity[integration_point_index]   = rProp[DENSITY_WATER];
                 mWaterViscosity[integration_point_index] = rProp[DYNAMIC_VISCOSITY];
            }
        }

        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    void UpdateMaterialProperties() 
    {
        const GeometryType& rGeom  = this->GetGeometry();
        const auto& r_N_container = GetGeometry().ShapeFunctionsValues(GetIntegrationMethod());

        for (unsigned int integration_point_index = 0;
             integration_point_index < GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
             ++integration_point_index) {

            const auto N = Vector{row(r_N_container, integration_point_index)};

            mWaterDensity[integration_point_index] =
                ThermalUtilities::CalculateWaterDensityOnIntegrationPoints(N, rGeom);

            mWaterViscosity[integration_point_index] =
                ThermalUtilities::CalculateWaterViscosityOnIntegrationPoints(N, rGeom);
        }
    }


    std::unique_ptr<GeoThermalLaw> CreateThermalLaw() const
    {
        const std::size_t number_of_dimensions = GetGeometry().LocalSpaceDimension();
        std::unique_ptr<GeoThermalLaw> law;
        
        if (GetProperties().Has(THERMAL_LAW_NAME)) {
            const std::string& rThermalLawName = GetProperties()[THERMAL_LAW_NAME];
            if (rThermalLawName == "GeoThermalDispersionLaw") {
                law = std::make_unique<GeoThermalDispersionLaw>(number_of_dimensions);
            } 
            else if (rThermalLawName == "GeoThermalFilterLaw") {
                law = std::make_unique<GeoThermalFilterLaw>();
            } 
            else {
                KRATOS_ERROR << "Undefined THERMAL_LAW_NAME! " << rThermalLawName << std::endl;
            }
        } 
        else {
            law = std::make_unique<GeoThermalDispersionLaw>(number_of_dimensions);
        }
        return law;
    }

    [[nodiscard]] DofsVectorType GetDofs() const
    {
        return Geo::DofUtilities::ExtractDofsFromNodes(GetGeometry(), TEMPERATURE);
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
