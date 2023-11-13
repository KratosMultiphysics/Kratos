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

    TransientThermalElement& operator=(TransientThermalElement const& rOther) = delete;

    TransientThermalElement(TransientThermalElement const& rOther) = delete;

    ~TransientThermalElement() override;

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
                    const ProcessInfo& rCurrentProcessInfo) const override;

    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;
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

    Vector CalculateIntegrationCoefficients(const Vector& detJContainer) const;

    BoundedMatrix<double, TNumNodes, TNumNodes> CalculateConductivityMatrix(const GeometryType::ShapeFunctionsGradientsType& rShapeFunctionGradients,
                                                                            const Vector& rIntegrationCoefficients,
                                                                            const ProcessInfo& rCurrentProcessInfo) const;

    BoundedMatrix<double, TNumNodes, TNumNodes> CalculateCapacityMatrix(const Vector& rIntegrationCoefficients) const;

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

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    void VerifyProperty(Kratos::Variable<double>& rVariable) const;
    void CheckDomainSize() const;
    void CheckSolutionStepsData(int rId, Kratos::Variable<double>& rVariable) const;

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
