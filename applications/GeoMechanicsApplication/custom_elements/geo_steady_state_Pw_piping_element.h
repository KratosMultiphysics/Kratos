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

#include "includes/element.h"
#include "includes/serializer.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoSteadyStatePwPipingElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoSteadyStatePwPipingElement);

    explicit GeoSteadyStatePwPipingElement(IndexType NewId = 0);

    GeoSteadyStatePwPipingElement(IndexType NewId, const GeometryType::Pointer& rpGeometry);

    GeoSteadyStatePwPipingElement(IndexType                      NewId,
                                  const GeometryType::Pointer&   rpGeometry,
                                  const PropertiesType::Pointer& rpProperties);

    Element::Pointer Create(IndexType               NewId,
                            const NodesArrayType&   rThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType               NewId,
                            GeometryType::Pointer   pGeometry,
                            PropertiesType::Pointer pProperties) const override;

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override;

    void CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                              VectorType&        rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    int Check(const ProcessInfo&) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) final;

    double CalculateEquilibriumPipeHeight(const PropertiesType& rProp, const GeometryType& rGeom, double);

    void CalculateOnIntegrationPoints(const Variable<bool>& rVariable,
                                      std::vector<bool>&    rValues,
                                      const ProcessInfo&    rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>&    rValues,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                      std::vector<array_1d<double, 3>>&    rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>&    rValues,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    using Element::CalculateOnIntegrationPoints;

    std::string Info() const override;

private:
    static double CalculateLength(const GeometryType& Geom);

    double CalculateHeadGradient(const PropertiesType& rProp, const GeometryType& rGeom);

    static void AddContributionsToLhsMatrix(MatrixType& rLeftHandSideMatrix,
                                            const BoundedMatrix<double, TNumNodes, TNumNodes>& rPermeabilityMatrix);

    void AddContributionsToRhsVector(VectorType& rRightHandSideVector,
                                     const BoundedMatrix<double, TNumNodes, TNumNodes>& rPermeabilityMatrix,
                                     const array_1d<double, TNumNodes>& rFluidBodyVector) const;

    Vector CalculateIntegrationCoefficients(const Vector& rDetJContainer) const;

    Matrix FillPermeabilityMatrix(double pipe_height) const;

    BoundedMatrix<double, TNumNodes, TNumNodes> CalculatePermeabilityMatrix(
        const GeometryType::ShapeFunctionsGradientsType& rShapeFunctionGradients,
        const Vector&                                    rIntegrationCoefficients) const;

    array_1d<double, TNumNodes> CalculateFluidBodyVector(const Matrix& rNContainer,
                                                         const GeometryType::ShapeFunctionsGradientsType& rShapeFunctionGradients,
                                                         const Vector& rIntegrationCoefficients) const;

    [[nodiscard]] DofsVectorType GetDofs() const;

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    bool mIsInitialized = false;
};
} // namespace Kratos
