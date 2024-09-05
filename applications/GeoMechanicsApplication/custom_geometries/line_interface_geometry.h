// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//                   Anne van de Graaf
//

#pragma once

#include "geometries/geometry.h"
#include "includes/node.h"

namespace Kratos
{

template <typename MidGeometryType>
class LineInterfaceGeometry : public Geometry<Node>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineInterfaceGeometry);

    using BaseType = Geometry<Node>;

    LineInterfaceGeometry() = default;

    explicit LineInterfaceGeometry(const PointsArrayType& rThisPoints)
        : LineInterfaceGeometry(0, rThisPoints)
    {
    }

    LineInterfaceGeometry(IndexType NewGeometryId, const PointsArrayType& rThisPoints)
        : BaseType(NewGeometryId, rThisPoints)
    {
        KRATOS_ERROR_IF_NOT((rThisPoints.size() == 4) || (rThisPoints.size() == 6))
            << "Number of nodes must be 2+2 or 3+3\n";

        mMidLineGeometry = std::make_unique<MidGeometryType>(CreatePointsOfMidLine());
        this->SetGeometryData(&mMidLineGeometry->GetGeometryData());
    }

    [[nodiscard]] BaseType::Pointer Create(const PointsArrayType& rThisPoints) const override
    {
        constexpr auto id = IndexType{0};
        return Create(id, rThisPoints);
    }

    [[nodiscard]] BaseType::Pointer Create(const IndexType NewGeometryId, const PointsArrayType& rThisPoints) const override
    {
        return std::make_shared<LineInterfaceGeometry>(NewGeometryId, rThisPoints);
    }

    [[nodiscard]] double ShapeFunctionValue(IndexType ShapeFunctionIndex,
                                            const CoordinatesArrayType& rLocalCoordinate) const override
    {
        return mMidLineGeometry->ShapeFunctionValue(ShapeFunctionIndex, rLocalCoordinate);
    }

    Vector& ShapeFunctionsValues(Vector& rResult, const CoordinatesArrayType& rLocalCoordinate) const override
    {
        return mMidLineGeometry->ShapeFunctionsValues(rResult, rLocalCoordinate);
    }

    Matrix& ShapeFunctionsLocalGradients(Matrix& rResult, const CoordinatesArrayType& rLocalCoordinate) const override
    {
        return mMidLineGeometry->ShapeFunctionsLocalGradients(rResult, rLocalCoordinate);
    }

    Matrix& Jacobian(Matrix& rResult, const CoordinatesArrayType& rLocalCoordinate) const override
    {
        return mMidLineGeometry->Jacobian(rResult, rLocalCoordinate);
    }

    [[nodiscard]] double DeterminantOfJacobian(const CoordinatesArrayType& rLocalCoordinate) const override
    {
        return mMidLineGeometry->DeterminantOfJacobian(rLocalCoordinate);
    }

    Matrix& InverseOfJacobian(Matrix& rResult, const CoordinatesArrayType& rLocalCoordinate) const override
    {
        KRATOS_ERROR << "Inverse of Jacobian is not implemented for the line interface geometry\n";
    }

    [[nodiscard]] double Length() const override { return mMidLineGeometry->Length(); }

    [[nodiscard]] double DomainSize() const override { return mMidLineGeometry->DomainSize(); }

    [[nodiscard]] std::string Info() const override
    {
        return "An interface geometry consisting of two sub-geometries with Info: " +
               mMidLineGeometry->Info();
    }

    CoordinatesArrayType& PointLocalCoordinates(CoordinatesArrayType& rResult,
                                                const CoordinatesArrayType& rGlobalCoordinate) const override
    {
        return mMidLineGeometry->PointLocalCoordinates(rResult, rGlobalCoordinate);
    }

    Matrix& PointsLocalCoordinates(Matrix& rResult) const override
    {
        return mMidLineGeometry->PointsLocalCoordinates(rResult);
    }

    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

    void PrintData(std::ostream& rOStream) const override { mMidLineGeometry->PrintData(rOStream); }

    [[nodiscard]] static std::string IntegrationSchemeFunctionalityNotImplementedMessage()
    {
        return "This Geometry type does not support functionality related to integration "
               "schemes.\n";
    }

    array_1d<double, 3> Normal(IndexType IntegrationPointIndex) const override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    array_1d<double, 3> Normal(IndexType IntegrationPointIndex, IntegrationMethod ThisMethod) const override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    array_1d<double, 3> UnitNormal(IndexType IntegrationPointIndex) const override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    array_1d<double, 3> UnitNormal(IndexType IntegrationPointIndex, IntegrationMethod ThisMethod) const override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    JacobiansType& Jacobian(JacobiansType& rResult, IntegrationMethod ThisMethod) const override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    JacobiansType& Jacobian(JacobiansType& rResult, IntegrationMethod ThisMethod, Matrix& DeltaPosition) const override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    Matrix& Jacobian(Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod) const override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    Matrix& Jacobian(Matrix&           rResult,
                     IndexType         IntegrationPointIndex,
                     IntegrationMethod ThisMethod,
                     const Matrix&     rDeltaPosition) const override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    Vector& DeterminantOfJacobian(Vector& rResult, IntegrationMethod ThisMethod) const override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    double DeterminantOfJacobian(IndexType IntegrationPointIndex, IntegrationMethod ThisMethod) const override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    JacobiansType& InverseOfJacobian(JacobiansType& rResult, IntegrationMethod ThisMethod) const override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    Matrix& InverseOfJacobian(Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod) const override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    void ShapeFunctionsIntegrationPointsGradients(ShapeFunctionsGradientsType& rResult,
                                                  IntegrationMethod ThisMethod) const override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    void ShapeFunctionsIntegrationPointsGradients(ShapeFunctionsGradientsType& rResult,
                                                  Vector&           rDeterminantsOfJacobian,
                                                  IntegrationMethod ThisMethod) const override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    void ShapeFunctionsIntegrationPointsGradients(ShapeFunctionsGradientsType& rResult,
                                                  Vector&           rDeterminantsOfJacobian,
                                                  IntegrationMethod ThisMethod,
                                                  Matrix& ShapeFunctionsIntegrationPointsValues) const override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    IntegrationInfo GetDefaultIntegrationInfo() const override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    void CreateIntegrationPoints(IntegrationPointsArrayType& rIntegrationPoints,
                                 IntegrationInfo&            rIntegrationInfo) const override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    void CreateQuadraturePointGeometries(GeometriesArrayType& rResultGeometries,
                                         IndexType            NumberOfShapeFunctionDerivatives,
                                         const IntegrationPointsArrayType& rIntegrationPoints,
                                         IntegrationInfo& rIntegrationInfo) override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    void CreateQuadraturePointGeometries(GeometriesArrayType& rResultGeometries,
                                         IndexType            NumberOfShapeFunctionDerivatives,
                                         IntegrationInfo&     rIntegrationInfo) override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

    void GlobalSpaceDerivatives(std::vector<CoordinatesArrayType>& rGlobalSpaceDerivatives,
                                IndexType                          IntegrationPointIndex,
                                const SizeType                     DerivativeOrder) const override
    {
        KRATOS_ERROR << IntegrationSchemeFunctionalityNotImplementedMessage();
    }

private:
    [[nodiscard]] PointerVector<Node> CreatePointsOfMidLine() const
    {
        const auto points                  = this->Points();
        const auto number_of_midline_nodes = std::size_t{points.size() / 2};

        auto result = PointerVector<Node>{};
        for (std::size_t i = 0; i < number_of_midline_nodes; ++i) {
            const auto mid_point = (points[i] + points[i + number_of_midline_nodes]) / 2;
            result.push_back(make_intrusive<Node>(points[i].Id(), mid_point));
        }
        return result;
    }

    std::unique_ptr<BaseType> mMidLineGeometry;
};

} // namespace Kratos
