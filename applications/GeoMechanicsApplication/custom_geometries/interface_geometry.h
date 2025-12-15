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

#include "custom_utilities/geometry_utilities.h"
#include "geometries/geometry.h"
#include "includes/node.h"

#include <algorithm>

namespace Kratos
{

template <typename MidGeometryType>
class InterfaceGeometry : public Geometry<Node>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceGeometry);

    using BaseType = Geometry<Node>;

    InterfaceGeometry() = default;

    explicit InterfaceGeometry(const PointsArrayType& rThisPoints)
        : InterfaceGeometry(0, rThisPoints)
    {
    }

    InterfaceGeometry(IndexType NewGeometryId, const PointsArrayType& rThisPoints)
        : BaseType(NewGeometryId, rThisPoints)
    {
        KRATOS_ERROR_IF_NOT((rThisPoints.size() == 4) || (rThisPoints.size() == 6 || rThisPoints.size() == 12) ||
                            (rThisPoints.size() == 8) || (rThisPoints.size() == 16))
            << "Number of nodes must be 2+2, 3+3, 6+6, 4+4 or 8+8\n";

        mpMidGeometry = std::make_shared<MidGeometryType>(CreatePointsOfMidGeometry());
        this->SetGeometryData(&mpMidGeometry->GetGeometryData());
    }

    [[nodiscard]] BaseType::Pointer Create(const PointsArrayType& rThisPoints) const override
    {
        constexpr auto id = IndexType{0};
        return Create(id, rThisPoints);
    }

    [[nodiscard]] BaseType::Pointer Create(const IndexType NewGeometryId, const PointsArrayType& rThisPoints) const override
    {
        return std::make_shared<InterfaceGeometry>(NewGeometryId, rThisPoints);
    }

    [[nodiscard]] GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return mpMidGeometry->GetGeometryFamily();
    }

    [[nodiscard]] GeometryData::KratosGeometryOrderType GetGeometryOrderType() const override
    {
        return mpMidGeometry->GetGeometryOrderType();
    }

    [[nodiscard]] double Area() const override { return mpMidGeometry->Area(); }

    [[nodiscard]] double ShapeFunctionValue(IndexType ShapeFunctionIndex,
                                            const CoordinatesArrayType& rLocalCoordinate) const override
    {
        return mpMidGeometry->ShapeFunctionValue(ShapeFunctionIndex, rLocalCoordinate);
    }

    Vector& ShapeFunctionsValues(Vector& rResult, const CoordinatesArrayType& rLocalCoordinate) const override
    {
        return mpMidGeometry->ShapeFunctionsValues(rResult, rLocalCoordinate);
    }

    Matrix& ShapeFunctionsLocalGradients(Matrix& rResult, const CoordinatesArrayType& rLocalCoordinate) const override
    {
        return mpMidGeometry->ShapeFunctionsLocalGradients(rResult, rLocalCoordinate);
    }

    Matrix& Jacobian(Matrix& rResult, const CoordinatesArrayType& rLocalCoordinate) const override
    {
        return mpMidGeometry->Jacobian(rResult, rLocalCoordinate);
    }

    [[nodiscard]] double DeterminantOfJacobian(const CoordinatesArrayType& rLocalCoordinate) const override
    {
        return mpMidGeometry->DeterminantOfJacobian(rLocalCoordinate);
    }

    Matrix& InverseOfJacobian(Matrix& rResult, const CoordinatesArrayType& rLocalCoordinate) const override
    {
        KRATOS_ERROR << "Inverse of Jacobian is not implemented for the interface geometry\n";
    }

    [[nodiscard]] double Length() const override { return mpMidGeometry->Length(); }

    [[nodiscard]] double DomainSize() const override { return mpMidGeometry->DomainSize(); }

    [[nodiscard]] std::string Info() const override
    {
        return "An interface geometry consisting of two sub-geometries with Info: " + mpMidGeometry->Info();
    }

    CoordinatesArrayType& PointLocalCoordinates(CoordinatesArrayType& rResult,
                                                const CoordinatesArrayType& rGlobalCoordinate) const override
    {
        return mpMidGeometry->PointLocalCoordinates(rResult, rGlobalCoordinate);
    }

    Matrix& PointsLocalCoordinates(Matrix& rResult) const override
    {
        return mpMidGeometry->PointsLocalCoordinates(rResult);
    }

    // The way for client code to access the mid-geometry
    const GeometryType::Pointer pGetGeometryPart(const IndexType) const override
    {
        return mpMidGeometry;
    }

    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

    void PrintData(std::ostream& rOStream) const override { mpMidGeometry->PrintData(rOStream); }

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

    GeometriesArrayType GenerateEdges() const override
    {
        KRATOS_ERROR_IF_NOT(mpMidGeometry->GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Linear)
            << "Edges can only be generated for line geometries. This is a surface interface "
               "geometry, which does not support edges.\n";

        return GenerateTwoSides();
    }

    GeometriesArrayType GenerateFaces() const override
    {
        KRATOS_ERROR_IF(mpMidGeometry->GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Linear)
            << "Faces can only be generated for surface geometries. This is a line "
               "interface geometry, which does not support faces.\n";

        return GenerateTwoSides();
    }

    GeometriesArrayType GenerateBoundariesEntities() const override
    {
        switch (mpMidGeometry->GetGeometryFamily()) {
            using enum GeometryData::KratosGeometryFamily;
        case Kratos_Linear:
            return this->GenerateEdges();
        case Kratos_Triangle:
        case Kratos_Quadrilateral:
            return this->GenerateFaces();
        default:
            KRATOS_ERROR << "Unsupported geometry type for generating boundaries entities\n";
        }
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
    [[nodiscard]] PointerVector<Node> CreatePointsOfMidGeometry() const
    {
        const auto points                       = this->Points();
        const auto number_of_mid_geometry_nodes = std::size_t{points.size() / 2};
        auto       result                       = PointerVector<Node>{number_of_mid_geometry_nodes};

        auto is_null = [](const auto& rNodePtr) { return rNodePtr == nullptr; };
        if (std::any_of(points.ptr_begin(), points.ptr_end(), is_null)) {
            // At least one point is not defined, so the points of the mid-geometry can't be
            // computed. As a result, all the mid-geometry points will be undefined.
            // This happens on element registration: the 'blue-print' element, creates a
            // geometry with a vector of null pointers to its nodes of the correct length.
            return result;
        }

        auto begin_of_first_side  = points.ptr_begin();
        auto begin_of_second_side = begin_of_first_side + number_of_mid_geometry_nodes;
        auto make_mid_point       = [](const auto& pPoint1, const auto& pPoint2) {
            return make_intrusive<Node>(pPoint1->Id(), Point{(*pPoint1 + *pPoint2) / 2});
        };
        std::transform(begin_of_first_side, begin_of_second_side, begin_of_second_side,
                       result.ptr_begin(), make_mid_point);
        return result;
    }

    [[nodiscard]] GeometriesArrayType GenerateTwoSides() const
    {
        const auto points = this->Points();

        // The first side is defined by the first half of the element nodes
        auto begin_of_second_side = points.ptr_begin() + (points.size() / 2);
        const auto nodes_of_first_side = PointerVector<Node>{points.ptr_begin(), begin_of_second_side};

        // The second side is defined by the second half of the element nodes. However, the
        // nodes must be traversed in opposite direction.
        auto nodes_of_second_side = PointerVector<Node>{begin_of_second_side, points.ptr_end()};
        GeometryUtilities::ReverseNodes(nodes_of_second_side, this->GetGeometryFamily(),
                                        this->GetGeometryOrderType());

        auto result = GeometriesArrayType{};
        result.push_back(std::make_shared<MidGeometryType>(nodes_of_first_side));
        result.push_back(std::make_shared<MidGeometryType>(nodes_of_second_side));
        return result;
    }

    [[nodiscard]] static std::string IntegrationSchemeFunctionalityNotImplementedMessage()
    {
        return "This Geometry type does not support functionality related to integration "
               "schemes.\n";
    }

    std::shared_ptr<BaseType> mpMidGeometry;
};

} // namespace Kratos
