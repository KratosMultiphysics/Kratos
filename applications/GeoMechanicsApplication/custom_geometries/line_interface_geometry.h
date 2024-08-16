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

class LineInterfaceGeometry : public Geometry<Node>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineInterfaceGeometry);

    LineInterfaceGeometry() = default;

    explicit LineInterfaceGeometry(const Geometry<Node>::PointsArrayType& rThisPoints)
        : LineInterfaceGeometry(0, rThisPoints)
    {
    }

    PointerVector<Node> CreatePointsOfMidLine()
    {
        const auto points                  = this->Points();
        auto       result                  = PointerVector<Node>{};
        const auto number_of_midline_nodes = std::size_t{points.size() / 2};

        for (std::size_t i = 0; i < number_of_midline_nodes; ++i) {
            auto mid_point = (points[i] + points[i + number_of_midline_nodes]) / 2;
            result.push_back(make_intrusive<Node>(i + 1, mid_point));
        }
        return result;
    }

    LineInterfaceGeometry(const IndexType NewGeometryId, const Geometry<Node>::PointsArrayType& rThisPoints)
        : Geometry<Node>(NewGeometryId, rThisPoints)
    {
        KRATOS_ERROR_IF_NOT((rThisPoints.size() == 4) || (rThisPoints.size() == 6))
            << "Number of nodes must be four or six\n";

        const auto new_points = CreatePointsOfMidLine();

        if (new_points.size() == 2) {
            mLineGeometry = std::make_unique<Line2D2<Node>>(new_points);
        } else {
            mLineGeometry = std::make_unique<Line2D3<Node>>(new_points);
        }
    }

    [[nodiscard]] Geometry<Node>::Pointer Create(const Geometry<Node>::PointsArrayType& rThisPoints) const override
    {
        constexpr auto id = IndexType{0};
        return Create(id, rThisPoints);
    }

    [[nodiscard]] Geometry<Node>::Pointer Create(const IndexType        NewGeometryId,
                                                 const PointsArrayType& rThisPoints) const override
    {
        return std::make_shared<LineInterfaceGeometry>(NewGeometryId, rThisPoints);
    }

    [[nodiscard]] double ShapeFunctionValue(IndexType                   ShapeFunctionIndex,
                                            const CoordinatesArrayType& rCoordinates) const override
    {
        return mLineGeometry->ShapeFunctionValue(ShapeFunctionIndex, rCoordinates);
    }

    Vector& ShapeFunctionsValues(Vector& rResult, const CoordinatesArrayType& rCoordinates) const override
    {
        return mLineGeometry->ShapeFunctionsValues(rResult, rCoordinates);
    }

    Matrix& ShapeFunctionsLocalGradients(Matrix& rResult, const CoordinatesArrayType& rPoint) const override
    {
        return mLineGeometry->ShapeFunctionsLocalGradients(rResult, rPoint);
    }

    Matrix& Jacobian(Matrix& rResult, const CoordinatesArrayType& rCoordinates) const override
    {
        return mLineGeometry->Jacobian(rResult, rCoordinates);
    }

    [[nodiscard]] double DeterminantOfJacobian(const CoordinatesArrayType& rPoint) const override
    {
        return mLineGeometry->DeterminantOfJacobian(rPoint);
    }

    Matrix& InverseOfJacobian(Matrix& rResult, const CoordinatesArrayType& rCoordinates) const override
    {
        KRATOS_ERROR << "Inverse of Jacobian is not implemented for the line interface geometry\n";
    }

    [[nodiscard]] double Length() const override { return mLineGeometry->Length(); }

    [[nodiscard]] double DomainSize() const override { return mLineGeometry->DomainSize(); }

    [[nodiscard]] std::string Info() const override
    {
        return "An interface geometry consisting of two sub-geometries with Info: " + mLineGeometry->Info();
    }

    CoordinatesArrayType& PointLocalCoordinates(CoordinatesArrayType&       rResult,
                                                const CoordinatesArrayType& rPoint) const override
    {
        return mLineGeometry->PointLocalCoordinates(rResult, rPoint);
    }

    Matrix& PointsLocalCoordinates(Matrix& rResult) const override
    {
        return mLineGeometry->PointsLocalCoordinates(rResult);
    }

    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

    void PrintData(std::ostream& rOStream) const override { mLineGeometry->PrintData(rOStream); }

private:
    std::unique_ptr<Geometry<Node>> mLineGeometry;
};

} // namespace Kratos
