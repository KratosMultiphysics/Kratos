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
//                   Marjan Fathian
//
#pragma once

#include "geo_mechanics_application_constants.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

class StressStatePolicy
{
public:
    virtual ~StressStatePolicy() = default;

    [[nodiscard]] virtual Matrix CalculateBMatrix(const Matrix&         rDN_DX,
                                                  const Vector&         rN,
                                                  const Geometry<Node>& rGeometry) const = 0;
    [[nodiscard]] virtual double CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                                                 double DetJ,
                                                                 const Geometry<Node>& rGeometry) const = 0;
    [[nodiscard]] virtual Vector CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const = 0;
    [[nodiscard]] virtual std::unique_ptr<StressStatePolicy> Clone() const          = 0;
    [[nodiscard]] virtual const Vector&                      GetVoigtVector() const = 0;
    [[nodiscard]] virtual SizeType                           GetVoigtSize() const   = 0;

    static constexpr std::size_t GetStressTensorSize(std::size_t Dimension)
    {
        return Dimension == N_DIM_3D ? STRESS_TENSOR_SIZE_3D : STRESS_TENSOR_SIZE_2D;
    }

protected:
    static const Vector VoigtVector2D;
    static const Vector VoigtVector3D;

    static constexpr SizeType GetVoigtSize2D() { return VOIGT_SIZE_2D_PLANE_STRAIN; }

    static constexpr SizeType GetVoigtSize3D() { return VOIGT_SIZE_3D; }

private:
    static Vector DefineVoigtVector(std::size_t Dimension);

    static constexpr std::size_t GetVoigtSize(std::size_t Dimension)
    {
        return Dimension == N_DIM_3D ? GetVoigtSize3D() : GetVoigtSize2D();
    }
};

} // namespace Kratos