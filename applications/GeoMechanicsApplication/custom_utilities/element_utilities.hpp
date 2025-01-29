// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

// Project includes
#include "includes/element.h"
#include "includes/kratos_export_api.h"
#include "utilities/math_utils.h"

// Application includes
#include "geo_aliases.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoElementUtilities
{
public:
    using IndexType    = std::size_t;
    using GeometryType = Geometry<Node>;

    template <unsigned int TDim, unsigned int TNumNodes>
    static inline void CalculateNuMatrix(BoundedMatrix<double, TDim, TDim * TNumNodes>& rNu,
                                         const Matrix&                                  NContainer,
                                         unsigned int                                   GPoint)
    {
        CalculateNuMatrix(TDim, TNumNodes, rNu, NContainer, GPoint);
    }

    template <typename MatrixType1, typename MatrixType2>
    static inline void CalculateNuMatrix(std::size_t        Dim,
                                         std::size_t        NumNodes,
                                         MatrixType1&       rNu,
                                         const MatrixType2& NContainer,
                                         unsigned int       GPoint)
    {
        for (unsigned int i = 0; i < Dim; ++i) {
            unsigned int index = i;
            for (unsigned int j = 0; j < NumNodes; ++j) {
                rNu(i, index) = NContainer(GPoint, j);
                index += Dim;
            }
        }
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    static inline void InterpolateVariableWithComponents(array_1d<double, TDim>& rVector,
                                                         const Matrix&           NContainer,
                                                         const array_1d<double, TDim * TNumNodes>& VariableWithComponents,
                                                         unsigned int GPoint)
    {
        noalias(rVector) = ZeroVector(TDim);

        unsigned int index = 0;
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            for (unsigned int j = 0; j < TDim; ++j) {
                rVector[j] += NContainer(GPoint, i) * VariableWithComponents[index];
                index++;
            }
        }
    }

    template <unsigned int TDof, unsigned int TNumNodes>
    static inline void InterpolateVariableWithComponents(Vector&       rVector,
                                                         const Matrix& NContainer,
                                                         const Vector& VariableWithComponents,
                                                         unsigned int  GPoint)
    {
        if (rVector.size() != TDof) rVector.resize(TDof, false);
        KRATOS_ERROR_IF(VariableWithComponents.size() != TDof * TNumNodes)
            << "Wrong size in InterpolateVariableWithComponents" << std::endl;

        noalias(rVector) = ZeroVector(TDof);

        unsigned int index = 0;
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            for (unsigned int j = 0; j < TDof; ++j) {
                rVector[j] += NContainer(GPoint, i) * VariableWithComponents[index++];
            }
        }
    }

    static void FillArray1dOutput(array_1d<double, 3>& rOutputValue, const array_1d<double, 2>& ComputedValue);

    static void FillArray1dOutput(array_1d<double, 3>& rOutputValue, const array_1d<double, 3>& ComputedValue);

    template <unsigned int TDim, unsigned int TNumNodes>
    static inline void GetNodalVariableVector(array_1d<double, TDim * TNumNodes>& rNodalVariableVector,
                                              const Element::GeometryType&         Geom,
                                              const Variable<array_1d<double, 3>>& Variable,
                                              IndexType SolutionStepIndex = 0)
    {
        array_1d<double, 3> NodalVariableAux;
        unsigned int        index = 0;
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            noalias(NodalVariableAux) = Geom[i].FastGetSolutionStepValue(Variable, SolutionStepIndex);
            for (unsigned int j = 0; j < TDim; ++j) {
                rNodalVariableVector[index] = NodalVariableAux[j];
                index++;
            }
        }
    }

    static void CheckPermeabilityProperties(const Element::PropertiesType& rProp,
                                                                 size_t Dimension);

    static void FillPermeabilityMatrix(BoundedMatrix<double, 1, 1>&   rPermeabilityMatrix,
                                       const Element::PropertiesType& Prop);

    static void FillPermeabilityMatrix(BoundedMatrix<double, 2, 2>&   rPermeabilityMatrix,
                                       const Element::PropertiesType& Prop);

    static void FillPermeabilityMatrix(BoundedMatrix<double, 3, 3>&   rPermeabilityMatrix,
                                       const Element::PropertiesType& Prop);

    static Matrix FillPermeabilityMatrix(const Element::PropertiesType& Prop, std::size_t Dimension);

    static void InvertMatrix2(BoundedMatrix<double, 2, 2>&       rInvertedMatrix,
                              const BoundedMatrix<double, 2, 2>& InputMatrix,
                              double&                            InputMatrixDet);

    static void InvertMatrix2(BoundedMatrix<double, 2, 2>&       rInvertedMatrix,
                              const BoundedMatrix<double, 2, 2>& InputMatrix);

    template <typename MatrixType1, typename MatrixType2>
    static inline void AssembleUUBlockMatrix(MatrixType1& rLeftHandSideMatrix, const MatrixType2& rUUBlockMatrix)
    {
        constexpr auto row_offset    = std::size_t{0};
        constexpr auto column_offset = row_offset;
        AddMatrixAtPosition(rUUBlockMatrix, rLeftHandSideMatrix, row_offset, column_offset);
    }

    template <typename MatrixType1, typename MatrixType2>
    static inline void AssembleUPBlockMatrix(MatrixType1& rLeftHandSideMatrix, const MatrixType2& rUPBlockMatrix)
    {
        constexpr auto row_offset    = std::size_t{0};
        const auto     column_offset = rLeftHandSideMatrix.size2() - rUPBlockMatrix.size2();
        AddMatrixAtPosition(rUPBlockMatrix, rLeftHandSideMatrix, row_offset, column_offset);
    }

    template <typename MatrixType1, typename MatrixType2>
    static inline void AssemblePUBlockMatrix(MatrixType1& rLeftHandSideMatrix, const MatrixType2& rPUBlockMatrix)
    {
        const auto     row_offset    = rLeftHandSideMatrix.size1() - rPUBlockMatrix.size1();
        constexpr auto column_offset = std::size_t{0};
        AddMatrixAtPosition(rPUBlockMatrix, rLeftHandSideMatrix, row_offset, column_offset);
    }

    template <typename MatrixType1, typename MatrixType2>
    static inline void AssemblePPBlockMatrix(MatrixType1& rLeftHandSideMatrix, const MatrixType2& rPPBlockMatrix)
    {
        const auto row_offset    = rLeftHandSideMatrix.size1() - rPPBlockMatrix.size1();
        const auto column_offset = row_offset;
        AddMatrixAtPosition(rPPBlockMatrix, rLeftHandSideMatrix, row_offset, column_offset);
    }

    template <typename VectorType1, typename VectorType2>
    static inline void AssembleUBlockVector(VectorType1& rRightHandSideVector, const VectorType2& rUBlockVector)
    {
        constexpr auto offset = std::size_t{0};
        AddVectorAtPosition(rUBlockVector, rRightHandSideVector, offset);
    }

    template <typename VectorType1, typename VectorType2>
    static inline void AssemblePBlockVector(VectorType1& rRightHandSideVector, const VectorType2& rPBlockVector)
    {
        const auto offset = rRightHandSideVector.size() - rPBlockVector.size();
        AddVectorAtPosition(rPBlockVector, rRightHandSideVector, offset);
    }

    static void CalculateNewtonCotesLocalShapeFunctionsGradients(BoundedMatrix<double, 2, 2>& DN_DeContainer);

    static void CalculateNewtonCotesLocalShapeFunctionsGradients(BoundedMatrix<double, 3, 3>& DN_DeContainer);

    static void CalculateNewtonCotesShapeFunctions(BoundedMatrix<double, 3, 3>& NContainer);

    static void CalculateEquallyDistributedPointsLineShapeFunctions3N(Matrix& NContainer);

    static void CalculateEquallyDistributedPointsLineGradientShapeFunctions3N(GeometryData::ShapeFunctionsGradientsType& DN_DeContainer);

    /**
     * Calculates the radius of axisymmetry
     * @param N: The Gauss Point shape function
     * @param Geom: The geometry studied
     * @return Radius: The radius of axisymmetry
     */

    static double CalculateRadius(const Vector& N, const GeometryType& Geom);

    static double CalculateAxisymmetricCircumference(const Vector& N, const GeometryType& Geom);

    static void CalculateExtrapolationMatrixTriangle(Matrix& rExtrapolationMatrix,
                                                     const GeometryData::IntegrationMethod& rIntegrationMethod);

    static void CalculateExtrapolationMatrixQuad(Matrix& rExtrapolationMatrix,
                                                 const GeometryData::IntegrationMethod& rIntegrationMethod);

    static void CalculateExtrapolationMatrixTetra(Matrix& rExtrapolationMatrix,
                                                  const GeometryData::IntegrationMethod& rIntegrationMethod);

    static void CalculateExtrapolationMatrixHexa(Matrix& rExtrapolationMatrix,
                                                 const GeometryData::IntegrationMethod& rIntegrationMethod);

    static Vector CalculateNodalHydraulicHeadFromWaterPressures(const GeometryType& rGeom,
                                                                const Properties&   rProp);

    static std::vector<Vector> EvaluateShapeFunctionsAtIntegrationPoints(const Geo::IntegrationPointVectorType& rIntegrationPoints,
                                                                         const Geometry<Node>& rGeometry);

private:
    template <typename VectorType1, typename VectorType2>
    static void AddVectorAtPosition(const VectorType1& rSourceVector, VectorType2& rDestinationVector, std::size_t Offset)
    {
        auto pos = std::begin(rDestinationVector) + Offset;
        std::transform(std::begin(rSourceVector), std::end(rSourceVector), pos, pos, std::plus<double>{});
    }

    template <typename MatrixType1, typename MatrixType2>
    static void AddMatrixAtPosition(const MatrixType1& rSourceMatrix,
                                    MatrixType2&       rDestinationMatrix,
                                    std::size_t        RowOffset,
                                    std::size_t        ColumnOffset)
    {
        for (auto i = std::size_t{0}; i < rSourceMatrix.size1(); ++i) {
            for (auto j = std::size_t{0}; j < rSourceMatrix.size2(); ++j) {
                rDestinationMatrix(i + RowOffset, j + ColumnOffset) += rSourceMatrix(i, j);
            }
        }
    }

    static int CheckPropertyExistsAndIsNotNegative(const Variable<double>&        rVariable,
                                                   const Element::PropertiesType& rProp);

}; /* Class GeoElementUtilities*/
} /* namespace Kratos.*/
