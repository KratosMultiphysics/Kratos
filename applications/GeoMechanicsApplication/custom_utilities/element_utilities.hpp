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
    static void CalculateNuMatrix(BoundedMatrix<double, TDim, TDim * TNumNodes>& rNu,
                                  const Matrix&                                  NContainer,
                                  unsigned int                                   GPoint)
    {
        CalculateNuMatrix(TDim, TNumNodes, rNu, NContainer, GPoint);
    }

    template <typename MatrixType1, typename MatrixType2>
    static void CalculateNuMatrix(std::size_t        Dim,
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
    static void InterpolateVariableWithComponents(array_1d<double, TDim>& rVector,
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
    static void InterpolateVariableWithComponents(Vector&       rVector,
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
    static void GetNodalVariableVector(array_1d<double, TDim * TNumNodes>&  rNodalVariableVector,
                                       const Element::GeometryType&         Geom,
                                       const Variable<array_1d<double, 3>>& Variable,
                                       IndexType                            SolutionStepIndex = 0)
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

    static void FillPermeabilityMatrix(BoundedMatrix<double, 1, 1>&   rPermeabilityMatrix,
                                       const Element::PropertiesType& Prop);

    static void FillPermeabilityMatrix(BoundedMatrix<double, 2, 2>&   rPermeabilityMatrix,
                                       const Element::PropertiesType& Prop);

    static void FillPermeabilityMatrix(BoundedMatrix<double, 3, 3>&   rPermeabilityMatrix,
                                       const Element::PropertiesType& Prop);

    static Matrix FillPermeabilityMatrix(const Element::PropertiesType& Prop, std::size_t Dimension);

    template <typename MatrixType1, typename MatrixType2>
    static void AssembleUUBlockMatrix(MatrixType1& rLeftHandSideMatrix, const MatrixType2& rUUBlockMatrix)
    {
        constexpr auto row_offset    = std::size_t{0};
        constexpr auto column_offset = row_offset;
        AddMatrixAtPosition(rUUBlockMatrix, rLeftHandSideMatrix, row_offset, column_offset);
    }

    template <typename MatrixType1, typename MatrixType2>
    static void AssembleUPBlockMatrix(MatrixType1& rLeftHandSideMatrix, const MatrixType2& rUPBlockMatrix)
    {
        constexpr auto row_offset    = std::size_t{0};
        const auto     column_offset = rLeftHandSideMatrix.size2() - rUPBlockMatrix.size2();
        AddMatrixAtPosition(rUPBlockMatrix, rLeftHandSideMatrix, row_offset, column_offset);
    }

    template <typename MatrixType1, typename MatrixType2>
    static void AssemblePUBlockMatrix(MatrixType1& rLeftHandSideMatrix, const MatrixType2& rPUBlockMatrix)
    {
        const auto     row_offset    = rLeftHandSideMatrix.size1() - rPUBlockMatrix.size1();
        constexpr auto column_offset = std::size_t{0};
        AddMatrixAtPosition(rPUBlockMatrix, rLeftHandSideMatrix, row_offset, column_offset);
    }

    template <typename MatrixType1, typename MatrixType2>
    static void AssemblePPBlockMatrix(MatrixType1& rLeftHandSideMatrix, const MatrixType2& rPPBlockMatrix)
    {
        const auto row_offset    = rLeftHandSideMatrix.size1() - rPPBlockMatrix.size1();
        const auto column_offset = row_offset;
        AddMatrixAtPosition(rPPBlockMatrix, rLeftHandSideMatrix, row_offset, column_offset);
    }

    template <typename VectorType1, typename VectorType2>
    static void AssembleUBlockVector(VectorType1& rRightHandSideVector, const VectorType2& rUBlockVector)
    {
        constexpr auto offset = std::size_t{0};
        AddVectorAtPosition(rUBlockVector, rRightHandSideVector, offset);
    }

    template <typename VectorType1, typename VectorType2>
    static void AssemblePBlockVector(VectorType1& rRightHandSideVector, const VectorType2& rPBlockVector)
    {
        const auto offset = rRightHandSideVector.size() - rPBlockVector.size();
        AddVectorAtPosition(rPBlockVector, rRightHandSideVector, offset);
    }

    template <typename MatrixType1, typename MatrixType2>
    static void AssignMatrixAtPosition(MatrixType1&       rDestinationMatrix,
                                       const MatrixType2& rSourceMatrix,
                                       std::size_t        RowOffset,
                                       std::size_t        ColumnOffset)
    {
        KRATOS_DEBUG_ERROR_IF(RowOffset + rSourceMatrix.size1() > rDestinationMatrix.size1())
            << "Can't assign submatrix: last row index (" << RowOffset + rSourceMatrix.size1()
            << ") exceeds the row size of the destination matrix (" << rDestinationMatrix.size1() << ")\n";
        KRATOS_DEBUG_ERROR_IF(ColumnOffset + rSourceMatrix.size2() > rDestinationMatrix.size2())
            << "Can't assign submatrix: last column index (" << ColumnOffset + rSourceMatrix.size2()
            << ") exceeds the column size of the destination matrix (" << rDestinationMatrix.size2() << ")\n";

        subrange(rDestinationMatrix, RowOffset, RowOffset + rSourceMatrix.size1(), ColumnOffset,
                 ColumnOffset + rSourceMatrix.size2()) = rSourceMatrix;
    }

    template <typename MatrixType1, typename MatrixType2>
    static void AssignUUBlockMatrix(MatrixType1& rDestinationMatrix, const MatrixType2& rUUBlockMatrix)
    {
        constexpr auto row_offset    = std::size_t{0};
        constexpr auto column_offset = std::size_t{0};
        AssignMatrixAtPosition(rDestinationMatrix, rUUBlockMatrix, row_offset, column_offset);
    }

    template <typename VectorType1, typename VectorType2>
    static void AssignVectorAtPosition(VectorType1& rDestinationVector, const VectorType2& rSourceVector, std::size_t Offset)
    {
        KRATOS_DEBUG_ERROR_IF(Offset + rSourceVector.size() > rDestinationVector.size())
            << "Can't assign subvector: last index (" << Offset + rSourceVector.size()
            << ") exceeds the size of the destination vector (" << rDestinationVector.size() << ")\n";

        subrange(rDestinationVector, Offset, Offset + rSourceVector.size()) = rSourceVector;
    }

    template <typename VectorType1, typename VectorType2>
    static void AssignUBlockVector(VectorType1& rDestinationVector, const VectorType2& rUBlockVector)
    {
        constexpr auto offset = std::size_t{0};
        AssignVectorAtPosition(rDestinationVector, rUBlockVector, offset);
    }

    /**
     * Calculates the radius of axisymmetry
     * @param rN: The Gauss Point shape function
     * @param rGeometry: The geometry studied
     * @return Radius: The radius of axisymmetry
     */
    static double CalculateRadius(const Vector& rN, const GeometryType& rGeometry);

    static double CalculateAxisymmetricCircumference(const Vector& rN, const GeometryType& rGeometry);

    static Vector CalculateNodalHydraulicHeadFromWaterPressures(const GeometryType& rGeom,
                                                                const Properties&   rProp);

    static std::size_t                     GetNumberOfIntegrationPointsOf(const Element& rElement);
    static Geo::IntegrationPointVectorType GetIntegrationPointsOf(const Element& rElement);
    static std::vector<Vector> EvaluateShapeFunctionsAtIntegrationPoints(const Geo::IntegrationPointVectorType& rIntegrationPoints,
                                                                         const Geometry<Node>& rGeometry);

    static Vector EvaluateDeterminantsOfJacobiansAtIntegrationPoints(const Geo::IntegrationPointVectorType& rIntegrationPoints,
                                                                     const Geometry<Node>& rGeometry);

    template <typename MatrixType1, typename MatrixType2>
    static void AddMatrixAtPosition(const MatrixType1& rSourceMatrix,
                                    MatrixType2&       rDestinationMatrix,
                                    const std::size_t  RowOffset,
                                    const std::size_t  ColumnOffset)
    {
        const std::size_t size1 = rSourceMatrix.size1();
        const std::size_t size2 = rSourceMatrix.size2();

        for (std::size_t i = 0; i < size1; ++i) {
            const std::size_t di = i + RowOffset;

            for (std::size_t j = 0; j < size2; ++j) {
                rDestinationMatrix(di, j + ColumnOffset) += rSourceMatrix(i, j);
            }
        }
    }

private:
    template <typename VectorType1, typename VectorType2>
    static void AddVectorAtPosition(const VectorType1& rSourceVector, VectorType2& rDestinationVector, std::size_t Offset)
    {
        auto pos = std::begin(rDestinationVector) + Offset;
        std::transform(std::begin(rSourceVector), std::end(rSourceVector), pos, pos, std::plus<double>{});
    }

}; /* Class GeoElementUtilities*/
} /* namespace Kratos.*/
