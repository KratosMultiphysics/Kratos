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

#if !defined(KRATOS_GEO_ELEMENT_UTILITIES )
#define  KRATOS_GEO_ELEMENT_UTILITIES

// Project includes
#include "utilities/math_utils.h"
#include "includes/element.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class GeoElementUtilities
{
public:
    using IndexType    = std::size_t;
    using GeometryType = Geometry<Node>;

    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void CalculateNuMatrix(BoundedMatrix<double,TDim,TDim*TNumNodes>& rNu,
                                         const Matrix& NContainer,
                                         unsigned int GPoint)
    {
        for (unsigned int i=0; i < TDim; ++i) {
            unsigned int index = i - TDim;
            for (unsigned int j=0; j < TNumNodes; ++j) {
                index += TDim;
                rNu(i, index) = NContainer(GPoint, j);
            }
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void InterpolateVariableWithComponents(array_1d<double, TDim>& rVector,
                                                         const Matrix& NContainer,
                                                         const array_1d<double, TDim*TNumNodes>& VariableWithComponents,
                                                         unsigned int GPoint)
    {
        noalias(rVector) = ZeroVector(TDim);

        unsigned int index = 0;
        for (unsigned int i=0; i<TNumNodes; ++i) {
            for (unsigned int j=0; j<TDim; ++j) {
                rVector[j] += NContainer(GPoint,i)*VariableWithComponents[index++];
            }
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template< unsigned int TDof, unsigned int TNumNodes >
    static inline void InterpolateVariableWithComponents(Vector& rVector,
                                                         const Matrix& NContainer,
                                                         const Vector& VariableWithComponents,
                                                         unsigned int GPoint)
    {
        if ( rVector.size() != TDof ) rVector.resize( TDof, false );
        KRATOS_ERROR_IF(VariableWithComponents.size() != TDof*TNumNodes) << "Wrong size in InterpolateVariableWithComponents" << std::endl;

        noalias(rVector) = ZeroVector(TDof);

        unsigned int index = 0;
        for (unsigned int i=0; i<TNumNodes; ++i) {
            for (unsigned int j=0; j<TDof; ++j) {
                rVector[j] += NContainer(GPoint,i)*VariableWithComponents[index++];
            }
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void FillArray1dOutput(array_1d<double,3>& rOutputValue,
                                         const array_1d<double,2>& ComputedValue)
    {
        rOutputValue[0] = ComputedValue[0];
        rOutputValue[1] = ComputedValue[1];
        rOutputValue[2] = 0.0;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void FillArray1dOutput(array_1d<double,3>& rOutputValue,
                                         const array_1d<double,3>& ComputedValue)
    {
        rOutputValue[0] = ComputedValue[0];
        rOutputValue[1] = ComputedValue[1];
        rOutputValue[2] = ComputedValue[2];
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void GetNodalVariableVector(array_1d<double,TDim*TNumNodes>& rNodalVariableVector,
                                              const Element::GeometryType& Geom,
                                              const Variable<array_1d<double,3>>& Variable,
                                              IndexType SolutionStepIndex = 0)
    {
        array_1d<double, 3> NodalVariableAux;
        unsigned int index = 0;
        for (unsigned int i=0; i < TNumNodes; ++i) {
            noalias(NodalVariableAux) = Geom[i].FastGetSolutionStepValue(Variable, SolutionStepIndex);
            for (unsigned int j=0; j < TDim; ++j) {
                rNodalVariableVector[index++] = NodalVariableAux[j];
            }
        }
    }

//----------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void GetNodalVariableMatrix(Matrix& rNodalVariableMatrix,
                                              const Element::GeometryType& rGeom,
                                              const Variable<array_1d<double,3>>& Variable,
                                              IndexType SolutionStepIndex = 0)
    {
        rNodalVariableMatrix.resize(TNumNodes , TDim, false);

        for ( IndexType node = 0; node < TNumNodes; ++node ) {
            const array_1d<double, 3>& NodalVariableAux = rGeom[node].FastGetSolutionStepValue(Variable, SolutionStepIndex);

            for ( IndexType iDim = 0; iDim < TDim; ++iDim )
                rNodalVariableMatrix(node, iDim) = NodalVariableAux[iDim];
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void FillPermeabilityMatrix(BoundedMatrix<double,2,2>& rPermeabilityMatrix,
                                              const Element::PropertiesType& Prop)
    {
        //2D
        rPermeabilityMatrix(0,0) = Prop[PERMEABILITY_XX];
        rPermeabilityMatrix(1,1) = Prop[PERMEABILITY_YY];

        rPermeabilityMatrix(0,1) = Prop[PERMEABILITY_XY];
        rPermeabilityMatrix(1,0) = rPermeabilityMatrix(0,1);
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void FillPermeabilityMatrix(BoundedMatrix<double,3,3>& rPermeabilityMatrix,
                                              const Element::PropertiesType& Prop)
    {
        //3D
        rPermeabilityMatrix(0,0) = Prop[PERMEABILITY_XX];
        rPermeabilityMatrix(1,1) = Prop[PERMEABILITY_YY];
        rPermeabilityMatrix(2,2) = Prop[PERMEABILITY_ZZ];

        rPermeabilityMatrix(0,1) = Prop[PERMEABILITY_XY];
        rPermeabilityMatrix(1,0) = rPermeabilityMatrix(0,1);

        rPermeabilityMatrix(1,2) = Prop[PERMEABILITY_YZ];
        rPermeabilityMatrix(2,1) = rPermeabilityMatrix(1,2);

        rPermeabilityMatrix(2,0) = Prop[PERMEABILITY_ZX];
        rPermeabilityMatrix(0,2) = rPermeabilityMatrix(2,0);
    }


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void InvertMatrix2(BoundedMatrix<double,2,2>& rInvertedMatrix,
                                     const BoundedMatrix<double,2,2>& InputMatrix,
                                     double &InputMatrixDet)
    {
        KRATOS_TRY

        const double numerical_limit = std::numeric_limits<double>::epsilon();

        InputMatrixDet = InputMatrix(0,0)*InputMatrix(1,1)-InputMatrix(0,1)*InputMatrix(1,0);

        if (InputMatrixDet < numerical_limit) {
            KRATOS_ERROR << "determinant zero or negative" << std::endl;
        }

        rInvertedMatrix(0,0) =  InputMatrix(1,1)/InputMatrixDet;
        rInvertedMatrix(0,1) = -InputMatrix(0,1)/InputMatrixDet;
        rInvertedMatrix(1,0) = -InputMatrix(1,0)/InputMatrixDet;
        rInvertedMatrix(1,1) =  InputMatrix(0,0)/InputMatrixDet;

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void InvertMatrix2(BoundedMatrix<double,2,2>& rInvertedMatrix,
                                     const BoundedMatrix<double,2,2>& InputMatrix)
    {
        KRATOS_TRY

        double InputMatrixDet;

        InvertMatrix2(rInvertedMatrix, InputMatrix, InputMatrixDet);

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template< unsigned int TDim>
    static inline void AssembleDensityMatrix(BoundedMatrix<double,TDim+1, TDim+1> &DensityMatrix,
                                             double Density)
    {
        for (unsigned int idim = 0; idim < TDim; ++idim) {
            for (unsigned int jdim = 0; jdim < TDim; ++jdim) {
                DensityMatrix(idim, jdim) = Density;
            }
        }
    }

    template <typename MatrixType>
    static inline void AssembleDensityMatrix(MatrixType& DensityMatrix, double Density)
    {
        for (unsigned int idim = 0; idim < DensityMatrix.size1(); ++idim) {
            for (unsigned int jdim = 0; jdim < DensityMatrix.size2(); ++jdim) {
                DensityMatrix(idim, jdim) = Density;
            }
        }
    }

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

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void CalculateNewtonCotesLocalShapeFunctionsGradients(BoundedMatrix<double,2,2>& DN_DeContainer)
    {
        //Line 2-noded
        // nodes:
        // 0------------1
        const unsigned int NumNodes = 2;
        const std::vector<double> Xi{-1.0, 1.0};

        noalias(DN_DeContainer) = ZeroMatrix(NumNodes,NumNodes);

        for (unsigned int integrationPoint = 0; integrationPoint < NumNodes; ++integrationPoint) {
            DN_DeContainer(integrationPoint,0) = - 0.5;
            DN_DeContainer(integrationPoint,1) =   0.5;
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void CalculateNewtonCotesLocalShapeFunctionsGradients(BoundedMatrix<double,3,3>& DN_DeContainer)
    {
        //Line 3-noded
        // nodes:
        // 0------2------1
        const unsigned int NumNodes = 3;
        const std::vector<double> Xi{-1.0, 0.0, 1.0};

        noalias(DN_DeContainer) = ZeroMatrix(NumNodes,NumNodes);

        for (unsigned int integrationPoint = 0; integrationPoint < NumNodes; ++integrationPoint) {
            DN_DeContainer(integrationPoint,0) =  Xi[integrationPoint] - 0.5;
            DN_DeContainer(integrationPoint,1) =  Xi[integrationPoint] + 0.5;
            DN_DeContainer(integrationPoint,2) = -Xi[integrationPoint] * 2.0;
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void CalculateNewtonCotesShapeFunctions(BoundedMatrix<double,3,3>& NContainer)
    {
        //Line 3-noded
        // nodes:
        // 0------2------1
        const unsigned int NumNodes = 3;

        for (unsigned int integrationPoint = 0; integrationPoint < NumNodes; ++integrationPoint) {
            for (unsigned int node = 0; node < NumNodes; ++node) {
                NContainer(integrationPoint,node) = (integrationPoint == node ? 1.0 : 0.0);
            }
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void CalculateEquallyDistributedPointsLineShapeFunctions3N(Matrix& NContainer)
    {
        //Line 3-noded
        // nodes:
        // 0------2------1
        const unsigned int NumNodes = 3;
        const std::vector<double> Xi{-1.0/3.0, 1.0/3.0};

        if ( NContainer.size1() != Xi.size() || NContainer.size2() != NumNodes) 
            NContainer.resize( Xi.size(), NumNodes, false );

        for (unsigned int integrationPoint = 0; integrationPoint < Xi.size(); ++integrationPoint) {
            const double& X = Xi[integrationPoint];
            NContainer(integrationPoint, 0) = -0.5 * (1.0 - X) * X;
            NContainer(integrationPoint, 1) =  0.5 * (1.0 + X) * X;
            NContainer(integrationPoint, 2) = (1.0 + X) * (1.0 - X);
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void CalculateEquallyDistributedPointsLineGradientShapeFunctions3N(GeometryData::ShapeFunctionsGradientsType& DN_DeContainer)
    {
        //Line 3-noded
        // nodes:
        // 0------2------1
        const unsigned int NumNodes = 3;
        const std::vector<double> Xi{-1.0/3.0, 1.0/3.0};

        if ( DN_DeContainer.size() != Xi.size()) 
            DN_DeContainer.resize( Xi.size());

        for (unsigned int integrationPoint = 0; integrationPoint < Xi.size(); ++integrationPoint) {
            if ( DN_DeContainer[integrationPoint].size1() != NumNodes || DN_DeContainer[integrationPoint].size2() != 1) 
                DN_DeContainer[integrationPoint].resize( NumNodes, 1);

            DN_DeContainer[integrationPoint](0,0) =  Xi[integrationPoint] - 0.5;
            DN_DeContainer[integrationPoint](1,0) =  Xi[integrationPoint] + 0.5;
            DN_DeContainer[integrationPoint](2,0) = -Xi[integrationPoint] * 2.0;
        }
    }

    /**
     * Calculates the radius of axisymmetry
     * @param N: The Gauss Point shape function
     * @param Geom: The geometry studied
     * @return Radius: The radius of axisymmetry
     */

    static inline double CalculateRadius(const Vector& N, const GeometryType& Geom)
    {
        double Radius = 0.0;

        for (unsigned int iNode = 0; iNode < Geom.size(); ++iNode) {
            // Displacement from the reference to the current configuration
            const array_1d<double, 3 >& CurrentPosition = Geom[iNode].Coordinates();
            Radius += CurrentPosition[0] * N[iNode];
        }

        return Radius;
    }

    static inline double CalculateAxisymmetricCircumference(const Vector& N, const GeometryType& Geom)
    {
        const double Radius = CalculateRadius(N, Geom);
        const double Circumference = 2.0 * Globals::Pi * Radius;
        return Circumference;
    }


    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void CalculateExtrapolationMatrixTriangle(Matrix& rExtrapolationMatrix, const GeometryData::IntegrationMethod& rIntegrationMethod)
    {
        /// The matrix contains the shape functions at each GP evaluated at each node.
        /// Rows: nodes
        /// Columns: GP

        //Triangle_2d_3
        //GI_GAUSS_2

        if (rIntegrationMethod == GeometryData::IntegrationMethod::GI_GAUSS_1)
        {
            if ((rExtrapolationMatrix.size1() != 3) || (rExtrapolationMatrix.size2() != 1))
            {
                rExtrapolationMatrix.resize(3, 1, false);
            }
                
            rExtrapolationMatrix(0, 0) = 1; 
            rExtrapolationMatrix(1, 0) = 1; 
            rExtrapolationMatrix(2, 0) = 1;
        }
        else if (rIntegrationMethod == GeometryData::IntegrationMethod::GI_GAUSS_2)
        {
            if ((rExtrapolationMatrix.size1() != 3) || (rExtrapolationMatrix.size2() != 3))
            {
                rExtrapolationMatrix.resize(3, 3, false);
            }
            rExtrapolationMatrix(0, 0) = 1.6666666666666666666; rExtrapolationMatrix(0, 1) = -0.33333333333333333333; rExtrapolationMatrix(0, 2) = -0.33333333333333333333;
            rExtrapolationMatrix(1, 0) = -0.33333333333333333333; rExtrapolationMatrix(1, 1) = 1.6666666666666666666; rExtrapolationMatrix(1, 2) = -0.33333333333333333333;
            rExtrapolationMatrix(2, 0) = -0.33333333333333333333; rExtrapolationMatrix(2, 1) = -0.33333333333333333333; rExtrapolationMatrix(2, 2) = 1.6666666666666666666;
        }
        else
        {
            KRATOS_ERROR << "Extrapolation matrix for triangle is only defined for IntegrationMethod GI_GAUSS_1 and GI_GAUSS_2"<< std::endl;
        }       

    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void CalculateExtrapolationMatrixQuad(Matrix& rExtrapolationMatrix, const GeometryData::IntegrationMethod& rIntegrationMethod)
    {
        if (rIntegrationMethod == GeometryData::IntegrationMethod::GI_GAUSS_1)
        {
            if ((rExtrapolationMatrix.size1() != 4) || (rExtrapolationMatrix.size2() != 1))
            {
                rExtrapolationMatrix.resize(4, 1, false);
            }

            rExtrapolationMatrix(0, 0) = 1;
            rExtrapolationMatrix(1, 0) = 1;
            rExtrapolationMatrix(2, 0) = 1;
            rExtrapolationMatrix(3, 0) = 1;
        }
        else if (rIntegrationMethod == GeometryData::IntegrationMethod::GI_GAUSS_2)
        {
            if ((rExtrapolationMatrix.size1() != 4) || (rExtrapolationMatrix.size2() != 4))
            {
                rExtrapolationMatrix.resize(4, 4, false);
            }
            //Quadrilateral_2d_4
            //GI_GAUSS_2

            rExtrapolationMatrix(0, 0) = 1.8660254037844386; rExtrapolationMatrix(0, 1) = -0.5; rExtrapolationMatrix(0, 2) = 0.13397459621556132; rExtrapolationMatrix(0, 3) = -0.5;
            rExtrapolationMatrix(1, 0) = -0.5; rExtrapolationMatrix(1, 1) = 1.8660254037844386; rExtrapolationMatrix(1, 2) = -0.5; rExtrapolationMatrix(1, 3) = 0.13397459621556132;
            rExtrapolationMatrix(2, 0) = 0.13397459621556132; rExtrapolationMatrix(2, 1) = -0.5; rExtrapolationMatrix(2, 2) = 1.8660254037844386; rExtrapolationMatrix(2, 3) = -0.5;
            rExtrapolationMatrix(3, 0) = -0.5; rExtrapolationMatrix(3, 1) = 0.13397459621556132; rExtrapolationMatrix(3, 2) = -0.5; rExtrapolationMatrix(3, 3) = 1.8660254037844386;
        }
        else
        {
            KRATOS_ERROR << "Extrapolation matrix for quad is only defined for IntegrationMethod GI_GAUSS_1 and GI_GAUSS_2" << std::endl;
        }
    }


    static inline void CalculateExtrapolationMatrixTetra(Matrix& rExtrapolationMatrix, const GeometryData::IntegrationMethod& rIntegrationMethod)
    {
        if (rIntegrationMethod == GeometryData::IntegrationMethod::GI_GAUSS_1)
        {
            if ((rExtrapolationMatrix.size1() != 4) || (rExtrapolationMatrix.size2() != 1))
            {
                rExtrapolationMatrix.resize(4, 1, false);
            }

            rExtrapolationMatrix(0, 0) = 1;
            rExtrapolationMatrix(1, 0) = 1;
            rExtrapolationMatrix(2, 0) = 1;
            rExtrapolationMatrix(3, 0) = 1;
        }

        else if (rIntegrationMethod == GeometryData::IntegrationMethod::GI_GAUSS_2)
        {
            if ((rExtrapolationMatrix.size1() != 4) || (rExtrapolationMatrix.size2() != 4))
            {
                rExtrapolationMatrix.resize(4, 4, false);
            }
            //Tetrahedra_3d_4
            //GI_GAUSS_2
            rExtrapolationMatrix(0, 0) = -0.309016988749894905; rExtrapolationMatrix(0, 1) = -0.3090169887498949046; rExtrapolationMatrix(0, 2) = -0.309016988749894905; rExtrapolationMatrix(0, 3) = 1.9270509662496847144;
            rExtrapolationMatrix(1, 0) = 1.9270509662496847144; rExtrapolationMatrix(1, 1) = -0.30901698874989490481; rExtrapolationMatrix(1, 2) = -0.3090169887498949049; rExtrapolationMatrix(1, 3) = -0.30901698874989490481;
            rExtrapolationMatrix(2, 0) = -0.30901698874989490473; rExtrapolationMatrix(2, 1) = 1.9270509662496847143; rExtrapolationMatrix(2, 2) = -0.3090169887498949049; rExtrapolationMatrix(2, 3) = -0.30901698874989490481;
            rExtrapolationMatrix(3, 0) = -0.3090169887498949048; rExtrapolationMatrix(3, 1) = -0.30901698874989490471; rExtrapolationMatrix(3, 2) = 1.9270509662496847143; rExtrapolationMatrix(3, 3) = -0.30901698874989490481;
        }
        else
        {
            KRATOS_ERROR << "Extrapolation matrix for tetrahedral is only defined for IntegrationMethod GI_GAUSS_1 and GI_GAUSS_2" << std::endl;
        }

    }


    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void CalculateExtrapolationMatrixHexa(Matrix& rExtrapolationMatrix, const GeometryData::IntegrationMethod& rIntegrationMethod)
    {
        if (rIntegrationMethod == GeometryData::IntegrationMethod::GI_GAUSS_1)
        {
            if ((rExtrapolationMatrix.size1() != 8) || (rExtrapolationMatrix.size2() != 1))
            {
                rExtrapolationMatrix.resize(8, 1, false);
            }
            rExtrapolationMatrix(0, 0) = 1;
            rExtrapolationMatrix(1, 0) = 1;
            rExtrapolationMatrix(2, 0) = 1;
            rExtrapolationMatrix(3, 0) = 1;
            rExtrapolationMatrix(4, 0) = 1;
            rExtrapolationMatrix(5, 0) = 1;
            rExtrapolationMatrix(6, 0) = 1;
            rExtrapolationMatrix(7, 0) = 1;
        }

        else if (rIntegrationMethod == GeometryData::IntegrationMethod::GI_GAUSS_2)
        {
            if ((rExtrapolationMatrix.size1() != 8) || (rExtrapolationMatrix.size2() != 8))
            {
                rExtrapolationMatrix.resize(8, 8, false);
            }
            //Hexahedra_3d_8
            //GI_GAUSS_2

            rExtrapolationMatrix(0, 0) = 2.549038105676658; rExtrapolationMatrix(0, 1) = -0.6830127018922192; rExtrapolationMatrix(0, 2) = 0.18301270189221927; rExtrapolationMatrix(0, 3) = -0.6830127018922192;
            rExtrapolationMatrix(0, 4) = -0.6830127018922192; rExtrapolationMatrix(0, 5) = 0.18301270189221927; rExtrapolationMatrix(0, 6) = -0.04903810567665795; rExtrapolationMatrix(0, 7) = 0.18301270189221927;

            rExtrapolationMatrix(1, 0) = -0.6830127018922192; rExtrapolationMatrix(1, 1) = 2.549038105676658; rExtrapolationMatrix(1, 2) = -0.6830127018922192; rExtrapolationMatrix(1, 3) = 0.18301270189221927;
            rExtrapolationMatrix(1, 4) = 0.18301270189221927; rExtrapolationMatrix(1, 5) = -0.6830127018922192; rExtrapolationMatrix(1, 6) = 0.18301270189221927; rExtrapolationMatrix(1, 7) = -0.04903810567665795;

            rExtrapolationMatrix(2, 0) = 0.18301270189221927; rExtrapolationMatrix(2, 1) = -0.6830127018922192; rExtrapolationMatrix(2, 2) = 2.549038105676658; rExtrapolationMatrix(2, 3) = -0.6830127018922192;
            rExtrapolationMatrix(2, 4) = -0.04903810567665795; rExtrapolationMatrix(2, 5) = 0.18301270189221927; rExtrapolationMatrix(2, 6) = -0.6830127018922192; rExtrapolationMatrix(2, 7) = 0.18301270189221927;

            rExtrapolationMatrix(3, 0) = -0.6830127018922192; rExtrapolationMatrix(3, 1) = 0.18301270189221927; rExtrapolationMatrix(3, 2) = -0.6830127018922192; rExtrapolationMatrix(3, 3) = 2.549038105676658;
            rExtrapolationMatrix(3, 4) = 0.18301270189221927; rExtrapolationMatrix(3, 5) = -0.04903810567665795; rExtrapolationMatrix(3, 6) = 0.18301270189221927; rExtrapolationMatrix(3, 7) = -0.6830127018922192;

            rExtrapolationMatrix(4, 0) = -0.6830127018922192; rExtrapolationMatrix(4, 1) = 0.18301270189221927; rExtrapolationMatrix(4, 2) = -0.04903810567665795; rExtrapolationMatrix(4, 3) = 0.18301270189221927;
            rExtrapolationMatrix(4, 4) = 2.549038105676658; rExtrapolationMatrix(4, 5) = -0.6830127018922192; rExtrapolationMatrix(4, 6) = 0.18301270189221927; rExtrapolationMatrix(4, 7) = -0.6830127018922192;

            rExtrapolationMatrix(5, 0) = 0.18301270189221927; rExtrapolationMatrix(5, 1) = -0.6830127018922192; rExtrapolationMatrix(5, 2) = 0.18301270189221927; rExtrapolationMatrix(5, 3) = -0.04903810567665795;
            rExtrapolationMatrix(5, 4) = -0.6830127018922192; rExtrapolationMatrix(5, 5) = 2.549038105676658; rExtrapolationMatrix(5, 6) = -0.6830127018922192; rExtrapolationMatrix(5, 7) = 0.18301270189221927;

            rExtrapolationMatrix(6, 0) = -0.04903810567665795; rExtrapolationMatrix(6, 1) = 0.18301270189221927; rExtrapolationMatrix(6, 2) = -0.6830127018922192; rExtrapolationMatrix(6, 3) = 0.18301270189221927;
            rExtrapolationMatrix(6, 4) = 0.18301270189221927; rExtrapolationMatrix(6, 5) = -0.6830127018922192; rExtrapolationMatrix(6, 6) = 2.549038105676658; rExtrapolationMatrix(6, 7) = -0.6830127018922192;

            rExtrapolationMatrix(7, 0) = 0.18301270189221927; rExtrapolationMatrix(7, 1) = -0.04903810567665795; rExtrapolationMatrix(7, 2) = 0.18301270189221927; rExtrapolationMatrix(7, 3) = -0.6830127018922192;
            rExtrapolationMatrix(7, 4) = -0.6830127018922192; rExtrapolationMatrix(7, 5) = 0.18301270189221927; rExtrapolationMatrix(7, 6) = -0.6830127018922192; rExtrapolationMatrix(7, 7) = 2.549038105676658;
        }
        else
        {
            KRATOS_ERROR << "Extrapolation matrix for hexahedral is only defined for IntegrationMethod GI_GAUSS_1 and GI_GAUSS_2" << std::endl;
        }
     }


    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static Vector CalculateNodalHydraulicHeadFromWaterPressures(const GeometryType& rGeom, const Properties& rProp)
    {
        const auto NumericalLimit = std::numeric_limits<double>::epsilon();
    	//Defining necessary variables

        Vector nodal_hydraulic_heads(rGeom.PointsNumber());
        for (unsigned int node = 0; node < rGeom.PointsNumber(); ++node) {
            array_1d<double, 3> node_volume_acceleration;
            noalias(node_volume_acceleration) = rGeom[node].FastGetSolutionStepValue(VOLUME_ACCELERATION, 0);
            const double g = norm_2(node_volume_acceleration);
            if (g > NumericalLimit) {
                const double FluidWeight = g * rProp[DENSITY_WATER];

                array_1d<double, 3> NodeCoordinates;
                noalias(NodeCoordinates) = rGeom[node].Coordinates();
                array_1d<double, 3> NodeVolumeAccelerationUnitVector;
                noalias(NodeVolumeAccelerationUnitVector) = node_volume_acceleration / g;

                const double WaterPressure = rGeom[node].FastGetSolutionStepValue(WATER_PRESSURE);
                nodal_hydraulic_heads[node] = -inner_prod(NodeCoordinates, NodeVolumeAccelerationUnitVector)
                    - PORE_PRESSURE_SIGN_FACTOR * WaterPressure / FluidWeight;
            }
            else {
                nodal_hydraulic_heads[node] = 0.0;
            }
        }
        return nodal_hydraulic_heads;
    }

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

}; /* Class GeoElementUtilities*/
} /* namespace Kratos.*/

#endif /* KRATOS_GEO_ELEMENT_UTILITIES defined */
