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

// System includes
//#include <cmath>

// Project includes
#include "utilities/math_utils.h"
#include "includes/element.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class GeoElementUtilities
{

typedef std::size_t IndexType;

public:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void CalculateNuMatrix(BoundedMatrix<double,TDim,TDim*TNumNodes>& rNu,
                                         const Matrix& NContainer,
                                         const unsigned int& GPoint)
    {
        for (unsigned int i=0; i < TDim; ++i)
        {
            unsigned int index = i - TDim;
            for (unsigned int j=0; j < TNumNodes; ++j)
            {
                index += TDim;
                rNu(i, index) = NContainer(GPoint, j);
            }
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void CalculateNuElementMatrix(BoundedMatrix<double, (TDim+1), TNumNodes*(TDim+1)>& rNut,
                                                const Matrix& NContainer,
                                                const unsigned int& GPoint)
    {
        const unsigned int offset = (TDim+1);

        for (unsigned int i=0; i < TDim; ++i)
        {
            unsigned int index = i - offset;
            for (unsigned int j=0; j < TNumNodes; ++j)
            {
                index += offset;
                rNut(i, index) = NContainer(GPoint, j);
            }
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void InterpolateVariableWithComponents(array_1d<double, TDim>& rVector,
                                                         const Matrix& NContainer,
                                                         const array_1d<double, TDim*TNumNodes>& VariableWithComponents,
                                                         const unsigned int& GPoint)
    {
        noalias(rVector) = ZeroVector(TDim);

        unsigned int index = 0;
        for (unsigned int i=0; i<TNumNodes; i++)
        {
            for (unsigned int j=0; j<TDim; j++)
            {
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
        for (unsigned int i=0; i < TNumNodes; i++)
        {
            noalias(NodalVariableAux) = Geom[i].FastGetSolutionStepValue(Variable, SolutionStepIndex);
            for (unsigned int j=0; j < TDim; j++)
            {
                rNodalVariableVector[index++] = NodalVariableAux[j];
            }
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void CalculatePermeabilityMatrix(BoundedMatrix<double,2,2>& rPermeabilityMatrix,
                                                   const Element::PropertiesType& Prop)
    {
        //2D
        rPermeabilityMatrix(0,0) = Prop[PERMEABILITY_XX];
        rPermeabilityMatrix(1,1) = Prop[PERMEABILITY_YY];

        rPermeabilityMatrix(0,1) = Prop[PERMEABILITY_XY];
        rPermeabilityMatrix(1,0) = rPermeabilityMatrix(0,1);
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void CalculatePermeabilityMatrix(BoundedMatrix<double,3,3>& rPermeabilityMatrix,
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
        KRATOS_TRY;

        const double numerical_limit = std::numeric_limits<double>::epsilon();

        InputMatrixDet = InputMatrix(0,0)*InputMatrix(1,1)-InputMatrix(0,1)*InputMatrix(1,0);

        if (InputMatrixDet < numerical_limit)
        {
            KRATOS_ERROR << "determinant zero or negative" << std::endl;
        }

        rInvertedMatrix(0,0) =  InputMatrix(1,1)/InputMatrixDet;
        rInvertedMatrix(0,1) = -InputMatrix(0,1)/InputMatrixDet;
        rInvertedMatrix(1,0) = -InputMatrix(1,0)/InputMatrixDet;
        rInvertedMatrix(1,1) =  InputMatrix(0,0)/InputMatrixDet;

        KRATOS_CATCH("");
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void InvertMatrix2(BoundedMatrix<double,2,2>& rInvertedMatrix,
                                     const BoundedMatrix<double,2,2>& InputMatrix)
    {
        KRATOS_TRY;

        double InputMatrixDet;

        InvertMatrix2(rInvertedMatrix, InputMatrix, InputMatrixDet);

        KRATOS_CATCH("");
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void AssembleUBlockMatrix(Matrix &rLeftHandSideMatrix,
                                            const BoundedMatrix<double,TDim*TNumNodes, TDim*TNumNodes> &UBlockMatrix)
    {
        unsigned int Global_i, Global_j, Local_i, Local_j;

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            Global_i = i * (TDim + 1);
            Local_i  = i * TDim;

            for (unsigned int j = 0; j < TNumNodes; j++)
            {
                Global_j = j * (TDim + 1);
                Local_j  = j * TDim;

                for (unsigned int idim = 0; idim < TDim; ++idim)
                {
                    for (unsigned int jdim = 0; jdim < TDim; ++jdim)
                    {
                        rLeftHandSideMatrix(Global_i+idim, Global_j+jdim) += UBlockMatrix(Local_i+idim, Local_j+jdim);
                    }
                }
            }
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void AssembleUBlockMatrix(Matrix &rLeftHandSideMatrix,
                                            const Matrix &UBlockMatrix,
                                            const double &TNumNodes,
                                            const double &TDim)
    {
        unsigned int Global_i, Global_j, Local_i, Local_j;

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            Global_i = i * (TDim + 1);
            Local_i  = i * TDim;

            for (unsigned int j = 0; j < TNumNodes; j++)
            {
                Global_j = j * (TDim + 1);
                Local_j  = j * TDim;

                for (unsigned int idim = 0; idim < TDim; ++idim)
                {
                    for (unsigned int jdim = 0; jdim < TDim; ++jdim)
                    {
                        rLeftHandSideMatrix(Global_i+idim, Global_j+jdim) += UBlockMatrix(Local_i+idim, Local_j+jdim);
                    }
                }
            }
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void AssembleUPBlockMatrix(Matrix& rLeftHandSideMatrix,
                                             const BoundedMatrix<double,TDim*TNumNodes,TNumNodes>& UPBlockMatrix)
    {
        unsigned int Global_i, Global_j, Local_i;

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            Global_i = i * (TDim + 1);
            Local_i = i * TDim;

            for (unsigned int j = 0; j < TNumNodes; j++)
            {
                Global_j = j * (TDim + 1) + TDim;

                for (unsigned int dim = 0; dim < TDim; ++dim)
                {
                    rLeftHandSideMatrix(Global_i + dim, Global_j)  += UPBlockMatrix(Local_i + dim, j);
                }
            }
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void AssemblePUBlockMatrix(Matrix& rLeftHandSideMatrix,
                                             const BoundedMatrix<double,TNumNodes,TNumNodes*TDim>& PUBlockMatrix)
    {
        unsigned int Global_i, Global_j, Local_j;

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            Global_i = i * (TDim + 1) + TDim;

            for (unsigned int j = 0; j < TNumNodes; j++)
            {
                Global_j = j * (TDim + 1);
                Local_j = j * TDim;

                for (unsigned int dim = 0; dim < TDim; ++dim)
                {
                    rLeftHandSideMatrix(Global_i, Global_j+dim) += PUBlockMatrix(i, Local_j+dim);
                }
            }
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void AssemblePBlockMatrix(Matrix& rLeftHandSideMatrix,
                                            const BoundedMatrix<double,TNumNodes,TNumNodes> &PBlockMatrix)
    {
        unsigned int Global_i, Global_j;

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            Global_i = i * (TDim + 1) + TDim;

            for (unsigned int j = 0; j < TNumNodes; j++)
            {
                Global_j = j * (TDim + 1) + TDim;

                rLeftHandSideMatrix(Global_i,Global_j) += PBlockMatrix(i,j);
            }
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void AssembleUBlockVector(Vector& rRightHandSideVector,
                                            const array_1d<double,TDim*TNumNodes>& UBlockVector)
    {
        unsigned int Global_i, Local_i;

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            Global_i = i * (TDim + 1);
            Local_i  = i * TDim;

            for (unsigned int dim = 0; dim < TDim; ++dim)
            {
                rRightHandSideVector[Global_i + dim] += UBlockVector[Local_i + dim];
            }
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void AssemblePBlockVector(Vector& rRightHandSideVector,
                                            const array_1d<double, TNumNodes> &PBlockVector)
    {
        unsigned int Global_i;

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            Global_i = i * (TDim + 1) + TDim;

            rRightHandSideVector[Global_i] += PBlockVector[i];
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void CalculateShapeFunctionsNodesGradients(BoundedMatrix<double,2,2>& DN_DXContainer)
    {
        //Line 2-noded
        const unsigned int NumNodes = 2;
        const std::vector<double> Xi{-1.0, 1.0};
        noalias(DN_DXContainer) = ZeroMatrix(NumNodes,NumNodes);

        for (unsigned int integrationPoint = 0; integrationPoint < NumNodes; ++integrationPoint)
        {
            DN_DXContainer(integrationPoint,0) = - 0.5;
            DN_DXContainer(integrationPoint,1) =   0.5;
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void CalculateShapeFunctionsNodesGradients(BoundedMatrix<double,3,3>& DN_DXContainer)
    {
        //Line 3-noded
        const unsigned int NumNodes = 3;
        const std::vector<double> Xi{-1.0, 0.0, 1.0};

        noalias(DN_DXContainer) = ZeroMatrix(NumNodes,NumNodes);

        for (unsigned int integrationPoint = 0; integrationPoint < NumNodes; ++integrationPoint)
        {
            DN_DXContainer(integrationPoint,0) =  Xi[integrationPoint] - 0.5;
            DN_DXContainer(integrationPoint,1) = -Xi[integrationPoint] * 2.0;
            DN_DXContainer(integrationPoint,2) =  Xi[integrationPoint] + 0.5;
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    /// CalculateExtrapolationMatrix
    /// The matrix contains the shape functions at each GP evaluated at each node.
    /// Rows: nodes
    /// Columns: GP

    static inline void Calculate2DExtrapolationMatrix(BoundedMatrix<double,3,3>& rExtrapolationMatrix)
    {
        //Triangle_2d_3
        //GI_GAUSS_2

        rExtrapolationMatrix(0,0) = 1.6666666666666666666; rExtrapolationMatrix(0,1) = -0.33333333333333333333; rExtrapolationMatrix(0,2) = -0.33333333333333333333;
        rExtrapolationMatrix(1,0) = -0.33333333333333333333; rExtrapolationMatrix(1,1) = 1.6666666666666666666; rExtrapolationMatrix(1,2) = -0.33333333333333333333;
        rExtrapolationMatrix(2,0) = -0.33333333333333333333; rExtrapolationMatrix(2,1) = -0.33333333333333333333; rExtrapolationMatrix(2,2) = 1.6666666666666666666;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void Calculate2DExtrapolationMatrix(BoundedMatrix<double,4,4>& rExtrapolationMatrix)
    {
        //Quadrilateral_2d_4
        //GI_GAUSS_2

        rExtrapolationMatrix(0,0) = 1.8660254037844386; rExtrapolationMatrix(0,1) = -0.5; rExtrapolationMatrix(0,2) = 0.13397459621556132; rExtrapolationMatrix(0,3) = -0.5;
        rExtrapolationMatrix(1,0) = -0.5; rExtrapolationMatrix(1,1) = 1.8660254037844386; rExtrapolationMatrix(1,2) = -0.5; rExtrapolationMatrix(1,3) = 0.13397459621556132;
        rExtrapolationMatrix(2,0) = 0.13397459621556132; rExtrapolationMatrix(2,1) = -0.5; rExtrapolationMatrix(2,2) = 1.8660254037844386; rExtrapolationMatrix(2,3) = -0.5;
        rExtrapolationMatrix(3,0) = -0.5; rExtrapolationMatrix(3,1) = 0.13397459621556132; rExtrapolationMatrix(3,2) = -0.5; rExtrapolationMatrix(3,3) = 1.8660254037844386;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void Calculate3DExtrapolationMatrix(BoundedMatrix<double,4,4>& rExtrapolationMatrix)
    {
        //Tetrahedra_3d_4
        //GI_GAUSS_2

        rExtrapolationMatrix(0,0) = -0.309016988749894905; rExtrapolationMatrix(0,1) = -0.3090169887498949046; rExtrapolationMatrix(0,2) = -0.309016988749894905; rExtrapolationMatrix(0,3) = 1.9270509662496847144;
        rExtrapolationMatrix(1,0) = 1.9270509662496847144; rExtrapolationMatrix(1,1) = -0.30901698874989490481; rExtrapolationMatrix(1,2) = -0.3090169887498949049; rExtrapolationMatrix(1,3) = -0.30901698874989490481;
        rExtrapolationMatrix(2,0) = -0.30901698874989490473; rExtrapolationMatrix(2,1) = 1.9270509662496847143; rExtrapolationMatrix(2,2) = -0.3090169887498949049; rExtrapolationMatrix(2,3) = -0.30901698874989490481;
        rExtrapolationMatrix(3,0) = -0.3090169887498949048; rExtrapolationMatrix(3,1) = -0.30901698874989490471; rExtrapolationMatrix(3,2) = 1.9270509662496847143; rExtrapolationMatrix(3,3) = -0.30901698874989490481;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void Calculate3DExtrapolationMatrix(BoundedMatrix<double,8,8>& rExtrapolationMatrix)
    {
        //Hexahedra_3d_8
        //GI_GAUSS_2

        rExtrapolationMatrix(0,0) = 2.549038105676658; rExtrapolationMatrix(0,1) = -0.6830127018922192; rExtrapolationMatrix(0,2) = 0.18301270189221927; rExtrapolationMatrix(0,3) = -0.6830127018922192;
        rExtrapolationMatrix(0,4) = -0.6830127018922192; rExtrapolationMatrix(0,5) = 0.18301270189221927; rExtrapolationMatrix(0,6) = -0.04903810567665795; rExtrapolationMatrix(0,7) = 0.18301270189221927;

        rExtrapolationMatrix(1,0) = -0.6830127018922192; rExtrapolationMatrix(1,1) = 2.549038105676658; rExtrapolationMatrix(1,2) = -0.6830127018922192; rExtrapolationMatrix(1,3) = 0.18301270189221927;
        rExtrapolationMatrix(1,4) = 0.18301270189221927; rExtrapolationMatrix(1,5) = -0.6830127018922192; rExtrapolationMatrix(1,6) = 0.18301270189221927; rExtrapolationMatrix(1,7) = -0.04903810567665795;

        rExtrapolationMatrix(2,0) = 0.18301270189221927; rExtrapolationMatrix(2,1) = -0.6830127018922192; rExtrapolationMatrix(2,2) = 2.549038105676658; rExtrapolationMatrix(2,3) = -0.6830127018922192;
        rExtrapolationMatrix(2,4) = -0.04903810567665795; rExtrapolationMatrix(2,5) = 0.18301270189221927; rExtrapolationMatrix(2,6) = -0.6830127018922192; rExtrapolationMatrix(2,7) = 0.18301270189221927;

        rExtrapolationMatrix(3,0) = -0.6830127018922192; rExtrapolationMatrix(3,1) = 0.18301270189221927; rExtrapolationMatrix(3,2) = -0.6830127018922192; rExtrapolationMatrix(3,3) = 2.549038105676658;
        rExtrapolationMatrix(3,4) = 0.18301270189221927; rExtrapolationMatrix(3,5) = -0.04903810567665795; rExtrapolationMatrix(3,6) = 0.18301270189221927; rExtrapolationMatrix(3,7) = -0.6830127018922192;

        rExtrapolationMatrix(4,0) = -0.6830127018922192; rExtrapolationMatrix(4,1) = 0.18301270189221927; rExtrapolationMatrix(4,2) = -0.04903810567665795; rExtrapolationMatrix(4,3) = 0.18301270189221927;
        rExtrapolationMatrix(4,4) = 2.549038105676658; rExtrapolationMatrix(4,5) = -0.6830127018922192; rExtrapolationMatrix(4,6) = 0.18301270189221927; rExtrapolationMatrix(4,7) = -0.6830127018922192;

        rExtrapolationMatrix(5,0) = 0.18301270189221927; rExtrapolationMatrix(5,1) = -0.6830127018922192; rExtrapolationMatrix(5,2) = 0.18301270189221927; rExtrapolationMatrix(5,3) = -0.04903810567665795;
        rExtrapolationMatrix(5,4) = -0.6830127018922192; rExtrapolationMatrix(5,5) = 2.549038105676658; rExtrapolationMatrix(5,6) = -0.6830127018922192; rExtrapolationMatrix(5,7) = 0.18301270189221927;

        rExtrapolationMatrix(6,0) = -0.04903810567665795; rExtrapolationMatrix(6,1) = 0.18301270189221927; rExtrapolationMatrix(6,2) = -0.6830127018922192; rExtrapolationMatrix(6,3) = 0.18301270189221927;
        rExtrapolationMatrix(6,4) = 0.18301270189221927; rExtrapolationMatrix(6,5) = -0.6830127018922192; rExtrapolationMatrix(6,6) = 2.549038105676658; rExtrapolationMatrix(6,7) = -0.6830127018922192;

        rExtrapolationMatrix(7,0) = 0.18301270189221927; rExtrapolationMatrix(7,1) = -0.04903810567665795; rExtrapolationMatrix(7,2) = 0.18301270189221927; rExtrapolationMatrix(7,3) = -0.6830127018922192;
        rExtrapolationMatrix(7,4) = -0.6830127018922192; rExtrapolationMatrix(7,5) = 0.18301270189221927; rExtrapolationMatrix(7,6) = -0.6830127018922192; rExtrapolationMatrix(7,7) = 2.549038105676658;
    }

}; /* Class GeoElementUtilities*/
} /* namespace Kratos.*/

#endif /* KRATOS_GEO_ELEMENT_UTILITIES defined */
