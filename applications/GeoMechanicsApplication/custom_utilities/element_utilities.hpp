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

    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;

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
    template< unsigned int TDim>
    static inline void AssembleDensityMatrix(BoundedMatrix<double,TDim+1, TDim+1> &DensityMatrix,
                                             const double &Density)
    {
        for (unsigned int idim = 0; idim < TDim; ++idim)
        {
            for (unsigned int jdim = 0; jdim < TDim; ++jdim)
            {
                DensityMatrix(idim, jdim) = Density;
            }
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void AssembleDensityMatrix(Matrix &DensityMatrix,
                                             const double &Density)
    {
        for (unsigned int idim = 0; idim < DensityMatrix.size1(); ++idim)
        {
            for (unsigned int jdim = 0; jdim < DensityMatrix.size2(); ++jdim)
            {
                DensityMatrix(idim, jdim) = Density;
            }
        }
    }


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void AssembleUBlockMatrix(Matrix &rLeftHandSideMatrix,
                                            const BoundedMatrix<double,TDim*TNumNodes, TDim*TNumNodes> &UBlockMatrix)
    {
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            const unsigned int Global_i = i * (TDim + 1);
            const unsigned int Local_i  = i * TDim;

            for (unsigned int j = 0; j < TNumNodes; j++)
            {
                const unsigned int Global_j = j * (TDim + 1);
                const unsigned int Local_j  = j * TDim;

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
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            const unsigned int Global_i = i * (TDim + 1);
            const unsigned int Local_i  = i * TDim;

            for (unsigned int j = 0; j < TNumNodes; j++)
            {
                const unsigned int Global_j = j * (TDim + 1);
                const unsigned int Local_j  = j * TDim;

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
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            const unsigned int Global_i = i * (TDim + 1);
            const unsigned int Local_i = i * TDim;

            for (unsigned int j = 0; j < TNumNodes; j++)
            {
                const unsigned int Global_j = j * (TDim + 1) + TDim;

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
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            const unsigned int Global_i = i * (TDim + 1) + TDim;

            for (unsigned int j = 0; j < TNumNodes; j++)
            {
                const unsigned int Global_j = j * (TDim + 1);
                const unsigned int Local_j = j * TDim;

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
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            const unsigned int Global_i = i * (TDim + 1) + TDim;

            for (unsigned int j = 0; j < TNumNodes; j++)
            {
                const unsigned int Global_j = j * (TDim + 1) + TDim;

                rLeftHandSideMatrix(Global_i,Global_j) += PBlockMatrix(i,j);
            }
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void AssembleUBlockVector(Vector& rRightHandSideVector,
                                            const array_1d<double,TDim*TNumNodes>& UBlockVector)
    {
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            const unsigned int Global_i = i * (TDim + 1);
            const unsigned int Local_i  = i * TDim;

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
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            const unsigned int Global_i = i * (TDim + 1) + TDim;

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

    /**
     * Calculates the radius of axisymmetry
     * @param N: The Gauss Point shape function
     * @param Geom: The geometry studied
     * @return Radius: The radius of axisymmetry
     */

    static inline double CalculateRadius(const Vector N, const GeometryType& Geom)
    {
        double Radius = 0.0;

        for (unsigned int iNode = 0; iNode < Geom.size(); iNode++)
        {
            // Displacement from the reference to the current configuration
            const array_1d<double, 3 >& CurrentPosition = Geom[iNode].Coordinates();
            Radius += CurrentPosition[0] * N[iNode];
        }

        return Radius;
    }


}; /* Class GeoElementUtilities*/
} /* namespace Kratos.*/

#endif /* KRATOS_GEO_ELEMENT_UTILITIES defined */
