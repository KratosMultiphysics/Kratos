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

#if !defined(KRATOS_CONDITION_UTILITIES )
#define  KRATOS_CONDITION_UTILITIES

// Project includes
#include "includes/element.h"
#include "geo_mechanics_application_variables.h"


namespace Kratos
{

class ConditionUtilities
{

public:

    //----------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void CalculateNuMatrix(BoundedMatrix<double, TDim, TDim*TNumNodes>& rNu,
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

    //----------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void InterpolateVariableWithComponents(array_1d<double,TDim>& rVector,
                                                         const Matrix& Ncontainer,
                                                         const array_1d<double,TDim*TNumNodes>& VariableWithComponents,
                                                         const unsigned int& GPoint)
    {
        noalias(rVector) = ZeroVector(TDim);

        unsigned int index = 0;
        for (unsigned int i=0; i<TNumNodes; i++)
        {
            for (unsigned int idim=0; idim<TDim; idim++)
            {
                rVector[idim] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            }
        }
    }

    //----------------------------------------------------------------------------------------
    static inline void GetDisplacementsVector(array_1d<double,4>& rDisplacementVector,
                                              const Element::GeometryType& Geom)
    {
        //Line_2d_2
        array_1d<double,3> DisplacementAux;
        unsigned int index = 0;
        for(unsigned int i=0; i<2; i++)
        {
            noalias(DisplacementAux) = Geom[i].FastGetSolutionStepValue(DISPLACEMENT);
            rDisplacementVector[index++] = DisplacementAux[0];
            rDisplacementVector[index++] = DisplacementAux[1];
        }
    }

    //----------------------------------------------------------------------------------------
    static inline void GetDisplacementsVector(array_1d<double,12>& rDisplacementVector,
                                              const Element::GeometryType& Geom)
    {
        //Quadrilateral_3d_4
        array_1d<double,3> DisplacementAux;
        unsigned int index = 0;
        for(unsigned int i=0; i<4; i++)
        {
            noalias(DisplacementAux) = Geom[i].FastGetSolutionStepValue(DISPLACEMENT);
            rDisplacementVector[index++] = DisplacementAux[0];
            rDisplacementVector[index++] = DisplacementAux[1];
            rDisplacementVector[index++] = DisplacementAux[2];
        }
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template<unsigned int TNumNodes >
    static inline void GetFaceLoadVector(array_1d<double,3*TNumNodes>& rFaceLoadVector,
                                         const Element::GeometryType& Geom)
    {

        // for 3D geometry
        const unsigned int TDim = 3;
        array_1d<double,3> FaceLoadAux;
        unsigned int index = 0;
        for (unsigned int i=0; i<TNumNodes; i++)
        {
            noalias(FaceLoadAux) = Geom[i].FastGetSolutionStepValue(SURFACE_LOAD);
            for (unsigned int idim=0; idim<TDim; idim++)
            {
                rFaceLoadVector[index++] = FaceLoadAux[idim];
            }
        }
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template<unsigned int TNumNodes >
    static inline void GetFaceLoadVector(array_1d<double,2*TNumNodes>& rFaceLoadVector,
                                         const Element::GeometryType& Geom)
    {

        // for 2D geometry
        const unsigned int TDim = 2;
        array_1d<double,3> FaceLoadAux;
        unsigned int index = 0;
        for (unsigned int i=0; i<TNumNodes; i++)
        {
            noalias(FaceLoadAux) = Geom[i].FastGetSolutionStepValue(LINE_LOAD);
            for (unsigned int idim=0; idim < TDim; idim++)
            {
                rFaceLoadVector[index++] = FaceLoadAux[idim];
            }
        }
    }

    //----------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void AssembleUBlockVector(Vector& rRightHandSideVector,
                                            const array_1d<double, TDim*TNumNodes>& UBlockVector)
    {
        unsigned int Global_i, Local_i;

        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            Global_i = i * (TDim + 1);
            Local_i  = i * TDim;
            for (unsigned int idim = 0; idim < TDim; ++idim)
            {
              rRightHandSideVector[Global_i + idim] += UBlockVector[Local_i + idim];
            }
        }
    }

    //----------------------------------------------------------------------------------------
    template< class TVectorType >
    static inline void AssemblePBlockVector(Vector& rRightHandSideVector,
                                            const TVectorType& PBlockVector,
                                            const unsigned int& Dim,
                                            const unsigned int& NumNodes)
    {
        unsigned int Global_i;

        for(unsigned int i = 0; i < NumNodes; i++)
        {
            Global_i = i * (Dim + 1) + Dim;

            rRightHandSideVector[Global_i] += PBlockVector[i];
        }
    }

    //----------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void AssembleUPMatrix(Matrix& rLeftHandSideMatrix,
                                        const BoundedMatrix<double,TDim*TNumNodes,TNumNodes>& UPBlockMatrix)
    {
        //Quadrilateral_3d_4
        unsigned int Global_i, Global_j, Local_i;

        for(unsigned int i = 0; i < TNumNodes; i++)
        {
            Global_i = i * (TDim + 1);
            Local_i = i * TDim;

            for(unsigned int j = 0; j < TNumNodes; j++)
            {
                Global_j = j * (TDim + 1) + TDim;
                for (unsigned int idim = 0; idim < TDim; ++idim)
                {
                   rLeftHandSideMatrix(Global_i+idim, Global_j) += UPBlockMatrix(Local_i+idim, j);
                }
            }
        }
    }

    //----------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void AssemblePUMatrix(Matrix& rLeftHandSideMatrix,
                                        const BoundedMatrix<double,TNumNodes,TDim*TNumNodes>& PUBlockMatrix)
    {
        //Quadrilateral_3d_4
        unsigned int Global_i, Global_j, Local_j;

        for(unsigned int i = 0; i < TNumNodes; i++)
        {
            Global_i = i * (TDim + 1) + TDim;

            for(unsigned int j = 0; j < TNumNodes; j++)
            {
                Global_j = j * (TDim + 1);
                Local_j = j * TDim;
                for (unsigned int idim = 0; idim < TDim; ++idim)
                {
                    rLeftHandSideMatrix(Global_i, Global_j+idim) += PUBlockMatrix(i, Local_j+idim);
                }
            }
        }
    }

}; /* Class ConditionUtilities*/
} /* namespace Kratos.*/

#endif /* KRATOS_CONDITION_UTILITIES defined */
